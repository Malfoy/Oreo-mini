extern crate rlimit;
use clap::Parser;
use bio::io::{fasta::{self, FastaRead}, fastq::{self, FastqRead}};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression as GzipCompression;
use zstd::stream::{Encoder, Decoder}; // Keep Decoder for input, Encoder for final output
use zstd::stream::raw::CParameter;
use num_cpus;
use std::sync::{Arc, Mutex};
use std::fs::{File, self};
use std::io::{self, BufReader, BufWriter, BufRead, Read, Write};
use std::time::Instant;
use rayon::prelude::*;
use rayon::{ThreadPoolBuilder, ThreadPool};
use std::path::{Path, PathBuf};
use indicatif::{ProgressBar, ProgressStyle};
use nthash::*;
use ahash::AHashMap;


/// Command-line arguments.
#[derive(Parser, Debug, Clone)]
struct Args {
    /// Input FASTA/FASTQ file (.gz or .zst compressed are supported)
    #[arg(short, long)]
    input: String,

    /// Output directory (will be created if needed)
    #[arg(short, long, default_value = "oreo-wdir")]
    output: String,

    /// Fingerprint bit-length P. The algorithm uses an array of length P and produces a P-bit fingerprint.
    /// This results in 2^P partition files. (Default: 8)
    #[arg(short, long, default_value = "8")]
    p: usize,

    /// k-mer length.
    #[arg(short, long, default_value = "21")]
    k: usize,

    // --- REMOVED: compression_level is no longer needed for intermediate files ---
    // #[arg(long, default_value = "-4")]
    // compression_level: i32,
    // --- END REMOVED ---

    /// Final compression algorithm ("zstd" or "gzip"). (Default: "zstd")
    #[arg(long, default_value = "zstd")]
    final_compression: String,

    /// Final compression level. For zstd this is typically in the range [1, 22]; for gzip in [0, 9].
    #[arg(long, default_value = "4")]
    final_compression_level: i32,

    /// Enable reverse complement sensitivity (uses canonical k-mers for fingerprinting and Bloom filter)
    #[arg(long, default_value_t = true)]
    rc_sensitivity: bool,

    /// Number of loops for reverse complement orientation update (0 disables)
    #[arg(long, default_value = "1")]
    rc_compression_loop: i32,

    /// Number of threads (0 uses all available cores)
    #[arg(short, long, default_value = "0")]
    thread: usize,

    // --- Optional Counting Filter ---
    /// Enable the counting filter strategy. If disabled, all subsampled k-mers are used for fingerprinting.
    #[arg(long, default_value_t = true)]
    use_counting_filter: bool,

    // --- Counting Filter Arguments (only used if --use-counting-filter is set) ---
    /// Total number of 2-bit counters (approximate size in bits = filter_counters * 2). Affects memory usage per shard.
    #[arg(long, default_value = "10000000000")]
    filter_counters: u64,

    // --- Subsampling Argument ---
    /// Number of trailing zeros required in k-mer hash for subsampling (0=no subsampling, 1=1/2, 2=1/4, etc.).
    #[arg(long, default_value = "3")] // Default to 1 (hash ends in 0)
    trailing_zeros: u32,
}


#[derive(Clone)]
enum Record {
    Fasta { id: String, seq: String },
    Fastq { id: String, seq: String, qual: String },
}

/// Open input file. Supports .gz and .zst.
fn open_input(path: &str) -> Box<dyn Read> {
    let file = File::open(path).expect("Cannot open input file");
    if path.ends_with(".gz") {
        Box::new(GzDecoder::new(file))
    } else if path.ends_with(".zst") {
        let decoder = Decoder::new(file).expect("Cannot create zstd decoder for input");
        Box::new(decoder)
    } else {
        Box::new(file)
    }
}

// --- Homopolymer Compression ---
fn homopolymer_compress(seq: &str) -> String {
    let mut compressed = String::with_capacity(seq.len());
    let mut last_char: Option<char> = None;
    for current_char in seq.chars() {
        if last_char.map(|c| c.to_ascii_uppercase()) != Some(current_char.to_ascii_uppercase()) {
            compressed.push(current_char);
            last_char = Some(current_char);
        }
    }
    compressed
}

// --- Simple Hash Function ---
#[inline]
fn simple_mix_hash(mut n: u64) -> u64 {
    const K: u64 = 0x517cc1b727220a95;
    n = n.wrapping_mul(K);
    n ^= n >> 32;
    n
}

// --- K-mer Subsampling Check ---
#[inline]
fn should_process_kmer(hash: u64, trailing_zeros: u32) -> bool {
    if trailing_zeros == 0 { true } else {
        let mixed_hash = simple_mix_hash(hash);
        let mask = (1u64 << trailing_zeros) - 1;
        (mixed_hash & mask) == 0
    }
}

// --- 2-bit K-mer Counting Filter Implementation ---
#[derive(Clone, Debug)]
struct KmerCounterFilter {
    counters: Vec<u64>, // Each u64 holds 32 2-bit counters
    num_counters: u64,  // Total number of logical counters
}

impl KmerCounterFilter {
    /// Creates a new KmerCounterFilter.
    pub fn new(num_counters: u64) -> Self {
        if num_counters == 0 {
            panic!("Number of counters must be greater than 0");
        }
        let num_u64 = (num_counters + 31) / 32;
        let counters = vec![0u64; num_u64 as usize];
        let actual_num_counters = num_u64 * 32;

        KmerCounterFilter {
            counters,
            num_counters: actual_num_counters,
        }
    }

    /// Calculates the index and bit offset for a counter.
    #[inline]
    fn get_counter_pos(&self, hash: u64) -> (usize, u32) {
        let counter_index = hash % self.num_counters;
        let vec_index = (counter_index / 32) as usize;
        let bit_offset = (counter_index % 32) * 2;
        (vec_index, bit_offset as u32)
    }

    /// Increments the counter for a given hash value (saturating at 3).
    #[inline]
    pub fn increment(&mut self, hash: u64) {
        let (vec_index, bit_offset) = self.get_counter_pos(hash);
        if vec_index < self.counters.len() {
            let current_val = (self.counters[vec_index] >> bit_offset) & 0b11;
            if current_val < 3 {
                let new_val = current_val + 1;
                self.counters[vec_index] &= !(0b11u64 << bit_offset);
                self.counters[vec_index] |= new_val << bit_offset;
            }
        } else {
            // This condition should ideally not be hit if num_counters is calculated correctly
            // eprintln!("Warning: Counter index out of bounds during increment: vec_idx={}, len={}", vec_index, self.counters.len());
        }
    }

    /// Checks if the counter for a hash value is "solid" (value is 3).
    #[inline]
    pub fn check_solid(&self, hash: u64) -> bool {
        let (vec_index, bit_offset) = self.get_counter_pos(hash);
        if vec_index < self.counters.len() {
            let current_val = (self.counters[vec_index] >> bit_offset) & 0b11;
            current_val == 3
        } else {
             // This condition should ideally not be hit
            // eprintln!("Warning: Counter index out of bounds during check_solid: vec_idx={}, len={}", vec_index, self.counters.len());
            false
        }
    }
}

/// 64-bit multipliers for rolling hash (different for first- and second-level).
const BASES: [u64; 9] = [0x9e3779b97f4a7c15, 0xc2b2ae3d27d4eb4f, 0x165667b19e3779f9, 0x27d4eb2f165667c5, 0xa0761d6478bd642f, 0xe7037ed1a0b428db, 0xbf58476d1ce4e5b9, 0x94d049bb133111eb, 0x2545f4914f6cdd1d];
const PBSTYLE: &str = "{prefix} [{bar:40.cyan/blue}] {pos}/{len} {elapsed} ETA: {eta}";
// --- Constants ---
const NUM_COUNTER_SHARDS: usize = 1024; // Number of shards for the counter filter
const READ_CHUNK_SIZE: usize = 1000;


// --- Counting Filter Pass Function ---

/// Performs the first pass to populate the KmerCounterFilter using homopolymer-compressed and subsampled k-mers.
/// Uses sharded filters to reduce lock contention. Processes reads in parallel chunks.
/// Returns the *sharded* counter filter.
fn run_counting_filter_pass(
    input_filename: &str,
    k: usize,
    args: &Args,
    pool: &ThreadPool,
) -> std::io::Result<Arc<Vec<Arc<Mutex<KmerCounterFilter>>>>> {

    // --- Sharded Counter Filter Initialization ---
    let num_counters_total = args.filter_counters;
    let trailing_zeros = args.trailing_zeros;

    if num_counters_total == 0 {
         return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "Filter counter number must be greater than 0"));
    }
    let counters_per_shard = (num_counters_total + NUM_COUNTER_SHARDS as u64 - 1) / NUM_COUNTER_SHARDS as u64;

    // Create counter shards vector
    let counter_shards_vec: Vec<Arc<Mutex<KmerCounterFilter>>> = (0..NUM_COUNTER_SHARDS)
        .map(|_i| { Arc::new(Mutex::new(KmerCounterFilter::new(counters_per_shard))) })
        .collect();
    let counter_shards = Arc::new(counter_shards_vec);

    let input_reader = open_input(input_filename);
    let mut buf_reader = BufReader::new(input_reader);
    let peek = buf_reader.fill_buf()?;
    let is_fastq = !peek.is_empty() && peek[0] == b'@';
    drop(buf_reader); // Close the peeking reader

    // --- Process Reads in Chunks ---
    let mut record_chunk = Vec::with_capacity(READ_CHUNK_SIZE);
    let input_reader_main = open_input(input_filename);
    let buf_reader_main = BufReader::new(input_reader_main);

    let rc_sensitivity = args.rc_sensitivity; // Local copy for closure

    if is_fastq {
        let mut reader = fastq::Reader::new(buf_reader_main); // Use mut reader
        loop {
            record_chunk.clear();
            for _ in 0..READ_CHUNK_SIZE {
                 let mut record = fastq::Record::new();
                 match reader.read(&mut record) {
                     Ok(_) => {
                         if record.is_empty() { continue; }
                         let id = record.id().to_owned();
                         let seq = String::from_utf8_lossy(record.seq()).into_owned();
                         let qual = String::from_utf8_lossy(record.qual()).into_owned();
                         record_chunk.push(Record::Fastq { id, seq, qual });
                     }
                     Err(e) => { eprintln!("Error reading FASTQ record: {}", e); break; } // Break on error too
                 }
            }

            if record_chunk.is_empty() { break; }

            pool.install(|| {
                record_chunk.par_iter().for_each(|record| {
                    process_record_for_counting(record, k, rc_sensitivity, trailing_zeros, &counter_shards);
                });
            });
         }

    } else { // FASTA
        let mut reader = fasta::Reader::new(buf_reader_main); // Use mut reader
         loop {
             record_chunk.clear();
             for _ in 0..READ_CHUNK_SIZE {
                 let mut record = fasta::Record::new();
                 match reader.read(&mut record) {
                     Ok(_) => {
                         if record.is_empty() { continue; }
                         let id = record.id().to_owned();
                         let seq = String::from_utf8_lossy(record.seq()).into_owned();
                         record_chunk.push(Record::Fasta { id, seq });
                     }
                     Err(e) => { eprintln!("Error reading FASTA record: {}", e); break; }
                 }
             }

            if record_chunk.is_empty() { break; }

            pool.install(|| {
                record_chunk.par_iter().for_each(|record| {
                    process_record_for_counting(record, k, rc_sensitivity, trailing_zeros, &counter_shards);
                });
            });
        }
    }

    Ok(counter_shards)
}

/// Helper function to process a single record for KmerCounterFilter population.
fn process_record_for_counting(
    record: &Record,
    k: usize,
    rc_sensitivity: bool,
    trailing_zeros: u32,
    counter_shards: &Arc<Vec<Arc<Mutex<KmerCounterFilter>>>>,
) {
    let seq = match record {
        Record::Fasta { ref seq, .. } => seq,
        Record::Fastq { ref seq, .. } => seq,
    };
    let compressed_seq = homopolymer_compress(seq);
    let seq_bytes = compressed_seq.as_bytes();
    if seq_bytes.len() < k { return; }

    let kmer_iterator: Box<dyn Iterator<Item = u64> + Send> = if rc_sensitivity {
        if let Ok(iter) = NtHashIterator::new(seq_bytes, k) { Box::new(iter) }
        else { return; }
    } else {
        if let Ok(iter) = NtHashForwardIterator::new(seq_bytes, k) { Box::new(iter) }
        else { return; }
    };

    for hash_val in kmer_iterator {
        if !should_process_kmer(hash_val, trailing_zeros) { continue; }

        let shard_idx = (hash_val as usize) % NUM_COUNTER_SHARDS;
        if shard_idx >= NUM_COUNTER_SHARDS { continue; }
        if let Ok(mut counter_guard) = counter_shards[shard_idx].try_lock() {
            counter_guard.increment(hash_val);
        }
    }
}

/// Computes the partition fingerprint using homopolymer-compressed and subsampled k-mers,
/// checking the appropriate shard of the sharded KmerCounterFilter if provided.
fn compute_partition(
    seq: &str, // Original sequence
    p: usize,
    k: usize,
    base: u64,
    counter_shards_opt: Option<&Arc<Vec<Arc<Mutex<KmerCounterFilter>>>>>,
    trailing_zeros: u32,
) -> u64 {
    let compressed_seq = homopolymer_compress(seq);
    let seq_bytes = compressed_seq.as_bytes();
    if seq_bytes.len() < k { return 0; }

    let index_bits = if p > 0 { (p as f64).log2().ceil() as usize } else { 0 };
    let index_bits = index_bits.min(64);
    let buckets_count = (1usize << index_bits).max(1);
    let mut buckets: Vec<Option<u64>> = vec![None; buckets_count];
    let value_bits = 64usize.saturating_sub(index_bits);
    let value_mask = if value_bits >= 64 { u64::MAX } else { (1u64 << value_bits) - 1 };

    if let Ok(iter) = NtHashForwardIterator::new(seq_bytes, k) {
        for kmer_hash in iter {
            if !should_process_kmer(kmer_hash, trailing_zeros) { continue; }
            let mut process_kmer = true;
            if let Some(counter_shards) = counter_shards_opt {
                let shard_idx = (kmer_hash as usize) % NUM_COUNTER_SHARDS;
                if shard_idx < counter_shards.len() {
                    if let Ok(counter_guard) = counter_shards[shard_idx].try_lock() {
                        process_kmer = counter_guard.check_solid(kmer_hash);
                    } else { process_kmer = false; }
                } else { process_kmer = false; }
            }

            if process_kmer {
                let hash = kmer_hash.wrapping_mul(base);
                let idx = if buckets_count > 1 { (hash >> (64 - index_bits)) as usize } else { 0 };
                let value = hash & value_mask;
                if idx < buckets_count {
                    buckets[idx] = Some(match buckets[idx] {
                        Some(current) if value >= current => current,
                        _ => value,
                    });
                }
            }
        }
    }

    let mut fingerprint = 0u64;
    for i in 0..p {
        let bucket_idx = i % buckets_count;
        let bit = buckets[bucket_idx].map_or(0, |v| (v >> (i / buckets_count)) & 1);
        fingerprint = (fingerprint << 1) | bit;
    }
    fingerprint
}

/// Computes the canonical partition fingerprint using homopolymer-compressed and subsampled k-mers,
/// checking the appropriate shard of the sharded KmerCounterFilter if provided.
fn compute_partition_canonique(
    seq: &str, // Original sequence
    p: usize,
    k: usize,
    base: u64,
    counter_shards_opt: Option<&Arc<Vec<Arc<Mutex<KmerCounterFilter>>>>>,
    trailing_zeros: u32,
) -> u64 {
    let compressed_seq = homopolymer_compress(seq);
    let seq_bytes = compressed_seq.as_bytes();
    if seq_bytes.len() < k { return 0; }

    let index_bits = if p > 0 { (p as f64).log2().ceil() as usize } else { 0 };
    let index_bits = index_bits.min(64);
    let buckets_count = (1usize << index_bits).max(1);
    let mut buckets: Vec<Option<u64>> = vec![None; buckets_count];
    let value_bits = 64usize.saturating_sub(index_bits);
    let value_mask = if value_bits >= 64 { u64::MAX } else { (1u64 << value_bits) - 1 };

    if let Ok(iter) = NtHashIterator::new(seq_bytes, k) { // Use canonical iterator
        for kmer_hash in iter {
            if !should_process_kmer(kmer_hash, trailing_zeros) { continue; }
            let mut process_kmer = true;
            if let Some(counter_shards) = counter_shards_opt {
                let shard_idx = (kmer_hash as usize) % NUM_COUNTER_SHARDS;
                 if shard_idx < counter_shards.len() {
                     if let Ok(counter_guard) = counter_shards[shard_idx].try_lock() {
                        process_kmer = counter_guard.check_solid(kmer_hash);
                    } else { process_kmer = false; }
                 } else { process_kmer = false; }
            }

            if process_kmer {
                let hash = kmer_hash.wrapping_mul(base);
                let idx = if buckets_count > 1 { (hash >> (64 - index_bits)) as usize } else { 0 };
                let value = hash & value_mask;
                 if idx < buckets_count {
                    buckets[idx] = Some(match buckets[idx] {
                        Some(current) if value >= current => current,
                        _ => value,
                    });
                }
            }
        }
    }

    let mut fingerprint = 0u64;
    for i in 0..p {
        let bucket_idx = i % buckets_count;
        let bit = buckets[bucket_idx].map_or(0, |v| (v >> (i / buckets_count)) & 1);
        fingerprint = (fingerprint << 1) | bit;
    }
    fingerprint
}

/// Generate a Gray code sequence for 'bits' bits.
fn generate_gray_code_order(bits: usize) -> Vec<usize> {
    let n = 1usize << bits;
    let mut order = Vec::with_capacity(n);
    for i in 0..n {
        order.push(i ^ (i >> 1));
    }
    order
}

fn reverse_complement(dna: &str) -> String {
    dna.chars()
        .rev()
        .map(|n| match n {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            'N' | 'n' => 'N',
            _ => n,
        })
        .collect()
}


// Helper function for RC update score calculation using *compressed* and *subsampled* sequences
fn calculate_scores_compressed(
    fwd_seq: &str, // Original forward sequence
    rev_seq: &str, // Original reverse sequence
    k: usize,
    kmer_counts: &AHashMap<u64, u64>, // Read-only access needed
    trailing_zeros: u32,
) -> (u64, u64) {
    let mut fwd_score = 0u64;
    let mut rev_score = 0u64;

    let compressed_fwd = homopolymer_compress(fwd_seq);
    let compressed_rev = homopolymer_compress(rev_seq);

    if compressed_fwd.len() >= k {
        if let Ok(fwd_iter) = NtHashForwardIterator::new(compressed_fwd.as_bytes(), k) {
            for fwd_kmer in fwd_iter {
                if should_process_kmer(fwd_kmer, trailing_zeros) {
                    fwd_score = fwd_score.saturating_add(*kmer_counts.get(&fwd_kmer).unwrap_or(&0));
                }
            }
        }
    }
    if compressed_rev.len() >= k {
         if let Ok(rev_iter) = NtHashForwardIterator::new(compressed_rev.as_bytes(), k) {
            for rev_kmer in rev_iter {
                 if should_process_kmer(rev_kmer, trailing_zeros) {
                     rev_score = rev_score.saturating_add(*kmer_counts.get(&rev_kmer).unwrap_or(&0));
                 }
            }
        }
    }
    (fwd_score, rev_score)
}


/// Processes a single partition file for RC orientation update.
/// Reads records, scores based on compressed/subsampled k-mers,
/// flips records if necessary, and writes back to a temporary file before replacing the original.
/// Assumes input partition file is UNCOMPRESSED. Writes temporary output UNCOMPRESSED.
fn process_partition_for_rc(
    partition_filename: &str,
    args: &Args,
    k: usize,
    pool: &ThreadPool,
) -> std::io::Result<()> {
    let rc_loops = args.rc_compression_loop;
    if rc_loops <= 0 { return Ok(()) };

    // --- Open UNCOMPRESSED input ---
    let input_file = match File::open(partition_filename) {
        Ok(f) => f,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(()),
        Err(e) => return Err(e),
    };
    let input_reader : Box<dyn Read> = Box::new(input_file);
    let buf_reader = BufReader::new(input_reader);

    let trailing_zeros = args.trailing_zeros;

    // Determine file type
    let mut peekable_reader = buf_reader.lines();
    let first_line = peekable_reader.next().transpose()?;
    let is_fastq = first_line.as_deref().map_or(false, |line| line.starts_with('@'));
    drop(peekable_reader);

    // --- Re-open UNCOMPRESSED for actual reading ---
    let input_file_main = File::open(partition_filename)?;
    let input_reader_main : Box<dyn Read> = Box::new(input_file_main);
    let buf_reader_main = BufReader::new(input_reader_main);

    let mut records: Vec<Record> = Vec::new();
    // Read all records from the partition file
    if is_fastq {
        let reader = fastq::Reader::new(buf_reader_main);
        for result in reader.records() {
            match result {
                Ok(record) => {
                    let id = record.id().to_owned();
                    let seq = String::from_utf8_lossy(record.seq()).into_owned();
                    let qual = String::from_utf8_lossy(record.qual()).into_owned();
                    records.push(Record::Fastq { id, seq, qual });
                }
                Err(e) => return Err(io::Error::new(io::ErrorKind::Other, e)),
            }
        }
    } else {
        let reader = fasta::Reader::new(buf_reader_main);
        for result in reader.records() {
             match result {
                 Ok(record) => {
                     let id = record.id().to_owned();
                     let seq = String::from_utf8_lossy(record.seq()).into_owned();
                     records.push(Record::Fasta { id, seq });
                 }
                 Err(e) => return Err(io::Error::new(io::ErrorKind::Other, e)),
             }
        }
    }

    if records.is_empty() {
        return Ok(());
    }

    let kmer_counts_arc = Arc::new(Mutex::new(AHashMap::<u64, u64>::new()));

    // Initial population (parallel)
    pool.install(|| {
        records.par_iter().for_each(|record| {
            let seq = match record { Record::Fastq { seq, .. } | Record::Fasta { seq, .. } => seq };
            let compressed_seq = homopolymer_compress(seq);
            let bytes = compressed_seq.as_bytes();
            if bytes.len() >= k {
                if let Ok(iter) = NtHashForwardIterator::new(bytes, k) {
                    let mut counts_guard = kmer_counts_arc.lock().unwrap();
                    for kmer_hash in iter {
                        if should_process_kmer(kmer_hash, trailing_zeros) {
                            *counts_guard.entry(kmer_hash).or_insert(0) += 1;
                        }
                    }
                }
            }
        });
    });


    for _loop_num in 0..rc_loops {
        // Score calculation (parallel)
        let flip_decisions: Vec<(usize, bool, String, Option<String>)> = pool.install(|| {
            records.par_iter().enumerate().map(|(idx, record)| {
                let kmer_counts_guard = kmer_counts_arc.lock().unwrap();
                let (should_flip, rev_seq_opt, rev_qual_opt) = match record {
                    Record::Fastq { seq, qual, .. } => {
                        let rev_seq = reverse_complement(seq);
                        let (fwd_score, rev_score) = calculate_scores_compressed(seq, &rev_seq, k, &kmer_counts_guard, trailing_zeros);
                        if rev_score > fwd_score { (true, Some(rev_seq), Some(qual.chars().rev().collect())) } else { (false, None, None) }
                    }
                    Record::Fasta { seq, .. } => {
                        let rev_seq = reverse_complement(seq);
                        let (fwd_score, rev_score) = calculate_scores_compressed(seq, &rev_seq, k, &kmer_counts_guard, trailing_zeros);
                        if rev_score > fwd_score { (true, Some(rev_seq), None) } else { (false, None, None) }
                    }
                };
                 drop(kmer_counts_guard);
                (idx, should_flip, rev_seq_opt.unwrap_or_default(), rev_qual_opt)
            }).filter(|(_, should_flip, _, _)| *should_flip).collect()
        });

        let flipped_count = flip_decisions.len();
        if flipped_count == 0 { break; }

        // Update dictionary and records (sequential)
        let mut kmer_counts_guard = kmer_counts_arc.lock().unwrap();
        for (idx, _should_flip, rev_seq, rev_qual_opt) in flip_decisions {
             match &mut records[idx] {
                 Record::Fastq { seq, qual, .. } => {
                     let old_compressed = homopolymer_compress(seq);
                     if old_compressed.len() >= k {
                         if let Ok(iter) = NtHashForwardIterator::new(old_compressed.as_bytes(), k) {
                             for kmer_hash in iter {
                                 if should_process_kmer(kmer_hash, trailing_zeros) {
                                    if let Some(count) = kmer_counts_guard.get_mut(&kmer_hash) {
                                        *count = count.saturating_sub(1);
                                    }
                                 }
                             }
                         }
                     }
                     *seq = rev_seq;
                     *qual = rev_qual_opt.expect("Missing rev_qual for FASTQ");
                     let new_compressed = homopolymer_compress(seq);
                      if new_compressed.len() >= k {
                         if let Ok(iter) = NtHashForwardIterator::new(new_compressed.as_bytes(), k) {
                             for kmer_hash in iter {
                                 if should_process_kmer(kmer_hash, trailing_zeros) {
                                     *kmer_counts_guard.entry(kmer_hash).or_insert(0) += 1;
                                 }
                             }
                         }
                     }
                 }
                 Record::Fasta { seq, .. } => {
                      let old_compressed = homopolymer_compress(seq);
                      if old_compressed.len() >= k {
                         if let Ok(iter) = NtHashForwardIterator::new(old_compressed.as_bytes(), k) {
                             for kmer_hash in iter {
                                 if should_process_kmer(kmer_hash, trailing_zeros) {
                                    if let Some(count) = kmer_counts_guard.get_mut(&kmer_hash) {
                                        *count = count.saturating_sub(1);
                                    }
                                 }
                             }
                         }
                     }
                      *seq = rev_seq;
                      let new_compressed = homopolymer_compress(seq);
                      if new_compressed.len() >= k {
                         if let Ok(iter) = NtHashForwardIterator::new(new_compressed.as_bytes(), k) {
                             for kmer_hash in iter {
                                 if should_process_kmer(kmer_hash, trailing_zeros) {
                                     *kmer_counts_guard.entry(kmer_hash).or_insert(0) += 1;
                                 }
                             }
                         }
                     }
                 }
             }
        }
        drop(kmer_counts_guard);
    }

    // --- Write updated records to an UNCOMPRESSED temporary file ---
    let temp_partition_filename = partition_filename.to_owned() + ".rc_temp";
    { // Scope for writer
        let temp_file = File::create(&temp_partition_filename)?;
        let mut writer = BufWriter::new(temp_file);

        if is_fastq {
             let mut fastq_writer = fastq::Writer::new(&mut writer);
             for record in &records {
                 if let Record::Fastq { id, seq, qual } = record {
                     let fastq_record = fastq::Record::with_attrs(id, None, seq.as_bytes(), qual.as_bytes());
                     fastq_writer.write_record(&fastq_record)?;
                 }
             }
             fastq_writer.flush()?;
        } else {
            let mut fasta_writer = fasta::Writer::new(&mut writer);
             for record in &records {
                 if let Record::Fasta { id, seq } = record {
                     let fasta_record = fasta::Record::with_attrs(id, None, seq.as_bytes());
                    fasta_writer.write_record(&fasta_record)?;
                 }
             }
            fasta_writer.flush()?;
        }
        writer.flush()?; // Flush BufWriter
    } // Writer scope ends, file is closed

    // Replace original partition file
    fs::rename(&temp_partition_filename, partition_filename)?;

    Ok(())
}



/// Creates UNCOMPRESSED bucket files. Fingerprint is calculated from compressed & subsampled sequence,
/// but the *original* sequence is written to the partition file.
fn create_bucket_files(
    filename_input: &str,
    filename_comp_base: &str, // Base name for partition files
    args: &Args,
    p: usize,
    k: usize,
    base: u64,
    counter_shards_opt: Option<Arc<Vec<Arc<Mutex<KmerCounterFilter>>>>>,
) -> std::io::Result<()> {

    let input = open_input(filename_input);
    let mut buf_reader = BufReader::new(input);
    let peek = buf_reader.fill_buf().map_err(|e| io::Error::new(io::ErrorKind::Other, format!("Error peeking input {}: {}", filename_input, e)))?;
    let is_fastq = !peek.is_empty() && peek[0] == b'@';
    drop(buf_reader); // Close peeking reader

    let rc_sensitivity = args.rc_sensitivity;
    let trailing_zeros = args.trailing_zeros;
    let size_array = 1usize << p;

    if let Some(parent_dir) = Path::new(filename_comp_base).parent() { fs::create_dir_all(parent_dir)?; }

    let mut writers: Vec<Arc<Mutex<BufWriter<File>>>> = Vec::with_capacity(size_array);
    for fp in 0..size_array {
        let filename_partition = get_filename_partition(filename_comp_base, fp, p);
        let file = File::create(&filename_partition).unwrap_or_else(|e| panic!("Cannot create partition file {}: {}", filename_partition, e));
        writers.push(Arc::new(Mutex::new(BufWriter::new(file))));
    }
    let files = Arc::new(writers);

    let num_workers = if args.thread > 0 { std::cmp::min(args.thread, num_cpus::get()) } else { num_cpus::get() };
    let pool = ThreadPoolBuilder::new().num_threads(num_workers).build().unwrap();

    // --- Process Reads in Chunks for Partitioning ---
    let mut record_chunk = Vec::with_capacity(READ_CHUNK_SIZE);
    let mut _total_reads_partitioned: u64 = 0;

    let input_reader_main = open_input(filename_input);
    let buf_reader_main = BufReader::new(input_reader_main);

    let counter_shards_ref = counter_shards_opt.clone();

    if is_fastq {
        let mut reader = fastq::Reader::new(buf_reader_main);
        loop {
            record_chunk.clear();
            for _ in 0..READ_CHUNK_SIZE {
                let mut record = fastq::Record::new();
                match reader.read(&mut record) {
                    Ok(_) => {
                        if record.is_empty() { continue; }
                        let id = record.id().to_owned();
                        let seq = String::from_utf8_lossy(record.seq()).into_owned();
                        let qual = String::from_utf8_lossy(record.qual()).into_owned();
                        record_chunk.push(Record::Fastq { id, seq, qual });
                    },
                    Err(e) => { eprintln!("Error reading FASTQ record in {}: {}", filename_input, e); break; }
                }
            }
            if record_chunk.is_empty() { break; }
            _total_reads_partitioned += record_chunk.len() as u64;

            let counter_shards_iter = counter_shards_ref.clone();
            pool.install(|| {
                record_chunk.par_iter().for_each(|record| {
                    let (id, original_seq, qual_opt) = match record {
                        Record::Fastq { id, seq, qual } => (id.clone(), seq.clone(), Some(qual.clone())),
                         _ => unreachable!(),
                    };
                    let counter_ref_inner = counter_shards_iter.as_ref();
                    let fp = if rc_sensitivity {
                        compute_partition_canonique(&original_seq, p, k, base, counter_ref_inner, trailing_zeros)
                    } else {
                        compute_partition(&original_seq, p, k, base, counter_ref_inner, trailing_zeros)
                    };
                    let fp_usize = fp as usize;
                    if fp_usize >= files.len() { eprintln!("Error: Fingerprint {} out of bounds ({}) for record {} in {}", fp, files.len(), id, filename_input); return; }

                     match files[fp_usize].lock() {
                         Ok(mut writer) => {
                             if let Some(qual) = qual_opt { // FASTQ
                                 if writeln!(writer, "@{}", id).is_err() { return; }
                                 if writeln!(writer, "{}", original_seq).is_err() { return; }
                                 if writeln!(writer, "+").is_err() { return; }
                                 if writeln!(writer, "{}", qual).is_err() { return; }
                             } else { unreachable!() }
                         },
                         Err(_) => { eprintln!("Partition writer mutex poisoned for fp={}", fp_usize); }
                     }
                });
            });
        }
    } else { // FASTA
        let mut reader = fasta::Reader::new(buf_reader_main);
        loop {
            record_chunk.clear();
             for _ in 0..READ_CHUNK_SIZE {
                 let mut record = fasta::Record::new();
                 match reader.read(&mut record) {
                    Ok(_) => {
                        if record.is_empty() { continue; }
                         let id = record.id().to_owned();
                         let seq = String::from_utf8_lossy(record.seq()).into_owned();
                         record_chunk.push(Record::Fasta { id, seq });
                     },
                    Err(e) => { eprintln!("Error reading FASTA record in {}: {}", filename_input, e); break; }
                 }
             }
             if record_chunk.is_empty() { break; }
            _total_reads_partitioned += record_chunk.len() as u64;

            let counter_shards_iter = counter_shards_ref.clone();
            pool.install(|| {
                 record_chunk.par_iter().for_each(|record| {
                    let (id, original_seq) = match record {
                        Record::Fasta { id, seq } => (id.clone(), seq.clone()),
                         _ => unreachable!(),
                    };

                    let counter_ref_inner = counter_shards_iter.as_ref();
                    let fp = if rc_sensitivity {
                        compute_partition_canonique(&original_seq, p, k, base, counter_ref_inner, trailing_zeros)
                    } else {
                        compute_partition(&original_seq, p, k, base, counter_ref_inner, trailing_zeros)
                    };
                    let fp_usize = fp as usize;
                     if fp_usize >= files.len() { eprintln!("Error: Fingerprint {} out of bounds ({}) for record {} in {}", fp, files.len(), id, filename_input); return; }

                     match files[fp_usize].lock() {
                        Ok(mut writer) => {
                             if writeln!(writer, ">{}", id).is_err() { return; }
                             if writeln!(writer, "{}", original_seq).is_err() { return; }
                         },
                         Err(_) => { eprintln!("Partition writer mutex poisoned for fp={}", fp_usize); }
                     }
                 });
            });
        }
    }

    // --- Flush and close UNCOMPRESSED partition writers ---
    for (idx, writer_mutex) in files.iter().enumerate() {
        match writer_mutex.lock() {
            Ok(mut writer_guard) => {
                if let Err(e) = writer_guard.flush() {
                    eprintln!("Flush error for partition writer {}: {}", idx, e);
                }
            },
            Err(_) => {
                 eprintln!("Partition writer mutex poisoned during close for index {}. Data might be lost.", idx);
            }
        }
    }

    Ok(())
}


/// Concatenates UNCOMPRESSED bucket files using Gray code order and cleans up partitions.
/// Handles final compression type (zstd/gzip).
fn concat_bucket_files(
    partition_base_filename: &str, // Base name used for creating partitions
    final_output_filename: &str,    // The actual final output file path
    args: &Args,
    p: usize,
    pool: &ThreadPool, // For parallel cleanup
) -> std::io::Result<()> {
    let gray_order = generate_gray_code_order(p);

    if let Some(parent_dir) = Path::new(final_output_filename).parent() { fs::create_dir_all(parent_dir)?; }

    let output_compression_level = args.final_compression_level;
    let use_gzip_final = args.final_compression.to_lowercase() == "gzip";

    let final_file = File::create(final_output_filename)?;

    if use_gzip_final {
        let gzip_level = GzipCompression::new(output_compression_level.clamp(0, 9) as u32);
        println!("Using Gzip level {} for final output: {}", output_compression_level.clamp(0, 9), final_output_filename);
        let encoder = GzEncoder::new(final_file, gzip_level);
        let mut buf_writer = BufWriter::new(encoder);

        for partition_idx in &gray_order {
            let part_filename = get_filename_partition(partition_base_filename, *partition_idx, p);
            match File::open(&part_filename) {
                Ok(file) => {
                    let mut reader = BufReader::new(file);
                    match std::io::copy(&mut reader, &mut buf_writer) {
                         Ok(_) => {  }
                         Err(_e) => { eprintln!("Error concatenating partition {} into {}", part_filename, final_output_filename); }
                     }
                }
                Err(ref e) if e.kind() == std::io::ErrorKind::NotFound => { /* Okay */ }
                Err(_e) => { eprintln!("Warning: Could not open partition file {} for concatenation", part_filename); }
            }
        }
        buf_writer.flush()?;
    } else {
        // Use zstd for final output
        let zstd_level = output_compression_level;
        println!("Using Zstd level {} for final output: {}", zstd_level, final_output_filename);
        let encoder_result = Encoder::new(final_file, zstd_level);
        let mut encoder = encoder_result.map_err(|e| io::Error::new(io::ErrorKind::Other, format!("Cannot create final zstd encoder for {}: {}", final_output_filename, e)))?;

        let num_threads = if args.thread > 0 { args.thread } else { num_cpus::get() };
        if num_threads > 0 {
            if encoder.set_parameter(CParameter::NbWorkers(num_threads as u32)).is_err() {
                 eprintln!("Warning: Failed to set {} threads for ZSTD compression.", num_threads);
            }
        }

        // Scope buf_writer so it's dropped before encoder.finish()
        {
            let mut buf_writer = BufWriter::new(&mut encoder);
            for partition_idx in &gray_order {
                let part_filename = get_filename_partition(partition_base_filename, *partition_idx, p);
                match File::open(&part_filename) {
                    Ok(file) => {
                        let mut reader = BufReader::new(file);
                        match std::io::copy(&mut reader, &mut buf_writer) {
                             Ok(_) => { }
                             Err(_e) => { eprintln!("Error concatenating partition {} into {}", part_filename, final_output_filename); }
                         }
                    }
                    Err(ref e) if e.kind() == std::io::ErrorKind::NotFound => { /* Okay */ }
                    Err(_e) => { eprintln!("Warning: Could not open partition file {} for concatenation", part_filename); }
                }
            }
            buf_writer.flush()?;
        } // buf_writer is dropped here, releasing borrow on encoder
        encoder.finish()?; // Finish ZstdEncoder
    }

    // --- Clean up partition files ---
    let cleanup_errors = Arc::new(Mutex::new(Vec::new()));
    pool.install(|| {
         (0..(1 << p)).into_par_iter().for_each(|i| {
            let filename_partition = get_filename_partition(partition_base_filename, i, p);
            if let Err(e) = fs::remove_file(&filename_partition) {
                if e.kind() != std::io::ErrorKind::NotFound {
                    if Path::new(&filename_partition).exists() {
                        eprintln!("Warning: could not remove partition {}: {}", filename_partition, e);
                        cleanup_errors.lock().expect("Cleanup mutex poisoned").push(filename_partition);
                    }
                }
            }
        });
    });
     let errors = cleanup_errors.lock().expect("Cleanup mutex poisoned");
     if !errors.is_empty() { eprintln!("Failed to remove {} partition files.", errors.len()); }

    Ok(())
}


/// Generates the filename for a specific partition based on a base name.
fn get_filename_partition(base_filename: &str, partition: usize, p: usize) -> String {
    let (stem, extension) = match base_filename.rfind('.') {
        Some(idx) if !base_filename[idx..].contains('/') && !base_filename[idx..].contains('\\') => {
             (&base_filename[..idx], &base_filename[idx..])
        },
        _ => (base_filename, ""),
    };
    let partition_id = format!("{:0width$b}", partition, width = p);
    format!("{}_{}{}", stem, partition_id, extension)
}

fn increase_file_limit(requested_limit: u64) -> io::Result<()> {
    let (soft_limit, hard_limit) = rlimit::getrlimit(rlimit::Resource::NOFILE)?;
    println!("Current file limits: soft={}, hard={}", soft_limit, hard_limit);

    let target_limit = std::cmp::min(requested_limit, hard_limit); // Don't exceed hard limit

    if target_limit > soft_limit {
        match rlimit::setrlimit(rlimit::Resource::NOFILE, target_limit, hard_limit) {
            Ok(_) => {
                let (new_soft, new_hard) = rlimit::getrlimit(rlimit::Resource::NOFILE)?;
                println!("Successfully set file limits: soft={}, hard={}", new_soft, new_hard);
            }
            Err(e) => {
                eprintln!("Warning: Failed to increase file limit to {}: {}", target_limit, e);
                eprintln!("Proceeding with current limit: {}", soft_limit);
                // Don't return error, just warn and continue
            }
        }
    } else {
        println!("Current soft limit ({}) is already sufficient.", soft_limit);
    }
    Ok(())
}


fn main() -> std::io::Result<()> {
    let args = Args::parse();
    println!("Running with arguments: {:?}", args);
    let start_total = Instant::now();
    fs::create_dir_all(&args.output).expect("Cannot create output directory");

    let k = args.k;
    let p = args.p;

    // --- Increase file limit ---
    let num_partitions = 1usize << p;
    // Request slightly more than needed for safety margin (e.g., + stdin/out/err + some buffer)
    let requested_limit = (num_partitions + 100) as u64;
    increase_file_limit(requested_limit)?;
    // --- End increase file limit ---


    let num_threads = if args.thread > 0 { args.thread } else { num_cpus::get() };
    println!("Using {} threads for processing.", num_threads);

    if p == 0 {
        eprintln!("Error: p value cannot be 0.");
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "p value cannot be 0"));
    }
    if p > 20 {
        eprintln!("Warning: p value {} is very large (>20), will create >1M partitions. This may exceed filesystem limits or performance.", p);
    } else if p > 16 {
        eprintln!("Warning: p value {} is large (>16), may lead to >65k partitions.", p);
    }


    let pool = ThreadPoolBuilder::new().num_threads(num_threads).build().unwrap();

    // --- Determine Final Output Filename and Intermediate Base ---
    let input_path = Path::new(&args.input);
    let input_stem = input_path.file_stem().unwrap_or_default().to_string_lossy();
    let p_string = p.to_string();
    let final_filename_stem = format!("{}_k{}_p{}", input_stem, k, p_string);
    let final_extension = match args.final_compression.to_lowercase().as_str() {
        "gzip" | "gz" => "gz",
        _ => "zst",
    };
    let final_output_filename = PathBuf::from(&args.output)
                                .join(format!("{}.{}", final_filename_stem, final_extension))
                                .to_string_lossy()
                                .into_owned();
    let intermediate_base_filename = PathBuf::from(&args.output)
                                .join(format!("{}.tmp_oreo", final_filename_stem))
                                .to_string_lossy()
                                .into_owned();


    println!("Input file: {}", args.input);
    println!("Working directory: {}", args.output);


    // --- Pipeline Steps ---

    // 1. Optional Counting Filter Pass
    let counter_shards_opt: Option<Arc<Vec<Arc<Mutex<KmerCounterFilter>>>>> = if args.use_counting_filter {
        println!("--- Phase 1: Counting Filter Population ---");
        let start_phase = Instant::now();
        let filters = run_counting_filter_pass(&args.input, k, &args, &pool)?;
        println!("Counting Filter Population finished in {:.3} seconds.", start_phase.elapsed().as_secs_f64());
        Some(filters)
    } else {
        println!("--- Phase 1: Counting Filter Population (Skipped) ---");
        None
    };

    // 2. Partitioning Pass (Single Level)
    println!("--- Phase 2: Partitioning (p={}) ---", p);
    let start_phase = Instant::now();
    let base_for_level = BASES[0];
    create_bucket_files(
        &args.input,
        &intermediate_base_filename,
        &args, p, k, base_for_level,
        counter_shards_opt.clone(),
    )?;
    println!("Partitioning finished in {:.3} seconds.", start_phase.elapsed().as_secs_f64());

    // 3. Release Counter Filter Memory
    drop(counter_shards_opt);

    // 4. Reverse Complementation Pass
     println!("--- Phase 3: Reverse Complement Update ---");
     let start_phase = Instant::now();
     if args.rc_compression_loop > 0 {
         let rc_bar = ProgressBar::new(num_partitions as u64);
         rc_bar.set_style(ProgressStyle::default_bar().template(PBSTYLE).unwrap().progress_chars("##-"));
         rc_bar.set_prefix(format!("RC Update (p={})", p));

         pool.install(|| {
             (0..num_partitions).into_par_iter().for_each(|i| {
                 let filename_partition = get_filename_partition(&intermediate_base_filename, i, p);
                 // No need to check exists here, process_partition_for_rc handles NotFound
                 if let Err(e) = process_partition_for_rc(&filename_partition, &args, k, &pool) {
                     eprintln!("Error during RC processing for partition {}: {}", filename_partition, e);
                 }
                 rc_bar.inc(1);
             });
         });
         rc_bar.finish_with_message("RC update finished.");
     } else {
         println!("Skipping RC update pass (rc_compression_loop = {}).", args.rc_compression_loop);
     }
     println!("RC Update finished in {:.3} seconds.", start_phase.elapsed().as_secs_f64());


    // 5. Concatenation and Cleanup Pass
    println!("--- Phase 4: Concatenation and Cleanup ---");
    let start_phase = Instant::now();
     concat_bucket_files(
         &intermediate_base_filename,
         &final_output_filename,
         &args, p, &pool
     )?;
     println!("Concatenation and Cleanup finished in {:.3} seconds.", start_phase.elapsed().as_secs_f64());


    // --- Final Statistics ---
    let final_meta = match fs::metadata(&final_output_filename) {
        Ok(meta) => meta,
        Err(e) => {
            eprintln!("Error getting metadata for final file {}: {}", final_output_filename, e);
            println!("\n--- Processing Finished (with errors) ---");
            println!("Total run time: {:.3} seconds", start_total.elapsed().as_secs_f64());
            println!("Final file {} may be incomplete or missing.", final_output_filename);
            return Err(e);
        }
    };
    let final_size = final_meta.len();
    let final_size_mb = final_size as f64 / (1024.0 * 1024.0);
    let elapsed_total = start_total.elapsed().as_secs_f64();

    println!("Final archive size: {:.3} MB", final_size_mb);
    println!("Total run time: {:.3} seconds", elapsed_total);
    Ok(())
}
