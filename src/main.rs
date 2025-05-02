use clap::Parser;
use bio::io::{fasta, fastq};
use flate2::read::GzDecoder;
use zstd::stream::{Encoder, Decoder};
use zstd::stream::raw::CParameter;
use num_cpus;
use std::sync::{Arc, Mutex};
use std::fs::{File, self};
use std::io::{self, BufReader, BufWriter, BufRead, Read, Write};
use std::time::Instant;
use rayon::prelude::*;
use rayon::{ThreadPoolBuilder, ThreadPool};
use std::path::{Path, PathBuf};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use nthash::*;
use ahash::AHashMap;

/// Command-line arguments.
#[derive(Parser, Debug)]
struct Args {
    /// Input FASTA/FASTQ file (.gz or .zst compressed are supported)
    #[arg(short, long)]
    input: String,

    /// Output directory (will be created if needed)
    #[arg(short, long, default_value = "oreo-wdir")]
    output: String,

    /// Fingerprint bit-length P. The algorithm uses an array of length P and produces a P-bit fingerprint.
    /// This results in 2^P partition files. (Default: 8)
    #[arg(short, long, value_delimiter = ',', use_value_delimiter = true, default_value = "8")]
    p: Vec<usize>,

    /// k-mer length.
    #[arg(short, long, default_value = "21")]
    k: usize,

    /// zstd compression level for writing intermediate partition files.
    #[arg(long, default_value = "1")]
    compression_level: i32,

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
    #[arg(long, default_value = "10000000")] // Default: 10M counters (~2.5 MiB) per shard
    filter_counters: u64,

    // --- Subsampling Argument ---
    /// Number of trailing zeros required in k-mer hash for subsampling (0=no subsampling, 1=1/2, 2=1/4, etc.).
    #[arg(long, default_value = "1")] // Default to 1 (hash ends in 0)
    trailing_zeros: u32,
}

/// A record (FASTA or FASTQ)
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
            eprintln!("Warning: Counter index out of bounds during increment: vec_idx={}, len={}", vec_index, self.counters.len());
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
            eprintln!("Warning: Counter index out of bounds during check_solid: vec_idx={}, len={}", vec_index, self.counters.len());
            false
        }
    }
}

/// 64-bit multipliers for rolling hash (different for first- and second-level).
const BASES: [u64; 9] = [0x9e3779b97f4a7c15, 0xc2b2ae3d27d4eb4f, 0x165667b19e3779f9, 0x27d4eb2f165667c5, 0xa0761d6478bd642f, 0xe7037ed1a0b428db, 0xbf58476d1ce4e5b9, 0x94d049bb133111eb, 0x2545f4914f6cdd1d];
const PBSTYLE: [&str; 9] = ["{prefix} [{bar:40.cyan/blue}] {pos}/{len} {elapsed} ETA: {eta}", "{prefix} [{bar:40.magenta/purple}] {pos}/{len} {elapsed} ETA: {eta}", "{prefix} [{bar:40.yellow/red}] {pos}/{len} {elapsed} ETA: {eta}", "{prefix} [{bar:40.green/black}] {pos}/{len} {elapsed} ETA: {eta}", "{prefix} [{bar:40.blue/white}] {pos}/{len} {elapsed} ETA: {eta}", "{prefix} [{bar:40.red/yellow}] {pos}/{len} {elapsed} ETA: {eta}", "{prefix} [{bar:40.white/green}] {pos}/{len} {elapsed} ETA: {eta}", "{prefix} [{bar:40.bright_blue/bright_black}] {pos}/{len} {elapsed} ETA: {eta}", "{prefix} [{bar:40.bright_cyan/bright_magenta}] {pos}/{len} {elapsed} ETA: {eta}"];
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
    println!("Starting K-mer counting filter pass (sharded, custom filter, homopolymer compression, subsampling)...");
    let start_pass = Instant::now();

    // --- Sharded Counter Filter Initialization ---
    let num_counters_total = args.filter_counters;
    let trailing_zeros = args.trailing_zeros;

    if num_counters_total == 0 {
         return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "Filter counter number must be greater than 0"));
    }
    let counters_per_shard = (num_counters_total + NUM_COUNTER_SHARDS as u64 - 1) / NUM_COUNTER_SHARDS as u64;
    let num_u64_per_shard = (counters_per_shard + 31) / 32;
    let shard_mem_mib = num_u64_per_shard * 8 / 1024 / 1024;
    println!( "Counter Filter Parameters: {} Shards, Counters/Shard={} (~{} MiB/shard), TrailingZeros={}", NUM_COUNTER_SHARDS, counters_per_shard, shard_mem_mib, trailing_zeros );
    println!( "Estimated Total Memory for Counter Filter: ~{} MiB", (NUM_COUNTER_SHARDS as u64) * shard_mem_mib );

    // Create counter shards vector
    let counter_shards_vec: Vec<Arc<Mutex<KmerCounterFilter>>> = (0..NUM_COUNTER_SHARDS)
        .map(|_i| { Arc::new(Mutex::new(KmerCounterFilter::new(counters_per_shard))) })
        .collect();
    let counter_shards = Arc::new(counter_shards_vec);

    let input_reader = open_input(input_filename);
    let mut buf_reader = BufReader::new(input_reader);
    let peek = buf_reader.fill_buf()?;
    let is_fastq = !peek.is_empty() && peek[0] == b'@';

    // --- Process Reads in Chunks ---
    let mut record_chunk = Vec::with_capacity(READ_CHUNK_SIZE);
    let mut total_reads_processed: u64 = 0;
    let input_reader_main = open_input(input_filename);
    let buf_reader_main = BufReader::new(input_reader_main);

    if is_fastq {
        let reader = fastq::Reader::new(buf_reader_main);
        let results: Vec<_> = reader.records().collect();
        for result_chunk in results.chunks(READ_CHUNK_SIZE) {
             record_chunk.clear();
             for result in result_chunk {
                 match result {
                     Ok(record) => {
                         let id = record.id().to_owned();
                         let seq = String::from_utf8_lossy(record.seq()).into_owned();
                         let qual = String::from_utf8_lossy(record.qual()).into_owned();
                         record_chunk.push(Record::Fastq { id, seq, qual });
                     }
                     Err(e) => { eprintln!("Error reading FASTQ record: {}", e); continue; }
                 }
             }
            if record_chunk.is_empty() { break; }
            total_reads_processed += record_chunk.len() as u64;

            pool.install(|| {
                record_chunk.par_iter().for_each(|record| {
                    process_record_for_counting(record, k, args.rc_sensitivity, trailing_zeros, &counter_shards);
                });
            });
        }
    } else { // FASTA
        let reader = fasta::Reader::new(buf_reader_main);
        let results: Vec<_> = reader.records().collect();
        for result_chunk in results.chunks(READ_CHUNK_SIZE) {
            record_chunk.clear();
            for result in result_chunk {
                match result {
                    Ok(record) => {
                        let id = record.id().to_owned();
                        let seq = String::from_utf8_lossy(record.seq()).into_owned();
                        record_chunk.push(Record::Fasta { id, seq });
                    }
                    Err(e) => { eprintln!("Error reading FASTA record: {}", e); continue; }
                }
            }
            if record_chunk.is_empty() { break; }
            total_reads_processed += record_chunk.len() as u64;

            pool.install(|| {
                record_chunk.par_iter().for_each(|record| {
                    process_record_for_counting(record, k, args.rc_sensitivity, trailing_zeros, &counter_shards);
                });
            });
        }
    }

    println!("Counting filter pass finished processing {} reads in {:.3} seconds.", total_reads_processed, start_pass.elapsed().as_secs_f64());

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

    let kmer_iterator: Box<dyn Iterator<Item = u64>> = if rc_sensitivity {
        if let Ok(iter) = NtHashIterator::new(seq_bytes, k) { Box::new(iter) }
        else { eprintln!("Warning: NtHashIterator failed (rc=true). Skipping record."); return; }
    } else {
        if let Ok(iter) = NtHashForwardIterator::new(seq_bytes, k) { Box::new(iter) }
        else { eprintln!("Warning: NtHashForwardIterator failed (rc=false). Skipping record."); return; }
    };

    for hash_val in kmer_iterator {
        if !should_process_kmer(hash_val, trailing_zeros) { continue; }

        let shard_idx = (hash_val as usize) % NUM_COUNTER_SHARDS;
        if shard_idx >= NUM_COUNTER_SHARDS {
             eprintln!("Warning: Calculated shard index {} out of bounds ({})", shard_idx, NUM_COUNTER_SHARDS);
             continue;
        }
        let counter_shard_mutex = &counter_shards[shard_idx];
        let mut counter_guard = counter_shard_mutex.lock().expect("Mutex counter shard poisoned");
        counter_guard.increment(hash_val);
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
                    let counter_guard = counter_shards[shard_idx].lock().expect("Mutex counter shard poisoned during check");
                    process_kmer = counter_guard.check_solid(kmer_hash);
                    drop(counter_guard);
                } else {
                    eprintln!("Warning: counter shard index out of bounds: {}", shard_idx);
                    process_kmer = false;
                }
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
                } else {
                     eprintln!("Warning: Calculated bucket index {} out of bounds ({})", idx, buckets_count);
                }
            }
        }
    } else {
         eprintln!("Warning: NtHashForwardIterator failed in compute_partition. Returning 0.");
         return 0;
    }

    let mut fingerprint = 0u64;
    for i in 0..p {
        let bucket_idx = i % buckets_count;
        let bit = buckets[bucket_idx].unwrap_or(0) & 1;
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

    if let Ok(iter) = NtHashIterator::new(seq_bytes, k) {
        for kmer_hash in iter {
            if !should_process_kmer(kmer_hash, trailing_zeros) { continue; }

            let mut process_kmer = true;
            if let Some(counter_shards) = counter_shards_opt {
                let shard_idx = (kmer_hash as usize) % NUM_COUNTER_SHARDS;
                 if shard_idx < counter_shards.len() {
                    let counter_guard = counter_shards[shard_idx].lock().expect("Mutex counter shard poisoned during check");
                    process_kmer = counter_guard.check_solid(kmer_hash);
                    drop(counter_guard);
                 } else {
                    eprintln!("Warning: counter shard index out of bounds: {}", shard_idx);
                    process_kmer = false;
                 }
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
                } else {
                     eprintln!("Warning: Calculated bucket index {} out of bounds ({}) in canonique", idx, buckets_count);
                }
            }
        }
    } else {
        eprintln!("Warning: NtHashIterator failed in compute_partition_canonique. Returning 0.");
        return 0;
    }

    let mut fingerprint = 0u64;
    for i in 0..p {
        let bucket_idx = i % buckets_count;
        let bit = buckets[bucket_idx].unwrap_or(0) & 1;
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
    kmer_counts: &AHashMap<u64, u64>,
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
        } else { eprintln!("Warning: NtHashForwardIterator failed in calculate_scores_compressed (fwd)."); }
    }
    if compressed_rev.len() >= k {
         if let Ok(rev_iter) = NtHashForwardIterator::new(compressed_rev.as_bytes(), k) {
            for rev_kmer in rev_iter {
                 if should_process_kmer(rev_kmer, trailing_zeros) {
                     // *** FIX: Get count from hash map ***
                     rev_score = rev_score.saturating_add(*kmer_counts.get(&rev_kmer).unwrap_or(&0));
                 }
            }
        } else { eprintln!("Warning: NtHashForwardIterator failed in calculate_scores_compressed (rev)."); }
    }
    (fwd_score, rev_score)
}



/// Processes a single partition file for RC orientation update.
/// Reads records, scores based on compressed/subsampled k-mers,
/// flips records if necessary, and writes back to a temporary file before replacing the original.
fn process_partition_for_rc(
    partition_filename: &str,
    args: &Args,
    k: usize,
    pool: &ThreadPool,
) -> std::io::Result<()> {
    let input_reader = open_input(partition_filename);
    let rc_loops = args.rc_compression_loop;
    let buf_reader = BufReader::new(input_reader);
    let compression_level = args.compression_level;
    let trailing_zeros = args.trailing_zeros;

    // Determine file type
    let mut peekable_reader = buf_reader.lines();
    let first_line = peekable_reader.next().transpose()?;
    let is_fastq = first_line.as_deref().map_or(false, |line| line.starts_with('@'));
    drop(peekable_reader);

    // Re-open for actual reading
    let input_reader_main = open_input(partition_filename);
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
                let kmer_counts_guard = kmer_counts_arc.lock().unwrap(); // Lock for read access
                let (should_flip, rev_seq_opt, rev_qual_opt) = match record {
                    Record::Fastq { seq, qual, .. } => {
                        let rev_seq = reverse_complement(seq);
                        let (fwd_score, rev_score) = calculate_scores_compressed(seq, &rev_seq, k, &kmer_counts_guard, trailing_zeros);
                        if fwd_score < rev_score { (true, Some(rev_seq), Some(qual.chars().rev().collect())) } else { (false, None, None) }
                    }
                    Record::Fasta { seq, .. } => {
                        let rev_seq = reverse_complement(seq);
                        let (fwd_score, rev_score) = calculate_scores_compressed(seq, &rev_seq, k, &kmer_counts_guard, trailing_zeros);
                        if fwd_score < rev_score { (true, Some(rev_seq), None) } else { (false, None, None) }
                    }
                };
                (idx, should_flip, rev_seq_opt.unwrap_or_default(), rev_qual_opt)
            }).filter(|(_, should_flip, _, _)| *should_flip).collect()
        });

        let flipped_count = flip_decisions.len();
        if flipped_count == 0 { break; }

        // Update dictionary and records (sequential)
        let mut kmer_counts_guard = kmer_counts_arc.lock().unwrap(); // Lock for updates

        for (idx, _should_flip, rev_seq, rev_qual_opt) in flip_decisions {
             match &mut records[idx] {
                 Record::Fastq { seq, qual, .. } => {
                     // Decrement old counts
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
                     // Update record
                     *seq = rev_seq;
                     *qual = rev_qual_opt.expect("Missing rev_qual for FASTQ");
                     // Increment new counts
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
                      // Decrement old counts
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
                      // Update record
                      *seq = rev_seq;
                      // Increment new counts
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

    // Write updated records to a temporary file
    let temp_partition_filename = partition_filename.to_owned() + ".rc_temp";
    {
        let temp_file = File::create(&temp_partition_filename)?;
        let encoder = Encoder::new(temp_file, compression_level)?;
        let mut writer = BufWriter::new(encoder);
        if is_fastq {
            for record in &records { if let Record::Fastq { id, seq, qual } = record { writeln!(writer, "@{}", id)?; writeln!(writer, "{}", seq)?; writeln!(writer, "+")?; writeln!(writer, "{}", qual)?; } }
        } else {
            for record in &records { if let Record::Fasta { id, seq } = record { writeln!(writer, ">{}", id)?; writeln!(writer, "{}", seq)?; } }
        }
        writer.flush()?;
        let encoder = writer.into_inner().map_err(|e| io::Error::new(io::ErrorKind::Other, format!("Writer error: {}", e)))?;
        encoder.finish()?;
    }

    // Replace original partition file
    fs::rename(&temp_partition_filename, partition_filename)?;

    Ok(())
}


// --- Modified Partitioning Functions ---

/// Creates bucket files. Fingerprint is calculated from compressed & subsampled sequence,
/// but the *original* sequence is written to the partition file.
fn create_bucket_files(
    filename_input: &str,
    filename_comp: &str, // Base name for partition files
    args: &Args,
    p: usize,
    k: usize,
    base: u64,
    counter_shards_opt: Option<Arc<Vec<Arc<Mutex<KmerCounterFilter>>>>>,
) -> std::io::Result<()> {
    let input = open_input(filename_input);
    let mut buf_reader = BufReader::new(input);
    let peek = buf_reader.fill_buf().expect("Error peeking input for partitioning");
    let is_fastq = !peek.is_empty() && peek[0] == b'@';
    let compression_level = args.compression_level;
    let rc_sensitivity = args.rc_sensitivity;
    let trailing_zeros = args.trailing_zeros;
    let size_array = 1usize << p;

    if let Some(parent_dir) = Path::new(filename_comp).parent() { fs::create_dir_all(parent_dir)?; }

    let mut writers: Vec<Arc<Mutex<BufWriter<Encoder<File>>>>> = Vec::with_capacity(size_array);
    for fp in 0..size_array {
        let filename_partition = get_filename_partition(filename_comp, fp, p);
        let file = File::create(&filename_partition).unwrap_or_else(|e| panic!("Cannot create partition file {}: {}", filename_partition, e));
        let encoder = Encoder::new(file, compression_level).expect("Cannot create zstd encoder for partition");
        writers.push(Arc::new(Mutex::new(BufWriter::new(encoder))));
    }
    let files = Arc::new(writers);

    let num_workers = if args.thread > 0 { std::cmp::min(args.thread, num_cpus::get()) } else { num_cpus::get() };

    // --- Process Reads in Chunks for Partitioning ---
    let mut record_chunk = Vec::with_capacity(READ_CHUNK_SIZE);
    let mut _total_reads_partitioned: u64 = 0;

    let input_reader_main = open_input(filename_input);
    let buf_reader_main = BufReader::new(input_reader_main);

    let pool = ThreadPoolBuilder::new().num_threads(num_workers).build().unwrap();

    if is_fastq {
        let reader = fastq::Reader::new(buf_reader_main);
        let results: Vec<_> = reader.records().collect();
        for result_chunk in results.chunks(READ_CHUNK_SIZE) {
             record_chunk.clear();
             for result in result_chunk {
                 match result {
                     Ok(record) => {
                         let id = record.id().to_owned();
                         let seq = String::from_utf8_lossy(record.seq()).into_owned();
                         let qual = String::from_utf8_lossy(record.qual()).into_owned();
                         record_chunk.push(Record::Fastq { id, seq, qual });
                     }
                     Err(e) => { eprintln!("Error reading FASTQ record: {}", e); continue; }
                 }
             }
            if record_chunk.is_empty() { break; }
            _total_reads_partitioned += record_chunk.len() as u64;

            pool.install(|| {
                record_chunk.par_iter().for_each(|record| {
                    let (id, original_seq, qual_opt) = match record {
                        Record::Fasta { id, seq } => (id.clone(), seq.clone(), None),
                        Record::Fastq { id, seq, qual } => (id.clone(), seq.clone(), Some(qual.clone())),
                    };
                    let counter_ref = counter_shards_opt.as_ref();
                    let fp = if rc_sensitivity {
                        compute_partition_canonique(&original_seq, p, k, base, counter_ref, trailing_zeros)
                    } else {
                        compute_partition(&original_seq, p, k, base, counter_ref, trailing_zeros)
                    };
                    let fp_usize = fp as usize;
                    if fp_usize >= files.len() { eprintln!("Error: Fingerprint {} out of bounds ({})", fp, files.len()); return; }
                    let writer_mutex = &files[fp_usize];
                    let mut writer = writer_mutex.lock().expect("Partition writer mutex poisoned");
                    if let Some(qual) = qual_opt { // FASTQ
                        if writeln!(writer, "@{}", id).is_err() { eprintln!("Write error (FASTQ ID)"); return; }
                        if writeln!(writer, "{}", original_seq).is_err() { eprintln!("Write error (FASTQ Seq)"); return; }
                        if writeln!(writer, "+").is_err() { eprintln!("Write error (FASTQ Plus)"); return; }
                        if writeln!(writer, "{}", qual).is_err() { eprintln!("Write error (FASTQ Qual)"); return; }
                    } else { // FASTA
                        if writeln!(writer, ">{}", id).is_err() { eprintln!("Write error (FASTA ID)"); return; }
                        if writeln!(writer, "{}", original_seq).is_err() { eprintln!("Write error (FASTA Seq)"); return; }
                    }
                });
            });
        }
    } else { // FASTA
        let reader = fasta::Reader::new(buf_reader_main);
        let results: Vec<_> = reader.records().collect();
        for result_chunk in results.chunks(READ_CHUNK_SIZE) {
            record_chunk.clear();
            for result in result_chunk {
                 match result {
                     Ok(record) => {
                         let id = record.id().to_owned();
                         let seq = String::from_utf8_lossy(record.seq()).into_owned();
                         record_chunk.push(Record::Fasta { id, seq });
                     }
                     Err(e) => { eprintln!("Error reading FASTA record: {}", e); continue; }
                 }
            }
            if record_chunk.is_empty() { break; }
            _total_reads_partitioned += record_chunk.len() as u64;

            pool.install(|| {
                 record_chunk.par_iter().for_each(|record| {
                    let (id, original_seq, _qual_opt) = match record {
                        Record::Fasta { id, seq } => (id.clone(), seq.clone(), None::<String>),
                        Record::Fastq { .. } => unreachable!(),
                    };
                    let counter_ref = counter_shards_opt.as_ref();
                    let fp = if rc_sensitivity {
                        compute_partition_canonique(&original_seq, p, k, base, counter_ref, trailing_zeros)
                    } else {
                        compute_partition(&original_seq, p, k, base, counter_ref, trailing_zeros)
                    };
                    let fp_usize = fp as usize;
                     if fp_usize >= files.len() { eprintln!("Error: Fingerprint {} out of bounds ({})", fp, files.len()); return; }
                    let writer_mutex = &files[fp_usize];
                    let mut writer = writer_mutex.lock().expect("Partition writer mutex poisoned");
                    if writeln!(writer, ">{}", id).is_err() { eprintln!("Write error (FASTA ID)"); return; }
                    if writeln!(writer, "{}", original_seq).is_err() { eprintln!("Write error (FASTA Seq)"); return; }
                 });
            });
        }
    }

    println!("Flushing and closing partition writers...");
    for (idx, writer_mutex) in files.iter().enumerate() {
        let mut writer_guard = writer_mutex.lock().expect("Partition writer mutex poisoned during close");
        let dummy_file = File::create("/dev/null").expect("Cannot open /dev/null");
        let dummy_encoder = Encoder::new(dummy_file, 0).expect("Failed to create dummy encoder");
        let mut writer = std::mem::replace(&mut *writer_guard, BufWriter::new(dummy_encoder));
        if let Err(e) = writer.flush() { eprintln!("Flush error for partition writer {}: {}", idx, e); }
        match writer.into_inner() {
            Ok(encoder) => { if let Err(e) = encoder.finish() { eprintln!("Error finishing encoder for partition {}: {}", idx, e); } }
            Err(e) => { eprintln!("Error retrieving inner encoder for partition {} (potential unflushed data): {}", idx, e); }
        }
    }
     println!("Finished closing partition writers.");

    Ok(())
}



/// Concatenates bucket files using Gray code order.
fn concat_bucket_files(filename:&str, args: &Args, p: usize, comp_level: i32) -> std::io::Result<()> {
    println!("Concatenating partitions for base: {}", filename);
    let gray_order = generate_gray_code_order(p);

    if let Some(parent_dir) = Path::new(filename).parent() { fs::create_dir_all(parent_dir)?; }

    let final_file = File::create(filename)?;
    let mut encoder = Encoder::new(final_file, comp_level).expect("Cannot create final zstd encoder");

    let num_threads = if args.thread > 0 { args.thread } else { num_cpus::get() };
    if num_threads > 0 {
        if encoder.set_parameter(CParameter::NbWorkers(num_threads as u32)).is_ok() {
             println!("Using {} threads for ZSTD compression of: {}", num_threads, filename);
        } else { eprintln!("Warning: Failed to set {} threads for ZSTD compression.", num_threads); }
    }

    let mut final_writer = BufWriter::new(encoder);
    let mut total_bytes_written: u64 = 0;

    for partition_idx in gray_order {
        let part_filename = get_filename_partition(filename, partition_idx, p);
        match File::open(&part_filename) {
            Ok(file) => {
                let mut decoder = match Decoder::new(BufReader::new(file)) {
                     Ok(d) => d,
                     Err(e) => { eprintln!("Error creating decoder for partition {}: {}. Skipping.", part_filename, e); continue; }
                };
                match std::io::copy(&mut decoder, &mut final_writer) {
                     Ok(bytes_copied) => { total_bytes_written += bytes_copied; }
                     Err(e) => { eprintln!("Error concatenating partition {}: {}", part_filename, e); }
                 }
            }
            Err(ref e) if e.kind() == std::io::ErrorKind::NotFound => {} // Okay if empty
            Err(e) => { eprintln!("Warning: Could not open partition file {}: {}", part_filename, e); }
        }
    }

    println!("Finished reading partitions. Flushing final writer for: {}", filename);
    final_writer.flush()?;
    match final_writer.into_inner() {
        Ok(encoder) => {
            match encoder.finish() {
                Ok(_) => println!("Successfully wrote {} bytes to {}", total_bytes_written, filename),
                Err(e) => eprintln!("Error finishing final zstd stream for {}: {}", filename, e),
            }
        }
        Err(e) => eprintln!("Error retrieving inner final zstd encoder for {}: {}", filename, e),
    }
    Ok(())
}


/// Generates the filename for a specific partition.
fn get_filename_partition(base_filename: &str, partition: usize, p: usize) -> String {
    let path = Path::new(base_filename);
    let parent = path.parent().unwrap_or_else(|| Path::new(""));
    let file_name_osstr = path.file_name().unwrap_or_default();
    let file_name_str = file_name_osstr.to_string_lossy();
    let (file_stem, extensions) = match file_name_str.find('.') {
        Some(idx) => (&file_name_str[..idx], &file_name_str[idx..]),
        None => (&file_name_str[..], ""),
    };
    let partition_id = format!("{:0width$b}", partition, width = p);
    let new_file_name = format!("{}_{}{}", file_stem, partition_id, extensions);
    parent.join(new_file_name).to_string_lossy().into_owned()
}


/// Recursive function to compute partitions. Accepts Option for sharded counter filter.
fn compute_all_file(
    filename_input: &str,          // Input for this level
    filename_output_base: &str,   // Base name for output/partitions
    args: &Args,
    k: usize,
    level: usize,
    counter_shards_opt: Option<Arc<Vec<Arc<Mutex<KmerCounterFilter>>>>>,
    multi: &MultiProgress,
    bars: &mut Vec<Option<ProgressBar>>,
    pool: &ThreadPool,
) -> std::io::Result<()> {
    let p = args.p[level];
    let base_for_level = BASES[level % BASES.len()];

    println!("Starting level {}, p={}, input='{}', output_base='{}'", level, p, filename_input, filename_output_base);

    // Pass the Option<Arc<Vec<...>>> for counter shards
    create_bucket_files(filename_input, filename_output_base, args, p, k, base_for_level, counter_shards_opt.clone())?;

    if level > 0 {
         println!("Removing intermediate input file: {}", filename_input);
         if let Err(e) = fs::remove_file(filename_input) {
            if Path::new(filename_input).exists() {
                eprintln!("Warning: could not remove intermediate input {}: {}", filename_input, e);
            }
         }
    }

    let num_partitions = 1usize << p;
    if level + 1 < args.p.len() { // Recurse
        println!("Level {}: Recursing into {} partitions...", level, num_partitions);
        let bar = if let Some(pb) = &bars[level] { pb.clone() }
                  else {
                      let pb = multi.add(ProgressBar::new(num_partitions as u64));
                      pb.set_prefix(format!("Level {} (p={})", level, p));
                      pb.set_style(ProgressStyle::default_bar().template(PBSTYLE[level % PBSTYLE.len()]).unwrap().progress_chars("##-"));
                      bars[level] = Some(pb.clone());
                      pb
                  };
        bar.reset();
        bar.set_length(num_partitions as u64);

        pool.install(|| {
            (0..num_partitions).into_par_iter().for_each(|i| {
                let mut local_bars = vec![None; args.p.len()];
                let filename_partition_input = get_filename_partition(filename_output_base, i, p);

                let should_recurse = match fs::metadata(&filename_partition_input) { Ok(meta) => meta.len() > 0, Err(_) => false };

                if should_recurse {
                    if let Err(e) = compute_all_file(
                        &filename_partition_input, &filename_partition_input, args, k, level + 1,
                        counter_shards_opt.clone(), multi, &mut local_bars, pool ) {
                        eprintln!("Error processing partition {} at level {}: {}", i, level + 1, e);
                    }
                }
                bar.inc(1);
            });
        });
        bar.finish_with_message(format!("Level {} (p={}) finished recursion.", level, p));

    } else { // Base case: Last level of partitioning
        println!("Level {}: Reached final partitioning level. Processing partitions for RC update...", level);
        let bar = multi.add(ProgressBar::new(num_partitions as u64));
        bar.set_prefix(format!("RC Update Level {} (p={})", level, p));
        bar.set_style(ProgressStyle::default_bar().template(PBSTYLE[level % PBSTYLE.len()]).unwrap().progress_chars("##-"));
        bar.set_length(num_partitions as u64);

        // Call process_partition_for_rc (no filter needed)
        if args.rc_compression_loop > 0 {
            pool.install(|| {
                (0..num_partitions).into_par_iter().for_each(|i| {
                    let filename_partition = get_filename_partition(filename_output_base, i, p);
                    if Path::new(&filename_partition).exists() {
                        if let Err(e) = process_partition_for_rc(&filename_partition, args, k, pool) {
                            eprintln!("Error during RC processing for partition {}: {}", filename_partition, e);
                        }
                    }
                    bar.inc(1);
                });
            });
            bar.finish_with_message(format!("Level {} (p={}) RC update finished.", level, p));
        } else {
             println!("Skipping RC update pass for level {} (rc_compression_loop = 0).", level);
             bar.finish_with_message(format!("Level {} (p={}) RC update skipped.", level, p));
        }
    }

    // Concatenate and Cleanup
    let concat_compression_level = if level == 0 { args.final_compression_level } else { args.compression_level };
    concat_bucket_files(filename_output_base, args, p, concat_compression_level)?;

    println!("Cleaning up partitions for level {} (base: {})", level, filename_output_base);
    let cleanup_errors = Arc::new(Mutex::new(Vec::new()));
    (0..num_partitions).into_par_iter().for_each(|i| {
        let filename_partition = get_filename_partition(filename_output_base, i, p);
        if let Err(e) = fs::remove_file(&filename_partition) {
            if Path::new(&filename_partition).exists() {
                 eprintln!("Warning: could not remove partition {}: {}", filename_partition, e);
                 cleanup_errors.lock().expect("Cleanup mutex poisoned").push(filename_partition);
            }
        }
    });
     let errors = cleanup_errors.lock().expect("Cleanup mutex poisoned");
     if !errors.is_empty() { eprintln!("Failed to remove {} partition files for level {}.", errors.len(), level); }

    Ok(())
}


fn main() -> std::io::Result<()> {
    let args = Args::parse();
    println!("Running with arguments: {:?}", args); // Print args for debugging
    let start_total = Instant::now();
    fs::create_dir_all(&args.output).expect("Cannot create output directory");

    let k = args.k; // Use default or provided k directly

    let num_threads = if args.thread > 0 { args.thread } else { num_cpus::get() };
    println!("Using {} threads for processing.", num_threads);

    if args.p.is_empty() {
        eprintln!("Error: No p values provided. Use -p argument.");
        return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "No p values provided"));
    }
     for (i, &p_val) in args.p.iter().enumerate() {
        if p_val == 0 {
             eprintln!("Error: p value at level {} cannot be 0.", i);
             return Err(std::io::Error::new(std::io::ErrorKind::InvalidInput, "p value cannot be 0"));
        }
        if p_val > 32 { eprintln!("Warning: p value {} at level {} is large, may lead to many partitions.", p_val, i); }
    }

    let pool = ThreadPoolBuilder::new().num_threads(num_threads).build().unwrap();

    let current_input_file = args.input.clone(); // Start with the original input

    // --- Pipeline Steps ---

    // 1. Optional Counting Filter Pass
    let counter_shards_opt: Option<Arc<Vec<Arc<Mutex<KmerCounterFilter>>>>> = if args.use_counting_filter {
        println!("Counting filter pass is enabled.");
        Some(run_counting_filter_pass(&current_input_file, k, &args, &pool)?)
    } else {
        println!("Counting filter pass is disabled.");
        None
    };

    // 2. Multi-level Partitioning Pass (RC pass is now inside this)
    println!("Starting multi-level partitioning...");
    let multi_pb = MultiProgress::new();
    let mut bars = vec![None; args.p.len()];

    let final_filename_stem = format!("final_k{}_p{}", k, args.p.iter().map(|v| v.to_string()).collect::<Vec<_>>().join("-"));
    let final_extension = match args.final_compression.to_lowercase().as_str() { "gzip" | "gz" => "gz", _ => "zst" };
    let final_filename = PathBuf::from(&args.output).join(format!("{}.{}", final_filename_stem, final_extension)).to_string_lossy().into_owned();
    println!("Final output will be: {}", final_filename);

    // Pass the Option<Arc<Vec<...>>> of counter shards
    compute_all_file(&current_input_file, &final_filename, &args, k, 0, counter_shards_opt, &multi_pb, &mut bars, &pool)?;

    // --- Final Statistics ---
    let final_meta = match fs::metadata(&final_filename) {
        Ok(meta) => meta,
        Err(e) => {
            eprintln!("Error getting metadata for final file {}: {}", final_filename, e);
            println!("\n--- Processing Finished (with errors) ---");
            println!("Total run time: {:.3} seconds", start_total.elapsed().as_secs_f64());
            println!("Final file {} may be incomplete or missing.", final_filename);
            return Err(e);
        }
    };
    let final_size = final_meta.len();
    let final_size_mb = final_size as f64 / (1024.0 * 1024.0);
    let elapsed_total = start_total.elapsed().as_secs_f64();

    println!("\n--- Processing Complete ---");
    println!("Final archive size: {:.3} MB", final_size_mb);
    println!("Total run time: {:.3} seconds", elapsed_total);
    println!("Final compressed file: {}", final_filename);
    // Updated parameter printout
    print!("Parameters: k={}, p=[{}], RC sensitive={}, RC loops={}",
        k, args.p.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(", "),
        args.rc_sensitivity, args.rc_compression_loop);
    if args.use_counting_filter {
        print!(", Counting Filter enabled (counters={}, shards={})",
            args.filter_counters, NUM_COUNTER_SHARDS);
    } else {
        print!(", Counting Filter disabled");
    }
    println!(", Subsampling Zeros={}", args.trailing_zeros);
    println!("Note: K-mer hashing used homopolymer-compressed sequences.");

    Ok(())
}
