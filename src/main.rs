use clap::Parser;
// use ahash::AHashMap;
use bio::io::{fasta, fastq};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression as GzCompression;
use zstd::stream::{Encoder, Decoder};
use zstd::stream::raw::CParameter;
use crossbeam_channel;
use num_cpus;
use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicU64, Ordering};
use std::fs::{File, self};
use std::io::{BufReader, BufWriter, BufRead, Read, Write};
use std::time::Instant;

use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

use std::path::{Path};

use indicatif::{MultiProgress, ProgressBar, ProgressStyle};


use nthash::*;

/// Command-line arguments.
#[derive(Parser)]
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

    /// k-mer length. If not provided, k is automatically chosen based on the first 1000 sequences.
    #[arg(short, long)]
    k: Option<usize>,

    /// zstd compression level for writing partition files.
    #[arg(short, long, default_value = "4")]
    compression_level: i32,

    /// Final compression algorithm ("zstd" or "gzip"). (Default: "zstd")
    #[arg(long, default_value = "zstd")]
    final_compression: String,

    /// Final compression level. For zstd this is typically in the range [1, 22]; for gzip in [0, 9].
    #[arg(long, default_value = "19")]
    final_compression_level: i32,

    // /// Enable second-level partitioning + reordering (using a different hash function, same p/k).
    // #[arg(long)]
    // second_level: bool,

    /// Enable reverse complement sensitivity
    #[arg(long)]
    rc_sensitivity: bool,

    #[arg(long, default_value = "1")]
    rc_compression_loop: i32,

    #[arg(short, long, default_value = "0")]
    thread: usize,
}

/// A record (FASTA or FASTQ)
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



// /// 64-bit multipliers for rolling hash (different for first- and second-level).
// const BASE1: u64 = 0x9e3779b97f4a7c15;
// const BASE2: u64 = 0xc2b2ae3d27d4eb4f;
const BASES: [u64; 9] = [
    0x9e3779b97f4a7c15, // golden ratio
    0xc2b2ae3d27d4eb4f, // murmur3 constant
    0x165667b19e3779f9, // splitmix64 constant
    0x27d4eb2f165667c5, // reversed Murmur3 + variation
    0xa0761d6478bd642f, // wyhash constant
    0xe7037ed1a0b428db, // wyhash v4 constant
    0xbf58476d1ce4e5b9, // splitmix64 constant
    0x94d049bb133111eb, // splitmix64 constant
    0x2545f4914f6cdd1d, // LCG constant from PCG
];

/// Computes the partition fingerprint for a sequence using the given BASE multiplier.
/// We maintain an array of length P. For each k-mer:
///   1. Compute a rolling hash.
///   2. Use top log₂(P) bits for a bucket index.
///   3. Keep the min hashed value in each bucket.
/// Then build the P-bit fingerprint from the LSBs of each bucket.
fn compute_partition(seq: &str, p: usize, k: usize, base: u64) -> u64 {
    let seq_bytes = seq.as_bytes();
    let index_bits = (p as f64).log2().ceil() as usize;
    let buckets_count = 1 << index_bits; // use an array of length p
    let mut buckets: Vec<Option<u64>> = vec![None; buckets_count];

    let iter = NtHashForwardIterator::new(seq_bytes, k).expect("NtHash problem");
    for it in iter {
        let hash = it.wrapping_mul(base);
        let idx = (hash >> (64 - index_bits)) as usize;
        let value = hash & ((1u64 << (64 - index_bits)) - 1);
        buckets[idx] = Some(match buckets[idx] {
            Some(current) if value > current => current,
            _ => value,
        });
    }
    let mut fingerprint = 0;
    for bucket in buckets {
        let bit = bucket.unwrap_or(0) & 1;
        fingerprint = (fingerprint << 1) | bit;
    }
    fingerprint & ((1 << p) -1 )
}

fn compute_partition_canonique(seq: &str, p: usize, k: usize, base: u64) -> u64 {
    let seq_bytes = seq.as_bytes();
    let index_bits = (p as f64).log2().ceil() as usize;
    let buckets_count = 1 << index_bits; // use an array of length p
    let mut buckets: Vec<Option<u64>> = vec![None; buckets_count];

    let iter = NtHashIterator::new(seq_bytes, k).expect("NtHash problem");
    for it in iter {
        let hash = it.wrapping_mul(base);
        let idx = (hash >> (64 - index_bits)) as usize;
        let value = hash & ((1u64 << (64 - index_bits)) - 1);
        buckets[idx] = Some(match buckets[idx] {
            Some(current) if value > current => current,
            _ => value,
        });
    }
    let mut fingerprint = 0;
    for bucket in buckets {
        let bit = bucket.unwrap_or(0) & 1;
        fingerprint = (fingerprint << 1) | bit;
    }
    fingerprint & ((1 << p) -1 )
}

/// Generate a Gray code sequence for 'bits' bits.
/// Returns a vector of length 2^(bits) with the Gray code order.
fn generate_gray_code_order(bits: usize) -> Vec<usize> {
    let n = 1 << bits;
    let mut order = Vec::with_capacity(n);
    for i in 0..n {
        order.push(i ^ (i >> 1));
    }
    order
}

/// Reads up to the first 1000 records from the input file (FASTA or FASTQ),
/// finds the longest read length L, and returns the smallest k such that 4^k >= 10 * L.
fn auto_detect_k(input: &str) -> std::io::Result<usize> {
    let input_reader = open_input(input);
    let mut buf_reader = BufReader::new(input_reader);
    let _ = buf_reader.fill_buf()?;
    let peek = buf_reader.fill_buf()?;
    let is_fastq = !peek.is_empty() && peek[0] == b'@';
    let mut max_len = 0;
    let mut count = 0;
    if is_fastq {
        let reader = fastq::Reader::new(buf_reader);
        for result in reader.records() {
            let record = result.map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
            let len = record.seq().len();
            if len > max_len {
                max_len = len;
            }
            count += 1;
            if count >= 1000 {
                break;
            }
        }
    } else {
        let reader = fasta::Reader::new(buf_reader);
        for result in reader.records() {
            let record = result.map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
            let len = record.seq().len();
            if len > max_len {
                max_len = len;
            }
            count += 1;
            if count >= 1000 {
                break;
            }
        }
    }
    let threshold = 10 * max_len;
    let k = (threshold as f64).ln() / 4f64.ln();
    Ok(k.ceil() as usize)
}

/// Concatenate all partition files in Gray code order into a final compressed file.
/// The final compression is chosen based on `final_comp`.
fn concatenate_partitions(
    out_dir: &str,
    p: usize,
    final_filename: &str,
    final_comp: &str,
    final_comp_level: i32,
) -> std::io::Result<()> {
    let gray_order = generate_gray_code_order(p);
    if final_comp.to_lowercase() == "gzip" {
        let final_file = File::create(final_filename)?;
        let gz_encoder = GzEncoder::new(final_file, GzCompression::new(final_comp_level as u32));
        let mut final_writer = BufWriter::new(gz_encoder);
        for partition in gray_order {
            let partition_id = format!("{:0width$b}", partition, width = p);
            let part_filename = format!("{}/partition_{}.fa.zst", out_dir, partition_id);
            if let Ok(file) = File::open(&part_filename) {
                let mut decoder = Decoder::new(file).expect("Cannot create zstd decoder for partition");
                let mut buffer = [0u8; 8192];
                loop {
                    let n = decoder.read(&mut buffer)?;
                    if n == 0 { break; }
                    final_writer.write_all(&buffer[..n])?;
                }
            }
        }
        final_writer.flush()?;
        let gz_encoder = final_writer.into_inner().expect("Error retrieving inner gzip encoder");
        gz_encoder.finish()?;
    } else {
        let final_file = File::create(final_filename)?;
        let mut encoder = Encoder::new(final_file, final_comp_level)
            .expect("Cannot create final zstd encoder");
        encoder
            .set_parameter(CParameter::NbWorkers(num_cpus::get() as u32))
            .expect("Failed to set number of threads");
        let mut final_writer = BufWriter::new(encoder);
        for partition in gray_order {
            let partition_id = format!("{:0width$b}", partition, width = p);
            let part_filename = format!("{}/partition_{}.fa.zst", out_dir, partition_id);
            if let Ok(file) = File::open(&part_filename) {
                let mut decoder = Decoder::new(file).expect("Cannot create zstd decoder for partition");
                let mut buffer = [0u8; 8192];
                loop {
                    let n = decoder.read(&mut buffer)?;
                    if n == 0 { break; }
                    final_writer.write_all(&buffer[..n])?;
                }
            }
        }
        final_writer.flush()?;
        match final_writer.into_inner() {
            Ok(encoder) => {
                if let Err(e) = encoder.finish() {
                    eprintln!("Error finishing final zstd stream: {}", e);
                }
            }
            Err(_) => eprintln!("Error retrieving inner final zstd encoder"),
        }
    }
    Ok(())
}

/// Remove partition files after final concatenation.
fn remove_partition_files(out_dir: &str, p: usize) {
    let num_partitions = 1 << p;
    for i in 0..num_partitions {
        let partition_id = format!("{:0width$b}", i, width = p);
        let filename = format!("{}/partition_{}.fa.zst", out_dir, partition_id);
        if let Err(e) = fs::remove_file(&filename) {
            eprintln!("Warning: could not remove {}: {}", filename, e);
        }
    }
}

fn reverse_complement(dna: &str) -> String {
    dna.chars()
        .rev() // On inverse la séquence
        .map(|n| match n {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => n, // Si un caractère inconnu est rencontré, on le laisse inchangé
        })
        .collect()
}
fn update_kmer_dict<F>(seq: &str, k: usize, dict: &mut Vec<u64>, f: &mut F)
where
    F: FnMut(&mut Vec<u64>, usize),
{
    let bytes = seq.as_bytes();
    if bytes.len() < k {
        return;
    }
    // Create an iterator over k-mer hashes.
    let iter = NtHashForwardIterator::new(bytes, k).expect("NtHash error");
    for kmer in iter {
        f(dict, kmer as usize);
    }
}
// --- The Updated update_rc_file Function ---

/// Process the file once and perform the reverse-complement update internally in `rc_loops` passes.
/// For FASTQ, the quality string is reversed if the record is flipped.
fn update_rc_file(input: &str, args: &Args, k: usize) -> std::io::Result<()> {
    let input_reader = open_input(input);
    let rc_loops = args.rc_compression_loop;
    let mut buf_reader = BufReader::new(input_reader);
    let compression_level = args.compression_level;

    // Determine whether the input is FASTQ (first byte '@') or FASTA.
    let peek = buf_reader.fill_buf()?;
    let is_fastq = !peek.is_empty() && peek[0] == b'@';

    // Read and convert all records from the input file into memory.
    let mut records: Vec<Record> = Vec::new();
    // Build the initial k-mer dictionary.
    let size_array : usize = 1 << k;
    let mut kmer_dict : Vec<u64> = vec![0; size_array];

    // Define update closures.
    let mut incr = |dict: &mut Vec<u64>, key: usize| {
        dict[key % (size_array as u64) as usize] += 1;
    };
    let mut decr = |dict: &mut Vec<u64>, key: usize| {
        dict[key % (size_array as u64) as usize] -= 1;
    };

    if is_fastq {
        let reader = fastq::Reader::new(buf_reader);
        for result in reader.records() {
            let rec = result.map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
            // Convert the bio::io::fastq::Record to our custom Record.
            let id = rec.id().to_owned();
            let seq = String::from_utf8_lossy(rec.seq()).into_owned();
            let qual = String::from_utf8_lossy(rec.qual()).into_owned();
            // Update dictionary with forward sequence.
            update_kmer_dict(&seq, k, &mut kmer_dict, &mut incr);
            records.push(Record::Fastq { id, seq, qual });
        }
    } else {
        let reader = fasta::Reader::new(buf_reader);
        for result in reader.records() {
            let rec = result.map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
            let id = rec.id().to_owned();
            let seq = String::from_utf8_lossy(rec.seq()).into_owned();
            update_kmer_dict(&seq, k, &mut kmer_dict, &mut incr);
            records.push(Record::Fasta { id, seq });
        }
    }

    // Now perform the reverse-complement update rc_loops times.
    for _ in 0..rc_loops {
        for record in records.iter_mut() {
            match record {
                Record::Fastq { ref mut seq, ref mut qual, .. } => {
                    // Compute the reverse complement once.
                    let rev_seq = reverse_complement(seq);
                    // Create rolling hash iterators using nthash.
                    let fwd_iter = NtHashForwardIterator::new(seq.as_bytes(), k)
                        .expect("Error creating forward k-mer iterator");
                    let rev_iter = NtHashForwardIterator::new(rev_seq.as_bytes(), k)
                        .expect("Error creating reverse k-mer iterator");
                    // Compute both scores in a single pass.
                    let (mut fwd_score, mut rev_score) = (0, 0);
                    for (fwd_kmer, rev_kmer) in fwd_iter.zip(rev_iter) {
                        fwd_score += kmer_dict[(fwd_kmer % (size_array as u64)) as usize];
                        rev_score += kmer_dict[(rev_kmer % (size_array as u64)) as usize];
                    }
                    // If the reverse complement yields a higher score, flip the record.
                    if fwd_score < rev_score {
                        // Remove the forward k-mer counts.
                        update_kmer_dict(seq, k, &mut kmer_dict, &mut decr);
                        // Update the sequence with its reverse complement.
                        *seq = rev_seq;
                        // For FASTQ, reverse the quality string.
                        *qual = qual.chars().rev().collect();
                        // Update the dictionary with the k-mer counts for the new orientation.
                        update_kmer_dict(seq, k, &mut kmer_dict, &mut incr);
                    }
                }
                Record::Fasta { ref mut seq, .. } => {
                    let rev_seq = reverse_complement(seq);
                    let fwd_iter = NtHashForwardIterator::new(seq.as_bytes(), k)
                        .expect("Error creating forward k-mer iterator");
                    let rev_iter = NtHashForwardIterator::new(rev_seq.as_bytes(), k)
                        .expect("Error creating reverse k-mer iterator");
                    let (mut fwd_score, mut rev_score) = (0, 0);
                    for (fwd_kmer, rev_kmer) in fwd_iter.zip(rev_iter) {
                        fwd_score += kmer_dict[(fwd_kmer % (size_array as u64)) as usize];
                        rev_score += kmer_dict[(rev_kmer % (size_array as u64)) as usize];
                    }
                    if fwd_score < rev_score {
                        update_kmer_dict(seq, k, &mut kmer_dict, &mut decr);
                        *seq = rev_seq;
                        update_kmer_dict(seq, k, &mut kmer_dict, &mut incr);
                    }
                }
            }
        }
    }

    // Write the final records out to a temporary file.
    let temp_path = input.to_owned() + ".temp";
    let file = File::create(&temp_path)
        .unwrap_or_else(|e| panic!("Cannot create {}: {}", temp_path, e));
    let encoder = Encoder::new(file, compression_level).expect("Cannot create zstd encoder for partition");
    let mut writer = BufWriter::new(encoder);

    if is_fastq {
        for record in &records {
            if let Record::Fastq { id, seq, qual } = record {
                writeln!(writer, "@{}", id)?;
                writeln!(writer, "{}", seq)?;
                writeln!(writer, "+")?;
                writeln!(writer, "{}", qual)?;
            }
        }
    } else {
        for record in &records {
            if let Record::Fasta { id, seq } = record {
                writeln!(writer, ">{}", id)?;
                writeln!(writer, "{}", seq)?;
            }
        }
    }
    writer.flush()?;
    let encoder = writer
        .into_inner()
        .unwrap_or_else(|_| panic!("Error retrieving inner encoder"));
    if let Err(e) = encoder.finish() {
        eprintln!("Error finishing zstd stream: {}", e);
    }
    // For FASTA files, rename the temporary file to overwrite the original.
    if !is_fastq {
        fs::rename(temp_path, input).expect("Error renaming temporary file");
    }
    Ok(())
}



fn create_bucket_files(filename_input:&str, filename_comp:&str, args: &Args, p: usize, k: usize, base: u64) -> std::io::Result<()> {
    let input = open_input(filename_input);
    let mut buf_reader = BufReader::new(input);
    let peek = buf_reader.fill_buf().expect("Error peeking input");
    let is_fastq = !peek.is_empty() && peek[0] == b'@';
    let compression_level = args.compression_level;
    // let p = args.p[0];
    let thread = args.thread;
    let rc_sensitivity = args.rc_sensitivity;
    let outdir = args.output.clone();
    let size_array = 1 << p;

    let mut writers: Vec<Arc<Mutex<BufWriter<Encoder<File>>>>> = Vec::with_capacity(size_array);

    for fp in 0..size_array {
        let filename_partition = get_filename_partition(filename_comp,fp,p);
        let file = File::create(&filename_partition)
            .unwrap_or_else(|e| panic!("Cannot create {}: {}", filename_partition, e));

        let encoder = Encoder::new(file, compression_level)
            .expect("Cannot create zstd encoder for partition");
        let writer = BufWriter::new(encoder);
        writers.push(Arc::new(Mutex::new(writer)));
    }

    let files = Arc::new(writers);

    let (record_tx, record_rx) = crossbeam_channel::bounded::<Record>(100);
    let num_workers = if thread > 0 {
        std::cmp::min(thread, num_cpus::get())
    } else {
        num_cpus::get()
    };
    let worker_handles: Vec<_> = (0..num_workers)
        .map(|_| {
            let record_rx = record_rx.clone();
            let p = p;
            let k = k;
            let files_clone = Arc::clone(&files);
            let rc_sensitivity = rc_sensitivity;
            std::thread::spawn(move || {
                while let Ok(record) = record_rx.recv() {
                    match record {
                        Record::Fasta { id, seq } => {
                            let fp = if rc_sensitivity {
                                compute_partition_canonique(&seq, p, k, base)
                            } else {
                                compute_partition(&seq, p, k, base)
                            };

                            let writer_mutex = &files_clone[fp as usize];
                            let mut writer = writer_mutex.lock().unwrap();
                            writeln!(writer, ">{}", id).expect("Write error");
                            writeln!(writer, "{}", seq).expect("Write error");
                        }
                        Record::Fastq { id, seq, qual } => {
                            let fp = if rc_sensitivity {
                                compute_partition_canonique(&seq, p, k, base)
                            } else {
                                compute_partition(&seq, p, k, base)
                            };

                            let writer_mutex = &files_clone[fp as usize];
                            let mut writer = writer_mutex.lock().unwrap();
                            writeln!(writer, "@{}", id).expect("Write error");
                            writeln!(writer, "{}", seq).expect("Write error");
                            writeln!(writer, "+").expect("Write error");
                            writeln!(writer, "{}", qual).expect("Write error");
                        }
                    }
                }
            })
        })
        .collect();



    if is_fastq {
        let reader = fastq::Reader::new(buf_reader);
        for result in reader.records() {
            let rec = result.map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
            let id = rec.id().to_owned();
            let seq = std::str::from_utf8(rec.seq()).unwrap_or("").to_owned();
            let qual = std::str::from_utf8(rec.qual()).unwrap_or("").to_owned();
            record_tx.send(Record::Fastq { id, seq, qual }).expect("Failed to send record");
        }
    } else {
        let reader = fasta::Reader::new(buf_reader);
        for result in reader.records() {
            let rec = result.map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
            let id = rec.id().to_owned();
            let seq = std::str::from_utf8(rec.seq()).unwrap_or("").to_owned();
            record_tx.send(Record::Fasta { id, seq }).expect("Failed to send record");
        }
    }
    drop(record_tx);

    for handle in worker_handles {
        handle.join().expect("Worker thread panicked");
    }

    for writer_mutex in &*files {
        let mut writer_guard = writer_mutex.lock().unwrap();

        let null_file = File::create("/dev/null").expect("Cannot open /dev/null");
        let null_encoder = Encoder::new(null_file, 0).expect("Failed to create dummy encoder");

        let mut writer = std::mem::replace(&mut *writer_guard, BufWriter::new(null_encoder));

        if let Err(e) = writer.flush() {
            eprintln!("Flush error: {}", e);
        }

        match writer.into_inner() {
            Ok(encoder) => {
                if let Err(e) = encoder.finish() {
                    eprintln!("Error finishing encoder: {}", e);
                }
            }
            Err(_) => {
                eprintln!("Error retrieving inner encoder");
            }
        }
    }

    Ok(())
}

fn update_level(filename:&str, args: &Args, k: usize, base: u64) -> std::io::Result<()> {
    let input = filename;
    let compression_level = args.compression_level;
    let p = args.p[0];
    let rc_sensitivity = args.rc_sensitivity;
    let input_reader = open_input(input);
    let mut buf_reader = BufReader::new(input_reader);

    // Determine whether the input is FASTQ (first byte '@') or FASTA.
    let peek = buf_reader.fill_buf()?;
    let is_fastq = !peek.is_empty() && peek[0] == b'@';

    // Read and convert all records from the input file into memory.
    let mut records: Vec<Record> = Vec::new();

    if is_fastq {
        let reader = fastq::Reader::new(buf_reader);
        for result in reader.records() {
            let rec = result.map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
            let id = rec.id().to_owned();
            let seq = String::from_utf8_lossy(rec.seq()).into_owned();
            let qual = String::from_utf8_lossy(rec.qual()).into_owned();
            records.push(Record::Fastq { id, seq, qual });
        }
    } else {
        let reader = fasta::Reader::new(buf_reader);
        for result in reader.records() {
            let rec = result.map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
            let id = rec.id().to_owned();
            let seq = String::from_utf8_lossy(rec.seq()).into_owned();
            records.push(Record::Fasta { id, seq });
        }
    }




    let mut subpartitions: Vec<Vec<Record>> = (0..(1 << p)).map(|_| Vec::new()).collect();
    for record in records {
        let (seq, _) = match &record {
            Record::Fasta { seq, .. } => (seq, false),
            Record::Fastq { seq, .. } => (seq, true),
        };
        let fp2;
        if rc_sensitivity {
            fp2 = compute_partition_canonique(seq, p, k, base) as usize;
        } else {
            fp2 = compute_partition(seq, p, k, base) as usize;
        }
        subpartitions[fp2].push(record);
    }
    let temp_path = input.to_owned() + ".temp";
    let file = File::create(&temp_path)
        .unwrap_or_else(|e| panic!("Cannot create {}: {}", temp_path, e));
    let encoder = Encoder::new(file, compression_level)
        .expect("Cannot create zstd encoder for partition");
    let mut writer = BufWriter::new(encoder);
    let gray_order2 = generate_gray_code_order(p);
    for idx in gray_order2 {
        for rec in &subpartitions[idx] {
            match rec {
                Record::Fasta { id, seq } => {
                    writeln!(writer, ">{}", id).expect("Write error");
                    writeln!(writer, "{}", seq).expect("Write error");
                }
                Record::Fastq { id, seq, qual } => {
                    writeln!(writer, "@{}", id).expect("Write error");
                    writeln!(writer, "{}", seq).expect("Write error");
                    writeln!(writer, "+").expect("Write error");
                    writeln!(writer, "{}", qual).expect("Write error");
                }
            }
        }
    }
    writer.flush()?;
    let encoder = writer
        .into_inner()
        .unwrap_or_else(|_| panic!("Error retrieving inner encoder"));
    if let Err(e) = encoder.finish() {
        eprintln!("Error finishing zstd stream: {}", e);
    }
    // For FASTA files, rename the temporary file to overwrite the original.
    if !is_fastq {
        fs::rename(temp_path, input).expect("Error renaming temporary file");
    }
    Ok(())
}






fn concat_bucket_files(filename:&str, args: &Args, p: usize, comp_level: i32) -> std::io::Result<()> {
    let gray_order = generate_gray_code_order(p);
    let final_file = File::create(filename)?;
    let mut encoder = Encoder::new(final_file, comp_level)
        .expect("Cannot create final zstd encoder");
    encoder
        .set_parameter(CParameter::NbWorkers(num_cpus::get() as u32))
        .expect("Failed to set number of threads");
    let mut final_writer = BufWriter::new(encoder);
    for partition in gray_order {
        let part_filename = get_filename_partition(filename,partition,p);
        if let Ok(file) = File::open(&part_filename) {
            let mut decoder = Decoder::new(file).expect("Cannot create zstd decoder for partition");
            let mut buffer = [0u8; 8192];
            loop {
                let n = decoder.read(&mut buffer)?;
                if n == 0 { break; }
                final_writer.write_all(&buffer[..n])?;
            }
        }
    }
    final_writer.flush()?;
    match final_writer.into_inner() {
        Ok(encoder) => {
            if let Err(e) = encoder.finish() {
                eprintln!("Error finishing final zstd stream: {}", e);
            }
        }
        Err(_) => eprintln!("Error retrieving inner final zstd encoder"),
    }
    Ok(())
}


fn get_filename_partition(filename: &str, partition: usize, p: usize) -> String {
    let path = Path::new(filename);

    // On récupère le parent (le dossier)
    let parent = path.parent().unwrap_or_else(|| Path::new(""));

    // On récupère le nom du fichier (ex: monfichier.fa.zst)
    let file_name = path.file_name().unwrap().to_string_lossy();

    // On cherche les extensions
    let parts: Vec<&str> = file_name.split('.').collect();

    let partition_id = format!("{:0width$b}", partition, width = p);

    if parts.len() > 1 {
        // On garde toutes les extensions
        let extensions = &parts[1..];
        // Le nom de base (sans extensions)
        let base = parts[0];

        // Ajoute le suffixe
        let new_base = format!("{}_{}", base, partition_id);
        let new_file_name = format!("{}.{}", new_base, extensions.join("."));

        parent.join(new_file_name).to_string_lossy().into_owned()
    } else {
        // Pas d'extension, on ajoute juste le suffixe
        let new_base = format!("{}_{}", file_name, partition_id);
        parent.join(new_base).to_string_lossy().into_owned()
    }
}


fn compute_all_file(filename_input:&str, filename:&str, args: &Args, k: usize, level: usize, multi: &MultiProgress, bars: &mut Vec<Option<ProgressBar>>) -> std::io::Result<()> {
    let p = args.p[level];
    let _ = create_bucket_files(filename_input, filename, &args, p, k, BASES[level]);
    if level > 0 {
        if let Err(e) = fs::remove_file(filename) {
            eprintln!("Warning: could not remove {}: {}", filename, e);
        }
    }



    if level + 1 < args.p.len()  {
        if bars[level].is_none() {
            let bar = multi.add(ProgressBar::new(1 << p));
            bar.set_prefix(format!("Niveau {}", level));
            bar.set_style(
                ProgressStyle::default_bar()
                    .template("{prefix} [{bar:40.cyan/blue}] {pos}/{len} {elapsed} ETA: {eta}")
                    .unwrap()
                    .progress_chars("##-"),
            );
            bars[level] = Some(bar);
        }
        let bar = bars[level].as_ref().unwrap().clone();
        bar.set_position(0);
        for i in 0..(1 << p) {
            let filename_partition = get_filename_partition(filename,i,p);
            let _ = compute_all_file(&filename_partition, &filename_partition, args, k, level+1, multi, bars);
            bar.inc(1);
        }
        // let pool = ThreadPoolBuilder::new()
        //     .num_threads(args.thread)
        //     .build()
        //     .unwrap();
        //
        // pool.install(|| {
        //     ( 0..(1 << p)).into_par_iter().for_each(|i| {
        //         let filename_partition = get_filename_partition(filename,i,p);
        //         let _ = compute_all_file(&filename_partition, &filename_partition, args, k, level+1);
        //         bar.inc(1);
        //
        //     });
        // });
        if level == 0 {
            bars[level].as_ref().unwrap().finish_with_message("Racine : terminé !");
        }
    } else {
        if args.rc_sensitivity {
            let pool = ThreadPoolBuilder::new()
                .num_threads(args.thread)
                .build()
                .unwrap();

            pool.install(|| {
                ( 0..(1 << p)).into_par_iter().for_each(|i| {
                    let filename_partition = get_filename_partition(filename,i,p);
                    let _ = update_rc_file(&filename_partition, args, k);

                });
            });
        }
    }
    if level == 0 {
        let _ = concat_bucket_files(filename, &args, p, args.final_compression_level);
    } else {
        let _ = concat_bucket_files(filename, &args, p, args.compression_level);
    }
    for i in 0..(1 << p) {
        let filename_partition = get_filename_partition(filename,i,p);
        if let Err(e) = fs::remove_file(filename_partition) {
            eprintln!("Warning: could not remove {}: {}", filename, e);
        }
    }

    Ok(())
}












fn main() -> std::io::Result<()> {
    let args = Args::parse();
    let start = Instant::now();
    fs::create_dir_all(&args.output).expect("Cannot create output directory");

    let k = match args.k {
        Some(val) => val,
        None => {
            let detected = auto_detect_k(&args.input)?;
            println!("Auto-detected k = {}", detected);
            detected
        }
    };

    println!("{:?}",args.p);


    let multi = MultiProgress::new();
    let mut bars = vec![None; args.p.len()];


    let final_filename = format!("{}/final.zst", args.output);
    let _ = compute_all_file(&args.input.clone(), &final_filename, &args, k, 0, &multi, &mut bars);



    // let _ = create_bucket_files(&args, k, BASE1);
    //
    //
    // //
    // //
    let nucleotide_count = Arc::new(AtomicU64::new(0));
    // let num_partitions = 1 << args.p[0];
    //
    //
    // let out_dir = args.output.clone();
    //
    //
    //
    // if args.second_level {
    //
    //     println!("Second Level Activate");
    //
    //     let pool = ThreadPoolBuilder::new()
    //     .num_threads(args.thread)
    //     .build()
    //     .unwrap();
    //
    //     pool.install(|| {
    //         (0..num_partitions).into_par_iter().for_each(|i| {
    //             let partition_id = format!("{:0width$b}", i, width = args.p[0]);
    //             let filename = format!("{}/partition_{}.fa.zst", out_dir, partition_id);
    //             let _ = update_level(&filename, &args, k, BASE2);
    //
    //         });
    //     });
    // }
    //
    // if args.rc_sensitivity {
    //
    //     println!("Sensitivity Activate");
    //
    //     let pool = ThreadPoolBuilder::new()
    //     .num_threads(args.thread)
    //     .build()
    //     .unwrap();
    //
    //     pool.install(|| {
    //         (0..num_partitions).into_par_iter().for_each(|i| {
    //             let partition_id = format!("{:0width$b}", i, width = args.p[0]);
    //             let filename = format!("{}/partition_{}.fa.zst", out_dir, partition_id);
    //             let _ = update_rc_file(&filename, &args, k);
    //
    //         });
    //     });
    // }
    //
    //
    // println!("Final Compression");
    //
    // let final_filename = format!("{}/final.zst", &args.output);
    // concatenate_partitions(
    //     &args.output,
    //     args.p[0],
    //     &final_filename,
    //     &args.final_compression,
    //     args.final_compression_level,
    // )?;

    // remove_partition_files(&args.output, args.p[0]);

    let total_nucleotides = nucleotide_count.load(Ordering::Relaxed);
    let final_meta = fs::metadata(&final_filename)?;
    let final_size = final_meta.len();
    let total_nucleotides_billions = total_nucleotides as f64 / 1_000_000_000.0;
    let final_size_mb = final_size as f64 / (1024.0 * 1024.0);
    let bits_per_nucleotide = if total_nucleotides > 0 {
        (final_size as f64 * 8.0) / (total_nucleotides as f64)
    } else {
        0.0
    };
    let elapsed = start.elapsed().as_secs_f64();
    let throughput = if elapsed > 0.0 {
        total_nucleotides_billions / elapsed
    } else {
        0.0
    };

    println!("Processing complete.");
    println!("Total nucleotides compressed: {:.3} billion", total_nucleotides_billions);
    println!("Final archive size: {:.3} MB", final_size_mb);
    println!("Effective bits per nucleotide: {:.3}", bits_per_nucleotide);
    println!("Run time: {:.3} seconds", elapsed);
    println!("Throughput: {:.3} million nucleotides/second", throughput);
    println!("Final compressed file: {}", final_filename);

    Ok(())
}
