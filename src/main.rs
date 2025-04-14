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
use std::sync::{Arc};
use std::sync::atomic::{AtomicU64, Ordering};
use std::fs::{File, self};
use std::io::{BufReader, BufWriter, BufRead, Read, Write};
use std::time::Instant;

use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

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
    #[arg(short, long, default_value = "8")]
    p: usize,

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

    /// Enable second-level partitioning + reordering (using a different hash function, same p/k).
    #[arg(long)]
    second_level: bool,

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

/// Convert a nucleotide to a 2-bit value (A:0, C:1, G:2, T:3)
// fn nt_to_val(b: u8) -> u64 {
//     match b {
//         b'A' | b'a' => 0,
//         b'C' | b'c' => 1,
//         b'G' | b'g' => 2,
//         b'T' | b't' => 3,
//         _ => 0,
//     }
// }

/// 64-bit multipliers for rolling hash (different for first- and second-level).
const BASE1: u64 = 0x9e3779b97f4a7c15;
const BASE2: u64 = 0xc2b2ae3d27d4eb4f;

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
fn update_rc_file(input: &str, k: usize, rc_loops: i32) -> std::io::Result<()> {
    let input_reader = open_input(input);
    let mut buf_reader = BufReader::new(input_reader);

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
    let encoder = Encoder::new(file, 4).expect("Cannot create zstd encoder for partition");
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

    let nucleotide_count = Arc::new(AtomicU64::new(0));
    let num_partitions = 1 << args.p;

    let mut partition_senders = Vec::with_capacity(num_partitions);
    let mut partition_receivers = Vec::with_capacity(num_partitions);
    for _ in 0..num_partitions {
        let (tx, rx) = crossbeam_channel::bounded::<Record>(100);
        partition_senders.push(tx);
        partition_receivers.push(rx);
    }

    let out_dir = args.output.clone();
    let mut writer_handles = Vec::with_capacity(num_partitions);
    for i in 0..num_partitions {
        let partition_id = format!("{:0width$b}", i, width = args.p);
        let filename = format!("{}/partition_{}.fa.zst", out_dir, partition_id);
        let rx = partition_receivers.remove(0);
        let compression_level = args.compression_level;
        let p = args.p;
        let k = k;
        let second_level = args.second_level;
        let handle = std::thread::spawn(move || {
            let mut partition_records = Vec::new();
            while let Ok(record) = rx.recv() {
                partition_records.push(record);
            }
            if second_level {
                let mut subpartitions: Vec<Vec<Record>> = (0..(1 << p)).map(|_| Vec::new()).collect();
                for record in partition_records {
                    let (seq, _) = match &record {
                        Record::Fasta { seq, .. } => (seq, false),
                        Record::Fastq { seq, .. } => (seq, true),
                    };
                    let fp2;
                    if args.rc_sensitivity {
                        fp2 = compute_partition_canonique(seq, p, k, BASE2) as usize;
                    } else {
                        fp2 = compute_partition(seq, p, k, BASE2) as usize;
                    }
                    subpartitions[fp2].push(record);
                }
                let file = File::create(&filename)
                    .unwrap_or_else(|e| panic!("Cannot create {}: {}", filename, e));
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
                writer.flush().expect("Flush error");
                match writer.into_inner() {
                    Ok(encoder) => {
                        if let Err(e) = encoder.finish() {
                            eprintln!("Error finishing partition zstd stream: {}", e);
                        }
                    }
                    Err(_) => eprintln!("Error retrieving inner encoder for partition"),
                }
            } else {
                let file = File::create(&filename)
                    .unwrap_or_else(|e| panic!("Cannot create {}: {}", filename, e));
                let encoder = Encoder::new(file, compression_level)
                    .expect("Cannot create zstd encoder for partition");
                let mut writer = BufWriter::new(encoder);
                for rec in partition_records {
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
                writer.flush().expect("Flush error");
                match writer.into_inner() {
                    Ok(encoder) => {
                        if let Err(e) = encoder.finish() {
                            eprintln!("Error finishing partition zstd stream: {}", e);
                        }
                    }
                    Err(_) => eprintln!("Error retrieving inner encoder for partition"),
                }
            }
        });
        writer_handles.push(handle);
    }

    let (record_tx, record_rx) = crossbeam_channel::bounded::<Record>(100);
    let partition_senders = Arc::new(partition_senders);
    let num_workers = num_cpus::get();
    let worker_handles: Vec<_> = (0..num_workers)
        .map(|_| {
            let record_rx = record_rx.clone();
            let partition_senders = Arc::clone(&partition_senders);
            let p = args.p;
            let k = k;
            let nucleotide_count = Arc::clone(&nucleotide_count);
            std::thread::spawn(move || {
                while let Ok(record) = record_rx.recv() {
                    let seq = match &record {
                        Record::Fasta { seq, .. } => seq,
                        Record::Fastq { seq, .. } => seq,
                    };
                    nucleotide_count.fetch_add(seq.len() as u64, Ordering::Relaxed);
                    // let fp = compute_partition(seq, p, k, BASE1);
                    let fp;
                    if args.rc_sensitivity {
                        fp = compute_partition_canonique(seq, p, k, BASE1);
                    } else {
                        fp = compute_partition(seq, p, k, BASE1);
                    }
                    partition_senders[fp as usize]
                        .send(record)
                        .expect("Failed to send record to partition");
                }
            })
        })
        .collect();

    let input = open_input(&args.input);
    let mut buf_reader = BufReader::new(input);
    let peek = buf_reader.fill_buf().expect("Error peeking input");
    let is_fastq = !peek.is_empty() && peek[0] == b'@';

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
    drop(partition_senders);

    for handle in writer_handles {
        handle.join().expect("Writer thread panicked");
    }



    if args.rc_sensitivity {

        println!("Sensitivity Activate");

        let pool = ThreadPoolBuilder::new()
        .num_threads(args.thread)
        .build()
        .unwrap();

        pool.install(|| {
            (0..num_partitions).into_par_iter().for_each(|i| {
                let partition_id = format!("{:0width$b}", i, width = args.p);
                let filename = format!("{}/partition_{}.fa.zst", out_dir, partition_id);
                let _ = update_rc_file(&filename, k,args.rc_compression_loop);

            });
        });
    }


    let final_filename = format!("{}/final.zst", &args.output);
    concatenate_partitions(
        &args.output,
        args.p,
        &final_filename,
        &args.final_compression,
        args.final_compression_level,
    )?;

    remove_partition_files(&args.output, args.p);

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
