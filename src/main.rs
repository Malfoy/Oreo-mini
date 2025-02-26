use clap::Parser;
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
fn nt_to_val(b: u8) -> u64 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

/// A 64-bit multiplier for our rolling hash (golden-ratio based)
const BASE: u64 = 0x9e3779b97f4a7c15;

/// Computes the partition fingerprint for a sequence.
///
/// We maintain an array of length P (where P is the fingerprint bit-length). For each k-mer:
///   1. Compute a 64-bit rolling hash (using BASE to spread the 2-bit nucleotide values).
///   2. Use the top log₂(P) bits as a bucket index (so if P=8, use 3 bits for an index in 0..7).
///   3. In that bucket, store the remaining (64 – log₂(P)) bits (keeping only the minimum seen).
/// Finally, build a fingerprint by concatenating the least-significant bit from each bucket.
/// The result is a P-bit fingerprint (range 0 … 2^P - 1).
fn compute_partition(seq: &str, p: usize, k: usize) -> u64 {
    let buckets_count = p; // use an array of length p
    let index_bits = (p as f64).log2() as usize; // bits used for bucket index
    let mut buckets: Vec<Option<u64>> = vec![None; buckets_count];
    let bytes = seq.as_bytes();
    if bytes.len() < k {
        return 0;
    }
    // Precompute power = BASE^(k-1) for rolling update.
    let mut power: u64 = 1;
    for _ in 0..(k - 1) {
        power = power.wrapping_mul(BASE);
    }
    // Compute hash for first k-mer.
    let mut hash: u64 = 0;
    for i in 0..k {
        let val = nt_to_val(bytes[i]);
        hash = hash.wrapping_mul(BASE).wrapping_add(val);
    }
    {
        let idx = (hash >> (64 - index_bits)) as usize;
        let value = hash & ((1u64 << (64 - index_bits)) - 1);
        buckets[idx] = Some(value);
    }
    // Rolling update for subsequent k-mers.
    for i in k..bytes.len() {
        let old_val = nt_to_val(bytes[i - k]);
        let new_val = nt_to_val(bytes[i]);
        hash = hash.wrapping_sub(old_val.wrapping_mul(power));
        hash = hash.wrapping_mul(BASE).wrapping_add(new_val);
        let idx = (hash >> (64 - index_bits)) as usize;
        let value = hash & ((1u64 << (64 - index_bits)) - 1);
        buckets[idx] = match buckets[idx] {
            Some(current) if value < current => Some(value),
            None => Some(value),
            other => other,
        };
    }
    // Build fingerprint by concatenating the LSB of each bucket.
    let mut fingerprint = 0;
    for bucket in buckets {
        let bit = bucket.unwrap_or(0) & 1;
        fingerprint = (fingerprint << 1) | bit;
    }
    fingerprint
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
    // Ensure data is buffered.
    let _ = buf_reader.fill_buf()?;
    // Detect format.
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

/// Concatenate all partition files (named by binary IDs) in Gray code order into a final compressed file.
/// Each partition file is decompressed on–the–fly and streamed into the final output.
/// The final compression algorithm is chosen based on `final_comp` ("zstd" or "gzip").
fn concatenate_partitions(
    out_dir: &str,
    p: usize,
    final_filename: &str,
    final_comp: &str,
    final_comp_level: i32,
) -> std::io::Result<()> {
    let gray_order = generate_gray_code_order(p);
    
    if final_comp.to_lowercase() == "gzip" {
        // Use gzip final compression.
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
        // Finish gzip compression.
        let gz_encoder = final_writer.into_inner().expect("Error retrieving inner gzip encoder");
        gz_encoder.finish()?;
    } else {
        // Use zstd final compression.
        let final_file = File::create(final_filename)?;
        let mut encoder = Encoder::new(final_file, final_comp_level)
            .expect("Cannot create final zstd encoder");
        // Enable parallel (multithreaded) compression.
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

fn main() -> std::io::Result<()> {
    let args = Args::parse();
    let start = Instant::now(); // Start runtime timer

    // Create output directory if needed.
    fs::create_dir_all(&args.output).expect("Cannot create output directory");
    
    // Determine k-mer length.
    let k = match args.k {
        Some(val) => val,
        None => {
            let detected = auto_detect_k(&args.input)?;
            println!("Auto-detected k = {}", detected);
            detected
        }
    };
    
    // Atomic counter for total nucleotides processed.
    let nucleotide_count = Arc::new(AtomicU64::new(0));
    
    let num_partitions = 1 << args.p; // Number of partition files.
    
    // Create channels for partitioning.
    let mut partition_senders = Vec::with_capacity(num_partitions);
    let mut partition_receivers = Vec::with_capacity(num_partitions);
    for _ in 0..num_partitions {
        let (tx, rx) = crossbeam_channel::bounded::<Record>(100);
        partition_senders.push(tx);
        partition_receivers.push(rx);
    }
    
    // Spawn writer threads—each writes one partition file.
    let out_dir = args.output.clone();
    let mut writer_handles = Vec::with_capacity(num_partitions);
    for i in 0..num_partitions {
        let partition_id = format!("{:0width$b}", i, width = args.p);
        let filename = format!("{}/partition_{}.fa.zst", out_dir, partition_id);
        let rx = partition_receivers.remove(0);
        let compression_level = args.compression_level;
        let handle = std::thread::spawn(move || {
            let file = File::create(&filename)
                .unwrap_or_else(|e| panic!("Cannot create {}: {}", filename, e));
            let encoder = Encoder::new(file, compression_level)
                .expect("Cannot create zstd encoder for partition");
            let mut writer = BufWriter::new(encoder);
            while let Ok(record) = rx.recv() {
                match record {
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
        });
        writer_handles.push(handle);
    }
    
    // Create a channel for sending records to worker threads.
    let (record_tx, record_rx) = crossbeam_channel::bounded::<Record>(100);
    
    // Spawn worker threads (one per available CPU).
    let num_workers = num_cpus::get();
    let partition_senders = Arc::new(partition_senders);
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
                    let fingerprint = compute_partition(seq, p, k);
                    partition_senders[fingerprint as usize]
                        .send(record)
                        .expect("Failed to send record to partition");
                }
            })
        })
        .collect();
    
    // Open input file and detect format (peeking without consuming).
    let input = open_input(&args.input);
    let mut buf_reader = BufReader::new(input);
    let peek = buf_reader.fill_buf().expect("Error peeking input");
    let is_fastq = !peek.is_empty() && peek[0] == b'@';
    
    if is_fastq {
        let reader = fastq::Reader::new(buf_reader);
        for result in reader.records() {
            let rec = result.expect("Error reading FASTQ record");
            let id = rec.id().to_owned();
            let seq = std::str::from_utf8(rec.seq()).unwrap_or("").to_owned();
            let qual = std::str::from_utf8(rec.qual()).unwrap_or("").to_owned();
            let record = Record::Fastq { id, seq, qual };
            record_tx.send(record).expect("Failed to send record");
        }
    } else {
        let reader = fasta::Reader::new(buf_reader);
        for result in reader.records() {
            let rec = result.expect("Error reading FASTA record");
            let id = rec.id().to_owned();
            let seq = std::str::from_utf8(rec.seq()).unwrap_or("").to_owned();
            let record = Record::Fasta { id, seq };
            record_tx.send(record).expect("Failed to send record");
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
    
    // Concatenate partition files (decompress on–the–fly) in Gray code order into final compressed output.
    let final_filename = format!("{}/final.zst", &args.output);
    concatenate_partitions(&args.output, args.p, &final_filename, &args.final_compression, args.final_compression_level)?;
    
    // Remove partition files.
    remove_partition_files(&args.output, args.p);
    
    // Retrieve statistics.
    let total_nucleotides = nucleotide_count.load(Ordering::Relaxed);
    let final_meta = fs::metadata(&final_filename)?;
    let final_size = final_meta.len(); // bytes
    let total_nucleotides_millions = total_nucleotides as f64 / 1_000_000.0;
    let final_size_mb = final_size as f64 / (1024.0 * 1024.0);
    let bits_per_nucleotide = if total_nucleotides > 0 {
        (final_size as f64 * 8.0) / (total_nucleotides as f64)
    } else {
        0.0
    };
    
    let elapsed = start.elapsed();
    let seconds = elapsed.as_secs_f64();
    let throughput = if seconds > 0.0 {
        total_nucleotides_millions / seconds
    } else {
        0.0
    };
    
    println!("Processing complete.");
    println!("Total nucleotides compressed: {:.3} million", total_nucleotides_millions);
    println!("Final archive size: {:.3} MB", final_size_mb);
    println!("Effective bits per nucleotide: {:.3}", bits_per_nucleotide);
    println!("Run time: {:.3} seconds", seconds);
    println!("Throughput: {:.3} million nucleotides/second", throughput);
    println!("Final compressed file: {}", final_filename);
    
    Ok(())
}
