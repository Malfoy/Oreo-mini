use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use clap::Parser;
use bio::io::{fasta, fastq};
use flate2::read::GzDecoder;
use zstd::stream::Encoder;
use tar::Builder;
use crossbeam_channel;
use std::sync::Arc;
use std::io::BufRead;


/// Command-line arguments.
#[derive(Parser)]
struct Args {
    /// Input FASTA/FASTQ file (can be gzipped if filename ends with .gz)
    #[arg(short, long)]
    input: String,

    /// Output directory (it will be created if needed)
    #[arg(short, long)]
    output: String,

    /// Number of bits (P) used to generate the partition integer.
    #[arg(short, long)]
    p: usize,

    /// k-mer length (k)
    #[arg(short, long)]
    k: usize,

    /// zstd compression level.
    #[arg(short, long, default_value = "3")]
    compression_level: i32,
}

/// Represents a sequence record (either FASTA or FASTQ).
enum Record {
    Fasta { id: String, seq: String },
    Fastq { id: String, seq: String, qual: String },
}

/// Open an input file, decompressing if the filename ends with ".gz".
fn open_input(path: &str) -> Box<dyn Read> {
    let file = File::open(path).expect("Cannot open input file");
    if path.ends_with(".gz") {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    }
}

/// Map a nucleotide to a 2-bit value (A->0, C->1, G->2, T->3). Other characters map to 0.
fn nt_to_val(b: u8) -> u64 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

/// A simple 64-bit mix function (based on splitmix64).
fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}

/// Update the bucket array for a given k-mer hash.
/// The most-significant P bits (after a modulo into buckets) are used.
fn update_bucket(mixed: u64, buckets: &mut [Option<u64>]) {
    // Use modulo p for bucket index.
    let idx = (mixed as usize) % buckets.len();
    match buckets[idx] {
        Some(current) if mixed < current => buckets[idx] = Some(mixed),
        None => buckets[idx] = Some(mixed),
        _ => {}
    }
}

const BASE: u64 = 0x9e3779b97f4a7c15;
fn compute_partition(seq: &str, p: usize, k: usize) -> u64 {
    // p must be a power of two.
    // Let L = log₂(p) be the number of bits needed to index p buckets.
    let L = (p as f64).log2() as usize;
    
    // Initialize an array of p buckets (each holds an Option<u64> for the minimal "rest-of-hash").
    let mut buckets: Vec<Option<u64>> = vec![None; p];
    
    let bytes = seq.as_bytes();
    if bytes.len() < k {
        return 0; // If the sequence is too short, return fingerprint 0.
    }
    
    // Precompute power = BASE^(k-1) mod 2^64 for the rolling update.
    let mut power: u64 = 1;
    for _ in 0..(k - 1) {
        power = power.wrapping_mul(BASE);
    }
    
    // Compute the hash for the first k-mer.
    let mut hash: u64 = 0;
    for i in 0..k {
        let val = nt_to_val(bytes[i]);
        hash = hash.wrapping_mul(BASE).wrapping_add(val);
    }
    {
        // Use the top L bits as bucket index.
        let idx = (hash >> (64 - L)) as usize;
        // The remaining bits:
        let value = hash & ((1u64 << (64 - L)) - 1);
        buckets[idx] = Some(value);
    }
    
    // Process subsequent k-mers using the rolling hash update.
    for i in k..bytes.len() {
        let old_val = nt_to_val(bytes[i - k]);
        let new_val = nt_to_val(bytes[i]);
        // Remove the contribution of the old nucleotide and add the new one.
        hash = hash.wrapping_sub(old_val.wrapping_mul(power));
        hash = hash.wrapping_mul(BASE).wrapping_add(new_val);
        let idx = (hash >> (64 - L)) as usize;
        let value = hash & ((1u64 << (64 - L)) - 1);
        buckets[idx] = match buckets[idx] {
            Some(current) if value < current => Some(value),
            None => Some(value),
            other => other, // keep the current minimum.
        };
    }
    
    // Build the final fingerprint by concatenating the least-significant bit from each bucket.
    let mut fingerprint = 0;
    for bucket in buckets {
        // If no value was stored for a bucket, use 0.
        let bit = bucket.unwrap_or(0) & 1;
        fingerprint = (fingerprint << 1) | bit;
    }
    fingerprint
}

fn main() {
    let args = Args::parse();

    // Create output directory if it does not exist.
    std::fs::create_dir_all(&args.output).expect("Cannot create output directory");

    let num_partitions = 1 << args.p;

    // Create partition channels (one per partition) for sending records to writer threads.
    let mut partition_senders = Vec::with_capacity(num_partitions);
    let mut partition_receivers = Vec::with_capacity(num_partitions);
    for _ in 0..num_partitions {
        let (tx, rx) = crossbeam_channel::bounded::<Record>(100);
        partition_senders.push(tx);
        partition_receivers.push(rx);
    }

    // Open input and determine format.
    let input = open_input(&args.input);
    let mut buf_reader = BufReader::new(input);
    let peek = buf_reader.fill_buf().expect("Error reading input");
    let is_fastq = !peek.is_empty() && peek[0] == b'@';
    let out_ext = if is_fastq { "fq" } else { "fa" };

    // Spawn writer threads—each owns its zstd-compressed output file.
    let mut writer_handles = Vec::with_capacity(num_partitions);
    for i in 0..num_partitions {
        // Take the corresponding receiver.
        let rx = partition_receivers.remove(0);
        let filename = format!("{}/partition_{}.{}", args.output, i, out_ext);
        let file = File::create(&filename)
            .unwrap_or_else(|e| panic!("Cannot create {}: {}", filename, e));
        let encoder = Encoder::new(file, args.compression_level)
            .expect("Cannot create zstd encoder");
        let handle = std::thread::spawn(move || {
            let mut writer = BufWriter::new(encoder);
            while let Ok(record) = rx.recv() {
                match record {
                    Record::Fasta { id, seq } => {
                        writeln!(writer, ">{}", id).expect("Write error");
                        writeln!(writer, "{}", seq).expect("Write error");
                    },
                    Record::Fastq { id, seq, qual } => {
                        writeln!(writer, "@{}", id).expect("Write error");
                        writeln!(writer, "{}", seq).expect("Write error");
                        writeln!(writer, "+").expect("Write error");
                        writeln!(writer, "{}", qual).expect("Write error");
                    },
                }
            }
            writer.flush().expect("Flush error");
            match writer.into_inner() {
                Ok(encoder) => {
                    if let Err(e) = encoder.finish() {
                        eprintln!("Error finishing zstd stream: {}", e);
                    }
                }
                Err(_) => eprintln!("Error retrieving inner zstd encoder"),
            }
        });
        writer_handles.push(handle);
    }

    // Create a channel for sending records from the input reader to worker threads.
    let (record_tx, record_rx) = crossbeam_channel::bounded::<Record>(100);

    // Spawn worker threads (number = available CPUs).
    let num_workers = num_cpus::get();
    let partition_senders = Arc::new(partition_senders);
    let worker_handles: Vec<_> = (0..num_workers)
        .map(|_| {
            let record_rx = record_rx.clone();
            let partition_senders = Arc::clone(&partition_senders);
            let p = args.p;
            let k = args.k;
            std::thread::spawn(move || {
                while let Ok(record) = record_rx.recv() {
                    let seq = match &record {
                        Record::Fasta { seq, .. } => seq,
                        Record::Fastq { seq, .. } => seq,
                    };
                    let partition = compute_partition(seq, p, k);
                    // println!("{}",partition);
                    // Send the record to its designated partition channel.
                    partition_senders[partition as usize]
                        .send(record)
                        .expect("Failed to send record to partition");
                }
            })
        })
        .collect();

    // Read the input file and send records into the worker channel.
    if is_fastq {
        let mut fastq_reader = fastq::Reader::new(buf_reader);
        for result in fastq_reader.records() {
            let rec = result.expect("Error reading FASTQ record");
            let id = rec.id().to_owned();
            let seq = std::str::from_utf8(rec.seq()).unwrap_or("").to_owned();
            let qual = std::str::from_utf8(rec.qual()).unwrap_or("").to_owned();
            let record = Record::Fastq { id, seq, qual };
            record_tx.send(record).expect("Failed to send record");
        }
    } else {
        let mut fasta_reader = fasta::Reader::new(buf_reader);
        for result in fasta_reader.records() {
            let rec = result.expect("Error reading FASTA record");
            let id = rec.id().to_owned();
            let seq = std::str::from_utf8(rec.seq()).unwrap_or("").to_owned();
            let record = Record::Fasta { id, seq };
            record_tx.send(record).expect("Failed to send record");
        }
    }
    // Signal no more records.
    drop(record_tx);

    // Wait for all worker threads to finish.
    for handle in worker_handles {
        handle.join().expect("Worker thread panicked");
    }
    // Dropping the Arc of partition_senders will eventually close all senders.
    drop(partition_senders);

    // Wait for all writer threads to finish.
    for handle in writer_handles {
        handle.join().expect("Writer thread panicked");
    }

    // Finally, create a tar archive containing all partition files.
    let tar_filename = format!("partitions.tar");
    let tar_file = File::create(&tar_filename).expect("Cannot create tar archive");
    let mut tar_builder = Builder::new(tar_file);
    tar_builder
        .append_dir_all("partitions", &args.output)
        .expect("Error creating tar archive");
    tar_builder.finish().expect("Error finishing tar archive");

    println!("Processing complete. Tar archive created at: {}", tar_filename);
}
