mod bucket_list;
mod buf_reader_entry;
mod data_bucket;
mod kmer;
mod kmer_backlog;
mod kmer_bits;
mod kmer_read;
mod kmer_with_position;
mod min_max_reads;
mod multi_bucket_reader;
mod read_grouper;
mod read_pair_kmer;

use std::{fs::File, io, process::exit};

use kmer::Kmer;
pub use kmer_bits::*;
pub use min_max_reads::MinMaxReads;
use multi_bucket_reader::BucketDataRead;
pub use read_grouper::ReadGrouper;

pub type ReadId = u32;

fn cwd() -> String {
    std::env::current_dir()
        .unwrap()
        .to_str()
        .unwrap()
        .to_string()
}

fn compare_kmers(mut args: std::env::Args) {
    let note = "Missing fasta file argument\nUsage: compare_kmers <kmer_file1> <kmer_file1>";
    let kmer_file_1 = args.next().expect(note);
    let kmer_file_2 = args.next().expect(note);

    let file1 = File::open(kmer_file_1).expect("Could not open file {kmer_file_1}");
    let file2 = File::open(kmer_file_2).expect("Could not open file {kmer_file_2}");
    let mut reader1 = io::BufReader::new(file1);
    let mut reader2 = io::BufReader::new(file2);
    // let path1 = Path::new(&kmer_file_1);
    // let path2 = Path::new(&kmer_file_2);
    // let mut reader1 = ReadGrouper::get_uncompressed_or_gz_reader(path1);
    // let mut reader2 = ReadGrouper::get_uncompressed_or_gz_reader(path2);
    let mut kmer1: Kmer<Kmer32> = Kmer::default();
    let mut kmer2: Kmer<Kmer32> = Kmer::default();
    let mut kmers_read_1 = 1;
    let mut kmers_read_2 = 1;
    let mut kmers_equal = 0;
    kmer1.read(&mut reader1).unwrap();
    kmer2.read(&mut reader2).unwrap();
    loop {
        #[allow(clippy::comparison_chain)]
        if kmer1 < kmer2 {
            match kmer1.read(&mut reader1) {
                Ok(_) => kmers_read_1 += 1,
                Err(_) => break,
            }
        } else if kmer1 > kmer2 {
            match kmer2.read(&mut reader2) {
                Ok(_) => kmers_read_2 += 1,
                Err(_) => break,
            }
        } else {
            // Kmers are equal
            match kmer1.read(&mut reader1) {
                Ok(_) => kmers_read_1 += 1,
                Err(_) => break,
            }
            match kmer2.read(&mut reader2) {
                Ok(_) => kmers_read_2 += 1,
                Err(_) => break,
            }
            kmers_equal += 1;
        }
    }
    println!("Found {kmers_equal} kmers in both files");
    println!(
        "Kmers read from file 1: {kmers_read_1}, {}% found",
        kmers_equal as f32 / kmers_read_1 as f32 * 100.0
    );
    println!(
        "Kmers read from file 2: {kmers_read_2}, {}% found",
        kmers_equal as f32 / kmers_read_2 as f32 * 100.0
    );
}

fn fasta2kmers(mut args: std::env::Args) {
    const DESIRED_NUMBER_OF_BUCKETS: u64 = 20;

    let fasta_file = args
        .next()
        .expect("Missing fasta file argument\nUsage: fasta2kmers <fasta_file> [-d <bucket_dir>]");

    let mut bucket_dir = cwd();
    let max_kmer_occurrence = Some(5);

    while let Some(param) = args.next() {
        match param.as_str() {
            "-d" => {
                bucket_dir = args.next().expect("Missing bucket dir argument");
            }
            _ => {
                eprintln!("Unknown parameter: {}", param);
                exit(1);
            }
        }
    }

    let _ = std::fs::create_dir_all(&bucket_dir); // Ignore result
    let mut rg: ReadGrouper<Kmer32> = ReadGrouper::new(&bucket_dir);

    let max_bucket_size = std::fs::metadata(&fasta_file).unwrap().len() / DESIRED_NUMBER_OF_BUCKETS;
    rg.set_max_bucket_size(max_bucket_size as usize);

    let bucket_list1 = rg.read_fasta_file(&fasta_file).unwrap();
    let bucket_list2 = rg.unique_kmers(&bucket_list1, max_kmer_occurrence).unwrap();
    println!("Found {} unique kmers", bucket_list2.number_of_reads());
    println!("Kmer index file is {}", bucket_list2.filenames()[0]);
}

fn fastq2kmers(mut args: std::env::Args) {
    const DESIRED_NUMBER_OF_BUCKETS: u64 = 20;

    let fastq_file = args
        .next()
        .expect("Missing fasta file argument\nUsage: fastq2kmers <fastq_file> [-d <bucket_dir>]");

    let mut bucket_dir = cwd();
    let max_kmer_occurrence = Some(5);

    while let Some(param) = args.next() {
        match param.as_str() {
            "-d" => {
                bucket_dir = args.next().expect("Missing bucket dir argument");
            }
            _ => {
                eprintln!("Unknown parameter: {}", param);
                exit(1);
            }
        }
    }

    let _ = std::fs::create_dir_all(&bucket_dir); // Ignore result
    let mut rg: ReadGrouper<Kmer32> = ReadGrouper::new(&bucket_dir);

    let max_bucket_size = std::fs::metadata(&fastq_file).unwrap().len() / DESIRED_NUMBER_OF_BUCKETS;
    rg.set_max_bucket_size(max_bucket_size as usize);

    let bucket_list1 = rg.read_fastq_file(&fastq_file).unwrap();
    let bucket_list2 = rg.unique_kmers(&bucket_list1, max_kmer_occurrence).unwrap();
    println!("Found {} unique kmers", bucket_list2.number_of_reads());
    println!("Kmer index file is {}", bucket_list2.filenames()[0]);
}

fn bam2kmers(mut args: std::env::Args) {
    const DESIRED_NUMBER_OF_BUCKETS: u64 = 100;

    let fasta_file = args
        .next()
        .expect("Missing fasta file argument\nUsage: bam2kmers <bam_file> [-d <bucket_dir>]");

    let mut bucket_dir = cwd();
    let max_kmer_occurrence = Some(5);

    while let Some(param) = args.next() {
        match param.as_str() {
            "-d" => {
                bucket_dir = args.next().expect("Missing bucket dir argument");
            }
            _ => {
                eprintln!("Unknown parameter: {}", param);
                exit(1);
            }
        }
    }

    let _ = std::fs::create_dir_all(&bucket_dir); // Ignore result
    let mut rg: ReadGrouper<Kmer32> = ReadGrouper::new(&bucket_dir);

    let max_bucket_size = std::fs::metadata(&fasta_file).unwrap().len() / DESIRED_NUMBER_OF_BUCKETS;
    rg.set_max_bucket_size(max_bucket_size as usize);

    let bucket_list1 = rg.read_bam_file_kmers(&fasta_file).unwrap();
    let bucket_list2 = rg.unique_kmers(&bucket_list1, max_kmer_occurrence).unwrap();
    println!("Found {} unique kmers", bucket_list2.number_of_reads());
    println!("Kmer index file is {}", bucket_list2.filenames()[0]);
}

fn main() {
    let mut args = std::env::args();

    let command = args.next().unwrap(); // Binary name
    let subcommand = match args.next() {
        Some(subcommand) => subcommand,
        None => {
            eprintln!("Usage: {command} <subcommand> [options]");
            exit(1);
        }
    };

    match subcommand.as_str() {
        "fasta2kmers" => fasta2kmers(args),
        "fastq2kmers" => fastq2kmers(args),
        "bam2kmers" => bam2kmers(args),
        "compare_kmers" => compare_kmers(args),
        _ => {
            eprintln!("Unknown subcommand: {}", subcommand);
            exit(1);
        }
    }
}

// fn _main_bam() {
//     let bam_file = "/Users/mm6/rust/read_grouper/SRR9217386.sorted.bam";
//     let bucket_dir = "/Users/mm6/rust/read_grouper/buckets";

//     let rg: ReadGrouper<Kmer16> = ReadGrouper::new(bucket_dir);
//     let bucket_list = rg.read_bam_file(bam_file).unwrap();
//     println!("Sample name: {}", bucket_list.sample_name());
//     println!("Number of reads: {}", bucket_list.number_of_reads());
//     println!("Files: {}", bucket_list.filenames().len());

//     let (bucket_list, _stats) = rg
//         .process_read_kmer_buckets(&bucket_list, &MinMaxReads::new(3, 50))
//         .unwrap();

//     rg.align_reads(&bucket_list, bam_file).unwrap();
// }

/*
\rm ./target/aarch64-apple-darwin/release/read_grouper ; \
RUSTFLAGS="-C target-cpu=native" cargo build --release  --target aarch64-apple-darwin ; \
time ./target/aarch64-apple-darwin/release/read_grouper
*/
