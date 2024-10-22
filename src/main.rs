mod bucket_list;
mod buf_reader_entry;
mod data_bucket;
mod kmer;
mod kmer_bits;
mod kmer_read;
mod min_max_reads;
mod multi_bucket_reader;
mod read_grouper;
mod read_pair_kmer;

pub use kmer_bits::*;
pub use min_max_reads::MinMaxReads;
pub use read_grouper::ReadGrouper;

pub type ReadId = u32;

fn main() {
    let bucket_dir = "/Users/mm6/rust/read_grouper/buckets";
    let bam_file = "/Users/mm6/rust/read_grouper/SRR9217386.sorted.bam";

    let rg: ReadGrouper<Kmer16> = ReadGrouper::new(bucket_dir);
    let bucket_list = rg.read_bam_file(bam_file).unwrap();
    println!("Sample name: {}", bucket_list.sample_name());
    println!("Number of reads: {}", bucket_list.number_of_reads());
    println!("Files: {}", bucket_list.filenames().len());

    let (_bucket_list, _stats) = rg
        .process_read_kmer_buckets(&bucket_list, &MinMaxReads::new(3, 50))
        .unwrap();

    // let mut keys = stats.keys().cloned().collect::<Vec<_>>();
    // keys.sort();
    // println!("reads_per_kmer\toccurrences");
    // for key in keys {
    //     println!("{}\t{}", key, stats[&key]);
    // }
}

/*
\rm ./target/aarch64-apple-darwin/release/read_grouper ; \
RUSTFLAGS="-C target-cpu=native" cargo build --release  --target aarch64-apple-darwin ; \
time ./target/aarch64-apple-darwin/release/read_grouper
*/
