mod bucket_list;
mod kmer;
mod kmer_bucket;
mod kmer_read;
mod read_grouper;

use read_grouper::ReadGrouper;

fn main() {
    let rg = ReadGrouper::new("/Users/mm6/rust/read_grouper/buckets");
    let bucket_list = rg
        .read_bam_file("/Users/mm6/rust/read_grouper/SRR9217386.sorted.bam")
        .unwrap();
    println!("Sample name: {}", bucket_list.sample_name());
    println!("Number of reads: {}", bucket_list.number_of_reads());
    println!("Files: {}", bucket_list.filenames().len());
}
