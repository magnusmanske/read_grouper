mod bucket_list;
mod buf_reader_kmer_read;
mod kmer;
mod kmer_bucket;
mod kmer_read;
mod min_max_reads;
mod read_grouper;

use min_max_reads::MinMaxReads;
use read_grouper::ReadGrouper;

pub type KmerBits = u32;
pub type ReadId = u32;

fn main() {
    let rg = ReadGrouper::new("/Users/mm6/rust/read_grouper/buckets");
    let bucket_list = rg
        .read_bam_file("/Users/mm6/rust/read_grouper/SRR9217386.sorted.bam")
        .unwrap();
    println!("Sample name: {}", bucket_list.sample_name());
    println!("Number of reads: {}", bucket_list.number_of_reads());
    println!("Files: {}", bucket_list.filenames().len());

    rg.process_buckets(&bucket_list, &MinMaxReads::new(3, 50))
        .unwrap();
}
