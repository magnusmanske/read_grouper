use crate::{
    bucket_list::BucketList, data_bucket::DataBucket, kmer::Kmer, kmer_read::KmerRead,
    min_max_reads::MinMaxReads, multi_buf_reader::MultiBufReader, read_pair_kmer::ReadPairKmer,
    KmerReverse, ReadId,
};
use anyhow::{anyhow, Result};
use bam::RecordReader;
use std::{collections::HashMap, path::Path};

const DEFAULT_MIN_BASE_QUALITY: u8 = 20;
const MAX_BUCKET_SIZE: usize = 1_000_000; // kmer-read-pairs

// type KmerBucket = DataBucket<KmerRead>;
// type ReadPairKmerBucket = DataBucket<ReadPairKmer>;

#[derive(Default, Debug)]
pub struct ReadGrouper<KmerBits> {
    bucket_dir: String,
    min_base_quality: u8,
    max_bucket_size: usize,
}

impl<KmerBits: KmerReverse> ReadGrouper<KmerBits> {
    pub fn new(bucket_dir: &str) -> Self {
        Self {
            bucket_dir: bucket_dir.to_string(),
            min_base_quality: DEFAULT_MIN_BASE_QUALITY,
            max_bucket_size: MAX_BUCKET_SIZE,
        }
    }

    pub fn read_bam_file(&self, file_path: &str) -> Result<BucketList> {
        let file_path = Path::new(file_path);
        let sample_name = Self::file_path_to_sample_name(file_path)?;
        let mut reader = bam::BamReader::from_path(file_path, 4)?;
        let mut record = bam::Record::new();
        let mut out_bucket: DataBucket<KmerRead<KmerBits>> = DataBucket::new(
            self.max_bucket_size,
            &self.bucket_dir,
            &sample_name,
            "pairs",
        );
        let mut read_number: ReadId = 0;

        loop {
            match reader.read_into(&mut record) {
                Ok(true) => {}
                Ok(false) => break,
                Err(e) => panic!("{}", e),
            }

            // Generate and process kmers
            let sequence = record.sequence().to_vec();
            let qualities = record.qualities().raw();
            let kmers: Kmer<KmerBits> =
                Kmer::kmers_from_record_incremental(&sequence, qualities, self.min_base_quality);
            for kmer in kmers {
                // out_bucket.add(KmerRead::new(Kmer::new(kmer), read_number));
                // TODO FIXME
            }
            read_number += 1;
        }

        // Write final bucket to disk
        let filenames = out_bucket.finish()?;

        // Create metadata to return
        let bucket_list = BucketList::new(sample_name, filenames, read_number);
        Ok(bucket_list)
    }

    fn process_kmer_grouped_reads(
        &self,
        kmer: &Kmer<KmerBits>,
        reads: &mut Vec<ReadId>,
        min_max: &MinMaxReads,
        bucket: &mut DataBucket<ReadPairKmer<KmerBits>>,
    ) {
        // Reads will be sorted already
        reads.dedup();
        if min_max.is_valid(reads.len()) {
            for read2_pos in 1..reads.len() {
                for read1_pos in 0..read2_pos {
                    let rpk: ReadPairKmer<KmerBits> =
                        ReadPairKmer::new(reads[read1_pos], reads[read2_pos], kmer);
                    bucket.add(rpk);
                }
            }
            // println!("{kmer}: {reads:?}");
        }
        reads.clear();
    }

    pub fn process_read_kmer_buckets(
        &self,
        bucket_list: &BucketList,
        min_max: &MinMaxReads,
    ) -> Result<(BucketList, HashMap<usize, usize>)> {
        let mut mbr: MultiBufReader<KmerRead> = MultiBufReader::new(bucket_list.filenames());

        let sample_name = bucket_list.sample_name().to_string();
        let mut out_bucket: DataBucket<ReadPairKmer<KmerBits>> = DataBucket::new(
            self.max_bucket_size,
            &self.bucket_dir,
            &sample_name,
            "read_pairs",
        );
        let mut stats = HashMap::new();
        let mut last_kmer: Kmer<KmerBits> = Kmer::new(0);
        let mut last_reads_ids = Vec::new();
        while !mbr.is_empty() {
            let kmer_read = match mbr.next() {
                Some(kmer_read) => kmer_read,
                None => break,
            };

            // Flush reads if new kmer
            if last_kmer != *kmer_read.kmer() {
                *stats.entry(last_reads_ids.len()).or_insert(0usize) += 1;
                self.process_kmer_grouped_reads(
                    &last_kmer,
                    &mut last_reads_ids,
                    min_max,
                    &mut out_bucket,
                );
                last_kmer.clone_from(kmer_read.kmer());
            }
            last_reads_ids.push(kmer_read.read_id());
        }

        *stats.entry(last_reads_ids.len()).or_insert(0) += 1;
        self.process_kmer_grouped_reads(&last_kmer, &mut last_reads_ids, min_max, &mut out_bucket);
        stats.remove(&0); // Remove 0-read group

        // Write final bucket to disk
        let filenames = out_bucket.finish()?;

        let bucket_list = BucketList::new(sample_name, filenames, 0);
        Ok((bucket_list, stats))
    }

    fn file_path_to_sample_name(file_path: &Path) -> Result<String> {
        Ok(file_path
            .file_stem()
            .ok_or_else(|| anyhow!("Could not generate sample name from filename [1]"))?
            .to_str()
            .ok_or_else(|| anyhow!("Could not generate sample name from filename [2]"))?
            .to_string())
    }
}
