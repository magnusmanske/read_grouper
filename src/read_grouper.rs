use crate::{
    bucket_list::BucketList, buf_reader_entry::BufReaderEntry, data_bucket::DataBucket, kmer::Kmer,
    kmer_read::KmerRead, min_max_reads::MinMaxReads, read_pair_kmer::ReadPairKmer, ReadId,
};
use anyhow::{anyhow, Result};
use bam::RecordReader;
use std::{collections::HashMap, path::Path, sync::Arc, thread};

const DEFAULT_MIN_BASE_QUALITY: u8 = 20;
const MAX_BUCKET_SIZE: usize = 1_000_000; // kmer-read-pairs

type KmerBucket = DataBucket<KmerRead>;
type ReadPairKmerBucket = DataBucket<ReadPairKmer>;

#[derive(Default, Debug)]
pub struct ReadGrouper {
    bucket_dir: String,
    min_base_quality: u8,
    max_bucket_size: usize,
}

impl ReadGrouper {
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
        let mut out_bucket = KmerBucket::new(
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
            let kmers = if true {
                // Faster
                Kmer::kmers_from_record_incremental(&sequence, qualities, self.min_base_quality)
            } else {
                // More kmers after bad bases/qualitiy scores
                Kmer::kmers_from_record_de_novo(&sequence, qualities, self.min_base_quality)
            };
            for kmer in kmers {
                out_bucket.add(KmerRead::new(Kmer::new(kmer), read_number));
            }

            // Flush bucket to disk and start a new one, if necessary
            if out_bucket.is_full() {
                let next_bucket = out_bucket.next_bucket();
                out_bucket.start_writing();
                let _ = thread::spawn(move || out_bucket.write_to_disk().unwrap());
                out_bucket = next_bucket;
            }

            read_number += 1;
        }

        // Write final bucket to disk
        let currently_writing = out_bucket.currently_writing().clone();
        let filenames = out_bucket.filenames().clone();
        out_bucket.start_writing();
        let _ = thread::spawn(move || out_bucket.write_to_disk().unwrap());

        Self::wait_for_zero_lock(currently_writing);

        // Create metadata to return
        let filenames = filenames.lock().unwrap().clone();
        let bucket_list = BucketList::new(sample_name, filenames, read_number);
        Ok(bucket_list)
    }

    fn process_kmer_grouped_reads(
        &self,
        kmer: &Kmer,
        reads: &mut Vec<ReadId>,
        min_max: &MinMaxReads,
        bucket: &mut DataBucket<ReadPairKmer>,
    ) {
        // Reads will be sorted already
        reads.dedup();
        if min_max.is_valid(reads.len()) {
            for read2_pos in 1..reads.len() {
                for read1_pos in 0..read2_pos {
                    let rpk = ReadPairKmer::new(reads[read1_pos], reads[read2_pos], kmer);
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
        type BufReaderKmerRead = BufReaderEntry<KmerRead>;
        let mut readers = Vec::new();
        for filename in bucket_list.filenames() {
            let reader = BufReaderKmerRead::new(filename)?;
            readers.push(reader);
        }

        let sample_name = bucket_list.sample_name().to_string();
        let mut out_bucket = ReadPairKmerBucket::new(
            self.max_bucket_size,
            &self.bucket_dir,
            &sample_name,
            "read_pairs",
        );
        let mut stats = HashMap::new();
        let mut last_kmer = Kmer::new(0);
        let mut last_reads_ids = Vec::new();
        while !readers.is_empty() {
            // Find the next kmer/read pair to process
            let min_key = readers
                .iter()
                .enumerate()
                .min_by_key(|&(_, value)| value)
                .map(|(key, _)| key);
            let min_key = match min_key {
                Some(min_key) => min_key,
                None => break,
            };
            let kmer_read = readers[min_key].last_entry_read();

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

                if out_bucket.is_full() {
                    let next_bucket = out_bucket.next_bucket();
                    out_bucket.start_writing();
                    let _ = thread::spawn(move || out_bucket.write_to_disk().unwrap());
                    out_bucket = next_bucket;
                }
            }
            last_reads_ids.push(kmer_read.read_id());

            // Remove reader if it has no more kmers
            if readers[min_key].read_next_kmer_failed() {
                readers.remove(min_key);
            }
        }

        *stats.entry(last_reads_ids.len()).or_insert(0) += 1;
        self.process_kmer_grouped_reads(&last_kmer, &mut last_reads_ids, min_max, &mut out_bucket);
        stats.remove(&0); // Remove 0-read group

        // Write final bucket to disk
        let currently_writing = out_bucket.currently_writing().clone();
        let filenames = out_bucket.filenames().clone();
        out_bucket.start_writing();
        let _ = thread::spawn(move || out_bucket.write_to_disk().unwrap());

        Self::wait_for_zero_lock(currently_writing);

        let filenames = filenames.lock().unwrap().clone();
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

    fn wait_for_zero_lock(counter_mutex: Arc<std::sync::Mutex<usize>>) {
        // Wait for all writing threads to finish
        loop {
            if *counter_mutex.lock().unwrap() == 0 {
                break;
            }
            thread::sleep(std::time::Duration::from_millis(10));
        }
    }
}
