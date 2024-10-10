use crate::{bucket_list::BucketList, kmer::Kmer, kmer_bucket::KmerBucket, kmer_read::KmerRead};
use anyhow::{anyhow, Result};
use bam::RecordReader;
use std::{
    path::Path,
    sync::{Arc, Mutex},
    thread,
};

const DEFAULT_MIN_BASE_QUALITY: u8 = 20;
const MAX_BUCKET_SIZE: usize = 1_000_000; // kmer-read-pairs

#[derive(Default, Debug)]
pub struct ReadGrouper {
    bucket_dir: String,
    min_base_quality: u8,
    max_bucket_size: usize,
    writing: Arc<Mutex<u32>>,
}

impl ReadGrouper {
    pub fn new(bucket_dir: &str) -> Self {
        Self {
            bucket_dir: bucket_dir.to_string(),
            min_base_quality: DEFAULT_MIN_BASE_QUALITY,
            max_bucket_size: MAX_BUCKET_SIZE,
            writing: Arc::new(Mutex::new(0)),
        }
    }

    pub fn read_bam_file(&self, file_path: &str) -> Result<BucketList> {
        let file_path = Path::new(file_path);
        let sample_name = file_path
            .file_stem()
            .ok_or_else(|| anyhow!("Could not generate sample name from filename [1]"))?
            .to_str()
            .ok_or_else(|| anyhow!("Could not generate sample name from filename [2]"))?
            .to_string();
        let mut reader = bam::BamReader::from_path(file_path, 4)?;
        let mut record = bam::Record::new();
        let mut kmer_bucket = KmerBucket::new(
            self.max_bucket_size,
            &self.bucket_dir,
            &sample_name,
            self.writing.clone(),
        );
        let filenames: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
        let mut read_number: u32 = 0;
        loop {
            match reader.read_into(&mut record) {
                Ok(true) => {}
                Ok(false) => break,
                Err(e) => panic!("{}", e),
            }

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
                let kmer = Kmer::new(kmer);
                let kmer_read = KmerRead::new(kmer, read_number);
                kmer_bucket.add(kmer_read);

                if kmer_bucket.is_full() {
                    let next_bucket = kmer_bucket.next_bucket();
                    let mut last_kmer_bucket = kmer_bucket;
                    let filenames_clone = filenames.clone();
                    let writing_clone = self.writing.clone();
                    *writing_clone.lock().unwrap() += 1;
                    let _ = thread::spawn(move || {
                        let result = last_kmer_bucket.write_to_disk();
                        let result = match result {
                            Ok(result) => result,
                            Err(e) => panic!("Error writing final bucket to disk: {}", e),
                        };
                        // TODO handle error
                        if let Some(filename) = result {
                            filenames_clone.lock().unwrap().push(filename);
                        }
                        *writing_clone.lock().unwrap() -= 1;
                    });
                    kmer_bucket = next_bucket;
                }
            }
            read_number += 1;
        }

        // Write final bucket to disk
        let filenames_clone = filenames.clone();
        let writing_clone = self.writing.clone();
        *writing_clone.lock().unwrap() += 1;
        let _ = thread::spawn(move || {
            let result = kmer_bucket.write_to_disk();
            let result = match result {
                Ok(result) => result,
                Err(e) => panic!("Error writing final bucket to disk: {}", e),
            };
            // TODO handle error
            if let Some(filename) = result {
                filenames_clone.lock().unwrap().push(filename);
            }
            *writing_clone.lock().unwrap() -= 1;
        });

        // Wait for all writing threads to finish
        loop {
            if *self.writing.lock().unwrap() == 0 {
                break;
            }
            // println!("Waiting for writing threads to finish");
            thread::sleep(std::time::Duration::from_millis(10));
        }

        println!("Processed {read_number} reads");
        let filenames = Arc::into_inner(filenames).unwrap().into_inner().unwrap(); // TODO error handling
        let bucket_list = BucketList::new(sample_name, filenames, read_number);
        Ok(bucket_list)
    }
}
