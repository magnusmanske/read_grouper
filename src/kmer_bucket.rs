use crate::kmer_read::KmerRead;
use anyhow::Result;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;
use std::sync::{Arc, Mutex};

#[derive(Debug)]
pub struct KmerBucket {
    kmer_reads: Vec<KmerRead>,
    bucket_id: u32,
    bucket_size: usize,
    bucket_dir: String,
    sample_name: String,
    writing: Arc<Mutex<u32>>,
}

impl KmerBucket {
    pub fn new(
        bucket_size: usize,
        bucket_dir: &str,
        sample_name: &str,
        writing: Arc<Mutex<u32>>,
    ) -> Self {
        Self {
            bucket_id: 0,
            bucket_size,
            bucket_dir: bucket_dir.to_string(),
            sample_name: sample_name.to_string(),
            kmer_reads: Vec::new(),
            writing,
        }
    }

    pub fn next_bucket(&self) -> Self {
        Self {
            bucket_id: self.bucket_id + 1,
            bucket_size: self.bucket_size,
            bucket_dir: self.bucket_dir.clone(),
            sample_name: self.sample_name.clone(),
            writing: self.writing.clone(),
            kmer_reads: Vec::new(),
        }
    }

    #[inline(always)]
    pub fn add(&mut self, kmer_read: KmerRead) {
        self.kmer_reads.push(kmer_read);
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.kmer_reads.len()
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.kmer_reads.is_empty()
    }

    #[inline(always)]
    pub fn is_full(&self) -> bool {
        self.len() >= self.bucket_size
    }

    pub fn write_to_disk(&mut self) -> Result<Option<String>> {
        if self.is_empty() {
            return Ok(None);
        }
        let filename = format!(
            "{}/{}_{}_{}.kmer_reads",
            self.bucket_dir, self.sample_name, self.bucket_size, self.bucket_id
        );
        if Path::new(&filename).exists() {
            return Ok(Some(filename));
        }
        self.kmer_reads.sort();
        let file = File::create(&filename)?;
        let mut file = BufWriter::new(file);
        for kmer_read in &self.kmer_reads {
            file.write_all(&kmer_read.kmer().to_le_bytes())?;
            file.write_all(&kmer_read.read_id().to_le_bytes())?;
        }
        Ok(Some(filename))
    }
}
