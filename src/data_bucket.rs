use anyhow::Result;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::sync::{Arc, Mutex};

pub trait BucketDataWrite {
    fn write(&self, buffer: &mut BufWriter<File>) -> Result<()>;
}

#[derive(Debug)]
pub struct DataBucket<T> {
    pairs: Vec<T>,
    bucket_id: usize,
    bucket_size: usize,
    bucket_dir: String,
    sample_name: String,
    ending: String,
    currently_writing: Arc<Mutex<usize>>,
    filenames: Arc<Mutex<Vec<String>>>,
}

impl<T: std::cmp::Ord + BucketDataWrite> DataBucket<T> {
    pub fn new(bucket_size: usize, bucket_dir: &str, sample_name: &str, ending: &str) -> Self {
        Self {
            bucket_id: 0,
            bucket_size,
            bucket_dir: bucket_dir.to_string(),
            sample_name: sample_name.to_string(),
            ending: ending.to_string(),
            pairs: Vec::with_capacity(bucket_size + 1),
            currently_writing: Arc::new(Mutex::new(0)),
            filenames: Arc::new(Mutex::new(Vec::new())),
        }
    }

    pub fn next_bucket(&self) -> Self {
        Self {
            bucket_id: self.bucket_id + 1,
            bucket_size: self.bucket_size,
            bucket_dir: self.bucket_dir.clone(),
            sample_name: self.sample_name.clone(),
            ending: self.ending.clone(),
            pairs: Vec::with_capacity(self.bucket_size + 1),
            currently_writing: self.currently_writing.clone(),
            filenames: self.filenames.clone(),
        }
    }

    #[inline(always)]
    pub fn filenames(&self) -> &Arc<Mutex<Vec<String>>> {
        &self.filenames
    }

    #[inline(always)]
    pub fn currently_writing(&self) -> &Arc<Mutex<usize>> {
        &self.currently_writing
    }

    #[inline(always)]
    pub fn start_writing(&self) {
        *self.currently_writing.lock().unwrap() += 1;
    }

    #[inline(always)]
    pub fn add(&mut self, pair: T) {
        self.pairs.push(pair);
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    #[inline(always)]
    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    #[inline(always)]
    pub fn is_full(&self) -> bool {
        self.len() >= self.bucket_size
    }

    fn set_filename(&self, filename: String) {
        // Add filename to list at bucket_id position
        let mut filenames = self.filenames.lock().unwrap();
        while filenames.len() <= self.bucket_id {
            filenames.push(String::new());
        }
        *filenames.get_mut(self.bucket_id).unwrap() = filename;
    }

    fn construct_filename(&self) -> String {
        format!(
            "{}/{}_{}_{}.{}",
            self.bucket_dir, self.sample_name, self.bucket_size, self.bucket_id, self.ending
        )
    }

    pub fn write_to_disk(&mut self) -> Result<()> {
        if self.is_empty() {
            *self.currently_writing.lock().unwrap() -= 1;
            return Ok(());
        }
        let filename = self.construct_filename();
        if Path::new(&filename).exists() {
            self.set_filename(filename);
            *self.currently_writing.lock().unwrap() -= 1;
            return Ok(());
        }
        self.pairs.sort();
        let file = File::create(&filename)?;
        let mut file = BufWriter::new(file);
        for pair in &self.pairs {
            pair.write(&mut file)?;
        }
        self.set_filename(filename);
        *self.currently_writing.lock().unwrap() -= 1;
        Ok(())
    }
}
