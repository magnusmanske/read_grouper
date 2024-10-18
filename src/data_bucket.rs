use anyhow::Result;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::{mem, thread};

pub trait BucketDataWrite {
    fn write(&self, buffer: &mut BufWriter<File>) -> Result<()>;
}

pub trait BucketDataRead {
    fn read(&mut self, file_buffer: &mut BufReader<File>) -> Result<()>;
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

impl<T: std::cmp::Ord + BucketDataWrite + Send + 'static> DataBucket<T> {
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

    #[inline(always)]
    pub fn add(&mut self, pair: T) {
        self.pairs.push(pair);
        self.flush();
    }

    pub fn finish(&mut self) -> Result<Vec<String>> {
        self.start_writing();
        self.write_to_disk()?; // Blocking main thread for final write
        self.wait_for_zero_lock();
        let filenames = self.filenames.lock().unwrap().clone();
        Ok(filenames)
    }

    #[inline(always)]
    fn flush(&mut self) {
        if !self.is_full() {
            return;
        }

        // Duplicate this bucket, with empty data
        let mut bucket_to_write = Self {
            bucket_id: self.bucket_id,
            bucket_size: self.bucket_size,
            bucket_dir: self.bucket_dir.clone(),
            sample_name: self.sample_name.clone(),
            ending: self.ending.clone(),
            pairs: Vec::with_capacity(self.bucket_size + 1),
            currently_writing: self.currently_writing.clone(),
            filenames: self.filenames.clone(),
        };

        // Swap data with the new bucket
        mem::swap(&mut bucket_to_write.pairs, &mut self.pairs);

        // Increase the bucket counter
        self.bucket_id += 1;

        // Start writing the bucket to disk in a new thread
        bucket_to_write.start_writing();
        let _ = thread::spawn(move || bucket_to_write.write_to_disk().unwrap());
    }

    fn wait_for_zero_lock(&self) {
        // Wait for all writing threads to finish
        loop {
            if *self.currently_writing.lock().unwrap() == 0 {
                break;
            }
            thread::sleep(std::time::Duration::from_millis(10));
        }
    }

    #[inline(always)]
    fn start_writing(&self) {
        *self.currently_writing.lock().unwrap() += 1;
    }

    #[inline(always)]
    fn end_writing(&self) {
        *self.currently_writing.lock().unwrap() -= 1;
    }

    #[inline(always)]
    fn is_full(&self) -> bool {
        self.pairs.len() >= self.bucket_size
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

    fn write_to_disk(&mut self) -> Result<()> {
        if self.pairs.is_empty() {
            self.end_writing();
            return Ok(());
        }
        let filename = self.construct_filename();
        if Path::new(&filename).exists() {
            self.set_filename(filename);
            self.end_writing();
            return Ok(());
        }
        self.pairs.sort();
        let file = File::create(&filename)?;
        let mut file = BufWriter::new(file);
        for pair in &self.pairs {
            pair.write(&mut file)?;
        }
        self.set_filename(filename);
        self.end_writing();
        Ok(())
    }
}
