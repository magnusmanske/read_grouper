use crate::multi_bucket_reader::BucketDataRead;
use anyhow::{anyhow, Result};
use std::cmp::Ordering;
use std::{fs::File, io::BufReader};

pub struct BufReaderEntry<T> {
    buffer: BufReader<File>,
    last_entry_read: T,
}

impl<T: std::cmp::Ord + BucketDataRead + Default> BufReaderEntry<T> {
    pub fn new(filename: &str) -> Result<Self> {
        let file = File::open(filename)?;
        let mut ret = Self {
            buffer: BufReader::new(file),
            last_entry_read: T::default(),
        };

        match ret.last_entry_read.read(&mut ret.buffer) {
            Ok(()) => {}
            Err(_e) => return Err(anyhow!("No entries found in buffer {filename}")),
        };

        Ok(ret)
    }

    #[inline(always)]
    pub fn last_entry_read(&self) -> &T {
        &self.last_entry_read
    }

    // Returns true if the file has reached the end
    #[inline(always)]
    pub fn read_next_entry_failed(&mut self) -> bool {
        match self.last_entry_read.read(&mut self.buffer) {
            Ok(()) => false,
            Err(_) => true,
        }
    }
}

impl<T: std::cmp::Ord + BucketDataRead + Default> Ord for BufReaderEntry<T> {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> Ordering {
        self.last_entry_read.cmp(&other.last_entry_read)
    }
}

impl<T: std::cmp::Ord + BucketDataRead + Default> PartialOrd for BufReaderEntry<T> {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T: std::cmp::Ord + BucketDataRead + Default> Eq for BufReaderEntry<T> {}

impl<T: std::cmp::Ord + BucketDataRead + Default> PartialEq for BufReaderEntry<T> {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.last_entry_read == other.last_entry_read
    }
}
