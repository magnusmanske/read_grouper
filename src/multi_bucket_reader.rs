use crate::buf_reader_entry::BufReaderEntry;
use anyhow::Result;
use std::{fs::File, io::BufReader};

pub trait BucketDataRead: std::cmp::Ord + Default + Clone {
    fn read(&mut self, file_buffer: &mut BufReader<File>) -> Result<()>;
}

pub struct MultiBucketReader<T> {
    readers: Vec<BufReaderEntry<T>>,
    files: Vec<String>,
    clean_after_read: bool,
}

impl<T: BucketDataRead> MultiBucketReader<T> {
    pub fn new(files: &[String]) -> Self {
        Self {
            readers: files
                .iter()
                .map(|f| BufReaderEntry::new(f).unwrap())
                .collect(),
            files: files.to_vec(),
            clean_after_read: false,
        }
    }

    /// If set to true, every file will be removed after the last entry is read.
    pub fn clean_after_read(mut self) -> Self {
        self.clean_after_read = true;
        self
    }
}

impl<T: BucketDataRead> Iterator for MultiBucketReader<T> {
    type Item = T;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        // Find the next entry to process
        let min_key = self
            .readers
            .iter()
            .enumerate()
            .min_by_key(|&(_, value)| value)
            .map(|(key, _)| key)?;
        let ret = (*self.readers[min_key].last_entry_read()).clone();

        // Load the next entry from the reader.
        // Remove reader if it has no more entries.
        if self.readers[min_key].read_next_entry_failed() {
            self.readers.remove(min_key);
            if self.clean_after_read {
                std::fs::remove_file(&self.files[min_key]).expect("Could not remove file");
            }
            self.files.remove(min_key);
        }

        Some(ret)
    }
}
