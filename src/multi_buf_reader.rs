use crate::{buf_reader_entry::BufReaderEntry, data_bucket::BucketDataRead};

pub struct MultiBufReader<T> {
    readers: Vec<BufReaderEntry<T>>,
}

impl<T: std::cmp::Ord + BucketDataRead + Default + Clone> MultiBufReader<T> {
    pub fn new(files: &[String]) -> Self {
        Self {
            readers: files
                .iter()
                .map(|f| BufReaderEntry::new(f).unwrap())
                .collect(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.readers.is_empty()
    }

    // TODO as iterator?
    pub fn next(&mut self) -> Option<T> {
        // Find the next entry to process
        let min_key = self
            .readers
            .iter()
            .enumerate()
            .min_by_key(|&(_, value)| value)
            .map(|(key, _)| key)?;
        let ret = (*self.readers[min_key].last_entry_read()).clone();

        // Remove reader if it has no more entries
        if self.readers[min_key].read_next_entry_failed() {
            self.readers.remove(min_key);
        }

        Some(ret)
    }
}
