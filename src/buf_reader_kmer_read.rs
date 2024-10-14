use crate::KmerBits;
use crate::{kmer::Kmer, kmer_read::KmerRead};
use anyhow::{anyhow, Result};
use std::cmp::Ordering;
use std::{
    fs::File,
    io::{BufReader, Read},
};

// #[derive(Default, Debug)]
pub struct BufReaderKmerRead {
    buffer: BufReader<File>,
    last_kmer_read: KmerRead,
}

impl BufReaderKmerRead {
    pub fn new(filename: &str) -> Result<Self> {
        let file = File::open(filename)?;
        let mut buffer = BufReader::new(file);
        let last_kmer_read = match Self::get_kmer_read_from_buffer(&mut buffer) {
            Some(kmer_read) => kmer_read,
            None => return Err(anyhow!("No kmer-reads found in buffer {filename}")),
        };

        Ok(Self {
            buffer,
            last_kmer_read,
        })
    }

    pub fn last_kmer_read(&self) -> &KmerRead {
        &self.last_kmer_read
    }

    // Returns true if the file has reached the end
    pub fn read_next_kmer_failed(&mut self) -> bool {
        match Self::get_kmer_read_from_buffer(&mut self.buffer) {
            Some(kmer_read) => {
                self.last_kmer_read = kmer_read;
                false
            }
            None => true,
        }
    }

    fn get_kmer_read_from_buffer(file_buffer: &mut BufReader<File>) -> Option<KmerRead> {
        let mut buffer = [0; 4];
        file_buffer.read_exact(&mut buffer[..]).ok()?;
        let kmer = KmerBits::from_le_bytes(buffer);
        file_buffer.read_exact(&mut buffer[..]).ok()?;
        let read_id = KmerBits::from_le_bytes(buffer);
        Some(KmerRead::new(Kmer::new(kmer), read_id))
    }
}

impl Ord for BufReaderKmerRead {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> Ordering {
        self.last_kmer_read.cmp(&other.last_kmer_read)
    }
}

impl PartialOrd for BufReaderKmerRead {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for BufReaderKmerRead {} // TODO?

impl PartialEq for BufReaderKmerRead {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.last_kmer_read == other.last_kmer_read
    }
}
