use crate::data_bucket::{BucketDataRead, BucketDataWrite};
use crate::KmerBits;
use crate::{kmer::Kmer, ReadId};
use anyhow::Result;
use std::io::{BufReader, Read, Write};
use std::{cmp::Ordering, fs::File, io::BufWriter};

/// A kmer paired with a read id.
#[derive(Debug, Default)]
pub struct KmerRead {
    kmer: Kmer,
    read_id: ReadId,
}

impl KmerRead {
    #[inline(always)]
    pub fn new(kmer: Kmer, read_id: ReadId) -> Self {
        Self { kmer, read_id }
    }

    #[inline(always)]
    pub fn kmer(&self) -> &Kmer {
        &self.kmer
    }

    #[inline(always)]
    pub fn read_id(&self) -> ReadId {
        self.read_id
    }
}

impl BucketDataWrite for KmerRead {
    #[inline(always)]
    fn write(&self, buffer: &mut BufWriter<File>) -> Result<()> {
        buffer.write_all(&self.kmer().to_le_bytes())?;
        buffer.write_all(&self.read_id().to_le_bytes())?;
        Ok(())
    }
}

impl BucketDataRead for KmerRead {
    #[inline(always)]
    fn read(&mut self, file_buffer: &mut BufReader<File>) -> Result<()> {
        let mut buffer = [0; 4];
        file_buffer.read_exact(&mut buffer[..])?;
        self.kmer = Kmer::new(KmerBits::from_le_bytes(buffer));
        file_buffer.read_exact(&mut buffer[..])?;
        self.read_id = KmerBits::from_le_bytes(buffer);
        Ok(())
    }
}

impl Ord for KmerRead {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> Ordering {
        self.kmer()
            .cmp(&other.kmer)
            .then_with(|| self.read_id.cmp(&other.read_id))
    }
}

impl PartialOrd for KmerRead {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for KmerRead {}

impl PartialEq for KmerRead {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.kmer == other.kmer && self.read_id == other.read_id
    }
}
