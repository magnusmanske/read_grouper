use crate::data_bucket::{BucketDataRead, BucketDataWrite};
use crate::KmerReverse;
use crate::{kmer::Kmer, ReadId};
use anyhow::Result;
use std::io::{BufReader, Read, Write};
use std::{cmp::Ordering, fs::File, io::BufWriter};

/// A kmer paired with a read id.
#[derive(Debug, Default, Clone)]
pub struct KmerRead<KmerBits> {
    kmer: Kmer<KmerBits>,
    read_id: ReadId,
}

impl<KmerBits: KmerReverse> KmerRead<KmerBits> {
    #[inline(always)]
    pub fn new(kmer: Kmer<KmerBits>, read_id: ReadId) -> Self {
        Self { kmer, read_id }
    }

    #[inline(always)]
    pub fn kmer(&self) -> &Kmer<KmerBits> {
        &self.kmer
    }

    #[inline(always)]
    pub fn read_id(&self) -> ReadId {
        self.read_id
    }
}

impl<KmerBits: KmerReverse> BucketDataWrite for KmerRead<KmerBits> {
    #[inline(always)]
    fn write(&self, buffer: &mut BufWriter<File>) -> Result<()> {
        self.kmer.write(buffer)?;
        buffer.write_all(&self.read_id().to_le_bytes())?;
        Ok(())
    }
}

impl<KmerBits: KmerReverse> BucketDataRead for KmerRead<KmerBits> {
    #[inline(always)]
    fn read(&mut self, file_buffer: &mut BufReader<File>) -> Result<()> {
        let mut buffer = [0; 4];
        self.kmer.read(file_buffer)?;
        file_buffer.read_exact(&mut buffer[..])?;
        self.read_id = ReadId::from_le_bytes(buffer);
        Ok(())
    }
}

impl<KmerBits: KmerReverse> Ord for KmerRead<KmerBits> {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> Ordering {
        self.kmer()
            .cmp(&other.kmer)
            .then_with(|| self.read_id.cmp(&other.read_id))
    }
}

impl<KmerBits: KmerReverse> PartialOrd for KmerRead<KmerBits> {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<KmerBits: KmerReverse> Eq for KmerRead<KmerBits> {}

impl<KmerBits: KmerReverse> PartialEq for KmerRead<KmerBits> {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.kmer == other.kmer && self.read_id == other.read_id
    }
}
