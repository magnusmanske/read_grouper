use crate::{
    data_bucket::{BucketDataRead, BucketDataWrite},
    kmer::Kmer,
    KmerBits, ReadId,
};
use anyhow::Result;
use std::{
    cmp::Ordering,
    fs::File,
    io::{BufReader, BufWriter, Read, Write},
};

#[derive(Debug)]
pub struct ReadPairKmer {
    read1: ReadId,
    read2: ReadId,
    kmer: Kmer,
}

impl ReadPairKmer {
    pub fn new(read1: ReadId, read2: ReadId, kmer: &Kmer) -> Self {
        Self {
            read1,
            read2,
            kmer: kmer.to_owned(),
        }
    }

    pub fn read1(&self) -> ReadId {
        self.read1
    }

    pub fn read2(&self) -> ReadId {
        self.read2
    }

    pub fn kmer(&self) -> Kmer {
        self.kmer.to_owned()
    }
}

impl BucketDataWrite for ReadPairKmer {
    #[inline(always)]
    fn write(&self, buffer: &mut BufWriter<File>) -> Result<()> {
        buffer.write_all(&self.read1().to_le_bytes())?;
        buffer.write_all(&self.read2().to_le_bytes())?;
        buffer.write_all(&self.kmer().to_le_bytes())?;
        Ok(())
    }
}

impl BucketDataRead for ReadPairKmer {
    #[inline(always)]
    fn read(&mut self, file_buffer: &mut BufReader<File>) -> Result<()> {
        let mut buffer = [0; 4];
        file_buffer.read_exact(&mut buffer[..])?;
        self.read1 = ReadId::from_le_bytes(buffer);
        file_buffer.read_exact(&mut buffer[..])?;
        self.read2 = ReadId::from_le_bytes(buffer);
        file_buffer.read_exact(&mut buffer[..])?;
        self.kmer = Kmer::new(KmerBits::from_le_bytes(buffer));
        Ok(())
    }
}

impl Ord for ReadPairKmer {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> Ordering {
        self.read1
            .cmp(&other.read1)
            .then_with(|| self.read2.cmp(&other.read2))
            .then_with(|| self.kmer.cmp(&other.kmer))
    }
}

impl PartialOrd for ReadPairKmer {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for ReadPairKmer {}

impl PartialEq for ReadPairKmer {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.kmer == other.kmer && self.read1 == other.read1 && self.read2 == other.read2
    }
}
