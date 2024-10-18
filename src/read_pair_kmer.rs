use crate::{
    data_bucket::{BucketDataRead, BucketDataWrite},
    kmer::Kmer,
    KmerReverse, ReadId,
};
use anyhow::Result;
use std::{
    cmp::Ordering,
    fs::File,
    io::{BufReader, BufWriter, Read, Write},
};

#[derive(Debug)]
pub struct ReadPairKmer<KmerBits> {
    read1: ReadId,
    read2: ReadId,
    kmer: Kmer<KmerBits>,
}

impl<KmerBits: KmerReverse> ReadPairKmer<KmerBits> {
    pub fn new(read1: ReadId, read2: ReadId, kmer: &Kmer<KmerBits>) -> Self {
        Self {
            read1,
            read2,
            kmer: kmer.to_owned(),
        }
    }
}

impl<KmerBits: KmerReverse> BucketDataWrite for ReadPairKmer<KmerBits> {
    #[inline(always)]
    fn write(&self, buffer: &mut BufWriter<File>) -> Result<()> {
        buffer.write_all(&self.read1.to_le_bytes())?;
        buffer.write_all(&self.read2.to_le_bytes())?;
        self.kmer.write(buffer)?;
        Ok(())
    }
}

impl<KmerBits: KmerReverse> BucketDataRead for ReadPairKmer<KmerBits> {
    #[inline(always)]
    fn read(&mut self, file_buffer: &mut BufReader<File>) -> Result<()> {
        let mut buffer = [0; 4];
        file_buffer.read_exact(&mut buffer[..])?;
        self.read1 = ReadId::from_le_bytes(buffer);
        file_buffer.read_exact(&mut buffer[..])?;
        self.read2 = ReadId::from_le_bytes(buffer);
        self.kmer.read(file_buffer)?;
        Ok(())
    }
}

impl<KmerBits: KmerReverse> Ord for ReadPairKmer<KmerBits> {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> Ordering {
        self.read1
            .cmp(&other.read1)
            .then_with(|| self.read2.cmp(&other.read2))
            .then_with(|| self.kmer.cmp(&other.kmer))
    }
}

impl<KmerBits: KmerReverse> PartialOrd for ReadPairKmer<KmerBits> {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<KmerBits: KmerReverse> Eq for ReadPairKmer<KmerBits> {}

impl<KmerBits: KmerReverse> PartialEq for ReadPairKmer<KmerBits> {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.kmer == other.kmer && self.read1 == other.read1 && self.read2 == other.read2
    }
}
