use std::cmp::Ordering;

use crate::kmer::Kmer;

#[derive(Debug)]
pub struct KmerRead {
    kmer: Kmer,
    read_id: u32,
}

impl KmerRead {
    #[inline(always)]
    pub fn new(kmer: Kmer, read_id: u32) -> Self {
        Self { kmer, read_id }
    }

    #[inline(always)]
    pub fn kmer(&self) -> &Kmer {
        &self.kmer
    }

    #[inline(always)]
    pub fn read_id(&self) -> u32 {
        self.read_id
    }
}

impl Ord for KmerRead {
    fn cmp(&self, other: &Self) -> Ordering {
        self.kmer()
            .cmp(&other.kmer)
            .then_with(|| self.read_id.cmp(&other.read_id))
    }
}

impl PartialOrd for KmerRead {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for KmerRead {}

impl PartialEq for KmerRead {
    fn eq(&self, other: &Self) -> bool {
        self.kmer == other.kmer && self.read_id == other.read_id
    }
}
