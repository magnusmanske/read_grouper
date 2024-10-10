use crate::{kmer::Kmer, ReadId};
use std::cmp::Ordering;

/// A kmer paired with a read id.
#[derive(Debug)]
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
