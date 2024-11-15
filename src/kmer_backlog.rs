use crate::{
    kmer_with_position::KmerWithPosition,
    read_grouper::{GenomicPosition, ScaffoldId},
    KmerReverse, ReadId,
};

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct KmerBacklog<KmerBits> {
    read1_id: ReadId,
    read1_scaffold: ScaffoldId,
    read1_position: GenomicPosition,
    read1_length: GenomicPosition,
    kwps: KmerWithPosition<KmerBits>,
}

impl<KmerBits: KmerReverse> KmerBacklog<KmerBits> {
    #[inline(always)]
    pub fn new(
        read1_id: ReadId,
        read1_scaffold: ScaffoldId,
        read1_position: GenomicPosition,
        read1_length: GenomicPosition,
        kwps: KmerWithPosition<KmerBits>,
    ) -> Self {
        Self {
            read1_id,
            read1_scaffold,
            read1_position,
            read1_length,
            kwps,
        }
    }

    pub fn kwps(&self) -> KmerWithPosition<KmerBits> {
        self.kwps.to_owned()
    }

    pub fn read1_id(&self) -> ReadId {
        self.read1_id
    }

    pub fn read1_position(&self) -> GenomicPosition {
        self.read1_position
    }

    pub fn read1_length(&self) -> GenomicPosition {
        self.read1_length
    }
}
