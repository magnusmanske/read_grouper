use std::fmt;

use crate::{read_grouper::GenomicPosition, KmerReverse};

pub type KmerPosition = u8; // TODO enough?

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct KmerWithPosition<KmerBits> {
    bits: KmerBits,
    position: KmerPosition,
    is_reverse_complement: bool,
}

impl<KmerBits: KmerReverse> KmerWithPosition<KmerBits> {
    #[inline(always)]
    pub fn from_pair(
        kmer: KmerBits,
        reverse_complement_kmer: KmerBits,
        position: KmerPosition,
    ) -> Self {
        match kmer <= reverse_complement_kmer {
            true => Self {
                bits: kmer,
                position,
                is_reverse_complement: false,
            },
            false => Self {
                bits: reverse_complement_kmer,
                position,
                is_reverse_complement: true,
            },
        }
    }

    #[inline(always)]
    pub fn genomic_position(&self, read_length: GenomicPosition) -> GenomicPosition {
        match self.is_reverse_complement {
            true => read_length - (self.position as GenomicPosition),
            false => self.position as GenomicPosition,
        }
    }

    #[inline(always)]
    pub fn kmer(&self) -> KmerBits {
        self.bits.to_owned()
    }

    // #[inline(always)]
    // pub fn position(&self) -> KmerPosition {
    //     self.position
    // }
}

impl<KmerBits: KmerReverse> fmt::Display for KmerWithPosition<KmerBits> {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        write!(
            fmt,
            "{}/{}/{}",
            self.bits, self.position, self.is_reverse_complement
        )?;
        Ok(())
    }
}
