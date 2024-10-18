use rayon::prelude::*;
use std::{cmp::Ordering, fmt};

use crate::KmerReverse;

#[derive(Clone, Debug, Default)]
pub struct Kmer<KmerBits> {
    bits: KmerBits,
}

impl<KmerBits: KmerReverse> Kmer<KmerBits> {
    #[inline(always)]
    pub fn new(kmer: KmerBits) -> Self {
        Self { bits: kmer }
    }

    #[inline(always)]
    pub fn to_le_bytes(&self) -> [u8] {
        self.0.to_le_bytes()
    }

    #[inline(always)]
    fn build_kmer_pair(
        sequence_bases: &[u8],
        quality_scores: &[u8],
        min_base_quality: u8,
    ) -> Option<(KmerBits, KmerBits)> {
        let mut kmer = KmerBits(0);
        for i in 0..KmerBits::bases() {
            if quality_scores[i] < min_base_quality {
                return None;
            }
            let base_forward: KmerBits = match sequence_bases[i] {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None,
            };
            kmer = (kmer << 2) | base_forward;
        }
        let reverse_kmer = kmer.reverse_complement(); // Reverse complement
        Some((kmer, reverse_kmer))
    }

    #[inline(always)]
    pub fn kmers_from_record_incremental(
        sequence: &[u8],
        quality_scores: &[u8],
        min_base_quality: u8,
    ) -> Vec<KmerBits> {
        let mut ret = Vec::with_capacity(sequence.len() - KmerBits::bases() + 1);

        // Generate first kmer
        let seq = &sequence[0..KmerBits::bases];
        let qual = &quality_scores[0..KmerBits::bases];
        let (mut kmer, mut reverse_complement_kmer) =
            match Self::build_kmer_pair(seq, qual, min_base_quality) {
                Some(x) => x,
                None => return ret,
            };
        ret.push(kmer.min(reverse_complement_kmer));

        for i in KmerBits::bases..(sequence.len()) {
            let base = sequence[i];
            if quality_scores[i] < min_base_quality {
                break; // Bad quality, abandon entire read
            }
            let base_forward: KmerBits = match base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => break, // Weird IUPAC letter, abandon entire read
            };
            kmer = (kmer << 2) | base_forward;
            reverse_complement_kmer = Self::reverse_by_two_bit_groups_u32(!kmer); // Reverse complement
            ret.push(kmer.min(reverse_complement_kmer));
        }
        ret.sort();
        ret.dedup();
        ret
    }
}

impl<KmerBits: KmerReverse> fmt::Display for Kmer<KmerBits> {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        for pos in 0..KmerBits::bases {
            let base = match (self.0 >> (2 * (KmerBits::bases - pos - 1))) & 0b11 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => unreachable!(),
            };
            write!(fmt, "{}", base as char)?;
        }
        Ok(())
    }
}

impl<KmerBits: KmerReverse> Ord for Kmer<KmerBits> {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
    }
}

impl<KmerBits: KmerReverse> PartialOrd for Kmer<KmerBits> {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<KmerBits: KmerReverse> Eq for Kmer<KmerBits> {}

impl<KmerBits: KmerReverse> PartialEq for Kmer<KmerBits> {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::*;

    #[test]
    fn test_build_kmer_pair() {
        let seq = b"ACGTACGTACGTGTAC";
        let qual = [
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
        ];
        let (kmer, reverse_complement): Kmer<Kmer16> =
            Kmer::build_kmer_pair(seq, &qual, 40).unwrap();
        assert_eq!(kmer, 0x1B1B1BB1);
        assert_eq!(
            reverse_complement,
            Kmer::reverse_by_two_bit_groups_u32(!kmer)
        );

        let qual = [
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39,
        ];
        assert_eq!(Kmer::build_kmer_pair(seq, &qual, 40), None);
    }

    #[test]
    fn test_kmers_from_record_de_novo() {
        let seq = b"ACGTACGTACGTGTACACGTACGTACGTGTAC";
        let qual = [
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
        ];
        let kmers = Kmer::_kmers_from_record_de_novo(seq, &qual, 40);
        assert_eq!(
            kmers,
            [
                296858043, 454761393, 454799643, 1187432172, 1819045572, 1819198572, 1858366572,
                2981214993, 2981826993
            ]
        );
    }

    #[test]
    fn test_kmers_from_record_incremental() {
        let seq = b"ACGTACGTACGTGTACACGTACGTACGTGTAC";
        let qual = [
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
        ];
        let kmers = Kmer::kmers_from_record_incremental(seq, &qual, 40);
        assert_eq!(
            kmers,
            [
                296858043, 454761393, 454799643, 1187432172, 1819045572, 1819198572, 1858366572,
                2981214993, 2981826993
            ]
        );
    }
}
