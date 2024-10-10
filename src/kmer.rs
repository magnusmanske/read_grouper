use rayon::prelude::*;
use std::cmp::Ordering;

use crate::KmerBits;

const BASES_PER_KMER: usize = 16;

#[derive(Clone, Debug)]
pub struct Kmer(KmerBits);

impl Kmer {
    #[inline(always)]
    pub fn new(kmer: KmerBits) -> Self {
        Self(kmer)
    }

    #[inline(always)]
    pub fn to_le_bytes(&self) -> [u8; 4] {
        self.0.to_le_bytes()
    }

    #[inline(always)]
    fn reverse_by_two_bit_groups_u32(value: KmerBits) -> KmerBits {
        // println!("Original : {value:032b}");
        let mut value = ((value & 0x33333333) << 2) | ((value & 0xCCCCCCCC) >> 2); // swap adjacent pairs

        // println!("Bit pairs: {value:032b}");
        value = ((value & 0x0F0F0F0F) << 4) | ((value & 0xF0F0F0F0) >> 4); // swap nibbles

        // println!("Nibbles  : {value:032b}");
        value = ((value & 0x00FF00FF) << 8) | ((value & 0xFF00FF00) >> 8); // swap bytes

        // println!("Bytes    : {value:032b}");
        value = ((value & 0x0000FFFF) << 16) | ((value & 0xFFFF0000) >> 16); // sawp words

        // println!("Words    : {value:032b}");
        value
    }

    #[inline(always)]
    fn build_kmer_min(seq: &[u8], quality_scores: &[u8], min_base_quality: u8) -> Option<KmerBits> {
        let mut kmer = 0;
        for i in 0..BASES_PER_KMER {
            let base = seq[i];
            if quality_scores[i] < min_base_quality {
                return None;
            }
            let base_forward: KmerBits = match base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None,
            };
            kmer = (kmer << 2) | base_forward;
        }
        let reverse_kmer = Self::reverse_by_two_bit_groups_u32(!kmer); // Reverse complement
        Some(kmer.min(reverse_kmer))
    }

    pub fn kmers_from_record_de_novo(
        sequence: &[u8],
        qualities: &[u8],
        min_base_quality: u8,
    ) -> Vec<KmerBits> {
        let range = 0..(sequence.len() - BASES_PER_KMER);
        let mut ret = range
            .into_par_iter()
            .filter_map(|start| {
                let seq = &sequence[start..start + BASES_PER_KMER];
                let qual = &qualities[start..start + BASES_PER_KMER];
                Self::build_kmer_min(seq, qual, min_base_quality)
            })
            .collect::<Vec<_>>();
        ret.sort();
        ret.dedup();
        ret
    }

    #[inline(always)]
    fn build_kmer_pair(
        sequence_bases: &[u8],
        quality_scores: &[u8],
        min_base_quality: u8,
    ) -> Option<(KmerBits, KmerBits)> {
        let mut kmer = 0;
        for i in 0..BASES_PER_KMER {
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
        let reverse_kmer = Self::reverse_by_two_bit_groups_u32(!kmer); // Reverse complement
        Some((kmer, reverse_kmer))
    }

    #[inline(always)]
    pub fn kmers_from_record_incremental(
        sequence: &[u8],
        quality_scores: &[u8],
        min_base_quality: u8,
    ) -> Vec<KmerBits> {
        let mut ret = Vec::with_capacity(sequence.len() - BASES_PER_KMER + 1);

        // Generate first kmer
        let seq = &sequence[0..BASES_PER_KMER];
        let qual = &quality_scores[0..BASES_PER_KMER];
        let (mut kmer, mut reverse_complement_kmer) =
            match Self::build_kmer_pair(seq, qual, min_base_quality) {
                Some(x) => x,
                None => return ret,
            };
        ret.push(kmer.min(reverse_complement_kmer));

        for i in BASES_PER_KMER..(sequence.len()) {
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

impl Ord for Kmer {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
    }
}

impl PartialOrd for Kmer {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for Kmer {}

impl PartialEq for Kmer {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_by_two_bit_groups_u32() {
        assert_eq!(Kmer::reverse_by_two_bit_groups_u32(3 << 30), 3);
        assert_eq!(Kmer::reverse_by_two_bit_groups_u32(3), 3 << 30);
        assert_eq!(Kmer::reverse_by_two_bit_groups_u32(3 << 2), 3 << 28);
        assert_eq!(Kmer::reverse_by_two_bit_groups_u32(3 << 18), 3 << 12);
        assert_eq!(Kmer::reverse_by_two_bit_groups_u32(3 << 16), 3 << 14);
        assert_eq!(Kmer::reverse_by_two_bit_groups_u32(3 << 14), 3 << 16);
        assert_eq!(Kmer::reverse_by_two_bit_groups_u32(3 << 12), 3 << 18);
        let value: KmerBits = 12345 | 6789 << 8 | 65432 << 16 | 23456 << 24;
        let reverse = Kmer::reverse_by_two_bit_groups_u32(value);
        let original = Kmer::reverse_by_two_bit_groups_u32(reverse);
        assert_eq!(value, original);
    }

    #[test]
    fn test_build_kmer_pair() {
        let seq = b"ACGTACGTACGTGTAC";
        let qual = [
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
        ];
        let (kmer, reverse_complement) = Kmer::build_kmer_pair(seq, &qual, 40).unwrap();
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
}
