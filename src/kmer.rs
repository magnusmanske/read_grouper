use crate::{data_bucket::BucketDataWrite, multi_bucket_reader::BucketDataRead, KmerReverse};
use anyhow::Result;
use std::{
    fmt,
    fs::File,
    io::{BufReader, BufWriter},
};

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct Kmer<KmerBits> {
    bits: KmerBits,
}

impl<KmerBits: KmerReverse> Kmer<KmerBits> {
    #[inline(always)]
    pub fn new(kmer: KmerBits) -> Self {
        Self { bits: kmer }
    }

    #[inline(always)]
    fn build_kmer_pair(
        sequence_bases: &[u8],
        quality_scores: &[u8],
        min_base_quality: u8,
    ) -> Option<(KmerBits, KmerBits)> {
        let mut kmer = KmerBits::default();
        for i in 0..KmerBits::bases() {
            if quality_scores[i] < min_base_quality {
                return None;
            }
            let base: u8 = match sequence_bases[i] {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None,
            };
            kmer.add_base(base);
        }
        let reverse_kmer = kmer.reverse_complement(); // Reverse complement
        Some((kmer, reverse_kmer))
    }

    #[inline(always)]
    pub fn kmers_from_record(
        sequence: &[u8],
        quality_scores: &[u8],
        min_base_quality: u8,
    ) -> Vec<KmerBits> {
        if sequence.len() + 1 < KmerBits::bases() {
            // Sequence is too short
            return Vec::new();
        }
        let mut ret = Vec::with_capacity(sequence.len() - KmerBits::bases() + 1);

        // Generate first kmer
        let seq = &sequence[0..KmerBits::bases()];
        let qual = &quality_scores[0..KmerBits::bases()];
        let (mut kmer, mut reverse_complement_kmer) =
            match Self::build_kmer_pair(seq, qual, min_base_quality) {
                Some(x) => x,
                None => return ret,
            };
        ret.push(kmer.to_owned().min(reverse_complement_kmer));

        for i in KmerBits::bases()..(sequence.len()) {
            let base = sequence[i];
            if quality_scores[i] < min_base_quality {
                break; // Bad quality, abandon entire read
            }
            let base: u8 = match base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => break, // Weird IUPAC letter, abandon entire read
            };
            kmer.add_base(base);
            reverse_complement_kmer = kmer.reverse_complement(); // Reverse complement
            ret.push(kmer.to_owned().min(reverse_complement_kmer));
        }
        ret.sort();
        ret.dedup();
        ret
    }
}

impl<KmerBits: KmerReverse> fmt::Display for Kmer<KmerBits> {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        write!(fmt, "{}", self.bits)?;
        Ok(())
    }
}

impl<KmerBits: KmerReverse> BucketDataWrite for Kmer<KmerBits> {
    #[inline(always)]
    fn write(&self, buffer: &mut BufWriter<File>) -> Result<()> {
        self.bits.write(buffer)
    }
}

impl<KmerBits: KmerReverse> BucketDataRead for Kmer<KmerBits> {
    #[inline(always)]
    fn read(&mut self, file_buffer: &mut BufReader<File>) -> Result<()> {
        self.bits.read(file_buffer)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Kmer16, Kmer32};

    type K16 = Kmer<Kmer16>;
    type K32 = Kmer<Kmer32>;

    #[test]
    fn test_build_kmer_pair16() {
        let seq = b"ACGTACGTACGTGTAC";
        let qual = [
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
        ];
        let (kmer, reverse_complement) = K16::build_kmer_pair(seq, &qual, 40).unwrap();
        assert_eq!(kmer, Kmer16(0x1B1B1BB1));
        assert_eq!(reverse_complement, kmer.reverse_complement());

        let qual = [
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39,
        ];
        assert_eq!(K16::build_kmer_pair(seq, &qual, 40), None);
    }

    #[test]
    fn test_build_kmer_pair32() {
        let seq = b"ACGTACGTACGTGTACACGTACGTACGTGTAC";
        let qual = [
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
        ];
        let (kmer, reverse_complement) = K32::build_kmer_pair(seq, &qual, 40).unwrap();
        assert_eq!(kmer, Kmer32(1953185310873164721));
        assert_eq!(reverse_complement, kmer.reverse_complement());

        let qual = [
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 40, 40, 40, 40, 40, 39,
        ];
        assert_eq!(K32::build_kmer_pair(seq, &qual, 40), None);
    }

    #[test]
    fn test_kmers_from_record_incremental16() {
        let seq = b"ACGTACGTACGTGTACACGTACGTACGTGTAC";
        let qual = [
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
        ];
        let kmers = K16::kmers_from_record(seq, &qual, 40);
        assert_eq!(
            kmers,
            [
                296858043, 454761393, 454799643, 1187432172, 1819045572, 1819198572, 1858366572,
                2981214993, 2981826993
            ]
            .into_iter()
            .map(Kmer16)
            .collect::<Vec<Kmer16>>()
        );
    }
}
