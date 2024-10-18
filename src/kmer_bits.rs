pub trait KmerReverse: num::Integer {
    fn reverse_complement(self) -> Self;
    fn bytes() -> usize;
    fn bases() -> usize {
        Self::bytes() * 4
    }
}

pub struct Kmer16(u32);
pub struct Kmer32(u64);

impl KmerReverse for Kmer16 {
    #[inline(always)]
    fn bytes() -> usize {
        4
    }

    #[inline(always)]
    fn reverse_complement(self) -> Self {
        // println!("Original : {value:032b}");
        let mut value = ((self.0 & 0x33333333) << 2) | ((self.0 & 0xCCCCCCCC) >> 2); // swap adjacent pairs
        value = ((value & 0x0F0F0F0F) << 4) | ((value & 0xF0F0F0F0) >> 4); // swap nibbles
        value = value.swap_bytes();
        Self(!value)
    }
}

impl KmerReverse for Kmer32 {
    #[inline(always)]
    fn bytes() -> usize {
        8
    }

    #[inline(always)]
    fn reverse_complement(self) -> Self {
        // println!("Original : {value:032b}");
        let mut value = ((self.0 & 0x3333333333333333) << 2) | ((self.0 & 0xCCCCCCCCCCCCCCCC) >> 2); // swap adjacent pairs
        value = ((value & 0x0F0F0F0F0F0F0F0F) << 4) | ((value & 0xF0F0F0F0F0F0F0F0) >> 4); // swap nibbles
        value = value.swap_bytes();
        Self(!value)
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn test_reverse_by_two_bit_groups_u32() {
//         assert_eq!(Kmer::reverse_by_two_bit_groups_u32(3 << 30), 3);
//         assert_eq!(Kmer::reverse_by_two_bit_groups_u32(3), 3 << 30);
//         assert_eq!(Kmer::reverse_by_two_bit_groups_u32(3 << 2), 3 << 28);
//         assert_eq!(Kmer::reverse_by_two_bit_groups_u32(3 << 18), 3 << 12);
//         assert_eq!(Kmer::reverse_by_two_bit_groups_u32(3 << 16), 3 << 14);
//         assert_eq!(Kmer::reverse_by_two_bit_groups_u32(3 << 14), 3 << 16);
//         assert_eq!(Kmer::reverse_by_two_bit_groups_u32(3 << 12), 3 << 18);
//         let value: KmerBits = 12345 | 6789 << 8 | 65432 << 16 | 23456 << 24;
//         let reverse = Kmer::reverse_by_two_bit_groups_u32(value);
//         let original = Kmer::reverse_by_two_bit_groups_u32(reverse);
//         assert_eq!(value, original);
//     }

//     #[test]
//     fn test_reverse_by_two_bit_groups_u64() {
//         assert_eq!(Kmer::reverse_by_two_bit_groups_u64(3 << 62), 3);
//         assert_eq!(Kmer::reverse_by_two_bit_groups_u64(3), 3 << 62);
//         assert_eq!(Kmer::reverse_by_two_bit_groups_u64(3 << 2), 3 << 60);
//         assert_eq!(Kmer::reverse_by_two_bit_groups_u64(3 << 16), 3 << 46);
//         assert_eq!(Kmer::reverse_by_two_bit_groups_u64(3 << 46), 3 << 16);
//         assert_eq!(Kmer::reverse_by_two_bit_groups_u64(3 << 34), 3 << 28);
//         assert_eq!(Kmer::reverse_by_two_bit_groups_u64(3 << 32), 3 << 30);
//         assert_eq!(Kmer::reverse_by_two_bit_groups_u64(3 << 30), 3 << 32);
//         assert_eq!(Kmer::reverse_by_two_bit_groups_u64(3 << 28), 3 << 34);
//         let value: u64 = 12345
//             | 6789 << 8
//             | 65432 << 16
//             | 23456 << 24
//             | 12345 << 32
//             | 6789 << 40
//             | 65432 << 48
//             | 23456 << 56;
//         let reverse = Kmer::reverse_by_two_bit_groups_u64(value);
//         let original = Kmer::reverse_by_two_bit_groups_u64(reverse);
//         assert_eq!(value, original);
//     }
// }
