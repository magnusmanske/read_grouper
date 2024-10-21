use crate::data_bucket::{BucketDataRead, BucketDataWrite};
use anyhow::Result;
use std::{
    fmt::{self, Display},
    fs::File,
    io::{BufReader, BufWriter, Read, Write},
    ops::{Shl, Shr},
};

pub trait KmerReverse:
    Default
    + Eq
    + PartialEq
    + Ord
    + PartialOrd
    + Clone
    + BucketDataWrite
    + BucketDataRead
    + Display
    + Shl
    + Shr
    + Send
    + 'static
{
    fn reverse_complement(&self) -> Self;
    fn bytes() -> usize;
    fn add_base(&mut self, base: u8);

    fn bases() -> usize {
        Self::bytes() * 4
    }
}

// #[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct Kmer16(pub u32);

// #[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct Kmer32(pub u64);

impl KmerReverse for Kmer16 {
    #[inline(always)]
    fn bytes() -> usize {
        4
    }

    #[inline(always)]
    fn add_base(&mut self, base: u8) {
        self.0 = (self.0 << 2) | (base as u32);
    }

    #[inline(always)]
    fn reverse_complement(&self) -> Self {
        // println!("Original : {value:032b}");
        let mut value = ((self.0 & 0x33333333) << 2) | ((self.0 & 0xCCCCCCCC) >> 2); // swap adjacent pairs
        value = ((value & 0x0F0F0F0F) << 4) | ((value & 0xF0F0F0F0) >> 4); // swap nibbles
        value = value.swap_bytes();
        Self(!value)
    }
}

impl Shl<Kmer16> for Kmer16 {
    type Output = Self;

    #[inline(always)]
    fn shl(self, Self(rhs): Self) -> Self::Output {
        let Self(lhs) = self;
        Self(lhs << rhs)
    }
}

impl Shr<Kmer16> for Kmer16 {
    type Output = Self;

    #[inline(always)]
    fn shr(self, Self(rhs): Self) -> Self::Output {
        let Self(lhs) = self;
        Self(lhs >> rhs)
    }
}

impl BucketDataWrite for Kmer16 {
    #[inline(always)]
    fn write(&self, buffer: &mut BufWriter<File>) -> Result<()> {
        buffer.write_all(&self.0.to_le_bytes())?;
        Ok(())
    }
}

impl BucketDataRead for Kmer16 {
    #[inline(always)]
    fn read(&mut self, file_buffer: &mut BufReader<File>) -> Result<()> {
        let mut buffer = [0; 4];
        file_buffer.read_exact(&mut buffer[..])?;
        self.0 = u32::from_le_bytes(buffer);
        Ok(())
    }
}

impl Display for Kmer16 {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        for pos in 0..Self::bases() {
            let base = match (self.0 >> (2 * (Self::bases() - pos - 1))) & 0b11 {
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

impl KmerReverse for Kmer32 {
    #[inline(always)]
    fn bytes() -> usize {
        8
    }

    #[inline(always)]
    fn add_base(&mut self, base: u8) {
        self.0 = (self.0 << 2) | (base as u64);
    }

    #[inline(always)]
    fn reverse_complement(&self) -> Self {
        // println!("Original : {:064b}", self.0);
        let mut value = ((self.0 & 0x3333333333333333) << 2) | ((self.0 & 0xCCCCCCCCCCCCCCCC) >> 2); // swap adjacent pairs
        value = ((value & 0x0F0F0F0F0F0F0F0F) << 4) | ((value & 0xF0F0F0F0F0F0F0F0) >> 4); // swap nibbles
        value = value.swap_bytes();
        Self(!value)
    }
}

impl BucketDataWrite for Kmer32 {
    #[inline(always)]
    fn write(&self, buffer: &mut BufWriter<File>) -> Result<()> {
        buffer.write_all(&self.0.to_le_bytes())?;
        Ok(())
    }
}

impl BucketDataRead for Kmer32 {
    #[inline(always)]
    fn read(&mut self, file_buffer: &mut BufReader<File>) -> Result<()> {
        let mut buffer = [0; 8];
        file_buffer.read_exact(&mut buffer[..])?;
        self.0 = u64::from_le_bytes(buffer);
        Ok(())
    }
}

impl Shl<Kmer32> for Kmer32 {
    type Output = Self;

    #[inline(always)]
    fn shl(self, Self(rhs): Self) -> Self::Output {
        let Self(lhs) = self;
        Self(lhs << rhs)
    }
}

impl Shr<Kmer32> for Kmer32 {
    type Output = Self;

    #[inline(always)]
    fn shr(self, Self(rhs): Self) -> Self::Output {
        let Self(lhs) = self;
        Self(lhs >> rhs)
    }
}

impl Display for Kmer32 {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        for pos in 0..Self::bases() {
            let base = match (self.0 >> (2 * (Self::bases() - pos - 1))) & 0b11 {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer16() {
        assert_eq!(Kmer16(3 << 30).reverse_complement(), Kmer16(!3));
        assert_eq!(Kmer16(3).reverse_complement(), Kmer16(!(3 << 30)));
        assert_eq!(Kmer16(3 << 2).reverse_complement(), Kmer16(!(3 << 28)));
        assert_eq!(Kmer16(3 << 18).reverse_complement(), Kmer16(!(3 << 12)));
        assert_eq!(Kmer16(3 << 16).reverse_complement(), Kmer16(!(3 << 14)));
        assert_eq!(Kmer16(3 << 14).reverse_complement(), Kmer16(!(3 << 16)));
        assert_eq!(Kmer16(3 << 12).reverse_complement(), Kmer16(!(3 << 18)));
        let value = Kmer16(12345 | 6789 << 8 | 65432 << 16 | 23456 << 24);
        let reverse = value.reverse_complement();
        let original = reverse.reverse_complement();
        assert_eq!(value, original);
    }

    #[test]
    fn test_kmer32() {
        assert_eq!(Kmer32(3 << 62).reverse_complement(), Kmer32(!3));
        assert_eq!(Kmer32(3).reverse_complement(), Kmer32(!(3 << 62)));
        assert_eq!(Kmer32(3 << 2).reverse_complement(), Kmer32(!(3 << 60)));
        assert_eq!(Kmer32(3 << 16).reverse_complement(), Kmer32(!(3 << 46)));
        assert_eq!(Kmer32(3 << 46).reverse_complement(), Kmer32(!(3 << 16)));
        assert_eq!(Kmer32(3 << 34).reverse_complement(), Kmer32(!(3 << 28)));
        assert_eq!(Kmer32(3 << 32).reverse_complement(), Kmer32(!(3 << 30)));
        assert_eq!(Kmer32(3 << 30).reverse_complement(), Kmer32(!(3 << 32)));
        assert_eq!(Kmer32(3 << 28).reverse_complement(), Kmer32(!(3 << 34)));
        let value = Kmer32(
            12345
                | 6789 << 8
                | 65432 << 16
                | 23456 << 24
                | 12345 << 32
                | 6789 << 40
                | 65432 << 48
                | 23456 << 56,
        );
        let reverse = value.reverse_complement();
        let original = reverse.reverse_complement();
        assert_eq!(value, original);
    }
}
