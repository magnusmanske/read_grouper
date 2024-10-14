/// The minimum and maximum number of reads to constitute a valid group.
/// Guaranteed to have min>=2.
#[derive(Debug)]
pub struct MinMaxReads {
    min: usize,
    max: usize,
}

impl MinMaxReads {
    #[inline(always)]
    pub fn new(min: usize, max: usize) -> Self {
        Self {
            min: min.max(2), // Ensure >1
            max,
        }
    }

    #[inline(always)]
    pub fn is_valid(&self, value: usize) -> bool {
        value >= self.min && value <= self.max
    }
}

impl Default for MinMaxReads {
    #[inline(always)]
    fn default() -> Self {
        Self {
            min: 2,
            max: usize::MAX,
        }
    }
}
