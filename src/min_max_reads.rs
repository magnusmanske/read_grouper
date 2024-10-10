/// The minimum and maximum nomber of reads to constitute a valid group.
/// Guaranteed to have min>=1.
#[derive(Debug)]
pub struct MinMaxReads {
    min: usize,
    max: usize,
}

impl MinMaxReads {
    pub fn new(min: usize, max: usize) -> Self {
        Self {
            min: min.max(1), // Ensure non-0
            max,
        }
    }

    pub fn is_valid(&self, value: usize) -> bool {
        value >= self.min && value <= self.max
    }
}

impl Default for MinMaxReads {
    fn default() -> Self {
        Self {
            min: 1,
            max: usize::MAX,
        }
    }
}
