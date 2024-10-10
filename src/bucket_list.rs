#[derive(Default, Debug)]
pub struct BucketList {
    filenames: Vec<String>,
    number_of_reads: u32,
    sample_name: String,
}

impl BucketList {
    pub fn new(sample_name: String, filenames: Vec<String>, number_of_reads: u32) -> Self {
        Self {
            filenames,
            number_of_reads,
            sample_name,
        }
    }

    pub fn filenames(&self) -> &Vec<String> {
        &self.filenames
    }

    pub fn number_of_reads(&self) -> u32 {
        self.number_of_reads
    }

    pub fn sample_name(&self) -> &str {
        &self.sample_name
    }
}
