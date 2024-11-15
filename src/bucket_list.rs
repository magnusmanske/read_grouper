use std::path::Path;

use crate::ReadId;
use anyhow::{anyhow, Result};

#[derive(Default, Debug)]
pub struct BucketList {
    filenames: Vec<String>,
    number_of_reads: ReadId,
    sample_name: String,
}

impl BucketList {
    pub fn new(sample_name: String, filenames: Vec<String>, number_of_reads: ReadId) -> Self {
        Self {
            filenames,
            number_of_reads,
            sample_name,
        }
    }

    pub fn filenames(&self) -> &Vec<String> {
        &self.filenames
    }

    pub fn number_of_reads(&self) -> ReadId {
        self.number_of_reads
    }

    pub fn sample_name(&self) -> &str {
        &self.sample_name
    }

    pub fn file_path_to_sample_name(file_path: &Path) -> Result<String> {
        Ok(file_path
            .file_stem()
            .ok_or_else(|| {
                anyhow!("Could not generate sample name from filename [1] {file_path:?}")
            })?
            .to_str()
            .ok_or_else(|| {
                anyhow!("Could not generate sample name from filename [2] {file_path:?}")
            })?
            .to_string())
    }
}
