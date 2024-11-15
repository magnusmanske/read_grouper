use crate::{
    bucket_list::BucketList,
    data_bucket::{BucketDataWrite, DataBucket},
    kmer::Kmer,
    kmer_backlog::KmerBacklog,
    kmer_read::KmerRead,
    kmer_with_position::KmerWithPosition,
    min_max_reads::MinMaxReads,
    multi_bucket_reader::MultiBucketReader,
    read_pair_kmer::ReadPairKmer,
    KmerReverse, ReadId,
};
use anyhow::{anyhow, Result};
use bam::RecordReader;
use flate2::read::GzDecoder;
use std::fs::File;
use std::{
    collections::HashMap,
    io::{self, BufRead, BufWriter},
    marker::PhantomData,
    path::Path,
};

const DEFAULT_MIN_BASE_QUALITY: u8 = 20;
const DEFAULT_MAX_BUCKET_SIZE: usize = 1_000_000; // kmer-read-pairs

pub type GenomicPosition = i32;
pub type ScaffoldId = u64;

#[derive(Default, Debug)]
pub struct ReadGrouper<KmerBits> {
    bucket_dir: String,
    min_base_quality: u8,
    max_bucket_size: usize,
    kmer_bits_dummy: PhantomData<KmerBits>,
}

impl<KmerBits: KmerReverse> ReadGrouper<KmerBits> {
    pub fn new(bucket_dir: &str) -> Self {
        Self {
            bucket_dir: bucket_dir.to_string(),
            min_base_quality: DEFAULT_MIN_BASE_QUALITY,
            max_bucket_size: DEFAULT_MAX_BUCKET_SIZE,
            kmer_bits_dummy: PhantomData,
        }
    }

    pub fn set_max_bucket_size(&mut self, max_bucket_size: usize) {
        self.max_bucket_size = max_bucket_size.max(DEFAULT_MAX_BUCKET_SIZE); // No less than default size
    }

    pub fn read_fasta_file(&self, file_path: &str) -> Result<BucketList> {
        let path = Path::new(file_path);
        let sample_name = BucketList::file_path_to_sample_name(path)?;
        let mut out_bucket: DataBucket<Kmer<KmerBits>> = DataBucket::new(
            self.max_bucket_size,
            &self.bucket_dir,
            &sample_name,
            "from_fasta",
        );

        let mut sequence_name = String::new();
        let mut sequence = String::new();
        let reader = Self::get_uncompressed_or_gz_reader(path)?;
        for line in reader.lines() {
            let line = line.unwrap();
            let line = line.trim();
            if let Some(end) = line.strip_prefix('>') {
                if !sequence_name.is_empty() {
                    self.store_sequence_kmers(&sequence, &mut out_bucket);
                }
                sequence_name = end.trim().to_string();
                sequence.clear();
            } else if !sequence_name.is_empty() {
                sequence += &line.to_ascii_uppercase();
            }
        }
        if !sequence_name.is_empty() {
            self.store_sequence_kmers(&sequence, &mut out_bucket);
        }

        // Write final bucket to disk
        let filenames = out_bucket.finish()?;

        // Create metadata to return
        let bucket_list = BucketList::new(sample_name, filenames, 0);
        Ok(bucket_list)
    }

    pub fn get_uncompressed_or_gz_reader(path: &Path) -> Result<Box<dyn BufRead>> {
        let file = File::open(path)?;
        if path.extension().is_some_and(|ext| ext == "gz") {
            // If the file is gzipped, use GzDecoder
            let decoder = GzDecoder::new(file);
            Ok(Box::new(io::BufReader::new(decoder)))
        } else {
            // Otherwise, read it normally
            Ok(Box::new(io::BufReader::new(file)))
        }
    }

    /// This function turns a fastq file into kmers and stores them in a bucket.
    /// Quality scores are used (like for BAM files).
    pub fn read_fastq_file(&self, file_path: &str) -> Result<BucketList> {
        let path = Path::new(file_path);
        let sample_name = BucketList::file_path_to_sample_name(path)?;
        let mut out_bucket: DataBucket<Kmer<KmerBits>> = DataBucket::new(
            self.max_bucket_size,
            &self.bucket_dir,
            &sample_name,
            "from_fasta",
        );

        let mut reader = Self::get_uncompressed_or_gz_reader(path)?;
        // let reader = io::BufReader::new(file);
        let mut read_number = 0;

        // let mut lines = Vec::new();
        // for line in reader.lines() {
        let mut line = String::new();
        loop {
            // lines.push(line?);
            // if lines.len() != 4 {
            //     continue;
            // }

            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            if !line.starts_with('@') {
                return Err(anyhow!("Expected line starting with '@'"));
            }

            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            let sequence = line.trim().as_bytes().to_owned();

            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            if line.trim() != "+" {
                return Err(anyhow!("Expected line starting with '+': {line}"));
            }

            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            let qualities = line.trim().as_bytes();

            let kmers: Vec<KmerBits> =
                Kmer::kmers_from_record(&sequence, qualities, self.min_base_quality);
            for kmer in kmers {
                out_bucket.add(Kmer::new(kmer));
            }
            read_number += 1;
            // lines.clear();
        }

        // Write final bucket to disk
        let filenames = out_bucket.finish()?;

        // Create metadata to return
        let bucket_list = BucketList::new(sample_name, filenames, read_number);
        Ok(bucket_list)
    }

    pub fn read_bam_file(&self, file_path: &str) -> Result<BucketList> {
        let file_path = Path::new(file_path);
        let sample_name = BucketList::file_path_to_sample_name(file_path)?;
        let mut reader = bam::BamReader::from_path(file_path, 4)?;
        let mut record = bam::Record::new();
        let mut out_bucket = DataBucket::new(
            self.max_bucket_size,
            &self.bucket_dir,
            &sample_name,
            "pairs",
        );
        let mut read_number: ReadId = 0;

        while reader.read_into(&mut record)? {
            // Generate and process kmers
            let sequence = record.sequence().to_vec();
            let qualities = record.qualities().raw();
            let kmers: Vec<KmerBits> =
                Kmer::kmers_from_record(&sequence, qualities, self.min_base_quality);
            for kmer in kmers {
                out_bucket.add(KmerRead::new(Kmer::new(kmer), read_number));
            }
            read_number += 1;
        }

        // Write final bucket to disk
        let filenames = out_bucket.finish()?;

        // Create metadata to return
        let bucket_list = BucketList::new(sample_name, filenames, read_number);
        Ok(bucket_list)
    }

    pub fn read_bam_file_kmers(&self, file_path: &str) -> Result<BucketList> {
        let file_path = Path::new(file_path);
        let sample_name = BucketList::file_path_to_sample_name(file_path)?;
        let mut reader = bam::BamReader::from_path(file_path, 4)?;
        let mut record = bam::Record::new();
        let mut out_bucket = DataBucket::new(
            self.max_bucket_size,
            &self.bucket_dir,
            &sample_name,
            "kmers",
        );
        let mut read_number: ReadId = 0;

        while reader.read_into(&mut record)? {
            // Generate and process kmers
            let sequence = record.sequence().to_vec();
            let qualities = record.qualities().raw();
            let kmers: Vec<KmerBits> =
                Kmer::kmers_from_record(&sequence, qualities, self.min_base_quality);
            for kmer in kmers {
                out_bucket.add(Kmer::new(kmer));
            }
            read_number += 1;
        }

        // Write final bucket to disk
        let filenames = out_bucket.finish()?;

        // Create metadata to return
        let bucket_list = BucketList::new(sample_name, filenames, read_number);
        Ok(bucket_list)
    }

    pub fn unique_kmers(
        &self,
        bucket_list: &BucketList,
        max_occurrence: Option<usize>,
    ) -> Result<BucketList> {
        let max_occurrence = max_occurrence.unwrap_or(usize::MAX);
        let sample_name = bucket_list.sample_name().to_string();

        let out_filename = format!("{}/{}.unique_kmers", self.bucket_dir, sample_name);
        let out_file = File::create(&out_filename)?;
        let mut writer = BufWriter::new(out_file);

        let mut unique_kmer_count = 0;
        let mbr = MultiBucketReader::new(bucket_list.filenames()).clean_after_read();
        let mut last_kmer: Kmer<KmerBits> = Kmer::default();
        let mut kmer_count = 0;
        let mut total_kmers = 0;
        mbr.for_each(|kmer: Kmer<KmerBits>| {
            total_kmers += 1;
            if kmer == last_kmer {
                kmer_count += 1;
            } else {
                if kmer_count <= max_occurrence {
                    last_kmer.write(&mut writer).unwrap();
                    unique_kmer_count += 1;
                }
                last_kmer = kmer;
                kmer_count = 1;
            }
        });
        if kmer_count <= max_occurrence {
            last_kmer.write(&mut writer).unwrap();
            unique_kmer_count += 1;
        }

        let filenames = vec![out_filename];
        let bucket_list = BucketList::new(sample_name, filenames, unique_kmer_count);
        Ok(bucket_list)
    }

    fn process_kmer_grouped_reads(
        &self,
        kmer: &Kmer<KmerBits>,
        reads: &mut Vec<ReadId>,
        min_max: &MinMaxReads,
        bucket: &mut DataBucket<ReadPairKmer<KmerBits>>,
    ) {
        // Reads will be sorted already
        reads.dedup();
        if min_max.is_valid(reads.len()) {
            for read2_pos in 1..reads.len() {
                for read1_pos in 0..read2_pos {
                    let rpk: ReadPairKmer<KmerBits> =
                        ReadPairKmer::new(reads[read1_pos], reads[read2_pos], kmer);
                    bucket.add(rpk);
                }
            }
        }
        reads.clear();
    }

    pub fn process_read_kmer_buckets(
        &self,
        bucket_list: &BucketList,
        min_max: &MinMaxReads,
    ) -> Result<(BucketList, HashMap<usize, usize>)> {
        let sample_name = bucket_list.sample_name().to_string();
        let mut out_bucket = DataBucket::new(
            self.max_bucket_size,
            &self.bucket_dir,
            &sample_name,
            "read_pairs",
        );
        let mut stats = HashMap::new();
        let mut last_kmer: Kmer<KmerBits> = Kmer::default();
        let mut last_reads_ids = Vec::new();

        let mbr = MultiBucketReader::new(bucket_list.filenames());
        mbr.for_each(|kmer_read: KmerRead<KmerBits>| {
            // Flush reads if new kmer
            if last_kmer != *kmer_read.kmer() {
                *stats.entry(last_reads_ids.len()).or_insert(0usize) += 1;
                self.process_kmer_grouped_reads(
                    &last_kmer,
                    &mut last_reads_ids,
                    min_max,
                    &mut out_bucket,
                );
                last_kmer.clone_from(kmer_read.kmer());
            }
            last_reads_ids.push(kmer_read.read_id());
        });

        *stats.entry(last_reads_ids.len()).or_insert(0) += 1;
        self.process_kmer_grouped_reads(&last_kmer, &mut last_reads_ids, min_max, &mut out_bucket);
        stats.remove(&0); // Remove 0-read group from stats

        // Write final bucket to disk
        let filenames = out_bucket.finish()?;

        let bucket_list = BucketList::new(sample_name, filenames, 0);
        Ok((bucket_list, stats))
    }

    pub fn align_reads(&self, bucket_list: &BucketList, bam_file: &str) -> Result<()> {
        let mut reader = bam::BamReader::from_path(bam_file, 4)?;
        let mut record = bam::Record::new();
        let mut read_number: ReadId = 0;
        let mut mbr: MultiBucketReader<ReadPairKmer<KmerBits>> =
            MultiBucketReader::new(bucket_list.filenames());
        let mut rpk = mbr
            .next()
            .expect("align_reads: No read pair kmer in buckets");

        let mut current_position: GenomicPosition = 10000;
        let mut current_scaffold: ScaffoldId = 1; // TODO CONST
        let mut backlog: HashMap<ReadId, Vec<KmerBacklog<KmerBits>>> = HashMap::new();

        while reader.read_into(&mut record)? {
            let current_read_number = read_number;
            read_number += 1;
            if current_read_number < rpk.read1() {
                let _ = backlog.remove(&current_read_number);
                continue;
            }

            let sequence = record.sequence().to_vec();
            let qualities = record.qualities().raw();

            let kmers_with_position: Vec<KmerWithPosition<KmerBits>> =
                Kmer::kmers_from_record_with_position(&sequence, qualities, self.min_base_quality);

            // Process backlog of read pair kmers matched to this read
            if let Some(kmer_backlog) = backlog.remove(&current_read_number) {
                self.process_backlog(
                    kmer_backlog,
                    &kmers_with_position,
                    current_read_number,
                    &mut current_scaffold,
                    &mut current_position,
                    sequence.len() as GenomicPosition,
                );
            } else {
                // No previous reads matches to this one, start a new scaffold
                current_scaffold += 1;
                current_position = 10000; // TODO CONST
            }

            // Add possible connections to backlog
            while current_read_number == rpk.read1() {
                Self::add2backlog(
                    &rpk,
                    &kmers_with_position,
                    current_read_number,
                    current_scaffold,
                    current_position,
                    sequence.len() as GenomicPosition,
                    &mut backlog,
                );

                // Next read-pair-kmer
                rpk = match mbr.next() {
                    Some(rpk) => rpk,
                    None => return Ok(()),
                }
            }
        }
        Ok(())
    }

    fn process_backlog(
        &self,
        backlog: Vec<KmerBacklog<KmerBits>>,
        kmers_with_position: &[KmerWithPosition<KmerBits>],
        current_read_number: ReadId,
        _current_scaffold: &mut ScaffoldId,
        _current_position: &mut GenomicPosition,
        read_length: GenomicPosition,
    ) {
        if current_read_number != 4 {
            return; // TODO TESTING FIXME
        }
        // println!("{backlog:?}");
        for kwps in backlog {
            let matches: Vec<_> = kmers_with_position
                .iter()
                .filter(|kwp| kwp.kmer() == kwps.kwps().kmer())
                .collect();
            for m in matches {
                let position_read1 = kwps.kwps().genomic_position(kwps.read1_length());
                let position_read2 = m.genomic_position(read_length);
                let new_position = kwps.read1_position() + position_read1 - position_read2;
                println!(
                    "{}: {position_read1} / {position_read2} => {new_position}",
                    kwps.read1_id()
                );
                // println!("{:?} {m:?}", &kwps);
            }
        }
    }

    fn add2backlog(
        rpk: &ReadPairKmer<KmerBits>,
        kmers_with_position: &[KmerWithPosition<KmerBits>],
        current_read_number: ReadId,
        current_scaffold: ScaffoldId,
        current_position: GenomicPosition,
        read_length: GenomicPosition,
        backlog: &mut HashMap<ReadId, Vec<KmerBacklog<KmerBits>>>,
    ) {
        let rpk_kmer = rpk.kmer().bits();

        // TODO find position of kmer in read1
        let kwps: Vec<_> = kmers_with_position
            .iter()
            .filter(|k| k.kmer() == rpk_kmer)
            .collect();

        // Paranoia
        if kwps.is_empty() {
            panic!("Did not find kmer again");
        }

        // Only process single matching kmer
        if kwps.len() == 1 {
            let log = KmerBacklog::new(
                current_read_number,
                current_scaffold,
                current_position,
                read_length,
                kwps[0].to_owned(),
            );
            backlog.entry(rpk.read2()).or_default().push(log);
        }
    }

    #[inline(always)]
    fn store_sequence_kmers(&self, sequence: &String, out_bucket: &mut DataBucket<Kmer<KmerBits>>) {
        let kmers: Vec<KmerBits> = Kmer::kmers_from_sequence(sequence.as_bytes());
        for kmer in kmers {
            out_bucket.add(Kmer::new(kmer));
        }
    }
}
