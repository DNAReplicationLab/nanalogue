#![deny(
    clippy::cast_possible_truncation,
    clippy::allow_attributes,
    clippy::allow_attributes_without_reason,
    missing_copy_implementations,
    missing_debug_implementations,
    missing_docs,
    trivial_casts,
    trivial_numeric_casts,
    unused_extern_crates,
    unused_import_braces,
    unused_qualifications,
    unused_results
)]
//! # Nanalogue Core (Nanalogue = Nucleic Acid Analogue)
//!
//! We process and calculate data associated with DNA molecules, their alignments to
//! reference genomes, modification information on them, and other miscellaneous
//! information.
//!
//! This file directly contains the functions associated with opening BAM files
//! and parsing modification information directly from BAM files. Other functions
//! in the crate are distributed over other files.

use bedrs::{Bed3, Coordinates};
use bio::alphabets::dna::revcomp;
use bio_types::sequence::SequenceRead;
use fibertools_rs::utils::basemods::{BaseMod, BaseMods};
use fibertools_rs::utils::bio_io::{convert_seq_uppercase, get_u8_tag};
use lazy_static::lazy_static;
use rand::random;
use regex::Regex;
use rust_htslib::{bam, bam::Read, bam::ext::BamRecordExtensions, bam::record::Aux, tpool};
use std::collections::HashSet;
use std::convert::TryFrom;
use std::fs::File;
use std::io::{BufRead, BufReader};

// Declare the modules.
pub mod analysis;
pub mod cli;
pub mod error;
pub mod file_utils;
pub mod read_utils;
pub mod simulate_mod_bam;
pub mod subcommands;
pub mod utils;

// Re-exports
pub use cli::{
    InputBam, InputModOptions, InputMods, InputRegionOptions, InputWindowing, OptionalTag,
    RequiredTag,
};
pub use error::Error;
pub use file_utils::{nanalogue_bam_reader, write_bam_denovo, write_fasta};
pub use read_utils::CurrRead;
pub use subcommands::{find_modified_reads, read_info, read_stats, reads_table, window_reads};
pub use utils::{
    Contains, DNARestrictive, F32AbsValAtMost1, F32Bw0and1, FilterByRefCoords, GenomicRegion,
    GetDNARestrictive, Intersects, ModChar, OrdPair, ReadState, ReadStates,
    RestrictModCalledStrand, ThresholdState, is_valid_dna_restrictive,
};

/// Extracts mod information from BAM record to the Fibertools-rs BaseMods Struct.
/// We are copying and modifying code from the fibertools-rs repository
/// (<https://github.com/fiberseq/fibertools-rs>).
///
/// Following is an example of reading and parsing modification data, with
/// no filters applied.
/// We are using an example mod BAM file which has very short reads and very few
/// modified positions, and just examining two reads in it below.
/// ```
/// use nanalogue_core::{Error, nanalogue_bam_reader, nanalogue_mm_ml_parser};
/// use rust_htslib::bam::Read;
/// use fibertools_rs::utils::basemods::{BaseMods, BaseMod};
/// use fibertools_rs::utils::bamranges::Ranges;
/// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
/// let mut count = 0;
/// for record in reader.records(){
///     let r = record?;
///     let Ok(BaseMods{base_mods: v}) = nanalogue_mm_ml_parser(&r, |&_| true, |&_| true, |&_, &_, &_| true, 0) else { unreachable!() };
///     match count {
///     0 => assert_eq!(v, vec![BaseMod{
///             modified_base: b'T',
///             strand: '+',
///             modification_type: 'T',
///             ranges: Ranges {
///                 starts: vec![Some(0), Some(3), Some(4), Some(7)],
///                 ends: vec![Some(1), Some(4), Some(5), Some(8)],
///                 lengths: vec![Some(1); 4],
///                 qual: vec![4, 7, 9, 6],
///                 reference_starts: vec![Some(9), Some(12), Some(13), Some(16)],
///                 reference_ends: vec![Some(9), Some(12), Some(13), Some(16)],
///                 reference_lengths: vec![Some(0); 4],
///                 seq_len: 8,
///                 reverse: false,
///             },
///             record_is_reverse: false,
///     }]),
///     2 => assert_eq!(v, vec![BaseMod{
///             modified_base: b'T',
///             strand: '+',
///             modification_type: 'T',
///             ranges: Ranges {
///                 starts: vec![Some(12), Some(13), Some(16), Some(19), Some(20)],
///                 ends: vec![Some(13), Some(14), Some(17), Some(20), Some(21)],
///                 lengths: vec![Some(1); 5],
///                 qual: vec![3, 3, 4, 3, 182],
///                 reference_starts: vec![Some(15), Some(16), Some(19), Some(22), Some(23)],
///                 reference_ends: vec![Some(15), Some(16), Some(19), Some(22), Some(23)],
///                 reference_lengths: vec![Some(0); 5],
///                 seq_len: 33,
///                 reverse: true,
///             },
///             record_is_reverse: true,
///     }]),
///     _ => {},
///     }
///     count = count + 1;
/// }
/// # Ok::<(), Error>(())
/// ```
pub fn nanalogue_mm_ml_parser<F, G, H>(
    record: &bam::Record,
    filter_mod_prob: F,
    filter_mod_pos: G,
    filter_mod_base_strand_tag: H,
    min_qual: u8,
) -> Result<BaseMods, Error>
where
    F: Fn(&u8) -> bool,
    G: Fn(&usize) -> bool,
    H: Fn(&u8, &char, &ModChar) -> bool,
{
    // regex for matching the MM tag
    lazy_static! {
        // MM:Z:([ACGTUN][-+]([A-Za-z]+|[0-9]+)[.?]?(,[0-9]+)*;)*
        static ref MM_RE: Regex =
            Regex::new(r"((([ACGTUN])([-+])([A-Za-z]+|[0-9]+)([.?]?))((,[0-9]+)*;)*)").unwrap();
    }
    // Array to store all the different modifications within the MM tag
    let mut rtn = vec![];

    let ml_tag = get_u8_tag(record, b"ML");

    let mut num_mods_seen = 0;

    // if there is an MM tag iterate over all the regex matches
    if let Ok(Aux::String(mm_text)) = record.aux(b"MM") {
        for cap in MM_RE.captures_iter(mm_text) {
            let mod_base = cap.get(3).map(|m| m.as_str().as_bytes()[0]).unwrap();
            let mod_strand = cap
                .get(4)
                .map_or("", |m| m.as_str())
                .chars()
                .next()
                .unwrap();

            // get modification type
            let modification_type: ModChar = cap
                .get(5)
                .map_or("", |m| m.as_str())
                .parse()
                .expect("error");

            let is_implicit = match cap.get(6).map_or("", |m| m.as_str()).as_bytes() {
                b"" => true,
                b"." => true,
                b"?" => false,
                _ => unreachable!(),
            };
            let mod_dists_str = cap.get(7).map_or("", |m| m.as_str());
            // parse the string containing distances between modifications into a vector of i64

            let mod_dists: Vec<usize> = mod_dists_str
                .trim_end_matches(';')
                .split(',')
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .map(|s| s.parse().unwrap())
                .collect();

            // base qualities; must reverse if rev comp.
            // NOTE this is always equal to number of bases in the sequence, otherwise
            // rust_htslib will throw an error, so we don't have to check this.
            let base_qual: Vec<u8> = match (min_qual, record.is_reverse()) {
                (0, _) => Vec::new(),
                (_, true) => record.qual().iter().rev().cloned().collect(),
                (_, false) => record.qual().to_vec(),
            };

            // get forward sequence bases from the bam record
            let forward_bases = if record.is_reverse() {
                revcomp(convert_seq_uppercase(record.seq().as_bytes()))
            } else {
                convert_seq_uppercase(record.seq().as_bytes())
            };

            // do we include bases with zero probabilities?
            let is_include_zero_prob = filter_mod_prob(&0);

            // find real positions in the forward sequence
            let mut cur_mod_idx: usize = 0;
            let mut dist_from_last_mod_base: usize = 0;

            // declare vectors with an approximate with_capacity
            let mod_data_len_approx = if is_implicit {
                // In implicit mode, there may be any number of bases
                // after the MM data is over, which must be assumed as unmodified.
                // So we cannot know the exact length of the data before actually
                // parsing it, and this is just a lower bound of the length.
                mod_dists.len() + mod_dists.iter().sum::<usize>()
            } else {
                mod_dists.len()
            };
            let mut modified_positions: Vec<i64> = Vec::with_capacity(mod_data_len_approx);
            let mut modified_probabilities: Vec<u8> = Vec::with_capacity(mod_data_len_approx);

            for (cur_seq_idx, &_) in forward_bases
                .iter()
                .enumerate()
                .filter(|&(_, &k)| mod_base == b'N' || k == mod_base)
            {
                if cur_mod_idx < mod_dists.len()
                    && dist_from_last_mod_base == mod_dists[cur_mod_idx]
                {
                    let prob = ml_tag
                        .get(cur_mod_idx + num_mods_seen)
                        .ok_or(Error::InvalidModProbs)?;
                    if filter_mod_prob(prob)
                        && filter_mod_pos(&cur_seq_idx)
                        && !(min_qual > 0 && base_qual[cur_seq_idx] < min_qual)
                    {
                        modified_positions.push(i64::try_from(cur_seq_idx).unwrap());
                        modified_probabilities.push(*prob);
                    }
                    dist_from_last_mod_base = 0;
                    cur_mod_idx += 1;
                } else {
                    if is_include_zero_prob
                        && is_implicit
                        && filter_mod_pos(&cur_seq_idx)
                        && !(min_qual > 0 && base_qual[cur_seq_idx] < min_qual)
                    {
                        modified_positions.push(i64::try_from(cur_seq_idx).unwrap());
                        modified_probabilities.push(0);
                    }
                    dist_from_last_mod_base += 1
                }
            }

            if cur_mod_idx == mod_dists.len() {
                num_mods_seen += cur_mod_idx;
            } else {
                return Err(Error::InvalidModCoords);
            }

            // if data matches filters, add to struct.
            if filter_mod_base_strand_tag(&mod_base, &mod_strand, &modification_type) {
                let mods = BaseMod::new(
                    record,
                    mod_base,
                    mod_strand,
                    modification_type.val(),
                    modified_positions,
                    modified_probabilities,
                );
                rtn.push(mods);
            }
        }
    }

    if num_mods_seen != ml_tag.len() {
        return Err(Error::InvalidModProbs);
    }

    // needed so I can compare methods
    rtn.sort();

    // Check for duplicate strand, modification_type combinations
    let mut seen_combinations = HashSet::new();
    for base_mod in &rtn {
        let combination = (base_mod.strand, base_mod.modification_type);
        if seen_combinations.contains(&combination) {
            return Err(Error::InvalidDuplicates(format!(
                "Duplicate strand '{}' and modification_type '{}' combination found",
                base_mod.strand, base_mod.modification_type
            )));
        }
        let _ = seen_combinations.insert(combination);
    }

    Ok(BaseMods { base_mods: rtn })
}

/// A global struct which contains BAM records for further usage.
/// NOTE: we don't derive many traits here as the RcRecords object
/// does not have many traits.
///
/// ```
/// use nanalogue_core::{Error, nanalogue_bam_reader, BamRcRecords, InputBam, InputMods};
/// use rust_htslib::bam::Read;
/// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
/// let BamRcRecords = BamRcRecords::new(&mut reader, &mut InputBam::default(),
///     &mut InputMods::default())?;
/// assert_eq!(BamRcRecords.header.tid(b"dummyI"), Some(0));
/// assert_eq!(BamRcRecords.header.tid(b"dummyII"), Some(1));
/// assert_eq!(BamRcRecords.header.tid(b"dummyIII"), Some(2));
/// # Ok::<(), Error>(())
/// ```
#[derive(Debug)]
pub struct BamRcRecords<'a> {
    /// RcRecords object output by rust htslib which we can iterate over
    pub rc_records: bam::RcRecords<'a, bam::Reader>,
    /// Header of the bam file
    pub header: bam::HeaderView,
}

impl<'a> BamRcRecords<'a> {
    /// Extracts RcRecords from a BAM Reader
    pub fn new(
        bam_reader: &'a mut bam::Reader,
        bam_opts: &mut InputBam,
        mod_opts: &mut impl InputRegionOptions,
    ) -> Result<Self, Error> {
        let header = bam_reader.header().clone();
        let tp = tpool::ThreadPool::new(bam_opts.threads.get())?;
        bam_reader.set_thread_pool(&tp)?;

        // Load read ID list if specified
        if let Some(file_path) = &bam_opts.read_id_list {
            let file = File::open(file_path).map_err(Error::InputOutputError)?;
            let reader = BufReader::new(file);
            let mut read_ids = HashSet::new();

            for line in reader.lines() {
                let line = line.map_err(Error::InputOutputError)?;
                let line = line.trim();
                if !line.is_empty() && !line.starts_with('#') {
                    let _ = read_ids.insert(line.to_string());
                }
            }
            bam_opts.read_id_set = Some(read_ids);
        }

        let rc_records = bam_reader.rc_records();
        if !(bam_opts.convert_region_to_bed3(header.clone())?) {
            return Err(Error::UnknownError);
        }
        if !(mod_opts.convert_region_to_bed3(header.clone())?) {
            return Err(Error::UnknownError);
        }
        Ok(BamRcRecords { rc_records, header })
    }
}

/// Trait that performs filtration
pub trait BamPreFilt {
    /// apply default filtration
    fn pre_filt(&self, _bam_opts: &InputBam) -> bool {
        todo!()
    }
    /// filtration by length
    fn filt_by_len(&self, _min_seq_len: u64, _include_zero_len: bool) -> bool {
        todo!()
    }
    /// filtration by alignment length
    fn filt_by_align_len(&self, _min_align_len: i64) -> bool {
        todo!()
    }
    /// filtration by read id
    fn filt_by_read_id(&self, _read_id: &str) -> bool {
        todo!()
    }
    /// filtration by read id set
    fn filt_by_read_id_set(&self, _read_id_set: &HashSet<String>) -> bool {
        todo!()
    }
    /// filtration using flags
    fn filt_by_bitwise_or_flags(&self, _states: &ReadStates) -> bool {
        todo!()
    }
    /// random filtration
    fn filt_random_subset(&self, _fraction: F32Bw0and1) -> bool {
        todo!()
    }
    /// filtration by mapq
    fn filt_by_mapq(&self, _min_mapq: u8, _exclude_mapq_unavail: bool) -> bool {
        todo!()
    }
    /// filtration by region
    fn filt_by_region(&self, _region: &Bed3<i32, u64>, _full_region: bool) -> bool {
        todo!()
    }
}

/// Trait that performs filtration on rust_htslib Record
///
/// Filter in action below, 100 bp read passed with a 20 bp filter
/// but fails with a 120 bp filter
/// ```
/// use nanalogue_core::{BamPreFilt, InputBam, Error};
/// use rust_htslib::bam::record;
/// use rust_htslib::bam::record::{Cigar, CigarString};
/// let mut bam_record = record::Record::new();
/// bam_record.set(&vec![b'r',b'e',b'a',b'd'], Some(&CigarString(vec![Cigar::Match(100)])),
///       &vec![ b'A' as u8; 100], &vec![255 as u8; 100]);
/// let mut bam_opts = InputBam::default();
/// bam_opts.min_seq_len = 20;
/// assert_eq!(bam_record.pre_filt(&bam_opts), true);
/// bam_opts.min_seq_len = 120;
/// assert_eq!(bam_record.pre_filt(&bam_opts), false);
/// # Ok::<(), Error>(())
/// ```
/// By default, zero-length records are excluded to avoid processing errors.
/// ```
/// # use nanalogue_core::{BamPreFilt, InputBam, Error};
/// # use rust_htslib::bam::record;
/// # use rust_htslib::bam::record::{Cigar, CigarString};
/// let mut bam_record = record::Record::new();
/// let mut bam_opts = InputBam::default();
/// bam_opts.min_seq_len = 20;
/// bam_opts.include_zero_len = false;  // default behavior
/// assert_eq!(bam_record.pre_filt(&bam_opts), false);  // excluded
/// # Ok::<(), Error>(())
/// ```
/// Zero-length records can be included if explicitly requested.
/// ```
/// # use nanalogue_core::{BamPreFilt, InputBam, Error};
/// # use rust_htslib::bam::record;
/// # use rust_htslib::bam::record::{Cigar, CigarString};
/// let mut bam_record = record::Record::new();
/// let mut bam_opts = InputBam::default();
/// bam_opts.min_seq_len = 20;
/// bam_opts.include_zero_len = true;
/// assert_eq!(bam_record.pre_filt(&bam_opts), true);  // included
/// # Ok::<(), Error>(())
/// ```
impl BamPreFilt for bam::Record {
    /// apply default filtration by read length
    fn pre_filt(&self, bam_opts: &InputBam) -> bool {
        self.filt_by_len(bam_opts.min_seq_len, bam_opts.include_zero_len)
            & self.filt_random_subset(bam_opts.sample_fraction)
            & self.filt_by_mapq(bam_opts.mapq_filter, bam_opts.exclude_mapq_unavail)
            & {
                if let Some(v) = &bam_opts.read_id {
                    self.filt_by_read_id(v.as_str())
                } else if let Some(read_id_set) = &bam_opts.read_id_set {
                    self.filt_by_read_id_set(read_id_set)
                } else {
                    true
                }
            }
            & {
                if let Some(v) = bam_opts.min_align_len {
                    self.filt_by_align_len(v)
                } else {
                    true
                }
            }
            & {
                if let Some(v) = &bam_opts.read_filter {
                    self.filt_by_bitwise_or_flags(v)
                } else {
                    true
                }
            }
            & {
                if let Some(v) = bam_opts.region_filter() {
                    self.filt_by_region(v, bam_opts.is_full_overlap())
                } else {
                    true
                }
            }
    }
    /// filtration by read length
    fn filt_by_len(&self, min_seq_len: u64, include_zero_len: bool) -> bool {
        // self.len() returns a usize which we convert to u64
        match (min_seq_len, self.len() as u64, include_zero_len) {
            (_, 0, false) => false, // Exclude zero-length by default
            (_, 0, true) => true,   // Include zero-length when explicitly requested
            (l_min, v, _) => v >= l_min,
        }
    }
    /// filtration by alignment length
    fn filt_by_align_len(&self, min_align_len: i64) -> bool {
        !self.is_unmapped() && (self.reference_end() - self.pos() >= min_align_len)
    }
    /// filtration by read id
    fn filt_by_read_id(&self, read_id: &str) -> bool {
        read_id.as_bytes() == self.name()
    }
    /// filtration by read id set
    fn filt_by_read_id_set(&self, read_id_set: &HashSet<String>) -> bool {
        if let Ok(name_str) = std::str::from_utf8(self.name()) {
            read_id_set.contains(name_str)
        } else {
            false
        }
    }
    /// filtration by flag list
    fn filt_by_bitwise_or_flags(&self, states: &ReadStates) -> bool {
        states.bam_flags().contains(&self.flags())
    }
    /// random filtration
    fn filt_random_subset(&self, fraction: F32Bw0and1) -> bool {
        match fraction.val() {
            1.0 => true,
            0.0 => false,
            v => {
                let random_number: f32 = random();
                random_number < v
            }
        }
    }
    /// filtration by mapq
    fn filt_by_mapq(&self, min_mapq: u8, exclude_mapq_unavail: bool) -> bool {
        !(exclude_mapq_unavail && self.mapq() == 255) && self.mapq() >= min_mapq
    }
    /// filtration by region
    fn filt_by_region(&self, region: &Bed3<i32, u64>, full_region: bool) -> bool {
        !self.is_unmapped() && (self.tid() == *region.chr()) && {
            let region_start = region.start();
            let region_end = region.end();
            (region_start == u64::MIN && region_end == u64::MAX) || {
                let start: u64 = self.pos().try_into().expect("no error");
                let end: u64 = self.reference_end().try_into().expect("no error");
                if full_region {
                    (start..end).contains(&region_start) && (start..(end + 1)).contains(&region_end)
                } else {
                    (start..end).intersects(&(region_start..region_end))
                }
            }
        }
    }
}

#[cfg(test)]
mod mod_parse_tests {
    use super::*;
    use fibertools_rs::utils::bamranges::Ranges;
    use fibertools_rs::utils::basemods::{BaseMod, BaseMods};
    use rust_htslib::bam::Read;

    /// Tests if Mod BAM modification parsing is alright,
    /// some of the test cases here may be a repeat of the doctest above.
    #[test]
    fn test_mod_bam_parsing_from_example_1_bam() -> Result<(), Error> {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        for (count, record) in reader.records().enumerate() {
            let r = record?;
            let BaseMods { base_mods: v } =
                nanalogue_mm_ml_parser(&r, |&_| true, |&_| true, |&_, &_, &_| true, 0).unwrap();
            match count {
                0 => assert_eq!(
                    v,
                    vec![BaseMod {
                        modified_base: b'T',
                        strand: '+',
                        modification_type: 'T',
                        ranges: Ranges {
                            starts: vec![Some(0), Some(3), Some(4), Some(7)],
                            ends: vec![Some(1), Some(4), Some(5), Some(8)],
                            lengths: vec![Some(1); 4],
                            qual: vec![4, 7, 9, 6],
                            reference_starts: vec![Some(9), Some(12), Some(13), Some(16)],
                            reference_ends: vec![Some(9), Some(12), Some(13), Some(16)],
                            reference_lengths: vec![Some(0); 4],
                            seq_len: 8,
                            reverse: false,
                        },
                        record_is_reverse: false,
                    }]
                ),
                1 => assert_eq!(
                    v,
                    vec![BaseMod {
                        modified_base: b'T',
                        strand: '+',
                        modification_type: 'T',
                        ranges: Ranges {
                            starts: vec![Some(3), Some(8), Some(27), Some(39), Some(47)],
                            ends: vec![Some(4), Some(9), Some(28), Some(40), Some(48)],
                            lengths: vec![Some(1); 5],
                            qual: vec![221, 242, 3, 47, 239],
                            reference_starts: vec![
                                Some(26),
                                Some(31),
                                Some(50),
                                Some(62),
                                Some(70)
                            ],
                            reference_ends: vec![Some(26), Some(31), Some(50), Some(62), Some(70)],
                            reference_lengths: vec![Some(0); 5],
                            seq_len: 48,
                            reverse: false,
                        },
                        record_is_reverse: false,
                    }]
                ),
                2 => assert_eq!(
                    v,
                    vec![BaseMod {
                        modified_base: b'T',
                        strand: '+',
                        modification_type: 'T',
                        ranges: Ranges {
                            starts: vec![Some(12), Some(13), Some(16), Some(19), Some(20)],
                            ends: vec![Some(13), Some(14), Some(17), Some(20), Some(21)],
                            lengths: vec![Some(1); 5],
                            qual: vec![3, 3, 4, 3, 182],
                            reference_starts: vec![
                                Some(15),
                                Some(16),
                                Some(19),
                                Some(22),
                                Some(23)
                            ],
                            reference_ends: vec![Some(15), Some(16), Some(19), Some(22), Some(23)],
                            reference_lengths: vec![Some(0); 5],
                            seq_len: 33,
                            reverse: true,
                        },
                        record_is_reverse: true,
                    }]
                ),
                3 => assert_eq!(
                    v,
                    vec![
                        BaseMod {
                            modified_base: b'G',
                            strand: '-',
                            modification_type: '\u{1C20}',
                            ranges: Ranges {
                                starts: vec![
                                    Some(28),
                                    Some(29),
                                    Some(30),
                                    Some(32),
                                    Some(43),
                                    Some(44)
                                ],
                                ends: vec![
                                    Some(29),
                                    Some(30),
                                    Some(31),
                                    Some(33),
                                    Some(44),
                                    Some(45)
                                ],
                                lengths: vec![Some(1); 6],
                                qual: vec![0, 0, 0, 0, 77, 0],
                                reference_starts: vec![None; 6],
                                reference_ends: vec![None; 6],
                                reference_lengths: vec![None; 6],
                                seq_len: 48,
                                reverse: false,
                            },
                            record_is_reverse: false,
                        },
                        BaseMod {
                            modified_base: b'T',
                            strand: '+',
                            modification_type: 'T',
                            ranges: Ranges {
                                starts: vec![Some(3), Some(8), Some(27), Some(39), Some(47)],
                                ends: vec![Some(4), Some(9), Some(28), Some(40), Some(48)],
                                lengths: vec![Some(1); 5],
                                qual: vec![221, 242, 0, 47, 239],
                                reference_starts: vec![None; 5],
                                reference_ends: vec![None; 5],
                                reference_lengths: vec![None; 5],
                                seq_len: 48,
                                reverse: false,
                            },
                            record_is_reverse: false,
                        }
                    ]
                ),
                _ => {}
            }
        }
        Ok(())
    }

    fn create_test_record_and_parse(mm_value: &str, ml_values: Vec<u8>) -> Result<BaseMods, Error> {
        let mut record = bam::Record::new();

        // Create a sequence - ATCG repeated
        let seq_bytes = b"ATCGATCGATCGATCGATCG";
        let qname = b"test_read";
        let cigar = bam::record::CigarString::from(vec![bam::record::Cigar::Match(20)]);
        let quals = vec![30u8; 20]; // Quality scores of 30 for all bases
        record.set(qname, Some(&cigar), seq_bytes, &quals);

        // Set MM tag
        record.push_aux(b"MM", Aux::String(mm_value))?;

        // Set ML tag (modification probabilities)
        record.push_aux(b"ML", Aux::ArrayU8((&ml_values).into()))?;

        // Call the parser and return result
        nanalogue_mm_ml_parser(
            &record,
            |&_| true,         // Accept all probabilities
            |&_| true,         // Accept all positions
            |&_, &_, &_| true, // Accept all base/strand/tag combinations
            0,                 // No quality threshold
        )
    }

    #[test]
    #[should_panic(expected = "InvalidDuplicates")]
    fn test_nanalogue_mm_ml_parser_detects_duplicates() {
        // Two T+ modifications with same strand and modification_type
        let mm_value = "T+T,0,3;T+T,1,2;";
        let ml_values = Vec::from([100u8, 200u8, 150u8, 180u8]);
        let _: BaseMods = create_test_record_and_parse(mm_value, ml_values).unwrap();
    }

    #[test]
    fn test_nanalogue_mm_ml_parser_accepts_unique_combinations() {
        // T+, C+, T- (all different combinations)
        let mm_value = "T+T,0,3;C+m,0,1;T-T,1,2;";
        let ml_values = Vec::from([100u8, 200u8, 150u8, 180u8, 120u8, 140u8]);
        let _: BaseMods = create_test_record_and_parse(mm_value, ml_values).unwrap();
    }

    #[test]
    #[should_panic(expected = "InvalidModCoords")]
    fn test_nanalogue_mm_ml_parser_detects_invalid_mod_coords_1() {
        // Invalid coordinates: seq does not have 11 As (4 + 1 + 5 + 1)
        let mm_value = "T+T,0,3;A-T,4,5;";
        let ml_values = Vec::from([100u8, 200u8, 150u8, 180u8]);
        let _: BaseMods = create_test_record_and_parse(mm_value, ml_values).unwrap();
    }

    #[test]
    #[should_panic(expected = "InvalidModCoords")]
    fn test_nanalogue_mm_ml_parser_detects_invalid_mod_coords_2() {
        // Invalid coords: there aren't 11 C in the sequence
        let mm_value = "T+T,0,3;C+m,4,5;T-T,1,2;";
        let ml_values = Vec::from([100u8, 200u8, 150u8, 180u8, 120u8, 140u8]);
        let _: BaseMods = create_test_record_and_parse(mm_value, ml_values).unwrap();
    }

    #[test]
    #[should_panic(expected = "InvalidModProbs")]
    fn test_nanalogue_mm_ml_parser_detects_ml_tag_longer_than_mm_data() {
        // ML tag has more values than modifications in MM tag
        let mm_value = "T+T,0,3;"; // 2 modifications
        let ml_values = Vec::from([100u8, 200u8, 150u8, 180u8]); // 4 values - too many!
        let _: BaseMods = create_test_record_and_parse(mm_value, ml_values).unwrap();
    }

    #[test]
    #[should_panic(expected = "InvalidModProbs")]
    fn test_nanalogue_mm_ml_parser_detects_ml_tag_shorter_than_mm_data() {
        // ML tag has fewer values than modifications in MM tag
        let mm_value = "T+T,0,3;"; // 2 modifications
        let ml_values = Vec::from([100u8]); // 1 value - too few!
        let _: BaseMods = create_test_record_and_parse(mm_value, ml_values).unwrap();
    }
}

#[cfg(test)]
mod zero_length_filtering_tests {
    use super::*;
    use rust_htslib::bam::Read;

    #[test]
    fn test_zero_length_filtering_with_example_2_zero_len_sam() -> Result<(), Error> {
        // Test with include_zero_len = false and min_seq_len = 1 (should get 1 read)
        let mut reader = nanalogue_bam_reader("examples/example_2_zero_len.sam")?;
        let bam_opts_exclude_zero: InputBam =
            serde_json::from_str(r#"{"min_seq_len": 1}"#).unwrap();

        let mut count_exclude_zero = 0;
        for record_result in reader.records() {
            let record = record_result?;
            if record.pre_filt(&bam_opts_exclude_zero) {
                count_exclude_zero += 1;
            }
        }

        // Test with include_zero_len = true and min_seq_len = 1 (should get 2 reads)
        let mut reader = nanalogue_bam_reader("examples/example_2_zero_len.sam")?;
        let bam_opts_include_zero: InputBam =
            serde_json::from_str(r#"{"min_seq_len": 1, "include_zero_len": true}"#).unwrap();

        let mut count_include_zero = 0;
        for record_result in reader.records() {
            let record = record_result?;
            if record.pre_filt(&bam_opts_include_zero) {
                count_include_zero += 1;
            }
        }

        assert_eq!(count_exclude_zero, 1);
        assert_eq!(count_include_zero, 2);

        Ok(())
    }
}

#[cfg(test)]
mod invalid_seq_length_tests {
    use super::*;

    #[test]
    fn test_invalid_seq_length_error_with_example_4() {
        let mut reader =
            nanalogue_bam_reader("examples/example_4_invalid_basequal_len.sam").unwrap();

        let mut cnt = 0;
        for record_result in reader.records() {
            cnt = cnt + 1;
            assert!(record_result.is_err());
        }
        assert_eq!(cnt, 1);
    }
}

#[cfg(test)]
mod base_qual_filtering_tests {
    use super::*;
    use fibertools_rs::utils::basemods::BaseMods;
    use rust_htslib::bam::Read;

    #[test]
    fn test_base_qual_filtering_with_example_5_valid_basequal_sam() -> Result<(), Error> {
        let mut reader = nanalogue_bam_reader("examples/example_5_valid_basequal.sam")?;

        for record_result in reader.records() {
            let record = record_result?;
            let BaseMods { base_mods: v } =
                nanalogue_mm_ml_parser(&record, |&_| true, |&_| true, |&_, &_, &_| true, 60)
                    .unwrap();

            assert_eq!(v.len(), 1);
            let base_mod = &v[0];
            assert_eq!(base_mod.modified_base, b'T');
            assert_eq!(base_mod.strand, '+');
            assert_eq!(base_mod.modification_type, 'T');
            assert_eq!(base_mod.ranges.qual, vec![7, 9]);
            assert_eq!(base_mod.ranges.starts, vec![Some(3), Some(4)]);
            assert_eq!(base_mod.ranges.ends, vec![Some(4), Some(5)]);
        }

        Ok(())
    }
}

#[cfg(test)]
mod bam_rc_record_tests {
    use super::*;
    use rand::random_range;
    use rust_htslib::bam::record;
    use rust_htslib::bam::record::{Cigar, CigarString};

    #[test]
    fn test_bam_rc_records() {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam").unwrap();
        let bam_rc_records = BamRcRecords::new(
            &mut reader,
            &mut InputBam::default(),
            &mut InputMods::default(),
        )
        .unwrap();
        assert_eq!(bam_rc_records.header.tid(b"dummyI"), Some(0));
        assert_eq!(bam_rc_records.header.tid(b"dummyII"), Some(1));
        assert_eq!(bam_rc_records.header.tid(b"dummyIII"), Some(2));
    }

    #[test]
    fn test_example_3_read_list_1() {
        // Test with example_3_subset_1 - should see 2 reads
        let mut reader = nanalogue_bam_reader("examples/example_3.bam").unwrap();
        let mut bam_opts: InputBam =
            serde_json::from_str(r#"{"read_id_list": "examples/example_3_subset_1"}"#).unwrap();

        let bam_rc_records =
            BamRcRecords::new(&mut reader, &mut bam_opts, &mut InputMods::default()).unwrap();

        // Count lines in the read ID file
        let file = File::open("examples/example_3_subset_1").unwrap();
        let reader_file = BufReader::new(file);
        let mut line_count = 0;
        for line in reader_file.lines() {
            let line = line.unwrap();
            let line = line.trim();
            if !line.is_empty() && !line.starts_with('#') {
                line_count += 1;
            }
        }
        assert_eq!(line_count, 2);

        let mut count = 0;
        for record_result in bam_rc_records.rc_records {
            let record = record_result.unwrap();
            if record.pre_filt(&bam_opts) {
                count += 1;
            }
        }
        assert_eq!(count, 2);
    }

    #[test]
    fn test_example_3_read_list_2() {
        // Test with example_3_subset_w_invalid - should see 2 reads
        // even though file contains 3 read IDs
        let mut reader = nanalogue_bam_reader("examples/example_3.bam").unwrap();
        let mut bam_opts: InputBam =
            serde_json::from_str(r#"{"read_id_list": "examples/example_3_subset_w_invalid"}"#)
                .unwrap();
        let bam_rc_records =
            BamRcRecords::new(&mut reader, &mut bam_opts, &mut InputMods::default()).unwrap();

        // Count lines in the read ID file
        let file = File::open("examples/example_3_subset_w_invalid").unwrap();
        let reader_file = BufReader::new(file);
        let mut line_count = 0;
        for line in reader_file.lines() {
            let line = line.unwrap();
            let line = line.trim();
            if !line.is_empty() && !line.starts_with('#') {
                line_count += 1;
            }
        }
        assert_eq!(line_count, 3);

        let mut count = 0;
        for record_result in bam_rc_records.rc_records {
            let record = record_result.unwrap();
            if record.pre_filt(&bam_opts) {
                count += 1;
            }
        }
        assert_eq!(count, 2);
    }

    #[test]
    fn test_random_retrieval() {
        let mut count_retained = 0;
        let bam_opts: InputBam =
            serde_json::from_str(r#"{"sample_fraction": 0.5, "include_zero_len": true}"#).unwrap();

        for _ in 0..10000 {
            let record = record::Record::new();
            if record.pre_filt(&bam_opts) {
                count_retained += 1;
            }
        }

        // 50% retention rate => 5000 reads, so we test if 4500-5500 reads come through
        // (this is quite lax actually)
        assert!(count_retained >= 4500 && count_retained <= 5500);
    }

    #[test]
    fn test_single_read_id_filtering() {
        let mut count_retained = 0;

        let bam_opts: InputBam =
            serde_json::from_str(r#"{"read_id": "read2", "include_zero_len": true}"#).unwrap();

        for i in 1..=1000 {
            let mut record = record::Record::new();
            let read_id = format!("read{}", i);
            record.set_qname(read_id.as_bytes());

            if record.pre_filt(&bam_opts) {
                count_retained += 1;
            }
        }

        assert_eq!(count_retained, 1);
    }

    #[test]
    fn test_seq_and_align_len_filtering() {
        let bam_opts_min_len: InputBam = serde_json::from_str(r#"{"min_seq_len": 5000}"#).unwrap();
        let bam_opts_min_align_len: InputBam =
            serde_json::from_str(r#"{"min_align_len": 2500}"#).unwrap();

        let mut count_retained = (0, 0);

        for _ in 1..=10000 {
            let mut record = record::Record::new();
            let seq_len: usize = random_range(2..=10000);
            let match_len = seq_len / 2;
            let hard_clip_len = seq_len - match_len;

            record.set(
                &vec![b'r', b'e', b'a', b'd'],
                Some(&CigarString(vec![
                    Cigar::Match(match_len as u32),
                    Cigar::HardClip(hard_clip_len as u32),
                ])),
                &vec![b'A'; seq_len],
                &vec![50; seq_len],
            );
            record.unset_flags();
            record.set_flags(loop {
                let random_state: ReadState = random();
                match random_state {
                    ReadState::Unmapped => continue,
                    v => break u16::from(v),
                }
            });
            if record.pre_filt(&bam_opts_min_align_len) {
                count_retained.0 += 1;
            }
            if record.pre_filt(&bam_opts_min_len) {
                count_retained.1 += 1;
            }
        }

        assert!(count_retained.0 >= 4500 && count_retained.0 <= 5500);
        assert!(count_retained.1 >= 4500 && count_retained.1 <= 5500);
    }

    #[test]
    fn test_filt_by_bitwise_or_flags() {
        // Create random subset of read states (ensure at least one)
        let selected_states = {
            let mut selected_states = Vec::new();
            let all_states = vec!["0", "4", "16", "256", "272", "2048", "2064"];
            for state in &all_states {
                if random::<f32>() < 0.5 {
                    selected_states.push(*state);
                }
            }
            if selected_states.is_empty() {
                selected_states.push(all_states[random_range(0..all_states.len())]);
            }
            selected_states
        };

        let states_string = selected_states.join(",");
        let bam_opts: InputBam = serde_json::from_str(
            format!("{{\"read_filter\": [{states_string}], \"include_zero_len\": true}}").as_str(),
        )
        .unwrap();

        let mut count_retained = 0;

        for _ in 1..=70000 {
            let mut record = record::Record::new();
            record.unset_flags();
            record.set_flags({
                let random_state: ReadState = random();
                u16::from(random_state)
            });
            if record.pre_filt(&bam_opts) {
                count_retained += 1;
            }
        }

        let expected_count = selected_states.len() * 10000;
        let tolerance = (expected_count as f32 * 0.2) as usize;
        let min_count = expected_count.saturating_sub(tolerance);
        let max_count = expected_count + tolerance;

        assert!(count_retained >= min_count && count_retained <= max_count);
    }

    #[test]
    fn test_filt_by_region() {
        let mut count_retained = (0, 0);

        for _ in 1..=10000 {
            let mut record = record::Record::new();

            // Draw four numbers from 0 to 10000
            let four_nums = [
                random_range(0..=10000),
                random_range(0..=10000),
                random_range(0..=10000),
                random_range(0..=10000u64),
            ];

            // First two for record
            let mut pos_nums = [four_nums[0], four_nums[1]];
            pos_nums.sort();
            let start_pos = pos_nums[0];
            let end_pos = pos_nums[1];

            // Use the last two to set up region, with a random contig chosen
            // from a set of two.
            let region_tid = if random::<bool>() { 0 } else { 1 };
            let mut region_nums = [four_nums[2], four_nums[3]];
            region_nums.sort();
            let region_start = region_nums[0];
            let region_end = region_nums[1];
            let bam_opts_no_full_region: InputBam = serde_json::from_str(
                format!("{{\"region_bed3\": [{region_tid}, {region_start}, {region_end}]}}")
                    .as_str(),
            )
            .unwrap();
            let bam_opts_full_region: InputBam = serde_json::from_str(
                format!("{{\"full_region\": true, \"region_bed3\": [{region_tid}, {region_start}, {region_end}]}}")
                    .as_str(),
            )
            .unwrap();

            // Calculate sequence length (ensure at least 1)
            let seq_len = if end_pos == start_pos {
                1
            } else {
                end_pos - start_pos
            };

            // Set sequence details first
            record.set(
                &vec![b'r', b'e', b'a', b'd'],
                Some(&CigarString(vec![Cigar::Match(seq_len as u32)])),
                &vec![b'A'; seq_len as usize],
                &vec![255; seq_len as usize],
            );

            // Set tid and position, contig chosen from a random set of two.
            let tid = if random::<bool>() { 0 } else { 1 };
            record.set_tid(tid);
            record.set_pos(start_pos as i64);

            // Verify reference_end calculation
            use rust_htslib::bam::ext::BamRecordExtensions;
            let expected_ref_end = start_pos as i64 + seq_len as i64;
            assert_eq!(record.reference_end(), expected_ref_end);

            if record.pre_filt(&bam_opts_no_full_region) {
                count_retained.0 += 1;
            }
            if record.pre_filt(&bam_opts_full_region) {
                count_retained.1 += 1;
            }
        }

        // the chance that two randomly chosen intervals on the interval [0, l]
        // intersect is 2/3, and then we have the added random element of one
        // of two contigs, so the probability is 2/3 * 1/2 = 1/3.
        let expected_count: usize = 10000 / 3; // ~3333
        let tolerance = (expected_count as f32 * 0.2) as usize;
        let min_count = expected_count.saturating_sub(tolerance);
        let max_count = expected_count + tolerance;
        assert!(count_retained.0 >= min_count && count_retained.0 <= max_count);

        // the chance that two randomly chosen intervals are such that one is
        // contained within the other is 1/3, and we are looking for the region
        // being contained completely within the read and NOT the read being contained
        // within the region, so the probability is 1/6. This is one fourth of the
        // 2/3 used above. So the same criterion will work with the second count quadrupled.
        assert!(4 * count_retained.1 >= min_count && 4 * count_retained.1 <= max_count);
    }
}
