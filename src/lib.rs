#![deny(
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
use fibertools_rs::utils::bio_io::get_u8_tag;
use lazy_static::lazy_static;
use regex::Regex;
use rust_htslib::{bam, bam::Read, bam::ext::BamRecordExtensions, bam::record::Aux, tpool};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::convert::TryFrom;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::str::FromStr;

// Declare the modules.
pub mod cli;
pub mod error;
pub mod read_utils;
pub mod subcommands;
pub mod utils;

// Re-exports
pub use cli::{
    InputBam, InputModOptions, InputMods, InputRegionOptions, InputWindowing, OptionalTag,
    RequiredTag,
};
pub use error::Error;
pub use read_utils::CurrRead;
pub use utils::{
    Contains, F32AbsValBelow1, F32Bw0and1, FilterByRefCoords, GenomicRegion, Intersects, ModChar,
    OrdPair, ReadState, RestrictModCalledStrand, ThresholdState,
};

/// Converts DNA bases to uppercase if needed, leaving other characters unchanged.
///
/// ```
/// use nanalogue_core::convert_seq_uppercase;
/// let x = convert_seq_uppercase(vec![b'a',b'b',b'c',b'g',b't',b'n',b'u',b'A',b'C',b'G',b'T',b'N']);
/// assert_eq!(x, vec![b'A',b'b',b'C',b'G',b'T',b'N',b'u',b'A',b'C',b'G',b'T',b'N']);
/// ```
pub fn convert_seq_uppercase(mut seq: Vec<u8>) -> Vec<u8> {
    // convert seq to uppercase, ignoring invalid characters
    for base in &mut seq {
        match *base {
            b'a' => *base = b'A',
            b'c' => *base = b'C',
            b'g' => *base = b'G',
            b't' => *base = b'T',
            b'n' => *base = b'N',
            _ => {}
        }
    }
    seq
}

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
                _ => return Err(Error::InvalidModNotation),
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

            // get forward sequence bases from the bam record
            let forward_bases = if record.is_reverse() {
                revcomp(convert_seq_uppercase(record.seq().as_bytes()))
            } else {
                convert_seq_uppercase(record.seq().as_bytes())
            };
            log::trace!(
                "mod_base: {}, mod_strand: {}, modification_type: {}, mod_dists: {:?}",
                mod_base as char,
                mod_strand,
                modification_type,
                mod_dists
            );

            // do we include bases with zero probabilities?
            let is_include_zero_prob = filter_mod_prob(&0);

            // base qualities
            let base_qual = record.qual();
            if min_qual > 0 && (base_qual.len() != forward_bases.len()) {
                continue;
            }

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
                    let prob = &ml_tag[cur_mod_idx + num_mods_seen];
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

            num_mods_seen += cur_mod_idx;

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

    // needed so I can compare methods
    rtn.sort();
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
            let file = File::open(file_path).map_err(|e| Error::InputOutputError(e))?;
            let reader = BufReader::new(file);
            let mut read_ids = HashSet::new();

            for line in reader.lines() {
                let line = line.map_err(|e| Error::InputOutputError(e))?;
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

/// Opens BAM file, also copied and edited from fiberseq repo.
///
/// ```
/// use nanalogue_core::{Error, nanalogue_bam_reader};
/// use rust_htslib::bam::Read;
/// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
/// // the above file should contain four reads, so we are checking
/// // if we load four records.
/// let mut count = 0;
/// for r in reader.records() {
///     count = count + 1;
/// }
/// assert_eq!(count, 4);
/// # Ok::<(), Error>(())
/// ```
pub fn nanalogue_bam_reader(bam_path: &str) -> Result<bam::Reader, Error> {
    if bam_path == "-" {
        Ok(bam::Reader::from_stdin()?)
    } else {
        Ok(bam::Reader::from_path(bam_path)?)
    }
}

/// Implements a collection-of-states of ReadState
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReadStates(Vec<u16>);

impl FromStr for ReadStates {
    type Err = Error;

    /// converts a comma-separated list of read states to ReadStates
    ///
    /// ```
    /// use nanalogue_core::{Error, ReadStates};
    /// use std::str::FromStr;
    /// let mut op = ReadStates::from_str("primary_forward")?;
    /// assert_eq!(op.bam_flags(), &[0]);
    /// let mut op = ReadStates::from_str("unmapped,secondary_reverse,supplementary_reverse")?;
    /// assert_eq!(op.bam_flags(), &[4, 256 + 16, 2048 + 16]);
    /// op = ReadStates::from_str("primary_reverse")?;
    /// assert_eq!(op.bam_flags(), &[16]);
    /// op = ReadStates::from_str("primary_reverse,primary_forward")?;
    /// assert_eq!(op.bam_flags(), &[16, 0]);
    /// # Ok::<(), Error>(())
    /// ```
    ///
    /// ```should_panic
    /// use nanalogue_core::{Error, ReadStates};
    /// use std::str::FromStr;
    /// let mut op = ReadStates::from_str("random")?;
    /// # Ok::<(), Error>(())
    /// ```
    fn from_str(s: &str) -> Result<ReadStates, Self::Err> {
        let mut states = Vec::<u16>::new();
        for part in s.split(",") {
            states.push(u16::try_from(ReadState::from_str(part)?)?);
        }
        Ok(ReadStates(states))
    }
}

impl ReadStates {
    /// Returns the flags contained within
    pub fn bam_flags(&self) -> &Vec<u16> {
        &self.0
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
                let random_number: f32 = rand::random();
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
mod tests {
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
            let Ok(BaseMods { base_mods: v }) =
                nanalogue_mm_ml_parser(&r, |&_| true, |&_| true, |&_, &_, &_| true, 0)
            else {
                unreachable!()
            };
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
}

#[cfg(test)]
mod zero_length_filtering_tests {
    use super::*;
    use rust_htslib::bam::Read;

    #[test]
    fn test_zero_length_filtering_with_example_2_zero_len_sam() -> Result<(), Error> {
        // Test with include_zero_len = false and min_seq_len = 1 (should get 1 read)
        let mut reader = nanalogue_bam_reader("examples/example_2_zero_len.sam")?;
        let mut bam_opts_exclude_zero = InputBam::default();
        bam_opts_exclude_zero.include_zero_len = false;
        bam_opts_exclude_zero.min_seq_len = 1;

        let mut count_exclude_zero = 0;
        for record_result in reader.records() {
            let record = record_result?;
            if record.pre_filt(&bam_opts_exclude_zero) {
                count_exclude_zero += 1;
            }
        }

        // Test with include_zero_len = true and min_seq_len = 1 (should get 2 reads)
        let mut reader = nanalogue_bam_reader("examples/example_2_zero_len.sam")?;
        let mut bam_opts_include_zero = InputBam::default();
        bam_opts_include_zero.include_zero_len = true;
        bam_opts_include_zero.min_seq_len = 1;

        let mut count_include_zero = 0;
        for record_result in reader.records() {
            let record = record_result?;
            if record.pre_filt(&bam_opts_include_zero) {
                count_include_zero += 1;
            }
        }

        assert_eq!(
            count_exclude_zero, 1,
            "Should get 1 read when include_zero_len is false with min_seq_len=1"
        );
        assert_eq!(
            count_include_zero, 2,
            "Should get 2 reads when include_zero_len is true, even with min_seq_len=1"
        );

        Ok(())
    }
}
