#![deny(missing_docs)]
//! # Nanalogue Core
//!
//! We process and calculate data associated with DNA molecules, their alignments to
//! reference genomes, modification information on them, and other miscellaneous
//! information.
//!
//! This file directly contains the functions associated with opening BAM files
//! and parsing modification information directly from BAM files. Other functions
//! in the crate are distributed over other files.

use bio::alphabets::dna::revcomp;
use bio_types::sequence::SequenceRead;
use fibertools_rs::utils::basemods::{BaseMod, BaseMods};
use fibertools_rs::utils::bio_io::get_u8_tag;
use lazy_static::lazy_static;
use regex::Regex;
use rust_htslib::{bam, bam::Read, bam::ext::BamRecordExtensions, bam::record::Aux, tpool};
use std::convert::TryFrom;

// Declare the modules.
pub mod cli;
pub mod error;
pub mod read_utils;
pub mod subcommands;
pub mod utils;

// Re-exports
pub use cli::{InputBam, InputWindowing, InputWindowingRestricted};
pub use error::Error;
pub use read_utils::{CurrRead, ReadState, ThresholdState};
pub use utils::{Contains, F32AbsValBelow1, F32Bw0and1, ModChar, OrdPair, RestrictModCalledStrand};

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
///     let BaseMods{base_mods: v} = nanalogue_mm_ml_parser(&r, |&_| true, |&_| true, |&_, &_, &_| true);
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
) -> BaseMods
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

            // get modification type and skip record if we don't find
            // mod of interest (if specified)
            let modification_type: ModChar = cap
                .get(5)
                .map_or("", |m| m.as_str())
                .parse()
                .expect("error");
            if !filter_mod_base_strand_tag(&mod_base, &mod_strand, &modification_type) {
                continue;
            }

            let is_implicit = match cap.get(6).map_or("", |m| m.as_str()).as_bytes() {
                b"" => true,
                b"." => true,
                b"?" => false,
                _ => panic!("unknown modification annotation!"),
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
                    if filter_mod_prob(prob) && filter_mod_pos(&cur_seq_idx) {
                        modified_positions.push(i64::try_from(cur_seq_idx).unwrap());
                        modified_probabilities.push(*prob);
                    }
                    dist_from_last_mod_base = 0;
                    cur_mod_idx += 1;
                } else {
                    if is_include_zero_prob && is_implicit && filter_mod_pos(&cur_seq_idx) {
                        modified_positions.push(i64::try_from(cur_seq_idx).unwrap());
                        modified_probabilities.push(0);
                    }
                    dist_from_last_mod_base += 1
                }
            }

            num_mods_seen += cur_mod_idx;

            // ensure equal lengths for either type of data
            assert_eq!(modified_positions.len(), modified_probabilities.len());

            // don't add empty basemods
            if modified_positions.is_empty() {
                continue;
            }
            // add to a struct
            let mods = BaseMod::new(
                record,
                mod_base,
                mod_strand,
                modification_type.get_val(),
                modified_positions,
                modified_probabilities,
            );
            rtn.push(mods);
        }
    } else {
        panic!("No MM tag found!")
    }

    // ensure that we have seen all mods
    assert_eq!(ml_tag.len(), num_mods_seen);

    // needed so I can compare methods
    rtn.sort();
    BaseMods { base_mods: rtn }
}

/// A global struct which contains BAM records for further usage.
/// NOTE: we don't derive many traits here as the RcRecords object
/// does not have many traits.
///
/// ```
/// use nanalogue_core::{Error, nanalogue_bam_reader, BamRcRecords, InputBam};
/// use rust_htslib::bam::Read;
/// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
/// let mut bam_opts = InputBam::default();
/// let BamRcRecords = BamRcRecords::new(&mut reader, &bam_opts)?;
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
    pub fn new(bam_reader: &'a mut bam::Reader, bam_opts: &InputBam) -> Result<Self, Error> {
        let header = bam_reader.header().clone();
        let tp = tpool::ThreadPool::new(bam_opts.threads.get())?;
        bam_reader.set_thread_pool(&tp)?;
        let rc_records = bam_reader.rc_records();
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

/// Trait that performs filtration
pub trait BamPreFilt {
    /// apply default filtration
    fn pre_filt(&self, _bam_opts: &InputBam) -> bool {
        todo!()
    }
    /// filtration by length
    fn filt_by_len(&self, _bam_opts: &InputBam) -> bool {
        todo!()
    }
    /// filtration by alignment length
    fn filt_by_align_len(&self, _bam_opts: &InputBam) -> bool {
        todo!()
    }
    /// filtration by read id
    fn filt_by_read_id(&self, _bam_opts: &InputBam) -> bool {
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
/// If zero-length records are expected, a specific input
/// flag is needed for the program to ignore them. Otherwise,
/// we get a panic.
/// ```should_panic
/// # use nanalogue_core::{BamPreFilt, InputBam, Error};
/// # use rust_htslib::bam::record;
/// # use rust_htslib::bam::record::{Cigar, CigarString};
/// let mut bam_record = record::Record::new();
/// let mut bam_opts = InputBam::default();
/// bam_opts.min_seq_len = 20;
/// bam_opts.exclude_zero_len = false;
/// assert_eq!(bam_record.pre_filt(&bam_opts), false);
/// # Ok::<(), Error>(())
/// ```
/// No panic when the filter is set.
/// ```
/// # use nanalogue_core::{BamPreFilt, InputBam, Error};
/// # use rust_htslib::bam::record;
/// # use rust_htslib::bam::record::{Cigar, CigarString};
/// let mut bam_record = record::Record::new();
/// let mut bam_opts = InputBam::default();
/// bam_opts.min_seq_len = 20;
/// bam_opts.exclude_zero_len = true;
/// assert_eq!(bam_record.pre_filt(&bam_opts), false);
/// # Ok::<(), Error>(())
/// ```
impl BamPreFilt for bam::Record {
    /// apply default filtration by read length
    fn pre_filt(&self, bam_opts: &InputBam) -> bool {
        self.filt_by_len(bam_opts)
            & self.filt_by_read_id(bam_opts)
            & self.filt_by_align_len(bam_opts)
    }
    /// filtration by read length
    fn filt_by_len(&self, bam_opts: &InputBam) -> bool {
        match self.len() as u64 {
            0 if !bam_opts.exclude_zero_len => panic!(
                "{}{}{}",
                "Cannot deal with 0 length seq in BAM file. ",
                "For instance, this could happen when some or all seq fields are set to '*', ",
                "although this is valid BAM. See the input options for how to avoid this error. "
            ),
            0 if bam_opts.exclude_zero_len => false,
            v => v >= bam_opts.min_seq_len,
        }
    }
    /// filtration by alignment length
    fn filt_by_align_len(&self, bam_opts: &InputBam) -> bool {
        match (&bam_opts.min_align_len, self.is_unmapped()) {
            (None, _) => true,
            (_, true) => false,
            (Some(v), false) => self.reference_end() - self.pos() >= *v,
        }
    }
    /// filtration by read id
    fn filt_by_read_id(&self, bam_opts: &InputBam) -> bool {
        match &bam_opts.read_id {
            Some(v) => v.as_bytes() == self.name(),
            None => true,
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
        let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
        let mut count = 0;
        for record in reader.records() {
            let r = record?;
            let BaseMods { base_mods: v } =
                nanalogue_mm_ml_parser(&r, |&_| true, |&_| true, |&_, &_, &_| true);
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
            count = count + 1;
        }
        Ok(())
    }
}
