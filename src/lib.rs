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
use rust_htslib::{bam, bam::Read, bam::record::Aux, tpool};
use std::convert::TryFrom;

// Declare the modules.
pub mod cli;
pub mod error;
pub mod read_utils;
pub mod subcommands;
pub mod utils;

// Re-exports
pub use cli::InputBam;
pub use error::Error;
pub use read_utils::{Contains, CurrRead, ReadState, ThresholdState};
pub use utils::{F32Bw0and1, ModChar, OrdPair};

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
pub fn nanalogue_mm_ml_parser<F, G>(
    record: &bam::Record,
    filter_mod_prob_pos: F,
    filter_mod_base_strand_tag: G,
) -> BaseMods
where
    F: Fn(&u8, &i64) -> bool,
    G: Fn(u8, char, ModChar) -> bool,
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
            if !filter_mod_base_strand_tag(mod_base, mod_strand, modification_type) {
                continue;
            }

            let _implicit = cap.get(6).map_or(".", |m| m.as_str()).as_bytes().first();
            let mod_dists_str = cap.get(7).map_or("", |m| m.as_str());
            // parse the string containing distances between modifications into a vector of i64

            let mod_dists: Vec<i64> = mod_dists_str
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
            // find real positions in the forward sequence
            let mut cur_mod_idx = 0;
            let mut cur_seq_idx = 0;
            let mut dist_from_last_mod_base = 0;
            let mut unfiltered_modified_positions: Vec<i64> = vec![0; mod_dists.len()];
            while cur_seq_idx < forward_bases.len() && cur_mod_idx < mod_dists.len() {
                let cur_base = forward_bases[cur_seq_idx];
                if (cur_base == mod_base || mod_base == b'N')
                    && dist_from_last_mod_base == mod_dists[cur_mod_idx]
                {
                    unfiltered_modified_positions[cur_mod_idx] =
                        i64::try_from(cur_seq_idx).unwrap();
                    dist_from_last_mod_base = 0;
                    cur_mod_idx += 1;
                } else if cur_base == mod_base {
                    dist_from_last_mod_base += 1
                }
                cur_seq_idx += 1;
            }
            // assert that we extract the same number of modifications as we have distances
            assert_eq!(
                cur_mod_idx,
                mod_dists.len(),
                "{:?} {}",
                String::from_utf8_lossy(record.qname()),
                record.is_reverse()
            );

            // check for the probability of modification.
            let num_mods_cur_end = num_mods_seen + unfiltered_modified_positions.len();
            let unfiltered_modified_probabilities = if num_mods_cur_end > ml_tag.len() {
                let needed_num_of_zeros = num_mods_cur_end - ml_tag.len();
                let mut to_add = vec![0; needed_num_of_zeros];
                let mut has = ml_tag[num_mods_seen..ml_tag.len()].to_vec();
                has.append(&mut to_add);
                log::warn!(
                    "{} {}",
                    "ML tag is too short for the number of modifications found in the MM tag.",
                    "Assuming an ML value of 0 after the first {num_mods_cur_end} modifications."
                );
                has
            } else {
                ml_tag[num_mods_seen..num_mods_cur_end].to_vec()
            };
            num_mods_seen = num_mods_cur_end;

            // must be true for filtering, and at this point
            assert_eq!(
                unfiltered_modified_positions.len(),
                unfiltered_modified_probabilities.len()
            );

            // Filter mods based on probabilities
            let (modified_probabilities, modified_positions): (Vec<u8>, Vec<i64>) =
                unfiltered_modified_probabilities
                    .iter()
                    .zip(unfiltered_modified_positions.iter())
                    .filter(|&(&ml, &mm)| filter_mod_prob_pos(&ml, &mm))
                    .unzip();

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
        log::trace!("No MM tag found");
    }

    if ml_tag.len() != num_mods_seen {
        log::warn!(
            "ML tag ({}) different number than MM tag ({}).",
            ml_tag.len(),
            num_mods_seen
        );
    }
    // needed so I can compare methods
    rtn.sort();
    BaseMods { base_mods: rtn }
}

/// A global struct which contains BAM records for further usage.
/// TODO: Filtering. If any filtering is done here, it is using non-base-mod
/// information in the BAM file such as alignment lengths or coordinates etc.
/// NOTE: we don't derive many traits here as the RcRecords object
/// does not have many traits.
#[derive(Debug)]
pub struct BamRcRecords<'a> {
    /// RcRecords object output by rust htslib which we can iterate over
    pub rc_records: bam::RcRecords<'a, bam::Reader>,
    /// Header of the bam file
    pub header: bam::Header,
    /// List of contig names
    pub contig_names: Vec<String>,
}

impl<'a> BamRcRecords<'a> {
    /// Extracts RcRecords from a BAM Reader
    pub fn new(bam_reader: &'a mut bam::Reader, bam_opts: &InputBam) -> Result<Self, Error> {
        let contig_names = bam_reader
            .header()
            .target_names()
            .iter()
            .map(|r| str::from_utf8(r).map(String::from))
            .collect::<Result<Vec<String>, _>>()?;
        let header = bam::Header::from_template(bam_reader.header());
        let tp = tpool::ThreadPool::new(bam_opts.threads.get())?;
        bam_reader.set_thread_pool(&tp)?;
        let rc_records = bam_reader.rc_records();
        Ok(BamRcRecords {
            rc_records,
            header,
            contig_names,
        })
    }
}

/// Opens BAM file, also copied and edited from fiberseq repo.
pub fn nanalogue_bam_reader(bam_path: &str) -> Result<bam::Reader, Error> {
    if bam_path == "-" {
        Ok(bam::Reader::from_stdin()?)
    } else {
        Ok(bam::Reader::from_path(bam_path)?)
    }
}

/// Trait that performs filtration on structs implementing SequenceRead
/// such as a rust_htslib Record
pub trait BamPreFilt: SequenceRead {
    /// apply some default filtration e.g. by read length
    fn pre_filt(&self, bam_opts: &InputBam) -> bool {
        self.len() as u64 >= bam_opts.min_seq_len
    }
}

impl BamPreFilt for bam::record::Record {}
