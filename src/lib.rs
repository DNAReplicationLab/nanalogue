use lazy_static::lazy_static;
use regex::Regex;
use rust_htslib::{bam, bam::record::Aux};
use std::convert::TryFrom;
use fibertools_rs::utils::basemods::{BaseMod, BaseMods};
use fibertools_rs::utils::bio_io::get_u8_tag;
use bio::alphabets::dna::revcomp;

// Declare the modules.
pub mod subcommands;
pub mod error;
pub mod read_utils;

// Re-export the error type
pub use error::Error;

// Re-export read utils
pub use read_utils::{ReadState, CurrRead};

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

pub fn process_mod_type(mod_type: &str) -> Result<char, Error> {
    // process the modification type, returning the first character if it is a letter,
    // or converting it to a character if it is a number
    let first_char = mod_type.chars().next().ok_or(Error::EmptyModType)?; 
    match first_char {
        'A' ..= 'Z' | 'a' ..= 'z' => Ok(first_char),
        '0' ..= '9' => char::from_u32(mod_type.parse()?).ok_or(Error::InvalidModType),
        _ => Err(Error::InvalidModType),
    }
}

// We are copying and modifying code from the fibertools-rs repository.
// https://github.com/fiberseq/fibertools-rs
pub fn nanalogue_mm_ml_parser(record: &bam::Record,
    min_ml_score: u8, mod_tag: Option<char>) -> BaseMods {
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
            let mod_strand = cap.get(4).map_or("", |m| m.as_str());

            // get modification type and skip record if we don't find
            // mod of interest (if specified)
            let modification_type = process_mod_type(cap.get(5).map_or("", |m| m.as_str())).unwrap();
            if let Some(v) = mod_tag {
                if v != modification_type {
                    continue;
                }
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
                    .filter(|&(&ml, &_mm)| ml >= min_ml_score)
                    .unzip();

            // don't add empty basemods
            if modified_positions.is_empty() {
                continue;
            }
            // add to a struct
            let mods = BaseMod::new(
                record,
                mod_base,
                mod_strand.chars().next().unwrap(),
                modification_type,
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

/// Opens BAM file, also copied and edited from fiberseq repo.
pub fn nanalogue_bam_reader(bam_path: &str) -> Result<bam::Reader, Error> {
    match bam_path {
        "-" => Ok(bam::Reader::from_stdin()?),
        s => Ok(bam::Reader::from_path(s)?),
    }
}
