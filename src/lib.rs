use bio::alphabets::dna::revcomp;
use fibertools_rs::utils::basemods::{BaseMod, BaseMods};
use fibertools_rs::utils::bio_io::*;
use lazy_static::lazy_static;
use regex::Regex;
use rust_htslib::{bam, bam::record::Aux};
use std::convert::TryFrom;
use std::fmt;
use std::collections::HashMap;

// Declare the modules.
pub mod subcommands;

// A read can exist in four states
pub enum ReadState {
    Primary,
    Secondary,
    Supplementary,
    Unmapped,
}

// Assume a read is primary unless otherwise stated,
// so we have three possible transitions
pub enum ReadTransition {
    PrimaryToSecondary,
    PrimaryToSupplementary,
    PrimaryToUnmapped,
}

pub struct CurrRead {
    state: ReadState,
    read_id: Option<String>,
    seq_len: Option<u64>,
    align_len: Option<u64>,
    mods: Option<Vec<BaseMod>>,
}

impl CurrRead {
    fn new() -> Self {
        Self {
            state: ReadState::Primary,
            read_id: None,
            seq_len: None,
            align_len: None,
            mods: None,
        }
    }
    fn transition(&mut self, transition: ReadTransition) {
        match (&self.state, transition) {
            (ReadState::Primary, ReadTransition::PrimaryToSecondary) => {
                self.state = ReadState::Secondary
            }
            (ReadState::Primary, ReadTransition::PrimaryToSupplementary) => {
                self.state = ReadState::Supplementary
            }
            (ReadState::Primary, ReadTransition::PrimaryToUnmapped) => {
                self.state = ReadState::Unmapped
            }
            _ => {
                eprintln!("Invalid state reached!");
                std::process::exit(1);
            }
        }
    }
    fn header_string() -> String {
        "read_id\tsequence_length_template\talign_length\talignment_type\tmod_count".to_string()
    }
    fn mod_count_per_mod(&self) -> Option<HashMap<char, u32>> {
        let mut output = HashMap::<char, u32>::new();
        match &self.mods {
            Some(v) => {
                if v.is_empty() {
                    None
                } else {
                    for k in v {
                        let mod_count = k.ranges.qual.len() as u32;
                        output.entry(k.modification_type).and_modify(|e| *e += mod_count ).or_insert(mod_count);
                    }
                    Some(output)
                }
            }
            None => None,
        }
    }
}

impl fmt::Display for CurrRead {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut output_string = String::from("");

        if let Some(v) = &self.read_id {
            output_string = output_string + v;
        } else {
            output_string += "NA";
        }

        if let Some(v) = self.seq_len {
            output_string = output_string + "\t" + &v.to_string();
        } else {
            output_string += "\tNA";
        }

        if let Some(v) = self.align_len {
            output_string = output_string + "\t" + &v.to_string();
        } else {
            output_string += "\tNA";
        }

        match self.state {
            ReadState::Primary => output_string = output_string + "\t" + "primary",
            ReadState::Secondary => output_string = output_string + "\t" + "secondary",
            ReadState::Supplementary => output_string = output_string + "\t" + "supplementary",
            ReadState::Unmapped => output_string = output_string + "\t" + "unmapped",
        }

        if let Some(v) = &self.mods {
            if !v.is_empty() {
                output_string += "\t";
                for k in v {
                    output_string += format!("{}{}{}:{};",
                        k.modified_base as char,
                        k.strand,
                        match k.modification_type {
                            'A'..='Z' | 'a'..='z' => k.modification_type.to_string(),
                            _ => format!("{}", k.modification_type as u32),
                        },
                        k.ranges.qual.len()
                    ).as_str();
                }
            } else {
                output_string += "\t0";
            }
        } else {
            output_string += "\t0";
        }

        write!(f, "{output_string}")
    }
}

fn convert_seq_uppercase(mut seq: Vec<u8>) -> Vec<u8> {
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

fn process_mod_type(mod_type: &str) -> Result<char, String> {
    // process the modification type, returning the first character if it is a letter,
    // or converting it to a character if it is a number
    let first_char = match mod_type.chars().next() {
        Some(c) => c,
        None => return Err(format!("Modification type is empty: {}", mod_type)),
    };
    match first_char {
        'A' ..= 'Z' | 'a' ..= 'z' => Ok(first_char),
        '0' ..= '9' => {
            let u: u32 = match mod_type.parse() {
                Ok(num) => num,
                Err(_) => return Err(format!("Invalid modification type: {}", mod_type)),
            };
            char::from_u32(u).ok_or_else(|| {
                format!("Invalid modification type: {} (not a valid Unicode character)", mod_type)
            })
        },
        _ => Err(format!("Invalid modification type: {}", mod_type)),
    }
}

// We are copying and modifying code from the fibertools-rs repository.
// https://github.com/fiberseq/fibertools-rs
pub fn nanalogue_mm_ml_parser(record: &bam::Record, min_ml_score: u8) -> BaseMods {
    // regex for matching the MM tag
    lazy_static! {
        // MM:Z:([ACGTUN][-+]([A-Za-z]+|[0-9]+)[.?]?(,[0-9]+)*;)*
        static ref MM_RE: Regex =
            Regex::new(r"((([ACGTUN])([-+])([A-Za-z]+|[0-9]+))[.?]?((,[0-9]+)*;)*)").unwrap();
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
            let modification_type = cap.get(5).map_or("", |m| m.as_str());
            let mod_dists_str = cap.get(6).map_or("", |m| m.as_str());
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
                process_mod_type(modification_type).unwrap(),
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
pub fn nanalogue_bam_reader(bam_path: &str) -> bam::Reader {
    match bam_path {
        "-" => bam::Reader::from_stdin().unwrap_or_else(|e| {
            eprintln!("Problem opening file, error: {e}");
            std::process::exit(1)
        }),
        s => bam::Reader::from_path(s).unwrap_or_else(|e| {
            eprintln!("Problem opening file, error: {e}");
            std::process::exit(1)
        }),
    }
}
