//! # Makes a table of information on reads
//!
//! This file contains routines used to calculate some information per
//! read such as sequence length, alignment length, type of alignment,
//! and modification counts, and displays them in table using the run
//! function. The routine reads both BAM and sequencing summary files
//! if provided, otherwise only reads the BAM file.

use crate::{CurrRead, Error, InputMods, ModChar, OptionalTag, ReadState, ThresholdState};
use bedrs::prelude::{Bed3, Coordinates};
use csv::ReaderBuilder;
use itertools::join;
use rust_htslib::bam;
use serde::{Deserialize, Serialize};
use std::{collections::HashMap, fmt, fmt::Write, fs::File, rc::Rc, str};

/// Write a vector as a CSV-formatted string
macro_rules! vec_csv {
    ( $b: expr ) => {
        join($b, ", ")
    };
}

/// Declare a custom type for ease of use
#[derive(Debug, Clone, PartialEq)]
struct ModCountTbl(HashMap<ModChar, u32>);

/// Implements a constructor
impl ModCountTbl {
    /// construct from a hashmap given as input
    fn new(value: HashMap<ModChar, u32>) -> Self {
        Self(value)
    }
    /// constructs using a blank hashmap
    fn blank() -> Self {
        Self(HashMap::<ModChar, u32>::new())
    }
}

/// Implements a display method
impl fmt::Display for ModCountTbl {
    /// display contents as key:value separated by semicolons.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut v: Vec<(ModChar, u32)> = self.0.clone().into_iter().collect();
        v.sort_by_key(|k| k.0);
        write!(
            f,
            "{}",
            if v.is_empty() {
                "NA".to_string()
            } else {
                join(v.into_iter().map(|k| format!("{}:{}", k.0, k.1)), ";")
            }
        )
    }
}

/// Enum showing state of a read composed from BAM and optionally
/// a sequencing summary file.
/// - only BAM file alignment information is available
/// - only sequencing summary file information is available
/// - both are available, in which case we retain only the sequence
///   length from the sequencing summary file and the alignment length
///   from the BAM file, discarding the sequence length from the BAM file.
#[derive(Debug, Clone, PartialEq)]
enum ReadInstance {
    OnlyAlign {
        align_len: Vec<u64>,
        seq_len: u64,
        mod_count: Vec<ModCountTbl>,
        read_state: Vec<ReadState>,
        seq: Vec<String>,
    },
    OnlyBc(u64),
    BothAlignBc {
        align_len: Vec<u64>,
        bc_len: u64,
        mod_count: Vec<ModCountTbl>,
        read_state: Vec<ReadState>,
        seq: Vec<String>,
    },
}

impl fmt::Display for ReadInstance {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ReadInstance::OnlyBc(u64) => write!(f, "{u64}"),
            ReadInstance::OnlyAlign {
                align_len: al,
                seq_len: sl,
                mod_count: mc,
                read_state: rs,
                seq: sq,
            }
            | ReadInstance::BothAlignBc {
                align_len: al,
                bc_len: sl,
                mod_count: mc,
                read_state: rs,
                seq: sq,
            } => write!(
                f,
                "{}\t{sl}\t{}{}{}",
                vec_csv!(al),
                vec_csv!(rs),
                match vec_csv!(mc) {
                    v if !v.is_empty() => format!("\t{v}"),
                    _ => String::new(),
                },
                match vec_csv!(sq) {
                    v if !v.is_empty() => format!("\t{v}"),
                    _ => String::new(),
                },
            ),
        }
    }
}

/// Implements a structure representing the state of a read w.r.t
/// sequencing summary information and BAM information
struct Read(ReadInstance);

/// Implements functions for Read
impl Read {
    /// initialization using basecalled length
    fn new_bc_len(l: u64) -> Self {
        Self(ReadInstance::OnlyBc(l))
    }

    /// initialization using alignment information
    fn new_align_len(
        align_len: u64,
        seq_len: u64,
        mod_count: Option<ModCountTbl>,
        read_state: ReadState,
        seq: Option<String>,
    ) -> Self {
        Self(ReadInstance::OnlyAlign {
            align_len: vec![align_len],
            seq_len,
            mod_count: mod_count.into_iter().collect(),
            read_state: vec![read_state],
            seq: seq.into_iter().collect(),
        })
    }

    /// adding alignment information to an entry
    fn add_align_len(
        &mut self,
        align_len: u64,
        mod_count: Option<ModCountTbl>,
        read_state: ReadState,
        seq: Option<String>,
    ) {
        match &mut self.0 {
            ReadInstance::OnlyAlign {
                align_len: al,
                mod_count: mc,
                read_state: rs,
                seq: sq,
                ..
            }
            | ReadInstance::BothAlignBc {
                align_len: al,
                mod_count: mc,
                read_state: rs,
                seq: sq,
                ..
            } => {
                al.push(align_len);
                rs.push(read_state);
                mc.extend(mod_count);
                sq.extend(seq);
            }
            ReadInstance::OnlyBc(bl) => {
                self.0 = ReadInstance::BothAlignBc {
                    align_len: vec![align_len],
                    bc_len: *bl,
                    mod_count: mod_count.into_iter().collect(),
                    read_state: vec![read_state],
                    seq: seq.into_iter().collect(),
                };
            }
        }
    }
}

/// Represents a record with the columns of interest from the TSV file.
#[derive(Debug, Serialize, Deserialize)]
struct TSVRecord {
    read_id: String,
    sequence_length_template: u64,
}

/// Opens a TSV file, extracts '`read_id`' and '`sequence_length_template`' columns,
/// and builds a `HashMap`.
fn process_tsv(file_path: &str) -> Result<HashMap<String, Read>, Error> {
    let mut data_map = HashMap::<String, Read>::new();

    match file_path {
        "" => {}
        fp => {
            let file = File::open(fp)?;
            let mut rdr = ReaderBuilder::new()
                .comment(Some(b'#'))
                .delimiter(b'\t')
                .from_reader(file);

            for result in rdr.deserialize() {
                let record: TSVRecord = result?;
                if data_map
                    .insert(
                        record.read_id,
                        Read::new_bc_len(record.sequence_length_template),
                    )
                    .is_some()
                {
                    return Err(Error::InvalidDuplicates(file_path.to_string()));
                }
            }
        }
    }

    Ok(data_map)
}

/// Processes a BAM file and optionally a sequencing summary file
/// to print a table of reads with alignment length, sequence length,
/// and optionally modification count per row.
///
/// # Errors
/// Returns an error if TSV processing, BAM record reading, sequence retrieval, or output writing fails.
pub fn run<W, D>(
    handle: &mut W,
    bam_records: D,
    mut mods: Option<InputMods<OptionalTag>>,
    seq_region: Option<Bed3<i32, u64>>,
    seq_summ_path: &str,
) -> Result<(), Error>
where
    W: std::io::Write,
    D: IntoIterator<Item = Result<Rc<bam::Record>, rust_htslib::errors::Error>>,
{
    // read TSV file and convert into hashmap
    let mut data_map = process_tsv(seq_summ_path)?;

    // set up a flag to check if sequencing summary file has data
    let is_seq_summ_data: bool = !data_map.is_empty();

    match &mut mods {
        None => {}
        Some(p) => match p.mod_prob_filter {
            ref mut v @ ThresholdState::GtEq(w) => *v = ThresholdState::GtEq(u8::max(128, w)),
            ref mut v @ ThresholdState::InvertGtEqLtEq(w) => *v = ThresholdState::Both((128, w)),
            ref mut v @ ThresholdState::Both((w, x)) => {
                *v = ThresholdState::Both((u8::max(128, w), x));
            }
        },
    }

    // Go record by record in the BAM file,
    // get the read id and the alignment length, and put it in the hashmap
    for r in bam_records {
        // read records
        let record = r?;

        // get information of current read
        let curr_read_state = CurrRead::default().try_from_only_alignment(&record)?;
        let qname = String::from(curr_read_state.read_id()?);
        let read_state = curr_read_state.read_state();
        let align_len = match curr_read_state.align_len() {
            Ok(v) => v,
            Err(Error::Unmapped) => 0,
            Err(_) => continue,
        };
        let Ok(seq_len) = curr_read_state.seq_len() else {
            continue;
        };

        // get sequence
        let sequence = match &seq_region {
            Some(v) if v.len() != 0 => Some(match curr_read_state.seq_on_ref_coords(&record, v) {
                Err(Error::UnavailableData) => Ok(String::from("*")),
                Err(e) => Err(e),
                Ok(w) => Ok(String::from_utf8(w)?),
            }?),
            Some(_) => Some(String::from_utf8(record.seq().as_bytes())?),
            None => None,
        };

        // get modification information
        let mod_count: Option<ModCountTbl> = match &mods {
            None => None,
            Some(v) => {
                match curr_read_state
                    .set_mod_data_restricted_options(&record, v)?
                    .base_count_per_mod()
                {
                    None => Some(ModCountTbl::blank()),
                    Some(v) => Some(ModCountTbl::new(v)),
                }
            }
        };

        // add data depending on whether an entry is already present
        // in the hashmap from the sequencing summary file
        let _: &mut _ = data_map
            .entry(qname)
            .and_modify(|entry| {
                entry.add_align_len(align_len, mod_count.clone(), read_state, sequence.clone());
            })
            .or_insert(Read::new_align_len(
                align_len, seq_len, mod_count, read_state, sequence,
            ));
    }

    // set up an output header string
    let mut output_header = String::new();
    if is_seq_summ_data {
        writeln!(output_header, "# seq summ file: {seq_summ_path}")?;
    }
    write!(
        output_header,
        "{}read_id\talign_length\tsequence_length_template\talignment_type{}{}",
        match &mods {
            Some(_) => "# mod counts using probability threshold of 0.5\n",
            None => "",
        },
        match &mods {
            Some(_) => "\tmod_count",
            None => "",
        },
        match &seq_region {
            Some(_) => "\tsequence",
            None => "",
        },
    )?;

    // print the output header
    writeln!(handle, "{output_header}")?;

    // print output tsv data
    // If both seq summ and BAM file are available, then the length in the seq
    // summ file takes precedence. If only the BAM file is available, then
    // we use the sequence length in the BAM file as the basecalled sequence length.
    #[expect(
        clippy::iter_over_hash_type,
        reason = "sorting would add unnecessary performance overhead; random iteration order is acceptable"
    )]
    for (key, val) in &data_map {
        match (&val, is_seq_summ_data) {
            (Read(ReadInstance::OnlyBc(_)), _) | (Read(ReadInstance::OnlyAlign { .. }), true) => {}
            (Read(ReadInstance::BothAlignBc { .. }), false) => {
                return Err(Error::InvalidState(
                    "invalid state while writing output".to_string(),
                ));
            }
            (Read(v), _) => writeln!(handle, "{key}\t{v}")?,
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nanalogue_bam_reader;
    use rust_htslib::bam::Read as BamRead;

    fn normalize_output(output: &str) -> Vec<String> {
        let lines: Vec<String> = output
            .lines()
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect();

        // Separate header and data lines
        let mut headers = Vec::new();
        let mut data = Vec::new();

        for line in lines {
            if line.starts_with('#') || line.starts_with("read_id") {
                headers.push(line);
            } else {
                data.push(line);
            }
        }

        // Sort data lines to handle non-deterministic ordering
        data.sort();

        // Combine headers and sorted data
        headers.extend(data);
        headers
    }

    fn run_read_table_test(mods: Option<InputMods<OptionalTag>>, expected_output_file: &str) {
        let mut reader = nanalogue_bam_reader("./examples/example_1.bam").expect("no error");
        let records: Vec<_> = reader.records().map(|r| r.map(Rc::new)).collect();

        let mut output = Vec::new();
        run(&mut output, records, mods, None, "").expect("no error");

        let actual_output = String::from_utf8(output).expect("Invalid UTF-8");
        let expected_output = std::fs::read_to_string(expected_output_file)
            .expect("Failed to read expected output file");

        let actual_lines = normalize_output(&actual_output);
        let expected_lines = normalize_output(&expected_output);

        assert_eq!(
            actual_lines,
            expected_lines,
            "\nActual output:\n{}\n\nExpected output:\n{}\n",
            actual_lines.join("\n"),
            expected_lines.join("\n")
        );
    }

    #[test]
    fn read_table_hide_mods() {
        run_read_table_test(None, "./examples/example_1_read_table_hide_mods");
    }

    #[test]
    fn read_table_show_mods() {
        run_read_table_test(
            Some(InputMods::<OptionalTag>::default()),
            "./examples/example_1_read_table_show_mods",
        );
    }
}
