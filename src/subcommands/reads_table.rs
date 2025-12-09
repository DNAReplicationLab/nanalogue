//! # Makes a table of information on reads
//!
//! This file contains routines used to calculate some information per
//! read such as sequence length, alignment length, type of alignment,
//! and modification counts, and displays them in table using the run
//! function. The routine reads both BAM and sequencing summary files
//! if provided, otherwise only reads the BAM file.

use crate::{
    CurrRead, Error, InputMods, ModChar, OptionalTag, ReadState, SeqCoordCalls, SeqDisplayOptions,
    ThresholdState,
};
use csv::ReaderBuilder;
use itertools::join;
use polars::prelude::*;
use rust_htslib::bam;
use serde::{Deserialize, Serialize};
use std::{collections::HashMap, fmt, fs::File, iter, rc::Rc, str};

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
}

/// Implements a display method
impl fmt::Display for ModCountTbl {
    /// display contents as key:value separated by semicolons.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut v: Vec<(ModChar, u32)> = self.0.clone().into_iter().collect();
        v.sort_by_key(|k| k.0);
        (if v.is_empty() {
            "NA".to_owned()
        } else {
            join(v.into_iter().map(|k| format!("{}:{}", k.0, k.1)), ";")
        })
        .fmt(f)
    }
}

/// Enum with instructions on how a base should be formatted
#[derive(Debug, Clone, PartialEq)]
enum BaseFmt {
    /// No special formatting needed
    No,
    /// Show as lowercase
    LowerCase,
    /// Show as highlight
    Highlight,
    /// Show as lowercase, highlight
    LowerCaseHighlight,
}

/// Enum showing state of a read composed from BAM and optionally
/// a sequencing summary file.
/// - only BAM file alignment information is available
/// - only sequencing summary file information is available
/// - both are available
#[derive(Debug, Clone, PartialEq)]
enum ReadInstance {
    /// Only information from the alignment file is available.
    /// NOTE that these are vectors as a read can have multiple BAM records.
    /// Sequencing length is not a vector as it is supposed to be unique.
    /// If multiple BAM records of one read have different sequencing lengths,
    /// we will have to make a choice.
    OnlyAlign {
        /// Alignment length from the BAM file
        align_len: Vec<u64>,
        /// Sequence length from the BAM file
        seq_len: u64,
        /// Count of various mods from the BAM file
        mod_count: Vec<ModCountTbl>,
        /// Alignment type (primary, secondary, reverse etc.) from the BAM file.
        read_state: Vec<ReadState>,
        /// Sequence from the BAM file.
        seq: Vec<Vec<(u8, BaseFmt)>>,
        /// Basecalling qualities from the BAM file.
        qual: Vec<Vec<u8>>,
    },
    /// Only basecalled info i.e. from the sequencing summary file is available
    OnlyBc(u64),
    /// Information from both sequencing summary and alignment file are available.
    /// NOTE that these are vectors as a read can have multiple BAM records,
    /// but the basecalled length is not a vector as a sequencing summary file
    /// has unique information per read. If the basecalled length from the
    /// sequencing summary file does not match that from the BAM file record(s),
    /// then we will have to make a choice.
    BothAlignBc {
        /// Alignment length from the BAM file
        align_len: Vec<u64>,
        /// Sequence length from the sequencing summary file
        bc_len: u64,
        /// Count of various mods from the BAM file
        mod_count: Vec<ModCountTbl>,
        /// Alignment type (primary, secondary, reverse etc.) from the BAM file.
        read_state: Vec<ReadState>,
        /// Sequence from the BAM file.
        seq: Vec<Vec<(u8, BaseFmt)>>,
        /// Basecalling qualities from the BAM file.
        qual: Vec<Vec<u8>>,
    },
}

impl fmt::Display for ReadInstance {
    #[expect(
        clippy::pattern_type_mismatch,
        reason = "expect no confusion from this code"
    )]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // we are fine with `unsafe` here as `rust-htslib` guarantees
        // that only 16 printable characters are allowed in the sequence.
        match self {
            ReadInstance::OnlyBc(u64) => format!("{u64}"),
            ReadInstance::OnlyAlign {
                align_len: al,
                seq_len: sl,
                mod_count: mc,
                read_state: rs,
                seq: sq,
                qual: q,
            }
            | ReadInstance::BothAlignBc {
                align_len: al,
                bc_len: sl,
                mod_count: mc,
                read_state: rs,
                seq: sq,
                qual: q,
            } => format!(
                "{}\t{sl}\t{}{}{}{}",
                vec_csv!(al),
                vec_csv!(rs),
                match vec_csv!(mc) {
                    v if !v.is_empty() => format!("\t{v}"),
                    _ => String::new(),
                },
                match vec_csv!(sq.iter().map(|v| unsafe {
                    String::from_utf8_unchecked(
                        v.clone()
                            .into_iter()
                            .map(|w| match w.1 {
                                BaseFmt::No => w.0,
                                BaseFmt::LowerCase => w.0.to_ascii_lowercase(),
                                BaseFmt::Highlight => b'Z',
                                BaseFmt::LowerCaseHighlight => b'z',
                            })
                            .collect::<Vec<u8>>(),
                    )
                })) {
                    v if !v.is_empty() => format!("\t{v}"),
                    _ => String::new(),
                },
                match vec_csv!(q.iter().map(|k| join(k, "."))) {
                    v if !v.is_empty() => format!("\t{v}"),
                    _ => String::new(),
                },
            ),
        }
        .fmt(f)
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
        seq: Option<Vec<(u8, BaseFmt)>>,
        qual: Option<Vec<u8>>,
    ) -> Self {
        Self(ReadInstance::OnlyAlign {
            align_len: vec![align_len],
            seq_len,
            mod_count: mod_count.into_iter().collect(),
            read_state: vec![read_state],
            seq: seq.into_iter().collect(),
            qual: qual.into_iter().collect(),
        })
    }

    /// adding alignment information to an entry
    #[expect(
        clippy::pattern_type_mismatch,
        reason = "expect no confusion from this code"
    )]
    fn add_align_len(
        &mut self,
        align_len: u64,
        seq_len: u64,
        mod_count: Option<ModCountTbl>,
        read_state: ReadState,
        seq: Option<Vec<(u8, BaseFmt)>>,
        qual: Option<Vec<u8>>,
    ) {
        match &mut self.0 {
            ReadInstance::OnlyAlign {
                align_len: al,
                mod_count: mc,
                read_state: rs,
                seq: sq,
                qual: q,
                seq_len: l,
            } => {
                al.push(align_len);
                rs.push(read_state);
                mc.extend(mod_count);
                sq.extend(seq);
                q.extend(qual);

                // BAM files are not supposed to have different sequence lengths for
                // the same read if there are multiple records corresponding to one read.
                // But, it is possible that some BAM files are in violation, or the sequence
                // is not stored in all reads. To account for this, we set our sequence length
                // to be the largest of all these records.
                if *l < seq_len {
                    *l = seq_len;
                }
            }
            ReadInstance::BothAlignBc {
                align_len: al,
                mod_count: mc,
                read_state: rs,
                seq: sq,
                qual: q,
                ..
            } => {
                // we ignore the seq length here, as we want the seq len from the
                // sequencing summary file to take precedence over the seq length
                // from the BAM file.
                al.push(align_len);
                rs.push(read_state);
                mc.extend(mod_count);
                sq.extend(seq);
                q.extend(qual);
            }
            ReadInstance::OnlyBc(bl) => {
                self.0 = ReadInstance::BothAlignBc {
                    align_len: vec![align_len],
                    bc_len: *bl,
                    mod_count: mod_count.into_iter().collect(),
                    read_state: vec![read_state],
                    seq: seq.into_iter().collect(),
                    qual: qual.into_iter().collect(),
                };
            }
        }
    }
}

/// Represents a record with the columns of interest from the TSV file.
#[derive(Debug, Serialize, Deserialize)]
struct TSVRecord {
    /// Read id of the molecule in the current row
    read_id: String,
    /// Sequence length of the molecule in the current row
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
                        record.read_id.clone(),
                        Read::new_bc_len(record.sequence_length_template),
                    )
                    .is_some()
                {
                    return Err(Error::InvalidDuplicates(format!(
                        "file: {file_path}, read: {0}",
                        record.read_id
                    )));
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
///
/// # Panics
/// If lists associated with sequence, modification, and/or base quality are malformed i.e.
/// not of correct lengths or missing data etc.
#[expect(clippy::too_many_lines, reason = "needs to process many options")]
pub fn run<W, D>(
    handle: &mut W,
    bam_records: D,
    mut mods: Option<InputMods<OptionalTag>>,
    seq_display: SeqDisplayOptions,
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

    match mods.as_mut() {
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
        let curr_read_state = {
            let (temp, is_non_zero_len) = match CurrRead::default().try_from_only_alignment(&record)
            {
                Ok(v) => (v, true),
                Err(Error::ZeroSeqLen(_)) => (
                    CurrRead::default().try_from_only_alignment_zero_seq_len(&record)?,
                    false,
                ),
                Err(e) => return Err(e),
            };
            let record_to_be_used = if mods.is_some() && is_non_zero_len {
                &record
            } else {
                // A new record contains no mod information.
                &bam::Record::new()
            };
            match mods.as_ref() {
                Some(v) => temp.set_mod_data_restricted_options(record_to_be_used, v)?,
                None => temp.set_mod_data(record_to_be_used, ThresholdState::GtEq(0), 0)?,
            }
        };

        let qname = String::from(curr_read_state.read_id());
        let read_state = curr_read_state.read_state();
        let align_len = match curr_read_state.align_len() {
            Ok(v) => v,
            Err(Error::Unmapped(_)) => 0,
            Err(_) => continue,
        };
        let Ok(seq_len) = curr_read_state.seq_len() else {
            continue;
        };

        // get sequence and basecalling qualities
        let (sequence, qualities) = {
            let seq = record.seq().as_bytes();
            let qual = record.qual().to_vec();
            match seq_display {
                SeqDisplayOptions::No => (None, None),
                SeqDisplayOptions::Full { show_base_qual } => (
                    Some({
                        if seq.is_empty() {
                            vec![(b'*', BaseFmt::No)]
                        } else {
                            seq.into_iter().zip(iter::repeat(BaseFmt::No)).collect()
                        }
                    }),
                    show_base_qual.then_some(if qual.is_empty() { vec![255u8] } else { qual }),
                ),
                SeqDisplayOptions::Region {
                    show_ins_lowercase,
                    show_base_qual,
                    region,
                    show_mod_z,
                } => {
                    let error_message = "no coordinate errors anticipated";
                    let (o_1, o_2): (Vec<(u8, BaseFmt)>, Vec<u8>) = {
                        let coord_map =
                            match curr_read_state.seq_coords_from_ref_coords(&record, &region) {
                                Err(Error::UnavailableData(_)) => vec![],
                                Err(e) => return Err(e),
                                Ok(x) => x,
                            };
                        if coord_map.is_empty() || seq.is_empty() {
                            (vec![(b'*', BaseFmt::No)], vec![255u8])
                        } else {
                            let mod_data =
                                match SeqCoordCalls::try_from(&curr_read_state.mod_data().0) {
                                    Err(Error::UnavailableData(_)) => vec![false; seq.len()],
                                    Err(e) => return Err(e),
                                    Ok(v) => v.collapse_mod_calls(),
                                };
                            coord_map
                                .into_iter()
                                .map(|y| {
                                    y.map_or(((b'.', BaseFmt::No), 255u8), |z| {
                                        (
                                            (
                                                *seq.get(z.1).expect(error_message),
                                                match (
                                                    show_ins_lowercase && !z.0,
                                                    show_mod_z
                                                        && *mod_data.get(z.1).expect(error_message),
                                                ) {
                                                    (true, true) => BaseFmt::LowerCaseHighlight,
                                                    (true, false) => BaseFmt::LowerCase,
                                                    (false, true) => BaseFmt::Highlight,
                                                    (false, false) => BaseFmt::No,
                                                },
                                            ),
                                            *qual.get(z.1).expect(error_message),
                                        )
                                    })
                                })
                                .collect()
                        }
                    };

                    (Some(o_1), show_base_qual.then_some(o_2))
                }
            }
        };

        // get modification information
        let mod_count = mods
            .as_ref()
            .map(|_| ModCountTbl::new(curr_read_state.base_count_per_mod()));

        // add data depending on whether an entry is already present
        // in the hashmap from the sequencing summary file
        let _: &mut _ = data_map
            .entry(qname)
            .and_modify(|entry| {
                entry.add_align_len(
                    align_len,
                    seq_len,
                    mod_count.clone(),
                    read_state,
                    sequence.clone(),
                    qualities.clone(),
                );
            })
            .or_insert(Read::new_align_len(
                align_len, seq_len, mod_count, read_state, sequence, qualities,
            ));
    }

    // print the output header
    writeln!(
        handle,
        "{}{}read_id\talign_length\tsequence_length_template\talignment_type{}{}{}",
        is_seq_summ_data
            .then_some(format!("# seq summ file: {seq_summ_path}\n"))
            .map_or(String::new(), |v| v),
        mods.clone()
            .map_or("", |_| "# mod-unmod threshold is 0.5\n"),
        mods.map_or("", |_| "\tmod_count"),
        match seq_display {
            SeqDisplayOptions::No => "",
            SeqDisplayOptions::Full { .. } | SeqDisplayOptions::Region { .. } => "\tsequence",
        },
        match seq_display {
            SeqDisplayOptions::No
            | SeqDisplayOptions::Full {
                show_base_qual: false,
            }
            | SeqDisplayOptions::Region {
                show_base_qual: false,
                ..
            } => "",
            SeqDisplayOptions::Full {
                show_base_qual: true,
            }
            | SeqDisplayOptions::Region {
                show_base_qual: true,
                ..
            } => "\tqualities",
        },
    )?;

    // print output tsv data
    // If both seq summ and BAM file are available, then the length in the seq
    // summ file takes precedence. If only the BAM file is available, then
    // we use the sequence length in the BAM file as the basecalled sequence length.
    #[expect(
        clippy::iter_over_hash_type,
        reason = "sorting would add unnecessary performance overhead; random iteration order is acceptable"
    )]
    for (key, val) in data_map {
        match (val, is_seq_summ_data) {
            (Read(ReadInstance::OnlyBc(_)), _) | (Read(ReadInstance::OnlyAlign { .. }), true) => {}
            (Read(ReadInstance::BothAlignBc { .. }), false) => {
                unreachable!("invalid state while writing output");
            }
            (Read(v), _) => writeln!(handle, "{key}\t{v}")?,
        }
    }

    Ok(())
}

/// Creates a `DataFrame` from read table data
///
/// This function calls [`run`] with a buffer handle, then parses the output into a Polars `DataFrame`.
/// Lines starting with '#' are removed as comments, and the remaining first line contains
/// tab-separated column names. Subsequent lines contain tab-separated data values.
///
/// The schema adapts to the columns present based on the options passed to [`run`]:
/// - Always present: `read_id`, `align_length`, `sequence_length_template`, `alignment_type`
/// - Optional: `mod_count` (if mods are requested)
/// - Optional: `sequence` (if `seq_display` is not `No`)
/// - Optional: `qualities` (if `seq_display` requests base quality)
///
/// # Errors
/// Returns an error if BAM record reading, output writing, or `DataFrame` construction fails.
///
#[expect(
    clippy::needless_pass_by_value,
    reason = "mods must be cloned to pass to run() which takes ownership and mutates it"
)]
pub fn run_df<D>(
    bam_records: D,
    mods: Option<InputMods<OptionalTag>>,
    seq_display: SeqDisplayOptions,
    seq_summ_path: &str,
) -> Result<DataFrame, Error>
where
    D: IntoIterator<Item = Result<Rc<bam::Record>, rust_htslib::errors::Error>>,
{
    // Create a buffer to capture output
    let mut buffer = Vec::new();

    // Call run with the buffer
    run(
        &mut buffer,
        bam_records,
        mods.clone(),
        seq_display,
        seq_summ_path,
    )?;

    // Convert buffer to string
    let output = String::from_utf8(buffer)?;

    // Remove comment lines (starting with '#') and find the header
    let lines: Vec<&str> = output
        .lines()
        .filter(|line| !line.starts_with('#'))
        .collect();

    if lines.is_empty() {
        return Err(Error::InvalidState("Output has no header line".to_string()));
    }

    // First non-comment line is the header
    let header_line = lines
        .first()
        .ok_or_else(|| Error::InvalidState("Output has no header line".to_string()))?;

    // Parse header to determine which columns are present
    let column_names: Vec<&str> = header_line.split('\t').collect();

    // Build schema dynamically based on detected columns
    let mut schema_fields = Vec::new();

    for col_name in &column_names {
        let field = match *col_name {
            "read_id" => Field::new("read_id".into(), DataType::String),
            "align_length" => Field::new("align_length".into(), DataType::String),
            "sequence_length_template" => {
                Field::new("sequence_length_template".into(), DataType::UInt64)
            }
            "alignment_type" => Field::new("alignment_type".into(), DataType::String),
            "mod_count" => Field::new("mod_count".into(), DataType::String),
            "sequence" => Field::new("sequence".into(), DataType::String),
            "qualities" => Field::new("qualities".into(), DataType::String),
            _ => {
                return Err(Error::InvalidState(format!(
                    "Unknown column name: {col_name}"
                )));
            }
        };
        schema_fields.push(field);
    }

    let schema = Schema::from_iter(schema_fields);

    // Reconstruct TSV without comment lines for parsing
    let tsv_without_comments = lines.join("\n");

    // Parse the TSV data with the schema
    let cursor = std::io::Cursor::new(tsv_without_comments.as_bytes());
    let df = CsvReadOptions::default()
        .with_has_header(true)
        .map_parse_options(|parse_options| parse_options.with_separator(b'\t'))
        .with_schema(Some(Arc::new(schema)))
        .into_reader_with_file_handle(cursor)
        .finish()?;

    Ok(df)
}

/// Sorts string that represents tabular data by first column (removing empty lines).
///
///
/// Useful in testing to compare expected and observed output.
/// Header must either start with '#' or with "`read_id`".
///
/// Example
///
/// ```
/// use nanalogue_core::reads_table::sort_output_lines;
/// let output = String::from("#comment\n\nread_id value\nbb 100\naa 50\ncc 75");
/// let expected_output = vec!["#comment","read_id value","aa 50","bb 100","cc 75"];
/// assert_eq!(expected_output, sort_output_lines(&output));
/// ```
#[must_use]
pub fn sort_output_lines(output: &str) -> Vec<String> {
    let (mut header, mut data): (Vec<String>, Vec<String>) = output
        .lines()
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .partition(|x| x.starts_with('#') | x.starts_with("read_id"));

    // our output table can return reads in any order, so we need to sort it
    // so that we can compare an expected and an observed output
    data.sort();
    header.append(&mut data);
    header
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nanalogue_bam_reader;
    use bedrs::Bed3;
    use rand::random;
    use rust_htslib::bam::Read as _;

    #[test]
    fn read_instance_only_bc_len_display() {
        assert_eq!("1000".to_owned(), ReadInstance::OnlyBc(1000u64).to_string());
    }

    fn run_read_table_test(
        bam_file: &str,
        mods: Option<InputMods<OptionalTag>>,
        seq_region: SeqDisplayOptions,
        seq_summ_file: Option<&str>,
        expected_output_file: &str,
    ) -> Result<(), Error> {
        let mut reader = nanalogue_bam_reader(bam_file).expect("no error");
        let records = reader.rc_records();

        let mut output = Vec::new();
        run(
            &mut output,
            records,
            mods,
            seq_region,
            seq_summ_file.unwrap_or(""),
        )?;

        let actual_output = String::from_utf8(output).expect("Invalid UTF-8");
        let expected_output = std::fs::read_to_string(expected_output_file)
            .expect("Failed to read expected output file");

        let actual_lines = sort_output_lines(&actual_output);
        let expected_lines = sort_output_lines(&expected_output);

        assert_eq!(
            actual_lines,
            expected_lines,
            "\nActual output:\n{}\n\nExpected output:\n{}\n",
            actual_lines.join("\n"),
            expected_lines.join("\n")
        );

        Ok(())
    }

    #[test]
    fn read_table_hide_mods() {
        run_read_table_test(
            "./examples/example_1.bam",
            None,
            SeqDisplayOptions::No,
            None,
            "./examples/example_1_read_table_hide_mods",
        )
        .expect("no error");
    }

    #[test]
    fn read_table_show_mods() {
        run_read_table_test(
            "./examples/example_1.bam",
            Some(InputMods::<OptionalTag>::default()),
            SeqDisplayOptions::No,
            None,
            "./examples/example_1_read_table_show_mods",
        )
        .expect("no error");
    }

    #[test]
    fn read_table_show_mods_but_one_record_no_mods() {
        run_read_table_test(
            "./examples/example_6.sam",
            Some(InputMods::<OptionalTag>::default()),
            SeqDisplayOptions::No,
            None,
            "./examples/example_6_read_table_show_mods",
        )
        .expect("no error");
    }

    #[test]
    fn read_table_show_mods_with_seq_summ() {
        run_read_table_test(
            "./examples/example_1.bam",
            Some(InputMods::<OptionalTag>::default()),
            SeqDisplayOptions::No,
            Some("./examples/example_1_sequencing_summary"),
            "./examples/example_1_read_table_show_mods_seq_summ",
        )
        .expect("no error");
    }

    #[test]
    #[should_panic(expected = "InvalidDuplicates")]
    fn read_table_show_mods_with_invalid_seq_summ() {
        // this sequencing summary file is invalid because a read id is repeated.
        run_read_table_test(
            "./examples/example_1.bam",
            Some(InputMods::<OptionalTag>::default()),
            SeqDisplayOptions::No,
            Some("./examples/example_1_invalid_sequencing_summary"),
            "./examples/example_1_read_table_show_mods_seq_summ",
        )
        .unwrap();
    }

    #[test]
    fn read_table_show_mods_seq_qual() {
        run_read_table_test(
            "./examples/example_5_valid_basequal.sam",
            Some(InputMods::<OptionalTag>::default()),
            SeqDisplayOptions::Full {
                show_base_qual: true,
            },
            None,
            "./examples/example_5_valid_basequal_read_table_show_mods",
        )
        .expect("no error");
    }

    #[test]
    fn read_table_show_mods_seq_qual_subset() {
        run_read_table_test(
            "./examples/example_5_valid_basequal.sam",
            Some(InputMods::<OptionalTag>::default()),
            SeqDisplayOptions::Region {
                region: Bed3::<i32, u64>::new(0, 10, 12),
                show_ins_lowercase: false,
                show_base_qual: true,
                show_mod_z: false,
            },
            None,
            "./examples/example_5_valid_basequal_read_table_show_mods_subset",
        )
        .expect("no error");
    }

    #[test]
    fn read_table_show_mods_seq_qual_subset_no_overlap() {
        // the region specified below does not overlap with the read in the sam file,
        // so we are testing read table outputs when sequence and basecalling quality
        // are requested in this scenario.
        run_read_table_test(
            "./examples/example_5_valid_basequal.sam",
            Some(InputMods::<OptionalTag>::default()),
            SeqDisplayOptions::Region {
                region: Bed3::<i32, u64>::new(1, 0, 2000),
                show_ins_lowercase: false,
                show_base_qual: true,
                show_mod_z: false,
            },
            None,
            "./examples/example_5_valid_basequal_read_table_show_mods_subset_no_overlap",
        )
        .expect("no error");
    }

    #[test]
    fn read_table_with_zero_seq_len_hide_mods() {
        run_read_table_test(
            "./examples/example_2_zero_len.sam",
            None,
            SeqDisplayOptions::No,
            None,
            "./examples/example_2_table_w_zero_hide_mods",
        )
        .expect("no error");
    }

    #[test]
    fn read_table_with_zero_seq_len_show_mods() {
        run_read_table_test(
            "./examples/example_2_zero_len.sam",
            Some(InputMods::<OptionalTag>::default()),
            SeqDisplayOptions::No,
            None,
            "./examples/example_2_table_w_zero_show_mods",
        )
        .expect("no error");
    }

    #[test]
    fn read_table_with_zero_seq_len_hide_mods_seq_summ() {
        run_read_table_test(
            "./examples/example_2_zero_len.sam",
            None,
            SeqDisplayOptions::No,
            Some("./examples/example_2_sequencing_summary"),
            "./examples/example_2_table_w_zero_hide_mods_seq_summ",
        )
        .expect("no error");
    }

    #[test]
    fn read_table_with_zero_seq_len_show_mods_seq_summ() {
        run_read_table_test(
            "./examples/example_2_zero_len.sam",
            Some(InputMods::<OptionalTag>::default()),
            SeqDisplayOptions::No,
            Some("./examples/example_2_sequencing_summary"),
            "./examples/example_2_table_w_zero_show_mods_seq_summ",
        )
        .expect("no error");
    }

    #[test]
    fn read_table_hide_mods_ins_lowercase_with_and_without() {
        run_read_table_test(
            "./examples/example_7.sam",
            None,
            SeqDisplayOptions::Region {
                region: Bed3::<i32, u64>::new(0, 0, 1000),
                show_ins_lowercase: true,
                show_base_qual: false,
                show_mod_z: false,
            },
            None,
            "./examples/example_7_table_hide_mods_ins_lowercase",
        )
        .expect("no error");

        run_read_table_test(
            "./examples/example_7.sam",
            None,
            SeqDisplayOptions::Region {
                region: Bed3::<i32, u64>::new(0, 0, 1000),
                show_ins_lowercase: false,
                show_base_qual: false,
                show_mod_z: false,
            },
            None,
            "./examples/example_7_table_hide_mods",
        )
        .expect("no error");
    }

    /// If a read has many sequence lengths (multiple records in the BAM file
    /// having mismatched lengths), and no sequencing summary file entry,
    /// test that we choose the largest length.
    #[test]
    #[expect(clippy::panic, reason = "panics are fine in tests")]
    #[expect(
        clippy::indexing_slicing,
        reason = "this is fine in tests, lists below are very short"
    )]
    fn maximum_length_replacement_test() {
        let lengths = [10, 20, 30, 40];
        let n = lengths.len();
        for count in 0..n {
            let random_state_1: ReadState = random();
            let mut read = Read::new_align_len(1, lengths[count], None, random_state_1, None, None);
            for l in 1..n {
                let random_state: ReadState = random();
                let indx = if count + l < n {
                    count + l
                } else {
                    count + l - n
                };
                read.add_align_len(1, lengths[indx], None, random_state, None, None);
            }
            match read.0 {
                ReadInstance::OnlyAlign { seq_len: 40, .. } => {}
                ReadInstance::OnlyAlign { seq_len, .. } => {
                    panic!("maximum length replacement test failed, we got {seq_len} instead of 40")
                }
                ReadInstance::BothAlignBc { .. } | ReadInstance::OnlyBc(_) => {
                    panic!("maximum length replacement test failed, wrong state")
                }
            }
        }
    }

    /// If a read has many sequence lengths (multiple records in the BAM file
    /// having mismatched lengths), and we have sequence length from the sequencing
    /// summary file, test that we retain the sequence length from the summary file
    /// no matter what.
    #[test]
    #[expect(clippy::panic, reason = "panics are fine in tests")]
    #[expect(
        clippy::indexing_slicing,
        reason = "this is fine in tests, lists below are very short"
    )]
    fn no_length_replacement_test() {
        let lengths = [10, 20, 30, 40];
        let n = lengths.len();
        for count in 0..n {
            let mut read = Read::new_bc_len(lengths[count]);
            for l in 1..n {
                let random_state: ReadState = random();
                let indx = if count + l < n {
                    count + l
                } else {
                    count + l - n
                };
                read.add_align_len(1, lengths[indx], None, random_state, None, None);
            }
            match read.0 {
                ReadInstance::BothAlignBc { bc_len, .. } if bc_len == lengths[count] => {}
                ReadInstance::BothAlignBc { bc_len, .. } => {
                    panic!(
                        "maximum length replacement test failed, we got {bc_len} instead of {}",
                        lengths[count]
                    )
                }
                ReadInstance::OnlyAlign { .. } | ReadInstance::OnlyBc(_) => {
                    panic!("maximum length replacement test failed, wrong state")
                }
            }
        }
    }

    /// Test that a record with zero sequence length outputs "*" for sequence
    /// and "255" for qualities when both sequence and qualities are requested.
    #[test]
    #[expect(
        clippy::indexing_slicing,
        reason = "this is fine in tests, parsing known output"
    )]
    fn zero_seq_len_record_with_sequence_and_qualities() {
        let mut reader = nanalogue_bam_reader("./examples/example_2_zero_len.sam")
            .expect("Failed to open example file");
        let mut bam_records = reader.rc_records();
        let first_record = bam_records.next().expect("No records in file");

        let mut output = Vec::new();
        run(
            &mut output,
            vec![first_record],
            None,
            SeqDisplayOptions::Full {
                show_base_qual: true,
            },
            "",
        )
        .expect("no error");

        let actual_output = String::from_utf8(output).expect("Invalid UTF-8");

        // Parse TSV using csv crate (similar to process_tsv)
        let mut rdr = ReaderBuilder::new()
            .comment(Some(b'#'))
            .delimiter(b'\t')
            .from_reader(actual_output.as_bytes());

        let headers = rdr.headers().expect("Failed to read headers");
        let seq_col_idx = headers
            .iter()
            .position(|h| h == "sequence")
            .expect("sequence column not found");
        let qual_col_idx = headers
            .iter()
            .position(|h| h == "qualities")
            .expect("qualities column not found");

        let mut csv_records = rdr.records();
        let record = csv_records
            .next()
            .expect("No data row found")
            .expect("Failed to parse row");

        assert_eq!(
            &record[seq_col_idx], "*",
            "Sequence column should contain '*'"
        );
        assert_eq!(
            &record[qual_col_idx], "255",
            "Qualities column should contain '255'"
        );
    }

    /// Test that a record with zero sequence length outputs "*" for sequence
    /// but no qualities column when only sequence is requested.
    #[test]
    #[expect(
        clippy::indexing_slicing,
        reason = "this is fine in tests, parsing known output"
    )]
    fn zero_seq_len_record_with_sequence_only() {
        let mut reader = nanalogue_bam_reader("./examples/example_2_zero_len.sam")
            .expect("Failed to open example file");
        let mut bam_records = reader.rc_records();
        let first_record = bam_records.next().expect("No records in file");

        let mut output = Vec::new();
        run(
            &mut output,
            vec![first_record],
            None,
            SeqDisplayOptions::Full {
                show_base_qual: false,
            },
            "",
        )
        .expect("no error");

        let actual_output = String::from_utf8(output).expect("Invalid UTF-8");

        // Parse TSV using csv crate (similar to process_tsv)
        let mut rdr = ReaderBuilder::new()
            .comment(Some(b'#'))
            .delimiter(b'\t')
            .from_reader(actual_output.as_bytes());

        let headers = rdr.headers().expect("Failed to read headers");
        let seq_col_idx = headers
            .iter()
            .position(|h| h == "sequence")
            .expect("sequence column not found");

        // Verify qualities column does not exist
        assert!(
            !headers.iter().any(|h| h == "qualities"),
            "Qualities column should not be present when show_base_qual is false"
        );

        let mut csv_records = rdr.records();
        let record = csv_records
            .next()
            .expect("No data row found")
            .expect("Failed to parse row");

        assert_eq!(
            &record[seq_col_idx], "*",
            "Sequence column should contain '*'"
        );
    }
}

#[cfg(test)]
mod stochastic_tests {
    use super::*;
    use crate::simulate_mod_bam::{
        ContigConfigBuilder, ReadConfigBuilder, SimulationConfigBuilder, TempBamSimulation,
    };
    use derive_builder::Builder;
    use rust_htslib::bam::Read as _;
    use std::ops::RangeInclusive;

    /// Helper to run reads table generation
    fn run_reads_table_generation(
        sim: &TempBamSimulation,
        mods: Option<InputMods<OptionalTag>>,
        seq_display: SeqDisplayOptions,
    ) -> Result<DataFrame, Error> {
        let mut bam_reader = bam::Reader::from_path(sim.bam_path())?;
        let bam_records = bam_reader.rc_records();
        run_df(bam_records, mods, seq_display, "")
    }

    /// Helper to keep track that all read states have been visited
    fn track_read_state_visits(val: &mut [bool; 7], state: &str) {
        match state {
            "unmapped" => val[0] = true,
            "primary_forward" => val[1] = true,
            "primary_reverse" => val[2] = true,
            "secondary_forward" => val[3] = true,
            "secondary_reverse" => val[4] = true,
            "supplementary_forward" => val[5] = true,
            "supplementary_reverse" => val[6] = true,
            &_ => unreachable!(),
        }
    }

    /// Helper struct used in the [`assert_expected_columns`] function
    /// to register what columns are expected.
    #[derive(Builder, Clone, Copy, Debug, Default, PartialEq)]
    #[builder(default, build_fn(error = "Error"))]
    struct Columns {
        /// Expect a `mod_count` column
        mod_count: bool,
        /// Expect a `sequence` column
        sequence: bool,
        /// Expect a `qualities` column
        qualities: bool,
    }

    /// Helper to assert dataframe has expected column headers
    fn assert_expected_columns(df: &DataFrame, columns: Columns) {
        let mut expected_columns = vec![
            "read_id",
            "align_length",
            "sequence_length_template",
            "alignment_type",
        ];
        if columns.mod_count {
            expected_columns.push("mod_count");
        }
        if columns.sequence {
            expected_columns.push("sequence");
        }
        if columns.qualities {
            expected_columns.push("qualities");
        }

        let actual_columns: Vec<String> = df
            .get_column_names()
            .iter()
            .map(ToString::to_string)
            .collect();
        assert_eq!(
            actual_columns, expected_columns,
            "DataFrame should have correct column headers"
        );
    }

    /// Helper struct used in the `check_seq_qual` function
    /// to register expected/observed counts of different non-upper-case A/C/G/T
    /// characters in a sequence
    #[derive(Builder, Clone, Copy, Debug, Default, PartialEq)]
    #[builder(default, build_fn(error = "Error"))]
    struct Counts {
        /// Number of occurences of '.'
        period: usize,
        /// Number of occurences of lowercase bases
        lowercase: usize,
        /// Number of occurences of 'Z'
        uppercase_z: usize,
        /// Number of occurences of 'z'
        lowercase_z: usize,
    }

    impl Counts {
        /// Increment appropriate count
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "we are ok with this in tests"
        )]
        fn increment(&mut self, value: char) -> Result<(), Error> {
            match value {
                '.' => self.period += 1,
                'A' | 'C' | 'G' | 'T' => {}
                'a' | 'c' | 'g' | 't' => {
                    self.lowercase += 1;
                }
                'z' => {
                    self.lowercase_z += 1;
                    self.lowercase += 1;
                }
                'Z' => {
                    self.uppercase_z += 1;
                }
                v => {
                    return Err(Error::InvalidSeq(format!(
                        "Invalid character {v} found in sequence!"
                    )));
                }
            }

            Ok(())
        }
    }

    /// Helper to check if sequence and quality columns are as desired
    fn check_seq_qual(
        seq_col: &ChunkedArray<StringType>,
        qual_col: &ChunkedArray<StringType>,
        seq_len: usize,
        expected_count: &Counts,
        qual_allowed: RangeInclusive<u8>,
    ) {
        let mut data_present = false;

        for k in seq_col.iter().zip(qual_col) {
            data_present = true;

            let seq = k.0.unwrap();
            assert_eq!(seq.len(), seq_len, "Sequence length mismatch");

            let qual =
                k.1.unwrap()
                    .split('.')
                    .map(|x| x.parse::<u8>().unwrap())
                    .collect::<Vec<u8>>();
            assert_eq!(qual.len(), seq_len, "Quality length mismatch");

            let mut count: Counts = Counts::default();

            for l in seq.chars().zip(qual) {
                match l {
                    ('.', 255u8) => count.increment('.').unwrap(),
                    (v, w) => {
                        count.increment(v).unwrap();
                        assert!(qual_allowed.contains(&w));
                    }
                }
            }
            assert!(
                count == *expected_count,
                "need correct number of . if . is expected"
            );
        }

        assert!(data_present, "Some data must be present in the table");
    }

    /// Simple test, produce some reads and see if we get expected statistics.
    /// Here `alignment length == sequence_length_template` and all read ids
    /// start with "0." as there is only one read group. The `alignment_type`
    /// is equally likely to be one of seven.
    #[test]
    fn run_df_simple() -> Result<(), Error> {
        // Create simulation config with no modifications
        let contig_config = ContigConfigBuilder::default()
            .number(2)
            .len_range((1000, 2000))
            .build()?;

        let read_config = ReadConfigBuilder::default()
            .number(100)
            .mapq_range((10, 20))
            .base_qual_range((30, 40))
            .len_range((0.5, 0.6));

        let sim_config = SimulationConfigBuilder::default()
            .contigs(contig_config)
            .reads(vec![read_config.build()?])
            .build()?;

        let sim = TempBamSimulation::new(sim_config)?;
        let df = run_reads_table_generation(&sim, None, SeqDisplayOptions::No)?;

        // Verify the dataframe is not empty and has expected columns
        assert!(df.height() > 0, "DataFrame should have some rows");
        assert_expected_columns(&df, ColumnsBuilder::default().build()?);

        let mut read_states = [false; 7];

        let read_id = df.column("read_id")?.str()?;
        assert!(
            read_id.iter().all(|k| k.unwrap().starts_with("0.")),
            "only one RG, must start with 0"
        );

        let align_length = df.column("align_length")?.str()?;
        let seq_length = df.column("sequence_length_template")?.u64()?;
        let alignment_type = df.column("alignment_type")?.str()?;

        for k in align_length.iter().zip(seq_length).zip(alignment_type) {
            let al = k.0.0.unwrap().parse::<u64>().unwrap();
            let sl = k.0.1.unwrap();
            let rs = k.1.unwrap();
            assert!(al == 0 || al == sl);
            if al == 0 {
                assert_eq!(rs, "unmapped");
            }
            track_read_state_visits(&mut read_states, rs);
            assert!((500..=1200).contains(&sl));
        }

        assert!(read_states.into_iter().all(|k| k));

        Ok(())
    }

    /// More complex, now sequence and alignment lengths are systematically
    /// different due to an insertion and a barcode.
    #[test]
    fn run_df_al_sl_different() -> Result<(), Error> {
        // Create simulation config with no modifications
        let contig_config = ContigConfigBuilder::default()
            .number(2)
            .len_range((1000, 2000))
            .build()?;

        let read_config = ReadConfigBuilder::default()
            .number(100)
            .mapq_range((10, 20))
            .base_qual_range((30, 40))
            .len_range((0.5, 0.6))
            .barcode("AATTGAA".into())
            .insert_middle("GGTT".into());

        let sim_config = SimulationConfigBuilder::default()
            .contigs(contig_config)
            .reads(vec![read_config.build()?])
            .build()?;

        let sim = TempBamSimulation::new(sim_config)?;
        let df = run_reads_table_generation(&sim, None, SeqDisplayOptions::No)?;

        // Verify the dataframe is not empty and has expected columns
        assert!(df.height() > 0, "DataFrame should have some rows");
        assert_expected_columns(&df, ColumnsBuilder::default().build()?);

        let mut read_states = [false; 7];

        let read_id = df.column("read_id")?.str()?;
        assert!(
            read_id.iter().all(|k| k.unwrap().starts_with("0.")),
            "only one RG, must start with 0"
        );

        let align_length = df.column("align_length")?.str()?;
        let seq_length = df.column("sequence_length_template")?.u64()?;
        let alignment_type = df.column("alignment_type")?.str()?;

        for k in align_length.iter().zip(seq_length).zip(alignment_type) {
            let al = k.0.0.unwrap().parse::<u64>().unwrap();
            let sl = k.0.1.unwrap();
            let rs = k.1.unwrap();
            assert!(al == 0 || al == sl - 18); // barcode is 7 bp, insertion is 4 bp, so 2 * 7 + 4
            if al == 0 {
                assert_eq!(rs, "unmapped");
            } else {
                assert!((500..=1200).contains(&al));
            }
            track_read_state_visits(&mut read_states, rs);
        }

        assert!(read_states.into_iter().all(|k| k));

        Ok(())
    }

    /// Test retrieval of sequence and qualities
    #[test]
    fn run_df_seq_qual_retrieval() -> Result<(), Error> {
        // Create simulation config with no modifications
        let contig_config = ContigConfigBuilder::default()
            .number(2)
            .len_range((1000, 1000))
            .build()?;

        let read_config = ReadConfigBuilder::default()
            .number(100)
            .mapq_range((10, 20))
            .base_qual_range((71, 79))
            .len_range((0.5, 0.5));

        let sim_config = SimulationConfigBuilder::default()
            .contigs(contig_config)
            .reads(vec![read_config.build()?])
            .build()?;

        let sim = TempBamSimulation::new(sim_config)?;
        let df = run_reads_table_generation(
            &sim,
            None,
            SeqDisplayOptions::Full {
                show_base_qual: true,
            },
        )?;

        // Verify the dataframe has expected columns
        assert_expected_columns(
            &df,
            ColumnsBuilder::default()
                .sequence(true)
                .qualities(true)
                .build()?,
        );

        let seq_col = df.column("sequence")?.str()?;
        let qual_col = df.column("qualities")?.str()?;

        check_seq_qual(seq_col, qual_col, 500, &Counts::default(), 71u8..=79u8);

        let df_2 = run_reads_table_generation(
            &sim,
            None,
            SeqDisplayOptions::Full {
                show_base_qual: false,
            },
        )?;

        // Verify the dataframe is not empty and has expected columns
        assert!(df_2.height() > 0, "DataFrame should have some rows");
        assert_expected_columns(&df_2, ColumnsBuilder::default().sequence(true).build()?);

        // extract seq but we won't be getting a column called qualities
        let seq_col_2 = df_2.column("sequence")?.str()?;
        drop(df_2.column("qualities").unwrap_err());

        for k in seq_col_2.iter() {
            assert_eq!(k.unwrap().len(), 500);
        }

        Ok(())
    }

    /// Test retrieval of sequence and qualities of the full
    /// sequence even with indels and barcode.
    #[test]
    fn run_df_seq_qual_retrieval_indels_barcode() -> Result<(), Error> {
        // Create simulation config with no modifications
        let contig_config = ContigConfigBuilder::default()
            .number(2)
            .len_range((1000, 1000))
            .build()?;

        let read_config = ReadConfigBuilder::default()
            .number(100)
            .mapq_range((10, 20))
            .base_qual_range((71, 79))
            .len_range((0.5, 0.5))
            .barcode("AATT".into())
            .delete((0.5, 0.6))
            .insert_middle("GGTTGG".into())
            .mismatch(0.1);

        let sim_config = SimulationConfigBuilder::default()
            .contigs(contig_config)
            .reads(vec![read_config.build()?])
            .build()?;

        let sim = TempBamSimulation::new(sim_config)?;
        let df = run_reads_table_generation(
            &sim,
            None,
            SeqDisplayOptions::Full {
                show_base_qual: true,
            },
        )?;

        // Verify the dataframe has expected columns
        assert_expected_columns(
            &df,
            ColumnsBuilder::default()
                .sequence(true)
                .qualities(true)
                .build()?,
        );

        let seq_col = df.column("sequence")?.str()?;
        let qual_col = df.column("qualities")?.str()?;

        check_seq_qual(seq_col, qual_col, 464, &Counts::default(), 71u8..=79u8);

        Ok(())
    }
}
