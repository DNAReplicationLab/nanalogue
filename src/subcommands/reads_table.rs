//! # Makes a table of information on reads
//!
//! This file contains routines used to calculate some information per
//! read such as sequence length, alignment length, type of alignment,
//! and modification counts, and displays them in table using the run
//! function. The routine reads both BAM and sequencing summary files
//! if provided, otherwise only reads the BAM file.

use crate::constants::shared::{MAX_RECORD_CAPACITY_BYTES, MAX_RECORDS};
use crate::{
    CurrRead, Error, InputMods, ModChar, OptionalTag, ReadState, SeqCoordCalls, SeqDisplayOptions,
    ThresholdState, assert_bounded_counter, assert_nonzero_counter, assert_record_data_capacity,
};
use polars::prelude::*;
use rust_htslib::bam;
use std::{collections::HashMap, fmt, fs::File, io::BufReader, io::Read as _, iter, rc::Rc, str};

/// Write an iterator as a separated string.
fn join_display<I>(items: I, separator: &str) -> Result<String, fmt::Error>
where
    I: IntoIterator,
    I::Item: fmt::Display,
{
    let mut output = String::new();
    for (index, item) in items.into_iter().enumerate() {
        if index > 0 {
            output.push_str(separator);
        }
        // we cap how many items can be joined this way
        if index
            > u32::MAX
                .try_into()
                .expect("no error expected on 32-bit platforms and above")
        {
            return Err(fmt::Error);
        }
        output.push_str(item.to_string().as_str());
    }
    Ok(output)
}

/// Write a vector as a CSV-formatted string
macro_rules! vec_csv {
    ( $b: expr ) => {
        join_display($b, ", ")?
    };
}

/// Read one sequencing summary line with a hard raw-byte cap.
#[expect(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "(1) line & buffer lengths are checked for smallness, \
(2) `bytes_to_consume` never exceeds buffer length and indexing only happens when it is > 0"
)]
fn read_seq_summ_line<R: std::io::BufRead>(
    reader: &mut R,
    line: &mut String,
) -> Result<usize, Error> {
    const MAX_SEQ_SUMM_RAW_LINE_BYTES: usize = 1002;

    line.clear();
    let mut total_bytes_read: usize = 0;

    loop {
        // With `std::io::BufReader`, this buffered slice is typically modest in
        // size rather than a huge chunk of memory, so using `fill_buf` here is
        // acceptable for our defensive line-length checks.
        let buffered = reader.fill_buf()?;
        if buffered.len() > 1000 * MAX_SEQ_SUMM_RAW_LINE_BYTES {
            return Err(Error::InvalidState(
                "sequencing summary tsv internal buffer unexpectedly large".to_owned(),
            ));
        }

        let (bytes_to_consume, is_newline_found) =
            match buffered.iter().position(|byte| *byte == b'\n') {
                Some(v) => (v + 1, true),
                None => (buffered.len(), false),
            };
        match bytes_to_consume {
            0 => break,
            v => {
                total_bytes_read += v;

                if total_bytes_read > MAX_SEQ_SUMM_RAW_LINE_BYTES {
                    return Err(Error::InvalidState(
                        "sequencing summary tsv line is too long (>1000 bytes)".to_owned(),
                    ));
                }

                let idx = if is_newline_found { v - 1 } else { v };
                line.push_str(str::from_utf8(&buffered[..idx])?);
                reader.consume(v);
                if is_newline_found {
                    break;
                }
            }
        }
    }

    if line.ends_with('\r') {
        let _: Option<char> = line.pop();
    }

    Ok(total_bytes_read)
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
            Ok("NA".to_owned())
        } else {
            join_display(v.into_iter().map(|k| format!("{}:{}", k.0, k.1)), ";")
        })?
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
        align_len: Vec<u32>,
        /// Sequence length from the BAM file
        seq_len: u32,
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
    OnlyBc(u32),
    /// Information from both sequencing summary and alignment file are available.
    /// NOTE that these are vectors as a read can have multiple BAM records,
    /// but the basecalled length is not a vector as a sequencing summary file
    /// has unique information per read. If the basecalled length from the
    /// sequencing summary file does not match that from the BAM file record(s),
    /// then we will have to make a choice.
    BothAlignBc {
        /// Alignment length from the BAM file
        align_len: Vec<u32>,
        /// Sequence length from the sequencing summary file
        bc_len: u32,
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
        // we are fine with `unsafe` here (`from_utf8_unchecked`) as `rust-htslib`
        // guarantees that only 16 printable characters are allowed in the sequence.
        match self {
            ReadInstance::OnlyBc(v) => format!("{v}"),
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
            } => {
                // assertions to enforce length constraints
                if (al.len() != rs.len())
                    || (!mc.is_empty() && mc.len() != al.len())
                    || (!sq.is_empty() && sq.len() != al.len())
                {
                    return Err(fmt::Error);
                }
                if !q.is_empty() {
                    if sq.len() != q.len() {
                        return Err(fmt::Error);
                    }
                    for k in sq.iter().zip(q) {
                        if k.0.len() != k.1.len() {
                            return Err(fmt::Error);
                        }
                    }
                }
                format!(
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
                    match vec_csv!(
                        q.iter()
                            .map(|k| join_display(k, "."))
                            .collect::<Result<Vec<String>, _>>()?
                    ) {
                        v if !v.is_empty() => format!("\t{v}"),
                        _ => String::new(),
                    },
                )
            }
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
    fn new_bc_len(l: u32) -> Self {
        Self(ReadInstance::OnlyBc(l))
    }

    /// initialization using alignment information
    fn new_align_len(
        align_len: u32,
        seq_len: u32,
        mod_count: Option<ModCountTbl>,
        read_state: ReadState,
        seq: Option<Vec<(u8, BaseFmt)>>,
        qual: Option<Vec<u8>>,
    ) -> Result<Self, Error> {
        let seq_cast: Vec<Vec<(u8, BaseFmt)>> = seq.into_iter().collect();
        let qual_cast: Vec<Vec<u8>> = qual.into_iter().collect();

        if qual_cast.is_empty() || (seq_cast.len() == qual_cast.len()) {
            Ok(Self(ReadInstance::OnlyAlign {
                align_len: vec![align_len],
                seq_len,
                mod_count: mod_count.into_iter().collect(),
                read_state: vec![read_state],
                seq: seq_cast,
                qual: qual_cast,
            }))
        } else {
            Err(Error::InvalidState(
                "seq and qual are of different lengths!".to_owned(),
            ))
        }
    }

    /// adding alignment information to an entry
    #[expect(
        clippy::pattern_type_mismatch,
        reason = "expect no confusion from this code"
    )]
    fn add_align_len(
        &mut self,
        align_len: u32,
        seq_len: u32,
        mod_count: Option<ModCountTbl>,
        read_state: ReadState,
        seq: Option<Vec<(u8, BaseFmt)>>,
        qual: Option<Vec<u8>>,
    ) -> Result<(), Error> {
        let seq_cast: Vec<Vec<(u8, BaseFmt)>> = seq.into_iter().collect();
        let qual_cast: Vec<Vec<u8>> = qual.into_iter().collect();

        if !qual_cast.is_empty() && (seq_cast.len() != qual_cast.len()) {
            return Err(Error::InvalidState(
                "seq and qual are of different lengths!".to_owned(),
            ));
        }
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
                sq.extend(seq_cast);
                q.extend(qual_cast);

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
                sq.extend(seq_cast);
                q.extend(qual_cast);
            }
            ReadInstance::OnlyBc(bl) => {
                self.0 = ReadInstance::BothAlignBc {
                    align_len: vec![align_len],
                    bc_len: *bl,
                    mod_count: mod_count.into_iter().collect(),
                    read_state: vec![read_state],
                    seq: seq_cast.into_iter().collect(),
                    qual: qual_cast.into_iter().collect(),
                };
            }
        }
        Ok(())
    }
}

/// Opens a TSV file, extracts '`read_id`' and '`sequence_length_template`' columns,
/// and builds a `HashMap`.
///
/// # Errors
/// Returns an error if TSV processing fails. Note that `#`-prefixed
/// lines are only treated as comments before the header row; after the
/// header they are parsed as data rows rather than skipped.
#[expect(
    clippy::too_many_lines,
    reason = "defensive checks for malformed sequencing summary input push the line count high"
)]
fn process_seq_summ(file_path: &str) -> Result<HashMap<String, Read>, Error> {
    const MAX_SEQ_SUMM_BYTES: u64 = 10u64 * 1024u64 * 1024u64 * 1024u64;

    let mut data_map = HashMap::<String, Read>::new();

    match file_path {
        "" => {}
        fp => {
            let file = File::open(fp)?;
            let mut reader = BufReader::new(file.take(MAX_SEQ_SUMM_BYTES + 1));
            let mut line = String::with_capacity(1002);

            let header = {
                let mut bounds_checker: u16 = 0;
                loop {
                    let bytes_read = read_seq_summ_line(&mut reader, &mut line)?;
                    if bytes_read == 0 {
                        if reader.get_ref().limit() == 0 {
                            return Err(Error::InvalidState(
                                "sequencing summary file too large (> 10 GB). Contact developer."
                                    .to_owned(),
                            ));
                        }
                        return Err(Error::InvalidState(
                            "sequencing summary tsv missing header row".to_owned(),
                        ));
                    }
                    if line.starts_with('#') {
                        bounds_checker = bounds_checker.saturating_add(1);
                        if bounds_checker > 1000 {
                            return Err(Error::InvalidState(
                                "sequencing summary tsv has >1000 comment lines before the header row"
                                    .to_owned(),
                            ));
                        }
                        continue;
                    }
                    if line.trim().is_empty() {
                        return Err(Error::InvalidState(
                            "sequencing summary tsv contains a whitespace-only line before the header row"
                                .to_owned(),
                        ));
                    }
                    break line.clone();
                }
            };
            let (read_id_idx, seq_len_idx): (usize, usize) = {
                let columns: Vec<&str> = header.split('\t').collect();
                let read_id_indexes: Vec<usize> = columns
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, col)| (*col == "read_id").then_some(idx))
                    .take(2)
                    .collect();
                let seq_len_indexes: Vec<usize> = columns
                    .iter()
                    .enumerate()
                    .filter_map(|(idx, col)| (*col == "sequence_length_template").then_some(idx))
                    .take(2)
                    .collect();
                let read_id_idx = match (read_id_indexes.first(), read_id_indexes.len()) {
                    (None, _) => {
                        return Err(Error::InvalidState(
                            "sequencing summary tsv missing `read_id` column".to_owned(),
                        ));
                    }
                    (Some(idx), 1) => *idx,
                    (Some(_), _) => {
                        return Err(Error::InvalidState(
                            "sequencing summary tsv contains multiple `read_id` columns".to_owned(),
                        ));
                    }
                };
                let seq_len_idx = match (seq_len_indexes.first(), seq_len_indexes.len()) {
                    (None, _) => {
                        return Err(Error::InvalidState(
                            "sequencing summary tsv missing `sequence_length_template` column"
                                .to_owned(),
                        ));
                    }
                    (Some(idx), 1) => *idx,
                    (Some(_), _) => {
                        return Err(Error::InvalidState(
                            "sequencing summary tsv contains multiple `sequence_length_template` columns"
                                .to_owned(),
                        ));
                    }
                };
                if read_id_idx == seq_len_idx {
                    return Err(Error::InvalidState(
                        "impossible state: read id and sequence length columns are the same"
                            .to_owned(),
                    ));
                }
                (read_id_idx, seq_len_idx)
            };

            let max_col_idx_of_interest = if read_id_idx < seq_len_idx {
                seq_len_idx
            } else {
                read_id_idx
            };

            loop {
                let bytes_read = read_seq_summ_line(&mut reader, &mut line)?;
                if bytes_read == 0 {
                    if reader.get_ref().limit() == 0 {
                        return Err(Error::InvalidState(
                            "sequencing summary file too large (> 10 GB). Contact developer."
                                .to_owned(),
                        ));
                    }
                    break;
                }
                if line.is_empty() {
                    continue;
                }
                let mut read_id_field: Option<&str> = None;
                let mut seq_len_field: Option<&str> = None;
                for (idx, field) in line.split('\t').enumerate() {
                    if idx == read_id_idx {
                        read_id_field = Some(field);
                    }
                    if idx == seq_len_idx {
                        seq_len_field = Some(field);
                    }
                    if idx == max_col_idx_of_interest {
                        break;
                    }
                }
                let read_id = read_id_field.ok_or_else(|| {
                    Error::InvalidState(format!(
                        "sequencing summary tsv row has no `read_id` field at column {read_id_idx}: `{line}`"
                    ))
                })?;
                if read_id.contains(['\'', '"', '`']) {
                    return Err(Error::InvalidReadID(format!(
                        "sequencing summary read_id contains unsupported quoting characters: `{read_id}`"
                    )));
                }
                // Long read ids take up disk space and memory. If we process a million
                // reads all with a 200 byte read id, we have to allocate a lot of memory.
                // So, we cap read names at 50 bytes.
                if read_id.len() > 50 {
                    return Err(Error::InvalidState(
                        "in sequencing summary tsv, read ids longer than 50 bytes encountered"
                            .to_owned(),
                    ));
                }
                let sequence_length_template: u32 = seq_len_field
                    .ok_or_else(|| {
                        Error::InvalidState(format!(
                            "sequencing summary tsv row has no `sequence_length_template` field at column {seq_len_idx}: `{line}`"
                        ))
                    })?
                    .parse()
                    .map_err(|err| {
                        Error::InvalidState(format!(
                            "sequencing summary tsv parse error in file `{file_path}` for read `{read_id}`: invalid `sequence_length_template` ({err})"
                        ))
                    })?;
                if data_map
                    .insert(
                        (*read_id).to_owned(),
                        Read::new_bc_len(sequence_length_template),
                    )
                    .is_some()
                {
                    return Err(Error::InvalidDuplicates(format!(
                        "file: {file_path}, read: {read_id}"
                    )));
                }
            }

            if data_map.is_empty() {
                return Err(Error::InvalidState(
                    "sequencing summary TSV did not contain any reads".to_owned(),
                ));
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
/// This command also errors on blank BAM files because otherwise sparse or empty output would be
/// ambiguous between "no rows matched" and "no reads were present".
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
    let mut data_map = process_seq_summ(seq_summ_path)?;

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

    let mut idx: u32 = 0;

    // Go record by record in the BAM file,
    // get the read id and the alignment length, and put it in the hashmap
    for r in bam_records {
        // read records
        let record = r?;
        assert_bounded_counter(&mut idx, MAX_RECORDS, "reads table")?;
        assert_record_data_capacity(
            record.inner().m_data,
            MAX_RECORD_CAPACITY_BYTES,
            "reads table",
        )?;

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
        let (sequence, qualities) = match seq_display {
            SeqDisplayOptions::No => (None, None),
            SeqDisplayOptions::Full { show_base_qual } => {
                let seq = record.seq().as_bytes();
                let sequence = if seq.is_empty() {
                    vec![(b'*', BaseFmt::No)]
                } else {
                    seq.into_iter().zip(iter::repeat(BaseFmt::No)).collect()
                };
                let qualities = show_base_qual.then(|| {
                    let qual = record.qual().to_vec();
                    if qual.is_empty() { vec![255u8] } else { qual }
                });
                if let Some(quality_values) = qualities.as_ref()
                    && sequence.len() != quality_values.len()
                {
                    return Err(Error::InvalidState(
                        "sequence and qualities are not of equal length".to_owned(),
                    ));
                }
                (Some(sequence), qualities)
            }
            SeqDisplayOptions::Region {
                show_ins_lowercase,
                show_base_qual,
                region,
                show_mod_z,
            } => {
                let error_message = "no coordinate errors anticipated";
                let seq = record.seq().as_bytes();
                let coord_map = match curr_read_state.seq_coords_from_ref_coords(&record, &region) {
                    Err(Error::UnavailableData(_)) => vec![],
                    Err(e) => return Err(e),
                    Ok(x) => x,
                };
                if coord_map.is_empty() || seq.is_empty() {
                    (
                        Some(vec![(b'*', BaseFmt::No)]),
                        show_base_qual.then_some(vec![255u8]),
                    )
                } else {
                    let qual = record.qual();
                    if seq.len() != qual.len() {
                        return Err(Error::InvalidState(
                            "sequence and qualities are not of equal length".to_owned(),
                        ));
                    }
                    let mod_data = match SeqCoordCalls::try_from(&curr_read_state.mod_data().0) {
                        Err(Error::UnavailableData(_)) => vec![false; seq.len()],
                        Err(e) => return Err(e),
                        Ok(v) => v.collapse_mod_calls(),
                    };
                    let (sequence, qualities): (Vec<(u8, BaseFmt)>, Vec<u8>) = coord_map
                        .into_iter()
                        .map(|y| {
                            y.map_or(((b'.', BaseFmt::No), 255u8), |z| {
                                (
                                    (
                                        *seq.get(z.1).expect(error_message),
                                        match (
                                            show_ins_lowercase && !z.0,
                                            show_mod_z && *mod_data.get(z.1).expect(error_message),
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
                        .collect();

                    (Some(sequence), show_base_qual.then_some(qualities))
                }
            }
        };

        // get modification information
        let mod_count = mods
            .as_ref()
            .map(|_| ModCountTbl::new(curr_read_state.base_count_per_mod()));

        // add data depending on whether an entry is already present
        // in the hashmap from the sequencing summary file
        let mut is_insert_error = false;
        let _: &mut _ = data_map
            .entry(qname)
            .and_modify(|entry| {
                if entry
                    .add_align_len(
                        align_len,
                        seq_len,
                        mod_count.clone(),
                        read_state,
                        sequence.clone(),
                        qualities.clone(),
                    )
                    .is_err()
                {
                    is_insert_error = true;
                }
            })
            .or_insert(Read::new_align_len(
                align_len, seq_len, mod_count, read_state, sequence, qualities,
            )?);
        if is_insert_error {
            return Err(Error::InvalidState(
                "invalid state encountered while populating sequence data".to_owned(),
            ));
        }
    }

    assert_nonzero_counter(idx, "records")?;

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

/// Removes comment lines from read-table output and returns the first
/// remaining line as the header along with the reconstructed TSV body.
fn strip_comment_lines(output: &str) -> (Option<&str>, String) {
    let lines: Vec<&str> = output
        .lines()
        .filter(|line| !line.starts_with('#'))
        .collect();
    let header_line = lines.first().copied();
    (header_line, lines.join("\n"))
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

    let (header_line_opt, tsv_without_comments) = strip_comment_lines(&output);

    let header_line = header_line_opt
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
                Field::new("sequence_length_template".into(), DataType::UInt32)
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
    use crate::{GenomicBed3, nanalogue_bam_reader};
    use rand::random;
    use rust_htslib::bam::Read as _;

    #[test]
    fn read_instance_only_bc_len_display() {
        assert_eq!("1000".to_owned(), ReadInstance::OnlyBc(1000u32).to_string());
    }

    #[test]
    fn strip_comment_lines_keeps_no_trailing_newline() {
        let output = "#comment\nread_id\talign_length\nread1\t42";
        let (header_line, tsv_without_comments) = strip_comment_lines(output);

        assert_eq!(header_line, Some("read_id\talign_length"));
        assert_eq!(tsv_without_comments, "read_id\talign_length\nread1\t42");
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
                region: GenomicBed3::new(0, 10, 12),
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
                region: GenomicBed3::new(1, 0, 2000),
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
                region: GenomicBed3::new(0, 0, 1000),
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
                region: GenomicBed3::new(0, 0, 1000),
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
            let mut read = Read::new_align_len(1, lengths[count], None, random_state_1, None, None)
                .expect("no error");
            for l in 1..n {
                let random_state: ReadState = random();
                let indx = if count + l < n {
                    count + l
                } else {
                    count + l - n
                };
                read.add_align_len(1, lengths[indx], None, random_state, None, None)
                    .expect("no error");
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
                read.add_align_len(1, lengths[indx], None, random_state, None, None)
                    .expect("no error");
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

        let lines: Vec<&str> = actual_output.lines().collect();
        let headers: Vec<&str> = lines
            .first()
            .expect("header line should be present")
            .split('\t')
            .collect();
        let seq_col_idx = headers
            .iter()
            .position(|h| *h == "sequence")
            .expect("sequence column not found");
        let qual_col_idx = headers
            .iter()
            .position(|h| *h == "qualities")
            .expect("qualities column not found");
        let record: Vec<&str> = lines
            .get(1)
            .expect("data row should be present")
            .split('\t')
            .collect();

        assert_eq!(
            record[seq_col_idx], "*",
            "Sequence column should contain '*'"
        );
        assert_eq!(
            record[qual_col_idx], "255",
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

        let lines: Vec<&str> = actual_output.lines().collect();
        let headers: Vec<&str> = lines
            .first()
            .expect("header line should be present")
            .split('\t')
            .collect();
        let seq_col_idx = headers
            .iter()
            .position(|h| *h == "sequence")
            .expect("sequence column not found");

        // Verify qualities column does not exist
        assert!(
            !headers.contains(&"qualities"),
            "Qualities column should not be present when show_base_qual is false"
        );

        let record: Vec<&str> = lines
            .get(1)
            .expect("data row should be present")
            .split('\t')
            .collect();

        assert_eq!(
            record[seq_col_idx], "*",
            "Sequence column should contain '*'"
        );
    }
}

#[cfg(test)]
mod sequencing_summary_tests {
    use super::*;
    use crate::uuid;
    use std::fs;
    use std::path::{Path, PathBuf};

    fn write_temp_seq_summary(contents: &str) -> PathBuf {
        let temp_path =
            std::env::temp_dir().join(format!("reads_table_seq_summary_{}.tsv", uuid::v4_random()));
        fs::write(&temp_path, contents).expect("should write temporary sequencing summary file");
        temp_path
    }

    fn remove_temp_file(path: &Path) {
        fs::remove_file(path).expect("should remove temporary sequencing summary file");
    }

    #[test]
    fn process_seq_summ_rejects_invalid_preheader_and_empty_input_cases() {
        let test_cases = [
            ("", "blank file should fail"),
            (" ", "whitespace-only file should fail"),
            ("\t", "tab-only file should fail"),
            ("\n\n", "multiple blank lines before header should fail"),
            (
                "# sequencing summary\n\n",
                "comment followed by blank line should fail",
            ),
            (
                "# sequencing summary\n   \n",
                "comment followed by whitespace-only line should fail",
            ),
            (
                "# sequencing summary\n\t\n",
                "comment followed by tab-only line should fail",
            ),
            (
                "# sequencing summary\n# still comments\n",
                "comments followed by eof should fail",
            ),
            (
                "# sequencing summary\n\n# still comments\nread_id\tsequence_length_template\nread-1\t1234\n",
                "blank line amidst pre-header comments should still fail even if valid data follows",
            ),
            (
                "# sequencing summary\nread_id\tsequence_length_template\n",
                "header without data rows should fail",
            ),
        ];

        for (contents, message) in test_cases {
            let temp_path = write_temp_seq_summary(contents);
            let result = process_seq_summ(
                temp_path
                    .to_str()
                    .expect("temporary file path should be valid UTF-8"),
            );
            remove_temp_file(&temp_path);
            assert!(result.is_err(), "{message}");
        }
    }

    #[test]
    fn process_seq_summ_accepts_header_with_data_row() {
        let temp_path = write_temp_seq_summary("read_id\tsequence_length_template\nread-1\t1234\n");
        let result = process_seq_summ(
            temp_path
                .to_str()
                .expect("temporary file path should be valid UTF-8"),
        );
        remove_temp_file(&temp_path);

        let data_map = result.expect("header plus data row should parse successfully");
        assert_eq!(data_map.len(), 1);
        assert!(matches!(
            data_map.get("read-1"),
            Some(Read(ReadInstance::OnlyBc(1234)))
        ));
    }

    #[test]
    fn process_seq_summ_allows_blank_lines_after_header() {
        let temp_path = write_temp_seq_summary(
            "read_id\tsequence_length_template\n\nread-1\t1234\n\nread-2\t5678\n",
        );
        let result = process_seq_summ(
            temp_path
                .to_str()
                .expect("temporary file path should be valid UTF-8"),
        );
        remove_temp_file(&temp_path);

        let data_map = result.expect("blank lines after header should be ignored");
        assert_eq!(data_map.len(), 2);
        assert!(matches!(
            data_map.get("read-1"),
            Some(Read(ReadInstance::OnlyBc(1234)))
        ));
        assert!(matches!(
            data_map.get("read-2"),
            Some(Read(ReadInstance::OnlyBc(5678)))
        ));
    }

    #[test]
    fn process_seq_summ_rejects_quoted_or_too_long_read_ids() {
        let failure_cases = [
            (
                "read_id\tsequence_length_template\n'read-1'\t1234\n",
                "single-quoted read_id should fail",
            ),
            (
                "read_id\tsequence_length_template\n\"read-1\"\t1234\n",
                "double-quoted read_id should fail",
            ),
            (
                "read_id\tsequence_length_template\n`read-1`\t1234\n",
                "backtick-quoted read_id should fail",
            ),
            (
                "read_id\tsequence_length_template\n123456789012345678901234567890123456789012345678901\t1234\n",
                "read_id longer than 50 bytes should fail",
            ),
        ];

        for (contents, message) in failure_cases {
            let temp_path = write_temp_seq_summary(contents);
            let result = process_seq_summ(
                temp_path
                    .to_str()
                    .expect("temporary file path should be valid UTF-8"),
            );
            remove_temp_file(&temp_path);
            assert!(result.is_err(), "{message}");
        }
    }

    #[test]
    // also checks that line starting with # in seq summ after header not ignored
    fn process_seq_summ_accepts_hash_prefixed_read_id() {
        let temp_path = write_temp_seq_summary("read_id\tsequence_length_template\n#read1\t1234\n");
        let result = process_seq_summ(
            temp_path
                .to_str()
                .expect("temporary file path should be valid UTF-8"),
        );
        remove_temp_file(&temp_path);

        let data_map = result.expect("read_id containing `#` should parse successfully");
        assert!(matches!(
            data_map.get("#read1"),
            Some(Read(ReadInstance::OnlyBc(1234)))
        ));
    }

    #[test]
    fn process_seq_summ_rejects_duplicate_read_id() {
        let temp_path = write_temp_seq_summary(
            "read_id\tsequence_length_template\nread-1\t1234\nread-1\t5678\n",
        );
        let result = process_seq_summ(
            temp_path
                .to_str()
                .expect("temporary file path should be valid UTF-8"),
        );
        remove_temp_file(&temp_path);

        assert!(result.is_err(), "duplicate read_id should fail");
    }

    #[test]
    fn process_seq_summ_requires_expected_columns() {
        let failure_cases = [
            (
                "foo\tbar\nvalue1\tvalue2\n",
                "unrelated columns should fail",
            ),
            ("read_id\nread-1\n", "header with only read_id should fail"),
            (
                "sequence_length_template\n1234\n",
                "header with only sequence_length_template should fail",
            ),
            (
                "read_id\tread_id\tsequence_length_template\nread-1\tread-1\t1234\n",
                "duplicate read_id columns should fail",
            ),
            (
                "read_id\tsequence_length_template\tsequence_length_template\nread-1\t1234\t1234\n",
                "duplicate sequence_length_template columns should fail",
            ),
            (
                "read_id\tread_id\tsequence_length_template\tsequence_length_template\nread-1\tread-1\t1234\t1234\n",
                "duplicate required columns should fail when both are repeated",
            ),
        ];

        for (contents, message) in failure_cases {
            let temp_path = write_temp_seq_summary(contents);
            let result = process_seq_summ(
                temp_path
                    .to_str()
                    .expect("temporary file path should be valid UTF-8"),
            );
            remove_temp_file(&temp_path);
            assert!(result.is_err(), "{message}");
        }
    }

    #[test]
    fn process_seq_summ_rejects_rows_with_missing_required_fields() {
        let failure_cases = [
            (
                "read_id\tsequence_length_template\nread-1\n",
                "row missing sequence_length_template should fail",
            ),
            (
                "sequence_length_template\tread_id\n1234\n",
                "row missing read_id should fail when columns are reordered",
            ),
        ];

        for (contents, message) in failure_cases {
            let temp_path = write_temp_seq_summary(contents);
            let result = process_seq_summ(
                temp_path
                    .to_str()
                    .expect("temporary file path should be valid UTF-8"),
            );
            remove_temp_file(&temp_path);
            assert!(result.is_err(), "{message}");
        }
    }

    #[test]
    fn process_seq_summ_accepts_additional_columns() {
        let temp_path = write_temp_seq_summary(
            "channel\tread_id\textra\tsequence_length_template\n1\tread-1\tfoo\t1234\n",
        );
        let result = process_seq_summ(
            temp_path
                .to_str()
                .expect("temporary file path should be valid UTF-8"),
        );
        remove_temp_file(&temp_path);

        let data_map = result.expect("additional columns should be ignored");
        assert_eq!(data_map.len(), 1);
        assert!(matches!(
            data_map.get("read-1"),
            Some(Read(ReadInstance::OnlyBc(1234)))
        ));
    }

    #[test]
    fn process_seq_summ_rejects_overlong_header_line() {
        let filler_columns = iter::repeat_n("extra_column", 82)
            .collect::<Vec<_>>()
            .join("\t");
        let contents = format!(
            "read_id\tsequence_length_template\t{filler_columns}\nread-1\t1234\t{}\n",
            iter::repeat_n("value", 82).collect::<Vec<_>>().join("\t")
        );
        let temp_path = write_temp_seq_summary(contents.as_str());
        let result = process_seq_summ(
            temp_path
                .to_str()
                .expect("temporary file path should be valid UTF-8"),
        );
        remove_temp_file(&temp_path);

        assert!(
            result.is_err(),
            "header longer than 1000 bytes should fail even when required columns are present"
        );
    }

    #[test]
    fn process_seq_summ_rejects_too_many_preheader_comments() {
        let comment_block =
            iter::repeat_n("# sequencing summary comment\n", 1001).collect::<String>();
        let contents = format!("{comment_block}read_id\tsequence_length_template\nread-1\t1234\n");
        let temp_path = write_temp_seq_summary(contents.as_str());
        let result = process_seq_summ(
            temp_path
                .to_str()
                .expect("temporary file path should be valid UTF-8"),
        );
        remove_temp_file(&temp_path);

        assert!(
            result.is_err(),
            "more than 1000 pre-header comment lines should fail"
        );
    }

    #[test]
    fn process_seq_summ_rejects_overlong_data_line() {
        let long_extra_value = iter::repeat_n("abcd", 250).collect::<String>();
        let contents = format!(
            "read_id\tsequence_length_template\tfoo\nread-1\t1234\tok\nread-2\t5678\t{long_extra_value}\nread-3\t9012\tok\n"
        );
        let temp_path = write_temp_seq_summary(contents.as_str());
        let result = process_seq_summ(
            temp_path
                .to_str()
                .expect("temporary file path should be valid UTF-8"),
        );
        remove_temp_file(&temp_path);

        assert!(
            result.is_err(),
            "nonblank data line longer than 1000 bytes should fail even when required fields are present"
        );
    }
}

#[cfg(test)]
mod stochastic_tests {
    use super::*;
    use crate::{
        GenomicBed3,
        simulate_mod_bam::{
            ContigConfigBuilder, ReadConfigBuilder, SimulationConfigBuilder, TempBamSimulation,
        },
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

        #[expect(
            clippy::needless_continue,
            clippy::redundant_else,
            reason = "I prefer it this way; I think this is more readable"
        )]
        for k in seq_col.iter().zip(qual_col) {
            data_present = true;

            let seq = k.0.unwrap();

            let qual =
                k.1.unwrap()
                    .split('.')
                    .map(|x| x.parse::<u8>().unwrap())
                    .collect::<Vec<u8>>();

            if seq == "*" && qual == vec![255u8] {
                continue;
            } else {
                assert_eq!(seq.len(), seq_len, "Sequence length mismatch");
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
                    "need correct number of bases and/or special characters!"
                );
            }
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
        let seq_length = df.column("sequence_length_template")?.u32()?;
        let alignment_type = df.column("alignment_type")?.str()?;

        for k in align_length.iter().zip(seq_length).zip(alignment_type) {
            let al = k.0.0.unwrap().parse::<u32>().unwrap();
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
        let seq_length = df.column("sequence_length_template")?.u32()?;
        let alignment_type = df.column("alignment_type")?.str()?;

        for k in align_length.iter().zip(seq_length).zip(alignment_type) {
            let al = k.0.0.unwrap().parse::<u32>().unwrap();
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

        // extract seq and mod_count but we won't be getting a column called qualities.
        // mod_count will all be NAs
        let df_2 = run_reads_table_generation(
            &sim,
            Some(InputMods::<OptionalTag>::default()),
            SeqDisplayOptions::Full {
                show_base_qual: false,
            },
        )?;

        // Verify the dataframe is not empty and has expected columns
        assert!(df_2.height() > 0, "DataFrame should have some rows");
        assert_expected_columns(
            &df_2,
            ColumnsBuilder::default()
                .sequence(true)
                .mod_count(true)
                .build()?,
        );

        let seq_col_2 = df_2.column("sequence")?.str()?;
        drop(df_2.column("qualities").unwrap_err());
        let mod_count = df_2.column("mod_count")?.str()?;

        assert!(
            seq_col_2.iter().all(|x| x.unwrap().len() == 500),
            "500 bp seq expected"
        );
        assert!(
            mod_count.iter().all(|x| x == Some("NA")),
            "no mod information expected"
        );

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

    /// Helper function to create a simulation with indels and barcodes for region testing.
    fn create_indels_barcode_simulation() -> Result<TempBamSimulation, Error> {
        let contig_config = ContigConfigBuilder::default()
            .number(1)
            .len_range((1000, 1000))
            .repeated_seq("ACGT".into())
            .build()?;

        let read_config = ReadConfigBuilder::default()
            .number(100)
            .mapq_range((10, 20))
            .base_qual_range((41, 49))
            .len_range((1.0, 1.0))
            .barcode("AATT".into())
            .delete((0.1, 0.2))
            .insert_middle("GGTTGG".into());

        let sim_config = SimulationConfigBuilder::default()
            .contigs(contig_config)
            .reads(vec![read_config.build()?])
            .build()?;

        TempBamSimulation::new(sim_config)
    }

    /// Helper function to test sequence retrieval from a specific region.
    fn test_region_retrieval(
        sim: &TempBamSimulation,
        seq_display_options: SeqDisplayOptions,
        expected_seq: &str,
        seq_len: usize,
        counts: &Counts,
    ) -> Result<(), Error> {
        // check that we get the correct kind of `SeqDisplayOptions`
        assert!(matches!(
            seq_display_options,
            SeqDisplayOptions::Region { .. }
        ));

        // first, do a check without retrieving mods
        let df = run_reads_table_generation(sim, None, seq_display_options)?;

        assert_expected_columns(
            &df,
            ColumnsBuilder::default()
                .sequence(true)
                .qualities(true)
                .build()?,
        );

        let seq_col = df.column("sequence")?.str()?;
        let qual_col = df.column("qualities")?.str()?;
        let alignment_type = df.column("alignment_type")?.str()?;

        for k in seq_col.iter().zip(alignment_type) {
            if k.1.unwrap() == "unmapped" {
                assert_eq!(k.0.unwrap(), "*");
            } else {
                assert_eq!(k.0.unwrap(), expected_seq);
            }
        }

        check_seq_qual(seq_col, qual_col, seq_len, counts, 41u8..=49u8);

        // do a check with retrieving mods. Must get N/As in the mod count column
        let df_2 = run_reads_table_generation(
            sim,
            Some(InputMods::<OptionalTag>::default()),
            seq_display_options,
        )?;

        assert_expected_columns(
            &df_2,
            ColumnsBuilder::default()
                .mod_count(true)
                .sequence(true)
                .qualities(true)
                .build()?,
        );

        let mod_count = df_2.column("mod_count")?.str()?;

        assert!(
            mod_count.iter().all(|k| k == Some("NA")),
            "should not get mod counts as no mod info available!"
        );

        Ok(())
    }

    #[test]
    fn region_0_to_10() -> Result<(), Error> {
        let sim = create_indels_barcode_simulation()?;
        test_region_retrieval(
            &sim,
            SeqDisplayOptions::Region {
                show_base_qual: true,
                show_ins_lowercase: true,
                show_mod_z: false,
                region: GenomicBed3::new(0, 0, 10),
            },
            "ACGTACGTAC",
            10,
            &Counts::default(),
        )
    }

    #[test]
    fn region_100_to_110() -> Result<(), Error> {
        let sim = create_indels_barcode_simulation()?;
        test_region_retrieval(
            &sim,
            SeqDisplayOptions::Region {
                show_base_qual: true,
                show_ins_lowercase: true,
                show_mod_z: false,
                region: GenomicBed3::new(0, 100, 110),
            },
            "..........",
            10,
            &CountsBuilder::default().period(10).build()?,
        )
    }

    #[test]
    fn region_195_to_205() -> Result<(), Error> {
        let sim = create_indels_barcode_simulation()?;
        test_region_retrieval(
            &sim,
            SeqDisplayOptions::Region {
                show_base_qual: true,
                show_ins_lowercase: true,
                show_mod_z: false,
                region: GenomicBed3::new(0, 195, 205),
            },
            ".....ACGTA",
            10,
            &CountsBuilder::default().period(5).build()?,
        )
    }

    #[test]
    fn region_495_to_505() -> Result<(), Error> {
        let sim = create_indels_barcode_simulation()?;
        test_region_retrieval(
            &sim,
            SeqDisplayOptions::Region {
                show_base_qual: true,
                show_ins_lowercase: true,
                show_mod_z: false,
                region: GenomicBed3::new(0, 495, 505),
            },
            "TACGTggttggACGTA",
            16,
            &CountsBuilder::default().lowercase(6).build()?,
        )
    }

    #[test]
    fn region_495_to_505_no_ins_lowercase() -> Result<(), Error> {
        let sim = create_indels_barcode_simulation()?;
        test_region_retrieval(
            &sim,
            SeqDisplayOptions::Region {
                show_base_qual: true,
                show_ins_lowercase: false,
                show_mod_z: false,
                region: GenomicBed3::new(0, 495, 505),
            },
            "TACGTGGTTGGACGTA",
            16,
            &CountsBuilder::default().build()?,
        )
    }

    #[test]
    fn region_990_to_1000() -> Result<(), Error> {
        let sim = create_indels_barcode_simulation()?;
        test_region_retrieval(
            &sim,
            SeqDisplayOptions::Region {
                show_base_qual: true,
                show_ins_lowercase: true,
                show_mod_z: false,
                region: GenomicBed3::new(0, 990, 1000),
            },
            "GTACGTACGT",
            10,
            &CountsBuilder::default().build()?,
        )
    }
}

#[cfg(test)]
mod stochastic_tests_with_mods {
    use super::*;
    use crate::{
        GenomicBed3, InputModsBuilder,
        simulate_mod_bam::{
            ContigConfigBuilder, ModConfigBuilder, ReadConfigBuilder, SimulationConfigBuilder,
            TempBamSimulation,
        },
    };
    use rust_htslib::bam::Read as _;

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

    /// Helper function to create a simulation with indels and barcodes for region testing.
    fn create_indels_barcode_simulation() -> Result<TempBamSimulation, Error> {
        let contig_config = ContigConfigBuilder::default()
            .number(1)
            .len_range((1000, 1000))
            .repeated_seq("ACGT".into())
            .build()?;

        let mod_config = ModConfigBuilder::default()
            .base('G')
            .is_strand_plus(true)
            .mod_code("g".into())
            .win(vec![1, 1])
            .mod_range(vec![(0.0, 0.0), (0.6, 0.6)])
            .build()?;

        let read_config = ReadConfigBuilder::default()
            .number(100)
            .mapq_range((10, 20))
            .base_qual_range((35, 35))
            .len_range((1.0, 1.0))
            .barcode("AATT".into())
            .delete((0.1, 0.2))
            .insert_middle("GGTTGG".into())
            .mods(vec![mod_config]);

        let sim_config = SimulationConfigBuilder::default()
            .contigs(contig_config)
            .reads(vec![read_config.build()?])
            .build()?;

        TempBamSimulation::new(sim_config)
    }

    /// Helper function to test sequence retrieval from a specific region.
    fn test_region_retrieval(
        sim: &TempBamSimulation,
        seq_display_options: SeqDisplayOptions,
        expected_seq_fwd: &str,
        expected_seq_rev: &str,
        qual_col_str: &str,
        mod_count_fwd_str: &str,
        mod_count_rev_str: &str,
    ) -> Result<(), Error> {
        // check that we get the correct kind of `SeqDisplayOptions`
        assert!(matches!(
            seq_display_options,
            SeqDisplayOptions::Region { .. }
        ));

        // first, do a check without retrieving mods
        let df = run_reads_table_generation(
            sim,
            Some(InputMods::<OptionalTag>::default()),
            seq_display_options,
        )?;
        let seq_col = df.column("sequence")?.str()?;
        let qual_col = df.column("qualities")?.str()?;
        let mod_count = df.column("mod_count")?.str()?;
        let alignment_type = df.column("alignment_type")?.str()?;

        for k in seq_col
            .into_iter()
            .zip(qual_col)
            .zip(mod_count)
            .zip(alignment_type)
            .map(|(((seq, qual), mod_count_value), alignment_type_value)| {
                (seq, qual, mod_count_value, alignment_type_value)
            })
        {
            match k.3.unwrap() {
                "unmapped" => {
                    assert_eq!(k.0.unwrap(), "*");
                    assert_eq!(k.1.unwrap(), "255");
                    assert_eq!(k.2.unwrap(), mod_count_fwd_str);
                }
                "primary_forward" | "secondary_forward" | "supplementary_forward" => {
                    assert_eq!(k.0.unwrap(), expected_seq_fwd);
                    assert_eq!(k.1.unwrap(), qual_col_str);
                    assert_eq!(k.2.unwrap(), mod_count_fwd_str);
                }
                "primary_reverse" | "secondary_reverse" | "supplementary_reverse" => {
                    assert_eq!(k.0.unwrap(), expected_seq_rev);
                    assert_eq!(k.1.unwrap(), qual_col_str);
                    assert_eq!(k.2.unwrap(), mod_count_rev_str);
                }
                _ => unreachable!("read states fall in the above 7 categories"),
            }
        }

        Ok(())
    }

    #[test]
    fn region_0_to_10() -> Result<(), Error> {
        let sim = create_indels_barcode_simulation()?;
        for k in [
            (false, "ACGTACGTAC", "ACGTACGTAC"),
            (true, "ACGTACZTAC", "ACGTAZGTAC"),
        ] {
            test_region_retrieval(
                &sim,
                SeqDisplayOptions::Region {
                    show_base_qual: true,
                    show_ins_lowercase: true,
                    show_mod_z: k.0,
                    region: GenomicBed3::new(0, 0, 10),
                },
                k.1,
                k.2,
                "35.35.35.35.35.35.35.35.35.35",
                "g:114",
                "g:112",
            )?;
        }

        Ok(())
    }

    #[test]
    fn region_100_to_110() -> Result<(), Error> {
        let sim = create_indels_barcode_simulation()?;
        for k in [false, true] {
            test_region_retrieval(
                &sim,
                SeqDisplayOptions::Region {
                    show_base_qual: true,
                    show_ins_lowercase: true,
                    show_mod_z: k,
                    region: GenomicBed3::new(0, 100, 110),
                },
                "..........",
                "..........",
                "255.255.255.255.255.255.255.255.255.255",
                "g:114",
                "g:112",
            )?;
        }

        Ok(())
    }

    #[test]
    fn region_195_to_205() -> Result<(), Error> {
        let sim = create_indels_barcode_simulation()?;
        for k in [
            (false, ".....ACGTA", ".....ACGTA"),
            (true, ".....ACZTA", ".....AZGTA"),
        ] {
            test_region_retrieval(
                &sim,
                SeqDisplayOptions::Region {
                    show_base_qual: true,
                    show_ins_lowercase: true,
                    show_mod_z: k.0,
                    region: GenomicBed3::new(0, 195, 205),
                },
                k.1,
                k.2,
                "255.255.255.255.255.35.35.35.35.35",
                "g:114",
                "g:112",
            )?;
        }

        Ok(())
    }

    #[test]
    fn region_495_to_505() -> Result<(), Error> {
        let sim = create_indels_barcode_simulation()?;
        for k in [
            (false, "TACGTggttggACGTA", "TACGTggttggACGTA"),
            (true, "TACZTgzttgzACGTA", "TAZGTggttggACGTA"),
        ] {
            test_region_retrieval(
                &sim,
                SeqDisplayOptions::Region {
                    show_base_qual: true,
                    show_ins_lowercase: true,
                    show_mod_z: k.0,
                    region: GenomicBed3::new(0, 495, 505),
                },
                k.1,
                k.2,
                "35.35.35.35.35.35.35.35.35.35.35.35.35.35.35.35",
                "g:114",
                "g:112",
            )?;
        }

        Ok(())
    }

    #[test]
    fn region_990_to_1000() -> Result<(), Error> {
        let sim = create_indels_barcode_simulation()?;
        for k in [
            (false, "GTACGTACGT", "GTACGTACGT"),
            (true, "GTACZTACGT", "GTAZGTACGT"),
        ] {
            test_region_retrieval(
                &sim,
                SeqDisplayOptions::Region {
                    show_base_qual: true,
                    show_ins_lowercase: true,
                    show_mod_z: k.0,
                    region: GenomicBed3::new(0, 990, 1000),
                },
                k.1,
                k.2,
                "35.35.35.35.35.35.35.35.35.35",
                "g:114",
                "g:112",
            )?;
        }

        Ok(())
    }

    /// Region-based modification counting
    /// This test verifies that when `InputMods::region_bed3` is set,
    /// mod counts reflect only modifications within that region
    #[test]
    fn region_based_mod_counting() -> Result<(), Error> {
        let sim = create_indels_barcode_simulation()?;

        // Test region 0-10
        let mods_for_region = InputModsBuilder::<OptionalTag>::default()
            .region_bed3(GenomicBed3::new(0, 0, 10))
            .build()?;

        let df = run_reads_table_generation(
            &sim,
            Some(mods_for_region),
            SeqDisplayOptions::Region {
                show_base_qual: true,
                show_ins_lowercase: true,
                show_mod_z: false,
                region: GenomicBed3::new(0, 0, 10),
            },
        )?;

        let mod_count_col = df.column("mod_count")?.str()?;
        let alignment_type = df.column("alignment_type")?.str()?;

        // Count mods in region 0-10 for forward and reverse reads
        for (mod_count, aln_type) in mod_count_col.into_iter().zip(alignment_type) {
            match aln_type.unwrap() {
                "unmapped" => {
                    // Unmapped reads show mod count of 0
                    assert_eq!(
                        mod_count.unwrap(),
                        "g:0",
                        "unmapped reads should have g:0 mod_count"
                    );
                }
                "primary_forward" | "secondary_forward" | "supplementary_forward" => {
                    let count_str = mod_count.unwrap();
                    // Parse the count (format is "g:N")
                    if let Some(count_part) = count_str.strip_prefix("g:") {
                        let count: u32 = count_part.parse().expect("valid number");
                        // Region 0-10 has 10 bases, with "ACGTACGTAC" pattern
                        assert_eq!(count, 1, "expected 1 'g' mod in region 0-10, got {count}");
                    } else {
                        unreachable!("unexpected mod_count format: {count_str}");
                    }
                }
                "primary_reverse" | "secondary_reverse" | "supplementary_reverse" => {
                    let count_str = mod_count.unwrap();
                    if let Some(count_part) = count_str.strip_prefix("g:") {
                        let count: u32 = count_part.parse().expect("valid number");
                        assert_eq!(count, 1, "expected 1 'g' mod in region 0-10, got {count}");
                    } else {
                        unreachable!("unexpected mod_count format: {count_str}");
                    }
                }
                _ => unreachable!("read states fall in the above 7 categories"),
            }
        }

        Ok(())
    }

    /// Test Multiple modifications on the same read
    /// This test creates reads with two different modification types
    /// and verifies both are counted correctly
    fn create_multi_mod_simulation() -> Result<TempBamSimulation, Error> {
        let contig_config = ContigConfigBuilder::default()
            .number(1)
            .len_range((1000, 1000))
            .repeated_seq("ACGT".into())
            .build()?;

        // First mod: 'm' on 'C' bases
        let mod_config_c = ModConfigBuilder::default()
            .base('C')
            .is_strand_plus(true)
            .mod_code("m".into())
            .win(vec![1])
            .mod_range(vec![(0.55, 0.55)])
            .build()?;

        // Second mod: 'a' on 'A' bases but on the opposite
        // strand to the basecalled strand
        let mod_config_a = ModConfigBuilder::default()
            .base('A')
            .is_strand_plus(false)
            .mod_code("a".into())
            .win(vec![1, 1])
            .mod_range(vec![(0.6, 0.6), (0.1, 0.1)])
            .build()?;

        let read_config = ReadConfigBuilder::default()
            .number(100)
            .mapq_range((10, 20))
            .base_qual_range((35, 35))
            .len_range((1.0, 1.0))
            .barcode("AATT".into())
            .delete((0.1, 0.2))
            .insert_middle("GGTTGG".into())
            .mods(vec![mod_config_c, mod_config_a]);

        let sim_config = SimulationConfigBuilder::default()
            .contigs(contig_config)
            .reads(vec![read_config.build()?])
            .build()?;

        TempBamSimulation::new(sim_config)
    }

    #[test]
    fn multiple_mod_count_on_same_read() -> Result<(), Error> {
        let sim = create_multi_mod_simulation()?;

        let df = run_reads_table_generation(
            &sim,
            Some(InputMods::<OptionalTag>::default()),
            SeqDisplayOptions::Full {
                show_base_qual: true,
            },
        )?;

        let mod_count_col = df.column("mod_count")?.str()?;
        let alignment_type = df.column("alignment_type")?.str()?;

        // Check that both 'a' and 'm' modifications are present
        for (mod_count, aln_type) in mod_count_col.into_iter().zip(alignment_type) {
            match aln_type.unwrap() {
                "unmapped" | "primary_forward" | "secondary_forward" | "supplementary_forward" => {
                    let count_str = mod_count.unwrap();
                    assert_eq!(
                        count_str, "a:115;m:225",
                        "incorrect mod counts in unmapped/forward: {count_str}"
                    );
                }
                "primary_reverse" | "secondary_reverse" | "supplementary_reverse" => {
                    // we have an insertion in the middle of the read,
                    // which will have mods on it if it is reverse complemented,
                    // as we assign mods to A and C.
                    let count_str = mod_count.unwrap();
                    assert_eq!(
                        count_str, "a:116;m:229",
                        "incorrect mod counts in reverse: {count_str}"
                    );
                }
                _ => unreachable!("read states fall in the above 7 categories"),
            }
        }

        Ok(())
    }

    #[test]
    fn multiple_mod_count_on_same_read_with_stronger_thresholding() -> Result<(), Error> {
        let sim = create_multi_mod_simulation()?;

        let df = run_reads_table_generation(
            &sim,
            Some(
                InputModsBuilder::<OptionalTag>::default()
                    .mod_prob_filter(ThresholdState::InvertGtEqLtEq((100u8, 150u8).try_into()?))
                    .build()?,
            ),
            SeqDisplayOptions::Full {
                show_base_qual: true,
            },
        )?;

        let mod_count_col = df.column("mod_count")?.str()?;
        let alignment_type = df.column("alignment_type")?.str()?;

        // Check that both 'a' and 'm' modifications are present
        for (mod_count, aln_type) in mod_count_col.into_iter().zip(alignment_type) {
            match aln_type.unwrap() {
                "unmapped" | "primary_forward" | "secondary_forward" | "supplementary_forward" => {
                    let count_str = mod_count.unwrap();
                    assert_eq!(
                        count_str, "a:115;m:0",
                        "incorrect mod counts in unmapped/forward: {count_str}"
                    );
                }
                "primary_reverse" | "secondary_reverse" | "supplementary_reverse" => {
                    // we have an insertion in the middle of the read,
                    // which will have mods on it if it is reverse complemented,
                    // as we assign mods to A and C.
                    let count_str = mod_count.unwrap();
                    assert_eq!(
                        count_str, "a:116;m:0",
                        "incorrect mod counts in reverse: {count_str}"
                    );
                }
                _ => unreachable!("read states fall in the above 7 categories"),
            }
        }

        Ok(())
    }

    /// Test multiple modifications with region filtering.
    /// Verifies that region filtering and retrieval works correctly when
    /// multiple mod types are present
    #[test]
    fn multiple_mods_with_region_filtering() -> Result<(), Error> {
        let sim = create_multi_mod_simulation()?;

        // Test a small region
        let mods_for_region = InputModsBuilder::<OptionalTag>::default()
            .region_bed3(GenomicBed3::new(0, 495, 505))
            .build()?;

        let df = run_reads_table_generation(
            &sim,
            Some(mods_for_region),
            SeqDisplayOptions::Region {
                show_base_qual: true,
                show_ins_lowercase: true,
                show_mod_z: true,
                region: GenomicBed3::new(0, 495, 505),
            },
        )?;

        let mod_count_col = df.column("mod_count")?.str()?;
        let alignment_type_col = df.column("alignment_type")?.str()?;
        let seq_col = df.column("sequence")?.str()?;

        // Verify both mod types are present but with region-limited counts and correct sequence.
        // The sequence here is "TACGTggttggACGTA"
        for ((mod_count, aln_type), seq) in mod_count_col
            .into_iter()
            .zip(alignment_type_col)
            .zip(seq_col)
        {
            match aln_type.unwrap() {
                "unmapped" => {
                    // Unmapped reads show mod count of 0 when region filtering is applied
                    assert_eq!(
                        mod_count.unwrap(),
                        "a:0;m:0",
                        "unmapped reads should have zero counts with region filtering"
                    );
                    assert_eq!(
                        seq,
                        Some("*"),
                        "unmapped reads should have a * for sequence"
                    );
                }
                "primary_forward" | "secondary_forward" | "supplementary_forward" => {
                    let count_str = mod_count.unwrap();
                    let seq_str = seq.unwrap();
                    assert_eq!(
                        count_str, "a:1;m:2",
                        "incorrect mod counts in unmapped/forward: {count_str}"
                    );
                    assert_eq!(
                        seq_str, "TAZGTggttggZZGTA",
                        "incorrect mod counts in unmapped/forward: {count_str}"
                    );
                }
                "primary_reverse" | "secondary_reverse" | "supplementary_reverse" => {
                    let count_str = mod_count.unwrap();
                    let seq_str = seq.unwrap();
                    assert_eq!(
                        count_str, "a:3;m:6",
                        "incorrect mod counts in reverse: {count_str}"
                    );
                    assert_eq!(
                        seq_str, "ZACZTzzztzzACZZA",
                        "incorrect mod counts in unmapped/forward: {count_str}"
                    );
                }
                _ => unreachable!("read states fall in the above 7 categories"),
            }
        }

        Ok(())
    }

    /// Test multiple modifications with region filtering and with
    /// very strict filters that would let nothing through and with
    /// inserts set to same case as the others.
    #[test]
    fn multiple_mods_with_region_filtering_strict_mod_prob_filter_no_ins_lowercase()
    -> Result<(), Error> {
        let sim = create_multi_mod_simulation()?;

        // Test a small region
        let mods_for_region = InputModsBuilder::<OptionalTag>::default()
            .region_bed3(GenomicBed3::new(0, 495, 505))
            .mod_prob_filter(ThresholdState::Both((250u8, (100u8, 150u8).try_into()?)))
            .build()?;

        let df = run_reads_table_generation(
            &sim,
            Some(mods_for_region),
            SeqDisplayOptions::Region {
                show_base_qual: false,
                show_ins_lowercase: false,
                show_mod_z: true,
                region: GenomicBed3::new(0, 495, 505),
            },
        )?;

        let mod_count_col = df.column("mod_count")?.str()?;
        let alignment_type_col = df.column("alignment_type")?.str()?;
        let seq_col = df.column("sequence")?.str()?;

        // Verify both mod types are present but with region-limited counts and correct sequence.
        // The sequence here is "TACGTggttggACGTA" but without lowercase for the insert.
        for ((mod_count, aln_type), seq) in mod_count_col
            .into_iter()
            .zip(alignment_type_col)
            .zip(seq_col)
        {
            match aln_type.unwrap() {
                "unmapped" => {
                    // Unmapped reads show mod count of 0 when region filtering is applied
                    assert_eq!(
                        mod_count.unwrap(),
                        "a:0;m:0",
                        "unmapped reads should have zero counts with region filtering"
                    );
                    assert_eq!(
                        seq,
                        Some("*"),
                        "unmapped reads should have a * for sequence"
                    );
                }
                "primary_forward" | "secondary_forward" | "supplementary_forward" => {
                    let count_str = mod_count.unwrap();
                    let seq_str = seq.unwrap();
                    assert_eq!(
                        count_str, "a:0;m:0",
                        "incorrect mod counts in unmapped/forward: {count_str}"
                    );
                    assert_eq!(
                        seq_str, "TACGTGGTTGGACGTA",
                        "incorrect mod counts in unmapped/forward: {count_str}"
                    );
                }
                "primary_reverse" | "secondary_reverse" | "supplementary_reverse" => {
                    let count_str = mod_count.unwrap();
                    let seq_str = seq.unwrap();
                    assert_eq!(
                        count_str, "a:0;m:0",
                        "incorrect mod counts in reverse: {count_str}"
                    );
                    assert_eq!(
                        seq_str, "TACGTGGTTGGACGTA",
                        "incorrect mod counts in unmapped/forward: {count_str}"
                    );
                }
                _ => unreachable!("read states fall in the above 7 categories"),
            }
        }

        Ok(())
    }
}
