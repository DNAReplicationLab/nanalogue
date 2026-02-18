//! # Window modification data on reads
//!
//! In this module, we window data along molecules, and then output
//! these windows

use crate::{
    AlignmentInfoBuilder, CurrRead, Error, F32AbsValAtMost1, InputMods, InputWindowing, ModChar,
    OptionalTag, ReadState,
};
use fibertools_rs::utils::basemods::BaseMod;
use polars::prelude::*;
use rust_htslib::bam::Record;
use serde::Serialize;
use std::rc::Rc;

/// A single modification type's windowed data
#[derive(Serialize)]
struct WindowedModTableEntry {
    /// The canonical base that is modified (e.g. 'C', 'T', 'A', 'G', 'N')
    base: char,
    /// Whether the modification strand is the plus strand
    is_strand_plus: bool,
    /// The modification code (e.g. 'm' for 5mC, 'T' for `BrdU`)
    mod_code: ModChar,
    /// Windowed data: `(win_start, win_end, win_val, mean_base_qual, ref_win_start, ref_win_end)`
    data: Vec<(i64, i64, F32AbsValAtMost1, u8, i64, i64)>,
}

/// A single read's windowed modification data for JSON output
#[derive(Serialize)]
struct WindowedReadEntry {
    /// The alignment type/state of this read
    alignment_type: ReadState,
    /// Alignment coordinates (absent for unmapped reads)
    #[serde(skip_serializing_if = "Option::is_none")]
    alignment: Option<crate::read_utils::AlignmentInfo>,
    /// Per-modification-type windowed data
    mod_table: Vec<WindowedModTableEntry>,
    /// The read identifier (QNAME)
    read_id: String,
    /// Basecalled sequence length
    seq_len: u64,
}

/// Computes windowed modification data for a single base modification type
///
/// Extracts annotation data from `base_mod`, slides a window of `win_size`
/// with step `slide_size`, and computes the window value, mean base quality,
/// and reference coordinate bounds for each window position.
#[expect(
    clippy::too_many_lines,
    reason = "windowing logic is simple but just takes up many lines"
)]
fn compute_windowed_mod_data<F>(
    base_mod: &BaseMod,
    base_qual: &[u8],
    win_size: usize,
    slide_size: usize,
    qname: &str,
    window_function: &F,
) -> Result<WindowedModTableEntry, Error>
where
    F: Fn(&[u8]) -> Result<F32AbsValAtMost1, Error>,
{
    // constant to mark windows with basecalled coordinates but no reference coordinates.
    const INVALID_REF_POS: i64 = -1;

    #[expect(
        clippy::type_complexity,
        reason = "I think a tuple of 5 `Vec` is fine if its readable"
    )]
    let (mod_data, starts, ends, ref_starts, ref_ends): (
        Vec<u8>,
        Vec<i64>,
        Vec<i64>,
        Vec<Option<i64>>,
        Vec<Option<i64>>,
    ) = base_mod
        .ranges
        .annotations
        .iter()
        .map(|k| (k.qual, k.start, k.end, k.reference_start, k.reference_end))
        .collect();
    let base = base_mod.modified_base as char;
    let mod_strand = base_mod.strand;
    let mod_type = ModChar::new(base_mod.modification_type);

    let mut windows: Vec<(i64, i64, F32AbsValAtMost1, u8, i64, i64)> = Vec::new();

    if let Some(v) = mod_data.len().checked_sub(win_size) {
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "a +1 on `ref_win_end`, no overflow as coords << 2^63, complex arithmetic in Q score avg"
        )]
        for window_idx in (0..=v).step_by(slide_size) {
            let win_val = match window_function(
                mod_data
                    .get(window_idx..)
                    .expect("window_idx <= v where v = len - win_size")
                    .get(0..win_size)
                    .expect("no error as we've checked data len >= win size"),
            ) {
                Ok(val) => val,
                Err(e) => {
                    eprintln!(
                        "Warning: Skipping {win_size} window starting at {qname}:{window_idx} due to error: {e}"
                    );
                    continue;
                }
            };
            // there is no way to trigger the errors below as we control how CurrRead is
            // populated quite strictly. Nevertheless, I am leaving these in for
            // future-proofing.
            let win_start = *starts.get(window_idx).expect("window_idx is valid");
            let win_end = *ends
                .get(window_idx..)
                .expect("window_idx <= v where v = len - win_size")
                .get(0..win_size)
                .expect("no error as we've checked data len >= win size")
                .last()
                .expect("no error as we've checked data len >= win size");

            let ref_win_start = ref_starts
                .get(window_idx..)
                .expect("window_idx <= v where v = len - win_size")
                .get(0..win_size)
                .expect("no error as we've checked data len >= win size")
                .iter()
                .flatten()
                .min()
                .copied()
                .unwrap_or(INVALID_REF_POS);
            let ref_win_end = ref_ends
                .get(window_idx..)
                .expect("window_idx <= v where v = len - win_size")
                .get(0..win_size)
                .expect("no error as we've checked data len >= win size")
                .iter()
                .flatten()
                .max()
                .copied()
                .map_or(INVALID_REF_POS, |x| x + 1);
            #[expect(
                clippy::cast_possible_truncation,
                clippy::cast_sign_loss,
                reason = "we are forced to do these due to the math itself, of taking a power, avg, and then log"
            )]
            let mean_base_qual = {
                let quals = base_qual
                    .get(usize::try_from(win_start)?..usize::try_from(win_end)?)
                    .expect("no error as `win_start`, `win_end` in range");
                if quals.is_empty() || quals.first() == Some(&255u8) {
                    // BAM format is such that all values are 255, or values are between
                    // 0 and 93. So if we see one 255, we can just return 255.
                    // Empty quals (win_start == win_end) also get 255 as a sentinel.
                    255u8
                } else {
                    // we do an average using the probability of errors,
                    // and not the Q scores directly.
                    let quals_min = quals.iter().min().expect("no error");
                    let data_size = f64::from(i32::try_from(win_end)? - i32::try_from(win_start)?);
                    let x = quals.iter().fold(0f64, |acc, x| {
                        acc + (10f64).powf(-0.1f64 * f64::from(x - quals_min))
                    });
                    quals_min.saturating_add((-10f64 * f64::log10(x / data_size)).round() as u8)
                }
            };
            windows.push((
                win_start,
                win_end,
                win_val,
                mean_base_qual,
                ref_win_start,
                ref_win_end,
            ));
        }
    }

    Ok(WindowedModTableEntry {
        base,
        is_strand_plus: mod_strand == '+',
        mod_code: mod_type,
        data: windows,
    })
}

/// Windowed modification data along molecules
///
/// # Examples
///
/// Windowing the first (mapped) and last (unmapped) records from `example_1.bam`
/// with a window size of 2 and step of 1 using [`threshold_and_mean`](crate::analysis::threshold_and_mean):
///
/// ```
/// use nanalogue_core::{window_reads, InputMods, InputWindowing, OptionalTag};
/// use nanalogue_core::analysis::threshold_and_mean;
/// use rust_htslib::bam::{self, Read as _};
/// use std::rc::Rc;
///
/// let mut bam_reader = bam::Reader::from_path("examples/example_1.bam").unwrap();
/// let records: Vec<Rc<bam::Record>> = bam_reader
///     .rc_records()
///     .collect::<Result<Vec<_>, _>>()
///     .unwrap();
/// let selected: Vec<Result<Rc<bam::Record>, rust_htslib::errors::Error>> = vec![
///     Ok(Rc::clone(&records[0])),
///     Ok(Rc::clone(records.last().unwrap())),
/// ];
///
/// let window_options: InputWindowing =
///     serde_json::from_str("{\"win\": 2, \"step\": 1}").unwrap();
/// let mods: InputMods<OptionalTag> = InputMods::default();
/// let mut output = Vec::new();
///
/// window_reads::run(&mut output, selected, window_options, &mods, |x| {
///     threshold_and_mean(x).map(Into::into)
/// })
/// .unwrap();
///
/// let output_str = String::from_utf8(output).unwrap();
/// let rows: Vec<Vec<&str>> = output_str
///     .lines()
///     .map(|line| line.split('\t').collect())
///     .collect();
///
/// // 1 header + 3 mapped windows + 9 unmapped windows = 13 rows
/// assert_eq!(rows.len(), 13);
///
/// // Header
/// assert_eq!(rows[0], ["#contig", "ref_win_start", "ref_win_end", "read_id", "win_val",
///     "strand", "base", "mod_strand", "mod_type", "win_start", "win_end", "basecall_qual"]);
///
/// // Mapped read (primary_forward on dummyI): 3 windows
/// //         contig  ref_s ref_e read_id                               val  str  b  ms mt  ws we bq
/// assert_eq!(rows[1], ["dummyI", "9",  "13", "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575", "0", "+", "T", "+", "T", "0", "4", "255"]);
/// assert_eq!(rows[2], ["dummyI", "12", "14", "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575", "0", "+", "T", "+", "T", "3", "5", "255"]);
/// assert_eq!(rows[3], ["dummyI", "13", "17", "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575", "0", "+", "T", "+", "T", "4", "8", "255"]);
///
/// // Unmapped read: 5 windows for G/-/7200, then 4 windows for T/+/T
/// //         contig ref_s ref_e read_id                               val    str b  ms mt     ws   we  bq
/// assert_eq!(rows[4],  [".", "-1", "-1", "a4f36092-b4d5-47a9-813e-c22c3b477a0c", "0",   ".", "G", "-", "7200", "28", "30", "255"]);
/// assert_eq!(rows[8],  [".", "-1", "-1", "a4f36092-b4d5-47a9-813e-c22c3b477a0c", "0",   ".", "G", "-", "7200", "43", "45", "255"]);
/// assert_eq!(rows[9],  [".", "-1", "-1", "a4f36092-b4d5-47a9-813e-c22c3b477a0c", "1",   ".", "T", "+", "T",    "3",  "9",  "255"]);
/// assert_eq!(rows[12], [".", "-1", "-1", "a4f36092-b4d5-47a9-813e-c22c3b477a0c", "0.5", ".", "T", "+", "T",    "39", "48", "255"]);
/// ```
///
/// # Errors
/// Returns an error if BAM record reading, or output writing fails.
///
pub fn run<W, F, D>(
    handle: &mut W,
    bam_records: D,
    window_options: InputWindowing,
    mods: &InputMods<OptionalTag>,
    window_function: F,
) -> Result<(), Error>
where
    W: std::io::Write,
    F: Fn(&[u8]) -> Result<F32AbsValAtMost1, Error>,
    D: IntoIterator<Item = Result<Rc<Record>, rust_htslib::errors::Error>>,
{
    // Get windowing parameters
    let win_size = window_options.win.get();
    let slide_size = window_options.step.get();

    // print header
    writeln!(
        handle,
        "#contig\tref_win_start\tref_win_end\tread_id\twin_val\tstrand\t\
base\tmod_strand\tmod_type\twin_start\twin_end\tbasecall_qual",
    )?;

    // Go record by record in the BAM file,
    for r in bam_records {
        // read records
        let record = r?;

        // set data in records
        let curr_read_state = CurrRead::default()
            .try_from_only_alignment(&record)?
            .set_mod_data_restricted_options(&record, mods)?;
        let qname = curr_read_state.read_id();
        let strand = curr_read_state.strand();
        let contig = if curr_read_state.read_state().is_unmapped() {
            "."
        } else {
            curr_read_state.contig_name()?
        };
        let base_qual = record.qual();

        // read and window modification data, then print the output
        for base_mod in &curr_read_state.mod_data().0.base_mods {
            let result = compute_windowed_mod_data(
                base_mod,
                base_qual,
                win_size,
                slide_size,
                qname,
                &window_function,
            )?;
            let mod_strand = if result.is_strand_plus { '+' } else { '-' };
            for &(win_start, win_end, win_val, mean_base_qual, ref_win_start, ref_win_end) in
                &result.data
            {
                writeln!(
                    handle,
                    "{contig}\t{ref_win_start}\t{ref_win_end}\t{qname}\t{win_val}\t{strand}\t\
                    {}\t{mod_strand}\t{}\t{win_start}\t{win_end}\t{mean_base_qual}",
                    result.base, result.mod_code,
                )?;
            }
        }
    }

    Ok(())
}

/// Creates a `DataFrame` from windowed modification data
///
/// This function calls [`run`] with a buffer handle, then parses the output into a Polars `DataFrame`.
/// The first line of output (after removing the leading '#') contains tab-separated column names,
/// and subsequent lines contain tab-separated data values.
///
/// # Errors
/// Returns an error if BAM record reading, output writing, or `DataFrame` construction fails.
///
pub fn run_df<F, D>(
    bam_records: D,
    window_options: InputWindowing,
    mods: &InputMods<OptionalTag>,
    window_function: F,
) -> Result<DataFrame, Error>
where
    F: Fn(&[u8]) -> Result<F32AbsValAtMost1, Error>,
    D: IntoIterator<Item = Result<Rc<Record>, rust_htslib::errors::Error>>,
{
    // Create a buffer to capture output
    let mut buffer = Vec::new();

    // Call run with the buffer
    run(
        &mut buffer,
        bam_records,
        window_options,
        mods,
        window_function,
    )?;

    // Convert buffer to string and remove leading '#' from header
    let output = String::from_utf8(buffer)?;
    let output_without_hash = output
        .strip_prefix('#')
        .ok_or_else(|| Error::InvalidState("Output does not start with '#'".to_string()))?;

    // Define schema based on column types
    let schema_fields = vec![
        Field::new("contig".into(), DataType::String),
        Field::new("ref_win_start".into(), DataType::Int64),
        Field::new("ref_win_end".into(), DataType::Int64),
        Field::new("read_id".into(), DataType::String),
        Field::new("win_val".into(), DataType::Float32),
        Field::new("strand".into(), DataType::String),
        Field::new("base".into(), DataType::String),
        Field::new("mod_strand".into(), DataType::String),
        Field::new("mod_type".into(), DataType::String),
        Field::new("win_start".into(), DataType::UInt64),
        Field::new("win_end".into(), DataType::UInt64),
        Field::new("basecall_qual".into(), DataType::UInt32),
    ];

    let schema = Schema::from_iter(schema_fields);

    // Parse the TSV data with the schema
    let cursor = std::io::Cursor::new(output_without_hash.as_bytes());
    let df = CsvReadOptions::default()
        .with_has_header(true)
        .map_parse_options(|parse_options| parse_options.with_separator(b'\t'))
        .with_schema(Some(Arc::new(schema)))
        .into_reader_with_file_handle(cursor)
        .finish()?;

    Ok(df)
}

/// Windowed modification data along molecules, output as JSON
///
/// Produces the same windowed data as [`run`] but serializes each read as a
/// JSON object with alignment info and a `mod_table` whose `data` entries
/// contain `[win_start, win_end, win_val, mean_base_qual, ref_win_start, ref_win_end]`.
///
/// # Examples
///
/// Windowing the first (mapped) and last (unmapped) records from `example_1.bam`
/// with a window size of 2 and step of 1 using [`threshold_and_mean`](crate::analysis::threshold_and_mean):
///
/// ```
/// use nanalogue_core::{window_reads, InputMods, InputWindowing, OptionalTag};
/// use nanalogue_core::analysis::threshold_and_mean;
/// use rust_htslib::bam::{self, Read as _};
/// use std::rc::Rc;
///
/// let mut bam_reader = bam::Reader::from_path("examples/example_1.bam").unwrap();
/// let records: Vec<Rc<bam::Record>> = bam_reader
///     .rc_records()
///     .collect::<Result<Vec<_>, _>>()
///     .unwrap();
/// let selected: Vec<Result<Rc<bam::Record>, rust_htslib::errors::Error>> = vec![
///     Ok(Rc::clone(&records[0])),
///     Ok(Rc::clone(records.last().unwrap())),
/// ];
///
/// let window_options: InputWindowing =
///     serde_json::from_str("{\"win\": 2, \"step\": 1}").unwrap();
/// let mods: InputMods<OptionalTag> = InputMods::default();
/// let mut output = Vec::new();
///
/// window_reads::run_json(&mut output, selected, window_options, &mods, |x| {
///     threshold_and_mean(x).map(Into::into)
/// })
/// .unwrap();
///
/// let output_str = String::from_utf8(output).unwrap();
/// let parsed: Vec<serde_json::Value> = serde_json::from_str(&output_str).unwrap();
/// assert_eq!(parsed.len(), 2);
///
/// // Mapped read (primary_forward on dummyI): 1 mod type, 3 windows
/// assert_eq!(parsed[0], serde_json::json!({
///     "alignment_type": "primary_forward",
///     "alignment": { "start": 9, "end": 17, "contig": "dummyI", "contig_id": 0 },
///     "read_id": "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575",
///     "seq_len": 8,
///     "mod_table": [{
///         "base": "T", "is_strand_plus": true, "mod_code": "T",
///         "data": [
///             [0, 4, 0.0, 255, 9,  13],
///             [3, 5, 0.0, 255, 12, 14],
///             [4, 8, 0.0, 255, 13, 17],
///         ]
///     }]
/// }));
///
/// // Unmapped read: 2 mod types (G/-/7200 with 5 windows, T/+/T with 4 windows)
/// assert_eq!(parsed[1], serde_json::json!({
///     "alignment_type": "unmapped",
///     "read_id": "a4f36092-b4d5-47a9-813e-c22c3b477a0c",
///     "seq_len": 48,
///     "mod_table": [
///         {
///             "base": "G", "is_strand_plus": false, "mod_code": "7200",
///             "data": [
///                 [28, 30, 0.0, 255, -1, -1],
///                 [29, 31, 0.0, 255, -1, -1],
///                 [30, 33, 0.0, 255, -1, -1],
///                 [32, 44, 0.0, 255, -1, -1],
///                 [43, 45, 0.0, 255, -1, -1],
///             ]
///         },
///         {
///             "base": "T", "is_strand_plus": true, "mod_code": "T",
///             "data": [
///                 [3,  9,  1.0, 255, -1, -1],
///                 [8,  28, 0.5, 255, -1, -1],
///                 [27, 40, 0.0, 255, -1, -1],
///                 [39, 48, 0.5, 255, -1, -1],
///             ]
///         }
///     ]
/// }));
/// ```
///
/// # Errors
/// Returns an error if BAM record reading, or output writing fails.
///
#[expect(
    clippy::missing_panics_doc,
    reason = "iterator `.next()` on repeat should not fail"
)]
pub fn run_json<W, F, D>(
    handle: &mut W,
    bam_records: D,
    window_options: InputWindowing,
    mods: &InputMods<OptionalTag>,
    window_function: F,
) -> Result<(), Error>
where
    W: std::io::Write,
    F: Fn(&[u8]) -> Result<F32AbsValAtMost1, Error>,
    D: IntoIterator<Item = Result<Rc<Record>, rust_htslib::errors::Error>>,
{
    // Get windowing parameters
    let win_size = window_options.win.get();
    let slide_size = window_options.step.get();

    let mut is_first_record_written = vec![false].into_iter().chain(std::iter::repeat(true));

    write!(handle, "[")?;

    // Go record by record in the BAM file,
    for r in bam_records {
        // read records
        let record = r?;

        // set data in records
        let curr_read_state = CurrRead::default()
            .try_from_only_alignment(&record)?
            .set_mod_data_restricted_options(&record, mods)?;
        let qname = curr_read_state.read_id().to_owned();
        let read_state = curr_read_state.read_state();
        let base_qual = record.qual();
        let seq_len = curr_read_state.seq_len()?;

        // Build alignment info for mapped reads
        let alignment = if read_state.is_unmapped() {
            None
        } else {
            let (contig_id, start) = curr_read_state.contig_id_and_start()?;
            let align_len = curr_read_state.align_len()?;
            let end = start
                .checked_add(align_len)
                .ok_or(Error::Arithmetic("alignment end overflow".to_owned()))?;
            Some(
                AlignmentInfoBuilder::default()
                    .start(start)
                    .end(end)
                    .contig(curr_read_state.contig_name()?.to_owned())
                    .contig_id(contig_id)
                    .build()?,
            )
        };

        // read and window modification data
        let mut mod_table: Vec<WindowedModTableEntry> = Vec::new();
        for base_mod in &curr_read_state.mod_data().0.base_mods {
            mod_table.push(compute_windowed_mod_data(
                base_mod,
                base_qual,
                win_size,
                slide_size,
                &qname,
                &window_function,
            )?);
        }

        let entry = WindowedReadEntry {
            alignment_type: read_state,
            alignment,
            mod_table,
            read_id: qname,
            seq_len,
        };

        if is_first_record_written.next().expect("no error") {
            writeln!(handle, ",")?;
        } else {
            writeln!(handle)?;
        }
        write!(handle, "{}", serde_json::to_string(&entry)?)?;
    }

    writeln!(handle, "\n]")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::F32Bw0and1;
    use crate::analysis::{threshold_and_mean, threshold_and_mean_and_thres_win};
    use rust_htslib::bam::{self, Read as _};

    /// Helper function to run `window_reads` tests with `threshold_and_mean_and_thres_win`
    ///
    /// This function encapsulates the common test setup and execution logic for `window_reads` tests
    /// that use a threshold value for filtering.
    fn run_window_reads_test_with_threshold(
        input_file: &str,
        threshold: Option<f32>,
        expected_output_file: &str,
    ) -> Result<(), Error> {
        // Set input, output, options
        let mut output = Vec::new();
        let mut bam_reader = bam::Reader::from_path(input_file)?;
        let bam_records = bam_reader.rc_records();
        let window_options: InputWindowing =
            serde_json::from_str("{\"win\": 2, \"step\": 1}").unwrap();
        let mods = InputMods::default();

        // Run the window_reads function with appropriate function based on threshold
        match threshold {
            None => {
                // Use threshold_and_mean when no threshold is specified
                run(&mut output, bam_records, window_options, &mods, |x| {
                    threshold_and_mean(x).map(Into::into)
                })?;
            }
            Some(thres_val) => {
                // Use threshold_and_mean_and_thres_win with specified threshold
                run(&mut output, bam_records, window_options, &mods, |x| {
                    threshold_and_mean_and_thres_win(x, F32Bw0and1::new(thres_val).unwrap())
                        .map(Into::into)
                })?;
            }
        }

        // Perform comparison
        let output_str = String::from_utf8(output)?;
        let expected_output = std::fs::read_to_string(expected_output_file)?;
        assert_eq!(output_str, expected_output);

        Ok(())
    }

    #[test]
    fn window_reads_example_1() -> Result<(), Error> {
        run_window_reads_test_with_threshold(
            "./examples/example_1.bam",
            None,
            "./examples/example_1_window_reads",
        )
    }

    #[test]
    fn window_reads_example_1_gt_0pt4() -> Result<(), Error> {
        run_window_reads_test_with_threshold(
            "./examples/example_1.bam",
            Some(0.4),
            "./examples/example_1_window_reads_gt_0pt4",
        )
    }

    #[test]
    fn window_reads_example_1_gt_0pt8() -> Result<(), Error> {
        run_window_reads_test_with_threshold(
            "./examples/example_1.bam",
            Some(0.8),
            "./examples/example_1_window_reads_gt_0pt8",
        )
    }

    #[test]
    fn window_reads_example_7() -> Result<(), Error> {
        run_window_reads_test_with_threshold(
            "./examples/example_7.sam",
            None,
            "./examples/example_7_window_reads",
        )
    }

    #[test]
    fn window_reads_json_example_7() -> Result<(), Error> {
        let mut output = Vec::new();
        let mut bam_reader = bam::Reader::from_path("./examples/example_7.sam")?;
        let bam_records = bam_reader.rc_records();
        let window_options: InputWindowing =
            serde_json::from_str("{\"win\": 2, \"step\": 1}").unwrap();
        let mods = InputMods::default();

        run_json(&mut output, bam_records, window_options, &mods, |x| {
            threshold_and_mean(x).map(Into::into)
        })?;

        let output_str = String::from_utf8(output)?;
        let parsed_output: serde_json::Value = serde_json::from_str(&output_str)?;
        let expected_output = std::fs::read_to_string("./examples/example_7_window_reads_json")?;
        let parsed_expected: serde_json::Value = serde_json::from_str(&expected_output)?;

        assert_eq!(
            parsed_output, parsed_expected,
            "Windowed JSON output should match expected file"
        );

        Ok(())
    }

    #[test]
    fn window_reads_json_example_1() -> Result<(), Error> {
        let mut output = Vec::new();
        let mut bam_reader = bam::Reader::from_path("./examples/example_1.bam")?;
        let bam_records = bam_reader.rc_records();
        let window_options: InputWindowing =
            serde_json::from_str("{\"win\": 2, \"step\": 1}").unwrap();
        let mods = InputMods::default();

        run_json(&mut output, bam_records, window_options, &mods, |x| {
            threshold_and_mean(x).map(Into::into)
        })?;

        let output_str = String::from_utf8(output)?;
        let parsed_output: serde_json::Value = serde_json::from_str(&output_str)?;
        let expected_output = std::fs::read_to_string("./examples/example_1_window_reads_json")?;
        let parsed_expected: serde_json::Value = serde_json::from_str(&expected_output)?;

        assert_eq!(
            parsed_output, parsed_expected,
            "Windowed JSON output should match expected file"
        );

        Ok(())
    }

    /// Helper to run `run_json` with `threshold_and_gradient` on a SAM/BAM file
    fn run_json_grad(input: &str, win: usize, step: usize) -> Result<String, Error> {
        let mut output = Vec::new();
        let mut bam_reader = bam::Reader::from_path(input)?;
        let bam_records = bam_reader.rc_records();
        let window_options: InputWindowing =
            serde_json::from_str(&format!("{{\"win\": {win}, \"step\": {step}}}"))?;
        let mods = InputMods::default();

        run_json(
            &mut output,
            bam_records,
            window_options,
            &mods,
            crate::analysis::threshold_and_gradient,
        )?;

        Ok(String::from_utf8(output)?)
    }

    /// Helper to compare JSON output with expected file (parsed comparison)
    fn assert_json_matches_file(actual_json: &str, expected_path: &str) {
        let parsed_actual: serde_json::Value =
            serde_json::from_str(actual_json).expect("parse actual JSON");
        let expected_str = std::fs::read_to_string(expected_path).expect("read expected JSON file");
        let parsed_expected: serde_json::Value =
            serde_json::from_str(&expected_str).expect("parse expected JSON");
        assert_eq!(
            parsed_actual, parsed_expected,
            "JSON output should match expected file {expected_path}"
        );
    }

    #[test]
    fn window_grad_json_example_10_win_10_step_1() -> Result<(), Error> {
        let output = run_json_grad("./examples/example_10.sam", 10, 1)?;
        assert_json_matches_file(&output, "./examples/example_10_win_grad_json_win_10_step_1");
        Ok(())
    }

    #[test]
    fn window_grad_json_example_10_win_20_step_2() -> Result<(), Error> {
        let output = run_json_grad("./examples/example_10.sam", 20, 2)?;
        assert_json_matches_file(&output, "./examples/example_10_win_grad_json_win_20_step_2");
        Ok(())
    }

    #[test]
    fn window_grad_json_example_11_win_10_step_1() -> Result<(), Error> {
        let output = run_json_grad("./examples/example_11.sam", 10, 1)?;
        assert_json_matches_file(&output, "./examples/example_11_win_grad_json_win_10_step_1");
        Ok(())
    }

    #[test]
    fn window_grad_json_example_11_win_20_step_2() -> Result<(), Error> {
        let output = run_json_grad("./examples/example_11.sam", 20, 2)?;
        assert_json_matches_file(&output, "./examples/example_11_win_grad_json_win_20_step_2");
        Ok(())
    }
}

#[cfg(test)]
mod stochastic_tests {
    use super::*;
    use crate::SimulationConfig;
    use crate::analysis;
    use crate::simulate_mod_bam::TempBamSimulation;
    use itertools::izip;
    use rust_htslib::bam::{self, Read as _};

    /// Helper to create a simulation from JSON config
    fn create_test_simulation(config_json: &str) -> Result<TempBamSimulation, Error> {
        let config: SimulationConfig = serde_json::from_str(config_json)?;
        TempBamSimulation::new(config)
    }

    /// Helper to run window analysis with `threshold_and_mean` aggregation function
    fn run_window_analysis_with_threshold(
        sim: &TempBamSimulation,
        win: usize,
        step: usize,
    ) -> Result<DataFrame, Error> {
        let mut bam_reader = bam::Reader::from_path(sim.bam_path())?;
        let bam_records = bam_reader.rc_records();

        let window_options: InputWindowing =
            serde_json::from_str(&format!("{{\"win\": {win}, \"step\": {step}}}"))?;
        let mods = InputMods::default();

        run_df(bam_records, window_options, &mods, |x| {
            analysis::threshold_and_mean(x).map(Into::into)
        })
    }

    /// Helper to assert dataframe has expected column headers
    fn assert_expected_columns(df: &DataFrame) {
        let expected_columns = vec![
            "contig",
            "ref_win_start",
            "ref_win_end",
            "read_id",
            "win_val",
            "strand",
            "base",
            "mod_strand",
            "mod_type",
            "win_start",
            "win_end",
            "basecall_qual",
        ];

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

    /// Test that `run_df` produces an empty dataframe for BAM files with no modification data
    ///
    /// This test creates a simulated BAM file without any modification information and verifies
    /// that `run_df` correctly returns an empty dataframe (no data rows, only column headers).
    #[test]
    fn run_df_empty_for_no_mods() -> Result<(), Error> {
        // Create simulation config with no modifications
        let config_json = r#"{
            "contigs": {
                "number": 2,
                "len_range": [100, 200]
            },
            "reads": [{
                "number": 100,
                "mapq_range": [10, 20],
                "base_qual_range": [10, 20],
                "len_range": [0.1, 0.8]
            }]
        }"#;

        let sim = create_test_simulation(config_json)?;
        let df = run_window_analysis_with_threshold(&sim, 2, 1)?;

        // Verify the dataframe is empty (no data rows)
        assert_eq!(
            df.height(),
            0,
            "DataFrame should have no rows for BAM with no modifications"
        );

        assert_expected_columns(&df);

        Ok(())
    }

    /// Test that `run_df` produces a non-empty dataframe with modification data
    ///
    /// This test creates a simulated BAM file with modification data and verifies that
    /// `run_df` correctly processes the modifications and returns a dataframe with data rows.
    #[test]
    fn run_df_with_mods() -> Result<(), Error> {
        // Create simulation config with modifications
        let config_json = r#"{
            "contigs": {
                "number": 4,
                "len_range": [100, 200]
            },
            "reads": [{
                "number": 1000,
                "mapq_range": [10, 20],
                "base_qual_range": [20, 30],
                "len_range": [0.1, 0.8],
                "mods": [{
                    "base": "T",
                    "is_strand_plus": true,
                    "mod_code": "T",
                    "win": [4],
                    "mod_range": [[0.1, 0.2]]
                }]
            }]
        }"#;

        let sim = create_test_simulation(config_json)?;
        let df = run_window_analysis_with_threshold(&sim, 2, 1)?;

        // Verify the dataframe is NOT empty (should have data rows with mods)
        assert!(
            df.height() > 0,
            "DataFrame should have rows when modifications are present"
        );

        assert_expected_columns(&df);

        let mod_qual = df.column("win_val")?.f32()?;
        assert!(mod_qual.iter().all(|x| x == Some(0.0)));

        let basecall_qual = df.column("basecall_qual")?.u32()?;
        assert!(
            basecall_qual
                .iter()
                .all(|x| (20..=30).contains(&x.unwrap()))
        );

        Ok(())
    }

    /// Test that `run_df` produces a non-empty dataframe with modification data,
    /// when reads that are 'noisy' i.e. not perfectly aligned are used.
    ///
    /// This test creates a simulated BAM file with modification data and verifies that
    /// `run_df` correctly processes the modifications and returns a dataframe with data rows.
    #[test]
    fn run_df_with_mods_and_non_perfectly_aligned_reads() -> Result<(), Error> {
        // Create simulation config with modifications
        let config_json = r#"{
            "contigs": {
                "number": 4,
                "len_range": [100000, 200000]
            },
            "reads": [{
                "number": 100,
                "mapq_range": [10, 20],
                "base_qual_range": [30, 40],
                "len_range": [0.1, 0.8],
                "delete": [0.5, 0.7],
                "insert_middle": "ATCGAATTGGAA",
                "mismatch": 0.2,
                "mods": [{
                    "base": "C",
                    "is_strand_plus": false,
                    "mod_code": "m",
                    "win": [4],
                    "mod_range": [[0.2, 0.8]]
                }]
            }]
        }"#;

        let sim = create_test_simulation(config_json)?;
        let df = run_window_analysis_with_threshold(&sim, 200, 100)?;

        // Verify the dataframe is NOT empty (should have data rows with mods)
        assert!(
            df.height() > 0,
            "DataFrame should have rows when modifications are present"
        );

        assert_expected_columns(&df);

        let mod_qual = df.column("win_val")?.f32()?;
        assert!(mod_qual.iter().all(|x| (0.2..=0.8).contains(&x.unwrap())));
        // when we window, we threshold. so a 20%-80% chance of mod on a base level
        // gets converted into a 0 or a 1 depending on whether the probability is
        // below or above 50%. Now, when we window over like 200 candidate bases,
        // the number of mod bases is 100 +- 7, so that's like a standard deviation
        // of 7% of the mean. So, most of the time we end up in the interval `(0.43..=0.57)`.
        // So `(0.2..=0.8)` is quite lax actually.

        let basecall_qual = df.column("basecall_qual")?.u32()?;
        assert!(
            basecall_qual
                .iter()
                .all(|x| (30..=40).contains(&x.unwrap()))
        );
        // as we are averaging 10^(-Q/10), the average basecalling quality will
        // be in a much tighter range around 30.. but I haven't calculated what
        // this is.. We are just using a very lax 30..=40 here.

        Ok(())
    }

    /// Test that `run_df` works as expected when we generate two types of reads,
    /// and that the statistics are as expected in the two groups of reads.
    #[test]
    #[expect(
        clippy::too_many_lines,
        clippy::arithmetic_side_effects,
        reason = "test with too many lines is ok; no chance of overflow due to small data len"
    )]
    fn run_df_with_two_types_of_mod_reads() -> Result<(), Error> {
        // Create simulation config with modifications
        let config_json = r#"{
            "contigs": {
                "number": 4,
                "len_range": [10000, 20000]
            },
            "reads": [
                {
                    "number": 100,
                    "mapq_range": [10, 20],
                    "base_qual_range": [30, 40],
                    "len_range": [0.1, 0.8],
                    "mods": [
                        {
                            "base": "C",
                            "is_strand_plus": false,
                            "mod_code": "m",
                            "win": [4],
                            "mod_range": [[0.2, 0.4]]
                        },
                        {
                            "base": "N",
                            "is_strand_plus": true,
                            "mod_code": "N",
                            "win": [4],
                            "mod_range": [[0.6, 0.8]]
                        }
                    ]
                },
                {
                    "number": 100,
                    "mapq_range": [10, 20],
                    "base_qual_range": [10, 20],
                    "len_range": [0.5, 0.6]
                }
            ]
        }"#;

        let sim = create_test_simulation(config_json)?;
        let df = run_window_analysis_with_threshold(&sim, 200, 100)?;

        // Verify the dataframe is NOT empty (should have data rows with mods)
        assert!(
            df.height() > 0,
            "DataFrame should have rows when modifications are present"
        );

        assert_expected_columns(&df);

        let mod_qual = df.column("win_val")?.f32()?;
        let basecall_qual = df.column("basecall_qual")?.u32()?;
        let read_id = df.column("read_id")?.str()?;
        let base = df.column("base")?.str()?;
        let mod_strand = df.column("mod_strand")?.str()?;
        let mod_type = df.column("mod_type")?.str()?;
        let ref_win_start = df.column("ref_win_start")?.i64()?;
        let ref_win_end = df.column("ref_win_end")?.i64()?;
        let win_start = df.column("win_start")?.u64()?;
        let win_end = df.column("win_end")?.u64()?;

        let mut previous_win_start: Option<u64> = None;
        let mut previous_win_end: Option<u64> = None;

        let mut sum_c_read_window_size: u64 = 0;
        let mut sum_c_ref_window_size: i64 = 0;
        let mut count_c_read_window_size: u32 = 0;
        let mut count_c_ref_window_size: u32 = 0;

        for k in izip!(
            mod_qual,
            basecall_qual,
            read_id,
            base,
            mod_strand,
            mod_type,
            ref_win_start,
            ref_win_end,
            win_start,
            win_end
        ) {
            assert!(
                k.2.unwrap().starts_with("0."),
                "only 1st read group comes through, 2nd read group has no mods"
            );
            assert!(
                (30..=40).contains(&k.1.unwrap()),
                "base call quals are 30 to 40"
            );
            if k.5 == Some("N") {
                assert_eq!(k.4.unwrap(), "+");
                assert_eq!(k.3.unwrap(), "N");
                assert_eq!(
                    k.9.unwrap() - k.8.unwrap(),
                    200,
                    "N mod should produce 200 bp windows"
                );
                assert_eq!(k.0, Some(1f32));

                let ref_st = k.6.unwrap();
                let ref_en = k.7.unwrap();

                if ref_en > -1 && ref_st > -1 {
                    assert_eq!(
                        ref_en - ref_st,
                        200,
                        "N mod should produce 200 bp windows on ref on mapped reads"
                    );
                }

                // if previous window data is available, and we are not at a transition from one
                // read to another, check if the windows have slid correctly
                match (previous_win_start, previous_win_end) {
                    (Some(s), Some(e)) if k.8.unwrap() != 0 => {
                        assert_eq!(k.9.unwrap() - e, 100, "100 bp sliding window on N mod");
                        assert_eq!(k.8.unwrap() - s, 100, "100 bp sliding window on N mod");
                    }
                    _ => {}
                }
                previous_win_start = k.8;
                previous_win_end = k.9;
            } else if k.5 == Some("m") {
                assert_eq!(k.4.unwrap(), "-");
                assert_eq!(k.3.unwrap(), "C");
                assert_eq!(k.0, Some(0f32));

                count_c_read_window_size += 1;
                sum_c_read_window_size += k.9.unwrap() - k.8.unwrap();

                let ref_st = k.6.unwrap();
                let ref_en = k.7.unwrap();
                if ref_en > -1 && ref_st > -1 {
                    sum_c_ref_window_size += k.7.unwrap() - k.6.unwrap();
                    count_c_ref_window_size += 1;
                }
            } else {
                unreachable!("Only N or m mods are present!");
            }
        }

        // Tolerate some spread around 800 (200 base windows with 25% chance of each base = 800).
        // I think the spread we have taken here is quite lax actually...
        assert!(
            (700..=900).contains(
                &sum_c_read_window_size
                    .checked_div(count_c_read_window_size.into())
                    .unwrap()
            )
        );
        assert!(
            (700..=900).contains(
                &sum_c_ref_window_size
                    .checked_div(count_c_ref_window_size.into())
                    .unwrap()
            )
        );

        Ok(())
    }

    /// Helper to run JSON window analysis with `threshold_and_mean` aggregation function
    fn run_json_window_analysis_with_threshold(
        sim: &TempBamSimulation,
        win: usize,
        step: usize,
    ) -> Result<Vec<serde_json::Value>, Error> {
        let mut bam_reader = bam::Reader::from_path(sim.bam_path())?;
        let bam_records = bam_reader.rc_records();

        let window_options: InputWindowing =
            serde_json::from_str(&format!("{{\"win\": {win}, \"step\": {step}}}"))?;
        let mods = InputMods::default();

        let mut output = Vec::new();
        run_json(&mut output, bam_records, window_options, &mods, |x| {
            analysis::threshold_and_mean(x).map(Into::into)
        })?;

        let output_str = String::from_utf8(output)?;
        let parsed: Vec<serde_json::Value> = serde_json::from_str(&output_str)?;
        Ok(parsed)
    }

    /// Test that `run_json` produces an empty array for BAM files with no modification data
    #[test]
    fn run_json_empty_for_no_mods() -> Result<(), Error> {
        let config_json = r#"{
            "contigs": {
                "number": 2,
                "len_range": [100, 200]
            },
            "reads": [{
                "number": 100,
                "mapq_range": [10, 20],
                "base_qual_range": [10, 20],
                "len_range": [0.1, 0.8]
            }]
        }"#;

        let sim = create_test_simulation(config_json)?;
        let entries = run_json_window_analysis_with_threshold(&sim, 2, 1)?;

        // One JSON record per BAM read
        assert_eq!(entries.len(), 100, "should have one JSON record per read");

        // Each read still appears in the JSON output, but mod_table.data should be empty
        for entry in &entries {
            let mod_table = entry["mod_table"].as_array().unwrap();
            for mod_entry in mod_table {
                let data = mod_entry["data"].as_array().unwrap();
                assert_eq!(
                    data.len(),
                    0,
                    "mod_table data should be empty for reads with no modifications"
                );
            }
        }

        Ok(())
    }

    /// Test that `run_json` produces entries with modification data
    #[test]
    fn run_json_with_mods() -> Result<(), Error> {
        let config_json = r#"{
            "contigs": {
                "number": 4,
                "len_range": [100, 200]
            },
            "reads": [{
                "number": 1000,
                "mapq_range": [10, 20],
                "base_qual_range": [20, 30],
                "len_range": [0.1, 0.8],
                "mods": [{
                    "base": "T",
                    "is_strand_plus": true,
                    "mod_code": "T",
                    "win": [4],
                    "mod_range": [[0.1, 0.2]]
                }]
            }]
        }"#;

        let sim = create_test_simulation(config_json)?;
        let entries = run_json_window_analysis_with_threshold(&sim, 2, 1)?;

        // One JSON record per BAM read
        assert_eq!(entries.len(), 1000, "should have one JSON record per read");

        let mut total_windows = 0usize;
        for entry in &entries {
            let mod_table = entry["mod_table"].as_array().unwrap();
            for mod_entry in mod_table {
                let data = mod_entry["data"].as_array().unwrap();
                for window in data {
                    total_windows += 1;
                    let win_val = window[2].as_f64().unwrap();
                    #[expect(
                        clippy::float_cmp,
                        reason = "exact 0.0 comparison is safe for thresholded-then-meaned low-probability mods"
                    )]
                    {
                        assert_eq!(
                            win_val, 0.0,
                            "win_val should be 0.0 for low mod probability"
                        );
                    }
                    let basecall_qual = window[3].as_u64().unwrap();
                    assert!(
                        (20..=30).contains(&basecall_qual),
                        "basecall_qual should be in range 20..=30"
                    );
                }
            }
        }

        assert!(
            total_windows > 0,
            "Should have at least one window across all reads"
        );

        Ok(())
    }

    /// Test that `run_json` works with non-perfectly aligned reads (deletions, insertions, mismatches)
    #[test]
    fn run_json_with_mods_and_non_perfectly_aligned_reads() -> Result<(), Error> {
        let config_json = r#"{
            "contigs": {
                "number": 4,
                "len_range": [100000, 200000]
            },
            "reads": [{
                "number": 100,
                "mapq_range": [10, 20],
                "base_qual_range": [30, 40],
                "len_range": [0.1, 0.8],
                "delete": [0.5, 0.7],
                "insert_middle": "ATCGAATTGGAA",
                "mismatch": 0.2,
                "mods": [{
                    "base": "C",
                    "is_strand_plus": false,
                    "mod_code": "m",
                    "win": [4],
                    "mod_range": [[0.2, 0.8]]
                }]
            }]
        }"#;

        let sim = create_test_simulation(config_json)?;
        let entries = run_json_window_analysis_with_threshold(&sim, 200, 100)?;

        // One JSON record per BAM read
        assert_eq!(entries.len(), 100, "should have one JSON record per read");

        let mut total_windows = 0usize;
        for entry in &entries {
            let mod_table = entry["mod_table"].as_array().unwrap();
            for mod_entry in mod_table {
                let data = mod_entry["data"].as_array().unwrap();
                for window in data {
                    total_windows += 1;
                    let win_val = window[2].as_f64().unwrap();
                    assert!(
                        (0.2..=0.8).contains(&win_val),
                        "win_val {win_val} should be in range 0.2..=0.8"
                    );
                    let basecall_qual = window[3].as_u64().unwrap();
                    assert!(
                        (30..=40).contains(&basecall_qual),
                        "basecall_qual should be in range 30..=40"
                    );
                }
            }
        }

        assert!(
            total_windows > 0,
            "Should have at least one window across all reads"
        );

        Ok(())
    }

    /// Test that `run_json` works with two types of mod reads and validates statistics per group
    #[test]
    #[expect(clippy::too_many_lines, reason = "test with too many lines is ok")]
    fn run_json_with_two_types_of_mod_reads() -> Result<(), Error> {
        let config_json = r#"{
            "contigs": {
                "number": 4,
                "len_range": [10000, 20000]
            },
            "reads": [
                {
                    "number": 100,
                    "mapq_range": [10, 20],
                    "base_qual_range": [30, 40],
                    "len_range": [0.1, 0.8],
                    "mods": [
                        {
                            "base": "C",
                            "is_strand_plus": false,
                            "mod_code": "m",
                            "win": [4],
                            "mod_range": [[0.2, 0.4]]
                        },
                        {
                            "base": "N",
                            "is_strand_plus": true,
                            "mod_code": "N",
                            "win": [4],
                            "mod_range": [[0.6, 0.8]]
                        }
                    ]
                },
                {
                    "number": 100,
                    "mapq_range": [10, 20],
                    "base_qual_range": [10, 20],
                    "len_range": [0.5, 0.6]
                }
            ]
        }"#;

        let sim = create_test_simulation(config_json)?;
        let entries = run_json_window_analysis_with_threshold(&sim, 200, 100)?;

        // One JSON record per BAM read (100 with mods + 100 without)
        assert_eq!(entries.len(), 200, "should have one JSON record per read");

        let mut previous_win_start: Option<i64> = None;
        let mut previous_win_end: Option<i64> = None;

        let mut sum_c_read_window_size: u64 = 0;
        let mut sum_c_ref_window_size: i64 = 0;
        let mut count_c_read_window_size: u32 = 0;
        let mut count_c_ref_window_size: u32 = 0;

        for entry in &entries {
            let read_id = entry["read_id"].as_str().unwrap();
            let mod_table = entry["mod_table"].as_array().unwrap();

            for mod_entry in mod_table {
                let base = mod_entry["base"].as_str().unwrap();
                let is_strand_plus = mod_entry["is_strand_plus"].as_bool().unwrap();
                let mod_code = &mod_entry["mod_code"];
                let data = mod_entry["data"].as_array().unwrap();

                for window in data {
                    let win_start = window[0].as_i64().unwrap();
                    let win_end = window[1].as_i64().unwrap();
                    let win_val = window[2].as_f64().unwrap();
                    let basecall_qual = window[3].as_u64().unwrap();
                    let ref_win_start = window[4].as_i64().unwrap();
                    let ref_win_end = window[5].as_i64().unwrap();

                    assert!(
                        read_id.starts_with("0."),
                        "only 1st read group comes through, 2nd read group has no mods"
                    );
                    assert!(
                        (30..=40).contains(&basecall_qual),
                        "base call quals are 30 to 40"
                    );

                    if mod_code == "N" {
                        assert!(is_strand_plus);
                        assert_eq!(base, "N");
                        assert_eq!(
                            win_end - win_start,
                            200,
                            "N mod should produce 200 bp windows"
                        );
                        #[expect(
                            clippy::float_cmp,
                            reason = "exact 1.0 comparison is safe for thresholded high-probability mods"
                        )]
                        {
                            assert_eq!(win_val, 1.0);
                        }

                        if ref_win_end > -1 && ref_win_start > -1 {
                            assert_eq!(
                                ref_win_end - ref_win_start,
                                200,
                                "N mod should produce 200 bp windows on ref on mapped reads"
                            );
                        }

                        // check sliding window consistency
                        match (previous_win_start, previous_win_end) {
                            (Some(s), Some(e)) if win_start != 0 => {
                                assert_eq!(win_end - e, 100, "100 bp sliding window on N mod");
                                assert_eq!(win_start - s, 100, "100 bp sliding window on N mod");
                            }
                            _ => {}
                        }
                        previous_win_start = Some(win_start);
                        previous_win_end = Some(win_end);
                    } else if mod_code == "m" {
                        assert!(!is_strand_plus);
                        assert_eq!(base, "C");
                        #[expect(
                            clippy::float_cmp,
                            reason = "exact 0.0 comparison is safe for thresholded low-probability mods"
                        )]
                        {
                            assert_eq!(win_val, 0.0);
                        }

                        count_c_read_window_size += 1;
                        #[expect(
                            clippy::cast_sign_loss,
                            reason = "window coordinates are non-negative and small"
                        )]
                        {
                            sum_c_read_window_size += (win_end - win_start) as u64;
                        }

                        if ref_win_end > -1 && ref_win_start > -1 {
                            sum_c_ref_window_size += ref_win_end - ref_win_start;
                            count_c_ref_window_size += 1;
                        }
                    } else {
                        unreachable!("Only N or m mods are present!");
                    }
                }
            }
        }

        // Tolerate some spread around 800 (200 base windows with 25% chance of each base = 800).
        assert!(
            (700..=900).contains(
                &sum_c_read_window_size
                    .checked_div(count_c_read_window_size.into())
                    .unwrap()
            )
        );
        assert!(
            (700..=900).contains(
                &sum_c_ref_window_size
                    .checked_div(count_c_ref_window_size.into())
                    .unwrap()
            )
        );

        Ok(())
    }

    /// Test that `run` (TSV) outputs only the header line when there are zero reads
    #[test]
    fn run_tsv_header_only_for_zero_reads() -> Result<(), Error> {
        let config_json = r#"{
            "contigs": {
                "number": 2,
                "len_range": [100, 200]
            },
            "reads": []
        }"#;

        let sim = create_test_simulation(config_json)?;
        let mut output = Vec::new();
        let mut bam_reader = bam::Reader::from_path(sim.bam_path())?;
        let bam_records = bam_reader.rc_records();
        let window_options: InputWindowing =
            serde_json::from_str("{\"win\": 2, \"step\": 1}").unwrap();
        let mods = InputMods::default();

        run(&mut output, bam_records, window_options, &mods, |x| {
            analysis::threshold_and_mean(x).map(Into::into)
        })?;

        let output_str = String::from_utf8(output)?;
        let lines: Vec<&str> = output_str.lines().collect();
        assert_eq!(lines.len(), 1, "should only have the header line");
        assert!(
            lines
                .first()
                .expect("already asserted len == 1")
                .starts_with("#contig"),
            "the single line should be the header"
        );

        Ok(())
    }

    /// Test that `run_json` outputs an empty JSON array when there are zero reads
    #[test]
    fn run_json_empty_array_for_zero_reads() -> Result<(), Error> {
        let config_json = r#"{
            "contigs": {
                "number": 2,
                "len_range": [100, 200]
            },
            "reads": []
        }"#;

        let sim = create_test_simulation(config_json)?;
        let entries = run_json_window_analysis_with_threshold(&sim, 2, 1)?;

        assert_eq!(
            entries.len(),
            0,
            "should have zero JSON records for zero reads"
        );

        Ok(())
    }
}
