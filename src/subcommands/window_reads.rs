//! # Window modification data on reads
//!
//! In this module, we window data along molecules, and then output
//! these windows

use crate::{CurrRead, Error, F32AbsValAtMost1, InputMods, InputWindowing, ModChar, OptionalTag};
use polars::prelude::*;
use rust_htslib::bam::Record;
use std::rc::Rc;

/// Windowed modification data along molecules
///
/// # Errors
/// Returns an error if BAM record reading, or output writing fails.
///
#[expect(
    clippy::missing_panics_doc,
    reason = "window slice indexing should not fail due to bounds checking"
)]
#[expect(
    clippy::too_many_lines,
    reason = "windowing logic is simple but just takes up many lines"
)]
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
    // constant to mark windows with basecalled coordinates but no reference coordinates.
    const INVALID_REF_POS: i64 = -1;

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
        #[expect(
            clippy::type_complexity,
            reason = "I think a tuple of 5 `Vec` is fine if its readable"
        )]
        for base_mod in &curr_read_state.mod_data().0.base_mods {
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
                    let win_start = starts.get(window_idx).expect("window_idx is valid");
                    let win_end = ends
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
                            .get(usize::try_from(*win_start)?..usize::try_from(*win_end)?)
                            .expect("no error as `win_start`, `win_end` in range");
                        if quals.first() == Some(&255u8) {
                            // BAM format is such that all values are 255, or values are between
                            // 0 and 93. So if we see one 255, we can just return 255
                            255u8
                        } else {
                            // we do an average using the probability of errors,
                            // and not the Q scores directly.
                            let quals_min = quals.iter().min().expect("no error");
                            let data_size =
                                f64::from(i32::try_from(*win_end)? - i32::try_from(*win_start)?);
                            let x = quals.iter().fold(0f64, |acc, x| {
                                acc + (10f64).powf(-0.1f64 * f64::from(x - quals_min))
                            });
                            quals_min + ((-10f64 * f64::log10(x / data_size)).round() as u8)
                        }
                    };
                    writeln!(
                        handle,
                        "{contig}\t{ref_win_start}\t{ref_win_end}\t{qname}\t{win_val}\t{strand}\t\
                    {base}\t{mod_strand}\t{mod_type}\t{win_start}\t{win_end}\t{mean_base_qual}",
                    )?;
                }
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
                if previous_win_start.is_some() && previous_win_end.is_some() && k.8.unwrap() != 0 {
                    assert_eq!(
                        k.9.unwrap() - previous_win_end.unwrap(),
                        100,
                        "100 bp sliding window on N mod"
                    );
                    assert_eq!(
                        k.8.unwrap() - previous_win_start.unwrap(),
                        100,
                        "100 bp sliding window on N mod"
                    );
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
}
