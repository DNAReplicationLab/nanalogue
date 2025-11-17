//! # Window modification data on reads
//!
//! In this module, we window data along molecules, and then output
//! these windows

use crate::{CurrRead, Error, F32AbsValAtMost1, InputMods, InputWindowing, ModChar, OptionalTag};
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

    // Go record by record in the BAM file,
    for r in bam_records {
        // read records
        let record = r?;

        // set data in records
        let curr_read_state = CurrRead::default()
            .try_from_only_alignment(&record)?
            .set_mod_data_restricted_options(&record, mods)?;
        let qname = curr_read_state.read_id()?;
        let strand = curr_read_state.strand();
        let contig = if curr_read_state.read_state().is_unmapped() {
            "."
        } else {
            curr_read_state.contig_name()?
        };

        // read and window modification data, then print the output
        for base_mod in &curr_read_state.mod_data().0.base_mods {
            let mod_data = &base_mod.ranges.qual;
            let starts = &base_mod.ranges.starts;
            let ends = &base_mod.ranges.ends;
            let ref_starts = &base_mod.ranges.reference_starts;
            let ref_ends = &base_mod.ranges.reference_ends;
            let base = base_mod.modified_base as char;
            let mod_strand = base_mod.strand;
            let mod_type = ModChar::new(base_mod.modification_type);
            if let Some(v) = mod_data.len().checked_sub(win_size) {
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
                    let win_start = starts
                        .get(window_idx)
                        .expect("window_idx is valid")
                        .ok_or_else(|| {
                            Error::InvalidState("Missing sequence start position".to_string())
                        })?;
                    let win_end = ends
                        .get(window_idx..)
                        .expect("window_idx <= v where v = len - win_size")
                        .get(0..win_size)
                        .expect("no error as we've checked data len >= win size")
                        .last()
                        .expect("no error as we've checked data len >= win size")
                        .ok_or_else(|| {
                            Error::InvalidState("Missing sequence end position".to_string())
                        })?;

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
                        .unwrap_or(INVALID_REF_POS);
                    writeln!(
                        handle,
                        "{contig}\t{ref_win_start}\t{ref_win_end}\t{qname}\t{win_val}\t{strand}\t\
{base}\t{mod_strand}\t{mod_type}\t{win_start}\t{win_end}"
                    )?;
                }
            }
        }
    }

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
        threshold: Option<f32>,
        expected_output_file: &str,
    ) -> Result<(), Error> {
        // Set input, output, options
        let mut output = Vec::new();
        let mut bam_reader = bam::Reader::from_path("./examples/example_1.bam")?;
        let bam_records = bam_reader.records().map(|r| r.map(Rc::new));
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
        run_window_reads_test_with_threshold(None, "./examples/example_1_window_reads")
    }

    #[test]
    fn window_reads_example_1_gt_0pt4() -> Result<(), Error> {
        run_window_reads_test_with_threshold(Some(0.4), "./examples/example_1_window_reads_gt_0pt4")
    }

    #[test]
    fn window_reads_example_1_gt_0pt8() -> Result<(), Error> {
        run_window_reads_test_with_threshold(Some(0.8), "./examples/example_1_window_reads_gt_0pt8")
    }
}
