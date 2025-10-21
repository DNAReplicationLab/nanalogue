//! # Window modification data on reads
//!
//! In this module, we window data along molecules, and then output
//! these windows

use crate::{
    CurrRead, Error, F32AbsValBelow1, InputMods, InputWindowing, ModChar, OptionalTag, ReadState,
};
use fibertools_rs::utils::basemods::BaseMods;
use rust_htslib::bam::Record;
use std::rc::Rc;

/// Windowed modification data along molecules
pub fn run<W, F, D>(
    handle: &mut W,
    bam_records: D,
    window_options: InputWindowing,
    mods: InputMods<OptionalTag>,
    window_function: F,
) -> Result<(), Error>
where
    W: std::io::Write,
    F: Fn(&[u8]) -> Result<F32AbsValBelow1, Error>,
    D: IntoIterator<Item = Result<Rc<Record>, rust_htslib::errors::Error>>,
{
    // Get windowing parameters
    let win_size: usize = window_options.win.get().try_into()?;
    let slide_size: usize = window_options.step.get().try_into()?;

    // constant to mark windows with basecalled coordinates but no reference coordinates.
    const INVALID_REF_POS: i64 = -1;

    // Go record by record in the BAM file,
    for r in bam_records {
        // read records
        let record = r?;

        // set data in records
        let curr_read_state = CurrRead::default()
            .try_from_only_alignment(&record)?
            .set_mod_data_restricted_options(&record, &mods)?;
        let qname = curr_read_state.read_id()?;
        let strand = curr_read_state.strand();
        let contig = match curr_read_state.read_state() {
            ReadState::Unmapped => ".",
            _ => curr_read_state.contig_name()?,
        };

        // read and window modification data, then print the output
        let (BaseMods { base_mods }, _) = curr_read_state.mod_data();
        for base_mod in base_mods {
            let mod_data = &base_mod.ranges.qual;
            let starts = &base_mod.ranges.starts;
            let ends = &base_mod.ranges.ends;
            let ref_starts = &base_mod.ranges.reference_starts;
            let ref_ends = &base_mod.ranges.reference_ends;
            let base = base_mod.modified_base as char;
            let mod_strand = base_mod.strand;
            let mod_type = ModChar::new(base_mod.modification_type);
            if win_size <= mod_data.len() {
                for window_idx in (0..=mod_data.len() - win_size).step_by(slide_size) {
                    let window_end_idx = window_idx + win_size;
                    let win_val = match window_function(&mod_data[window_idx..window_end_idx]) {
                        Ok(val) => val,
                        Err(e) => {
                            eprintln!(
                                "Warning: Skipping window at {}:{}-{} due to error: {}",
                                qname, window_idx, window_end_idx, e
                            );
                            continue;
                        }
                    };
                    // there is no way to trigger the errors below as we control how CurrRead is
                    // populated quite strictly. Nevertheless, I am leaving these in for
                    // future-proofing.
                    let win_start = starts[window_idx].ok_or_else(|| {
                        Error::InvalidState("Missing sequence start position".to_string())
                    })?;
                    let win_end = ends[window_end_idx - 1].ok_or_else(|| {
                        Error::InvalidState("Missing sequence end position".to_string())
                    })?;

                    let ref_win_start = ref_starts[window_idx..window_end_idx]
                        .iter()
                        .flatten()
                        .min()
                        .copied()
                        .unwrap_or(INVALID_REF_POS);
                    let ref_win_end = ref_ends[window_idx..window_end_idx]
                        .iter()
                        .flatten()
                        .max()
                        .copied()
                        .unwrap_or(INVALID_REF_POS);
                    writeln!(
                        handle,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        contig,
                        ref_win_start,
                        ref_win_end,
                        qname,
                        win_val,
                        strand,
                        base,
                        mod_strand,
                        mod_type,
                        win_start,
                        win_end
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
    use rust_htslib::bam::{self, Read};
    use std::rc::Rc;

    /// Helper function to run window_reads tests with threshold_and_mean_and_thres_win
    ///
    /// This function encapsulates the common test setup and execution logic for window_reads tests
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
                run(&mut output, bam_records, window_options, mods, |x| {
                    threshold_and_mean(x).map(|y| y.into())
                })?;
            }
            Some(thres_val) => {
                // Use threshold_and_mean_and_thres_win with specified threshold
                let threshold = F32Bw0and1::new(thres_val).unwrap();
                run(&mut output, bam_records, window_options, mods, |x| {
                    threshold_and_mean_and_thres_win(x, threshold).map(|y| y.into())
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
    fn test_window_reads_example_1() -> Result<(), Error> {
        run_window_reads_test_with_threshold(None, "./examples/example_1_window_reads")
    }

    #[test]
    fn test_window_reads_example_1_gt_0pt4() -> Result<(), Error> {
        run_window_reads_test_with_threshold(Some(0.4), "./examples/example_1_window_reads_gt_0pt4")
    }

    #[test]
    fn test_window_reads_example_1_gt_0pt8() -> Result<(), Error> {
        run_window_reads_test_with_threshold(Some(0.8), "./examples/example_1_window_reads_gt_0pt8")
    }
}
