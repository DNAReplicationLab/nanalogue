//! # Find modified reads using criteria on windowed mod data
//!
//! In this module, we window data along molecules, and then use
//! filtration criteria on these windows using user-supplied parameters
//! and output these reads.

use crate::{CurrRead, Error, F32Bw0and1, InputMods, InputWindowing, RequiredTag};
use rust_htslib::bam::Record;
use std::rc::Rc;

/// Finds read ids of molecules that fit filtration criteria on
/// windowed modification data.
pub fn run<W, F, G, D>(
    handle: &mut W,
    bam_records: D,
    window_options: InputWindowing,
    mod_options: InputMods<RequiredTag>,
    window_function: F,
    window_filter: G,
) -> Result<(), Error>
where
    W: std::io::Write,
    F: Fn(&[u8]) -> Result<F32Bw0and1, Error>,
    G: Fn(&Vec<F32Bw0and1>) -> bool,
    D: IntoIterator<Item = Result<Rc<Record>, rust_htslib::errors::Error>>,
{
    // Go record by record in the BAM file,
    for r in bam_records {
        // read records
        let record = r?;
        let curr_read_state = CurrRead::default().try_from_only_alignment(&record)?;
        let read_id = String::from(curr_read_state.read_id()?);
        // apply our windowing function and then the windowing filter
        if match curr_read_state
            .set_mod_data_restricted_options(&record, &mod_options)?
            .windowed_mod_data_restricted(&window_function, window_options, mod_options.tag())?
        {
            v if !v.is_empty() => window_filter(&v),
            _ => false,
        } {
            writeln!(handle, "{}", read_id)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        InputMods, InputWindowing, ModChar, RequiredTag, ThresholdState, nanalogue_bam_reader,
    };
    use rust_htslib::bam::Read as BamRead;
    use std::num::NonZeroU32;

    /// Helper function that runs find_modified_reads with specified parameters
    /// and returns the output as a vector of read IDs
    fn run_find_modified_reads_test<F, G>(
        bam_path: &str,
        window_size: u32,
        step_size: u32,
        tag: char,
        window_function: F,
        window_filter: G,
    ) -> Result<Vec<String>, Error>
    where
        F: Fn(&[u8]) -> Result<F32Bw0and1, Error>,
        G: Fn(&Vec<F32Bw0and1>) -> bool,
    {
        // Open the BAM file and collect records
        let mut reader = nanalogue_bam_reader(bam_path)?;
        let records: Vec<_> = reader.records().map(|r| r.map(Rc::new)).collect();

        // Set up windowing options
        let window_options = InputWindowing {
            win: NonZeroU32::new(window_size).unwrap(),
            step: NonZeroU32::new(step_size).unwrap(),
        };

        // Set up mod options
        let mod_options = InputMods::<RequiredTag> {
            tag: RequiredTag {
                tag: ModChar::new(tag),
            },
            mod_strand: None,
            mod_prob_filter: ThresholdState::default(),
            trim_read_ends_mod: 0,
            base_qual_filter_mod: 0,
            mod_region: None,
            region_bed3: None,
        };

        // Run the function
        let mut output = Vec::new();
        run(
            &mut output,
            records,
            window_options,
            mod_options,
            window_function,
            window_filter,
        )?;

        // Parse output into vector of read IDs
        let output_str = String::from_utf8(output).expect("Invalid UTF-8 output");
        let mut read_ids: Vec<String> = output_str.lines().map(|s| s.to_string()).collect();
        read_ids.sort();

        Ok(read_ids)
    }

    /// Helper function that computes the mean of modification probabilities
    fn mean_window_function(mod_data: &[u8]) -> Result<F32Bw0and1, Error> {
        if mod_data.is_empty() {
            return F32Bw0and1::new(0.0);
        }
        let sum: f32 = mod_data.iter().map(|&x| x as f32).sum();
        let mean = sum / (mod_data.len() as f32);
        F32Bw0and1::new(mean / 255.0)
    }

    /// Tests the find_modified_reads::run function with a simple filter
    /// that looks for reads with average modification probabilities above a threshold
    #[test]
    fn test_find_modified_reads_basic() -> Result<(), Error> {
        // Define window filter: keep reads with at least one window above threshold
        let window_filter =
            |windows: &Vec<F32Bw0and1>| -> bool { windows.iter().any(|&w| w.val() > 0.7) };

        let read_ids = run_find_modified_reads_test(
            "./examples/example_1.bam",
            2,
            1,
            'T',
            mean_window_function,
            window_filter,
        )?;

        // Based on example_1.bam, we expect certain reads to pass the filter
        // Read 1 has high quality mods, should pass
        // this read is repeated, so we should get two copies of this read.
        assert_eq!(read_ids, vec!["a4f36092-b4d5-47a9-813e-c22c3b477a0c"; 2]);

        Ok(())
    }

    /// Tests find_modified_reads with a filter that rejects all reads
    #[test]
    fn test_find_modified_reads_reject_all() -> Result<(), Error> {
        // Filter that rejects everything
        let window_filter = |_windows: &Vec<F32Bw0and1>| -> bool { false };

        let read_ids = run_find_modified_reads_test(
            "./examples/example_1.bam",
            2,
            1,
            'T',
            mean_window_function,
            window_filter,
        )?;

        assert!(
            read_ids.is_empty(),
            "Expected no reads to pass the reject-all filter"
        );

        Ok(())
    }

    /// Tests find_modified_reads with a filter that accepts all reads
    #[test]
    fn test_find_modified_reads_accept_all() -> Result<(), Error> {
        // Filter that accepts everything with non-empty windows
        let window_filter = |windows: &Vec<F32Bw0and1>| -> bool { !windows.is_empty() };

        let read_ids = run_find_modified_reads_test(
            "./examples/example_1.bam",
            2,
            1,
            'T',
            mean_window_function,
            window_filter,
        )?;

        // check that all reads come through
        assert_eq!(
            read_ids,
            vec![
                "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575",
                "a4f36092-b4d5-47a9-813e-c22c3b477a0c",
                "a4f36092-b4d5-47a9-813e-c22c3b477a0c",
                "fffffff1-10d2-49cb-8ca3-e8d48979001b"
            ]
        );

        Ok(())
    }
}
