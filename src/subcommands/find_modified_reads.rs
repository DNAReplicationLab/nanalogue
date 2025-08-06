//! # Find modified reads using criteria on windowed mod data
//!
//! In this module, we window data along molecules, and then use
//! filtration criteria on these windows using user-supplied parameters
//! and output these reads.

use crate::{CurrRead, Error, F32Bw0and1, ModChar, ThresholdState};
use rust_htslib::bam::Record;
use std::num::NonZeroU32;
use std::rc::Rc;

/// Finds read ids of molecules that fit filtration criteria on
/// windowed modification data.
pub fn run<F, G, D>(
    bam_records: D,
    tag: ModChar,
    win: NonZeroU32,
    step: NonZeroU32,
    window_function: F,
    window_filter: G,
) -> Result<bool, Error>
where
    F: Fn(&[u8]) -> Result<F32Bw0and1, Error>,
    G: Fn(Vec<F32Bw0and1>) -> bool,
    D: IntoIterator<Item = Result<Rc<Record>, rust_htslib::errors::Error>>,
{
    // prepare output string
    let mut output_string = String::from("");

    // Go record by record in the BAM file,
    for r in bam_records {
        // read records
        let mut curr_read_state = CurrRead::default();
        let record = r?;
        curr_read_state.set_read_id(&record)?;

        // set the modified read state
        curr_read_state.set_mod_data_one_tag(&record, ThresholdState::GtEq(0), tag);

        // apply our windowing function and then the windowing filter
        if match curr_read_state.windowed_mod_data(
            &window_function,
            win.get().try_into()?,
            step.get().try_into()?,
            tag,
        )? {
            Some(v) => window_filter(v),
            None => false,
        } {
            output_string = output_string + curr_read_state.read_id()? + "\n";
        }
    }

    if !output_string.is_empty() {
        print!("{output_string}");
    }

    Ok(true)
}
