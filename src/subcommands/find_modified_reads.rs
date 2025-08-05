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
pub fn run<F, D>(
    bam_records: D,
    tag: ModChar,
    win: NonZeroU32,
    step: NonZeroU32,
    dens_filter: F,
    invert: bool,
) -> Result<bool, Error>
where
    F: Fn(&F32Bw0and1) -> bool,
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

        // catch if one window meets our criterion,
        // and react accordingly using invert's state
        if match curr_read_state.windowed_mod_data(
            win.get().try_into()?,
            step.get().try_into()?,
            tag,
            ThresholdState::GtEq(128),
        )? {
            Some(v) => !(v.iter().any(|k| !dens_filter(k)) ^ invert),
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
