//! # Find modified reads using criteria on windowed mod data
//!
//! In this module, we window data along molecules, and then use
//! filtration criteria on these windows using user-supplied parameters
//! and output these reads.

use crate::{CurrRead, Error, F32Bw0and1, InputWindowingRestricted};
use rust_htslib::bam::Record;
use std::rc::Rc;

/// Finds read ids of molecules that fit filtration criteria on
/// windowed modification data.
pub fn run<W, F, G, D>(
    handle: &mut W,
    bam_records: D,
    window_options: InputWindowingRestricted,
    window_function: F,
    window_filter: G,
) -> Result<bool, Error>
where
    W: std::io::Write,
    F: Fn(&[u8]) -> Result<F32Bw0and1, Error>,
    G: Fn(&Vec<F32Bw0and1>) -> bool,
    D: IntoIterator<Item = Result<Rc<Record>, rust_htslib::errors::Error>>,
{
    let trim_end_bp = window_options.trim_read_ends;
    let mut curr_read_state = CurrRead::default();

    // Go record by record in the BAM file,
    for r in bam_records {
        // read records
        let record = r?;
        curr_read_state.reset();
        curr_read_state.set_read_id(&record)?;
        let seq_len: usize = curr_read_state.set_seq_len(&record)?.try_into().unwrap();

        // set the modified read state
        match (&trim_end_bp, &window_options.mod_strand) {
            (0, &Some(v)) => curr_read_state.set_mod_data_restricted(
                &record,
                window_options.mod_prob_filter,
                |&_| true,
                |_, &s, &t| t == window_options.tag && s == char::from(v),
            ),
            (&w, &Some(v)) => curr_read_state.set_mod_data_restricted(
                &record,
                window_options.mod_prob_filter,
                |&x| x >= w && x <= seq_len - w,
                |_, &s, &t| t == window_options.tag && s == char::from(v),
            ),
            (0, None) => curr_read_state.set_mod_data_restricted(
                &record,
                window_options.mod_prob_filter,
                |&_| true,
                |_, _, &t| t == window_options.tag,
            ),
            (&w, None) => curr_read_state.set_mod_data_restricted(
                &record,
                window_options.mod_prob_filter,
                |&x| x >= w && x <= seq_len - w,
                |_, _, &t| t == window_options.tag,
            ),
        }

        // apply our windowing function and then the windowing filter
        if match curr_read_state.windowed_mod_data_restricted(
            &window_function,
            window_options.win_params.win.get().try_into()?,
            window_options.win_params.step.get().try_into()?,
            window_options.tag,
        )? {
            v if !v.is_empty() => window_filter(&v),
            _ => false,
        } {
            writeln!(handle, "{}", curr_read_state.read_id()?)?;
        }
    }

    Ok(true)
}
