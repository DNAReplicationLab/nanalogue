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
) -> Result<bool, Error>
where
    W: std::io::Write,
    F: Fn(&[u8]) -> Result<F32Bw0and1, Error>,
    G: Fn(&Vec<F32Bw0and1>) -> bool,
    D: IntoIterator<Item = Result<Rc<Record>, rust_htslib::errors::Error>>,
{
    let mut curr_read_state = CurrRead::default();

    // Go record by record in the BAM file,
    for r in bam_records {
        // read records
        let record = r?;
        curr_read_state.reset();
        curr_read_state.set_read_state(&record)?;
        curr_read_state.set_seq_len(&record)?;
        curr_read_state.set_read_id(&record)?;
        match curr_read_state.set_mod_data_restricted_options(&record, &mod_options) {
            Ok(_) | Err(Error::NoModInfo) => {}
            Err(e) => return Err(e),
        };
        // apply our windowing function and then the windowing filter
        if match curr_read_state.windowed_mod_data_restricted(
            &window_function,
            window_options.win.get().try_into()?,
            window_options.step.get().try_into()?,
            mod_options.tag(),
        )? {
            v if !v.is_empty() => window_filter(&v),
            _ => false,
        } {
            writeln!(handle, "{}", curr_read_state.read_id()?)?;
        }
    }

    Ok(true)
}
