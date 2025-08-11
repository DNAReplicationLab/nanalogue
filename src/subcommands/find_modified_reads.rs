//! # Find modified reads using criteria on windowed mod data
//!
//! In this module, we window data along molecules, and then use
//! filtration criteria on these windows using user-supplied parameters
//! and output these reads.

use crate::{CurrRead, Error, F32Bw0and1, InputWindowing};
use rust_htslib::bam::Record;
use std::io::{self, Write};
use std::rc::Rc;

/// Finds read ids of molecules that fit filtration criteria on
/// windowed modification data.
pub fn run<F, G, D>(
    bam_records: D,
    window_options: InputWindowing,
    window_function: F,
    window_filter: G,
) -> Result<bool, Error>
where
    F: Fn(&[u8]) -> Result<F32Bw0and1, Error>,
    G: Fn(&Vec<F32Bw0and1>) -> bool,
    D: IntoIterator<Item = Result<Rc<Record>, rust_htslib::errors::Error>>,
{
    // This apparently helps writing to the terminal faster,
    // according to https://rust-cli.github.io/book/tutorial/output.html
    let stdout = io::stdout();
    let mut handle = io::BufWriter::new(stdout);

    // Go record by record in the BAM file,
    for r in bam_records {
        // read records
        let mut curr_read_state = CurrRead::default();
        let record = r?;
        curr_read_state.set_read_id(&record)?;

        // set the modified read state
        if let Some(v) = window_options.mod_strand {
            curr_read_state.set_mod_data_restrictive(
                &record,
                window_options.mod_prob_filter,
                |_, s, t| t == window_options.tag && s == char::from(v),
            );
        } else {
            curr_read_state.set_mod_data_restrictive(
                &record,
                window_options.mod_prob_filter,
                |_, _, t| t == window_options.tag,
            );
        }

        // trim ends of reads if requested
        match window_options.trim_read_ends {
            0 => {}
            v => curr_read_state.filter_starts_at_read_ends(v),
        }

        // apply our windowing function and then the windowing filter
        if match curr_read_state.windowed_mod_data(
            &window_function,
            window_options.win.get().try_into()?,
            window_options.step.get().try_into()?,
            window_options.tag,
        )? {
            Some(v) => window_filter(&v),
            None => false,
        } {
            writeln!(handle, "{}", curr_read_state.read_id()?)?;
        }
    }

    Ok(true)
}
