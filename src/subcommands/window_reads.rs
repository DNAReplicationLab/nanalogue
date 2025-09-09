//! # Window modification data on reads
//!
//! In this module, we window data along molecules, and then output
//! these windows

use crate::{CurrRead, Error, F32AbsValBelow1, InputWindowing, ModChar, ReadState, ThresholdState};
use fibertools_rs::utils::basemods::BaseMods;
use rust_htslib::bam::Record;
use std::rc::Rc;

/// Windowed modification data along molecules
pub fn run<W, F, D>(
    handle: &mut W,
    bam_records: D,
    window_options: InputWindowing,
    window_function: F,
) -> Result<bool, Error>
where
    W: std::io::Write,
    F: Fn(&[u8]) -> Result<F32AbsValBelow1, Error>,
    D: IntoIterator<Item = Result<Rc<Record>, rust_htslib::errors::Error>>,
{
    // Get windowing parameters
    let win_size: usize = window_options.win.get().try_into()?;
    let slide_size: usize = window_options.step.get().try_into()?;

    let mut curr_read_state = CurrRead::default();

    // Go record by record in the BAM file,
    for r in bam_records {
        // read records
        let record = r?;
        curr_read_state.reset();

        // set data in records
        curr_read_state.set_read_state(&record)?;
        curr_read_state.set_mod_data(&record, ThresholdState::GtEq(0), 0);
        let qname = curr_read_state.set_read_id(&record)?.to_string();
        let strand = curr_read_state.strand()?;
        let contig = match curr_read_state.read_state() {
            ReadState::Unknown => Err(Error::InvalidState(
                "unclear why we are in an invalid state!".to_string(),
            )),
            ReadState::Unmapped => Ok(".".to_string()),
            _ => Ok(String::from(curr_read_state.set_contig_name(&record)?)),
        }?;

        // read and window modification data, then print the output
        if let Ok((BaseMods { base_mods: v }, _)) = curr_read_state.mod_data() {
            for k in v {
                let mod_data = &k.ranges.qual;
                let start = &k.ranges.starts;
                let end = &k.ranges.ends;
                let ref_start = &k.ranges.reference_starts;
                let ref_end = &k.ranges.reference_ends;
                let base = k.modified_base as char;
                let mod_strand = &k.strand;
                let mod_type = ModChar::new(k.modification_type);
                if win_size <= mod_data.len() {
                    for l in (0..=mod_data.len() - win_size).step_by(slide_size) {
                        let win_val = window_function(&mod_data[l..l + win_size])?;
                        let win_start = match start[l] {
                            Some(w) => Ok(w),
                            None => Err(Error::InvalidState(
                                "unable to retrieve mod data!".to_string(),
                            )),
                        }?;
                        let win_end = match end[l + win_size - 1] {
                            Some(w) => Ok(w),
                            None => Err(Error::InvalidState(
                                "unable to retrieve mod data!".to_string(),
                            )),
                        }?;
                        let ref_win_start = match ref_start[l..l + win_size].iter().flatten().min()
                        {
                            Some(w) => *w,
                            None => -1,
                        };
                        let ref_win_end = match ref_end[l..l + win_size].iter().flatten().max() {
                            Some(w) => *w,
                            None => -1,
                        };
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
    }

    Ok(true)
}
