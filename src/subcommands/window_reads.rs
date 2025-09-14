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
) -> Result<bool, Error>
where
    W: std::io::Write,
    F: Fn(&[u8]) -> Result<F32AbsValBelow1, Error>,
    D: IntoIterator<Item = Result<Rc<Record>, rust_htslib::errors::Error>>,
{
    // Get windowing parameters
    let win_size: usize = window_options.win.get().try_into()?;
    let slide_size: usize = window_options.step.get().try_into()?;

    // Go record by record in the BAM file,
    for r in bam_records {
        // read records
        let record = r?;

        // set data in records
        let curr_read_state = CurrRead::default()
            .try_from_only_alignment(&record)?
            .set_mod_data_restricted_options(&record, &mods)?;
        let qname = curr_read_state.read_id()?.to_string();
        let strand = curr_read_state.strand();
        let contig = match curr_read_state.read_state() {
            ReadState::Unmapped => ".".to_string(),
            _ => curr_read_state.contig_name()?.to_string(),
        };

        // read and window modification data, then print the output
        let (BaseMods { base_mods: v }, _) = curr_read_state.mod_data();
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
                    let win_start = start[l].ok_or(Error::InvalidState(
                        "unable to retrieve mod data!".to_string(),
                    ))?;
                    let win_end = end[l + win_size - 1].ok_or(Error::InvalidState(
                        "unable to retrieve mod data!".to_string(),
                    ))?;
                    let ref_win_start = *(ref_start[l..l + win_size]
                        .iter()
                        .flatten()
                        .min()
                        .unwrap_or(&-1));
                    let ref_win_end = *(ref_end[l..l + win_size]
                        .iter()
                        .flatten()
                        .max()
                        .unwrap_or(&-1));
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

    Ok(true)
}
