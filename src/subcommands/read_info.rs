//! # Gets information on one read id
//!
//! This module retrieves information about reads
//! from a BAM file and writes it as a JSON to the standard output.
use crate::{CurrRead, Error};
use rust_htslib::bam;
use std::iter;
use std::rc::Rc;

/// Gets information on reads and prints it to standard output
/// in a JSON format.
pub fn run<W, D>(handle: &mut W, bam_records: D) -> Result<bool, Error>
where
    W: std::io::Write,
    D: IntoIterator<Item = Result<Rc<bam::Record>, rust_htslib::errors::Error>>,
{
    let mut is_first_record_written = vec![false].into_iter().chain(iter::repeat(true));

    write!(handle, "[")?;

    // Go record by record in the BAM file, and print entries
    for k in bam_records {
        let record = k?;

        match is_first_record_written.next().expect("no error") {
            false => writeln!(handle)?,
            true => writeln!(handle, ",")?,
        }

        let curr_read = CurrRead::try_from(record)?;
        write!(handle, "{}", curr_read)?;
    }

    match is_first_record_written.next().expect("no error") {
        false => writeln!(handle, "]")?,
        true => writeln!(handle, "\n]")?,
    };
    Ok(true)
}
