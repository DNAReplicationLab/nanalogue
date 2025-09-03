//! # Gets information on one read id
//!
//! This module retrieves information about one read id
//! from a BAM file and writes it as a JSON to the standard
//! output. If more than one entry has the same read id,
//! then both are output.
use crate::{CurrRead, Error};
use rust_htslib::bam;
use std::rc::Rc;

/// Gets information on one read id and prints it to standard output
/// in a JSON format.
pub fn run<W, D>(handle: &mut W, bam_records: D, read_id: &str) -> Result<bool, Error>
where
    W: std::io::Write,
    D: IntoIterator<Item = Result<Rc<bam::Record>, rust_htslib::errors::Error>>,
{
    // convert read id into bytes
    let read_id_bytes = read_id.as_bytes();

    // Go record by record in the BAM file,
    // and collect entries that match our read id
    for k in bam_records
        .into_iter()
        .filter(|r| match r {
            Ok(v) => v.qname() == read_id_bytes,
            Err(_) => true,
        })
        .collect::<Result<Vec<Rc<bam::Record>>, _>>()?
    {
        let record = CurrRead::try_from(k)?;
        writeln!(handle, "{}", record)?;
    }
    Ok(true)
}
