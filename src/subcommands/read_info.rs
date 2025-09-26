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

// Tests follow

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nanalogue_bam_reader;
    use rust_htslib::bam::{self, Read};
    use serde_json::Value;
    use std::rc::Rc;

    #[test]
    fn test_run_with_example_2_zero_len() -> Result<(), Error> {

        // Collect records and filter out zero-length sequences (like the main program does)
        let mut reader = nanalogue_bam_reader("./examples/example_2_zero_len.sam")?;
        let records: Vec<Result<Rc<bam::Record>, rust_htslib::errors::Error>> = reader.records()
            .map(|r| r.map(Rc::new))
            .filter(|r| {
                // Filter out records that would cause InvalidSeqLength (zero-length sequences)
                match r {
                    Ok(record) => record.seq_len() > 0,
                    Err(_) => true, // Keep errors to let run() handle them
                }
            })
            .collect();

        // Gets an output from the function and compares with expected
        let mut output_buffer = Vec::new();
        assert!(run(&mut output_buffer, records.into_iter())?);
        let output_json = String::from_utf8(output_buffer)?;
        let parsed: Value = serde_json::from_str(&output_json)?;
        let expected = serde_json::json!([
            {
                "read_id": "read2",
                "sequence_length": 48,
                "contig": "dummyIII",
                "reference_start": 23,
                "reference_end": 71,
                "alignment_length": 48,
                "alignment_type": "primary_forward",
                "mod_count": "NA"
            }
        ]);
        assert_eq!(parsed, expected);

        Ok(())
    }
}
