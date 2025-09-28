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

    writeln!(handle, "\n]")?;
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
        let records: Vec<Result<Rc<bam::Record>, rust_htslib::errors::Error>> = reader
            .records()
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

    #[test]
    fn test_run_with_unmapped_filter() -> Result<(), Error> {
        // Collect records and filter to only include unmapped reads (like --read-filter unmapped)
        let mut reader = nanalogue_bam_reader("./examples/example_1.bam")?;
        let records: Vec<Result<Rc<bam::Record>, rust_htslib::errors::Error>> = reader
            .records()
            .map(|r| r.map(Rc::new))
            .filter(|r| {
                // Filter to only unmapped reads (flag 4)
                match r {
                    Ok(record) => record.flags() == 4, // Unmapped flag
                    Err(_) => true,                    // Keep errors to let run() handle them
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
                "read_id": "a4f36092-b4d5-47a9-813e-c22c3b477a0c",
                "sequence_length": 48,
                "alignment_type": "unmapped",
                "mod_count": "G-7200:0;T+T:3;(probabilities >= 0.5, PHRED base qual >= 0)"
            }
        ]);
        assert_eq!(parsed, expected);

        Ok(())
    }

    #[test]
    fn test_run_with_region_filter_dummy_i() -> Result<(), Error> {
        // Collect records and filter to only include those in the dummyI region (like --region dummyI)
        let mut reader = nanalogue_bam_reader("./examples/example_1.bam")?;
        let records: Vec<Result<Rc<bam::Record>, rust_htslib::errors::Error>> = reader
            .records()
            .map(|r| r.map(Rc::new))
            .filter(|r| {
                // Filter to only records in dummyI contig (TID 0)
                match r {
                    Ok(record) => !record.is_unmapped() && record.tid() == 0, // dummyI is TID 0
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
                "read_id": "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575",
                "sequence_length": 8,
                "contig": "dummyI",
                "reference_start": 9,
                "reference_end": 17,
                "alignment_length": 8,
                "alignment_type": "primary_forward",
                "mod_count": "T+T:0;(probabilities >= 0.5, PHRED base qual >= 0)"
            }
        ]);
        assert_eq!(parsed, expected);

        Ok(())
    }

    #[test]
    fn test_run_with_example_6() -> Result<(), Error> {
        // Collect records from example_6.sam
        let mut reader = nanalogue_bam_reader("./examples/example_6.sam")?;
        let records: Vec<Result<Rc<bam::Record>, rust_htslib::errors::Error>> =
            reader.records().map(|r| r.map(Rc::new)).collect();

        // Gets an output from the function and compares with expected
        let mut output_buffer = Vec::new();
        assert!(run(&mut output_buffer, records.into_iter())?);
        let output_json = String::from_utf8(output_buffer)?;
        let parsed: Value = serde_json::from_str(&output_json)?;
        let expected = serde_json::json!([
            {
                "read_id": "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575",
                "sequence_length": 8,
                "contig": "dummyI",
                "reference_start": 9,
                "reference_end": 17,
                "alignment_length": 8,
                "alignment_type": "primary_forward",
                "mod_count": "NA"
            },
            {
                "read_id": "fffffff1-10d2-49cb-8ca3-e8d48979001b",
                "sequence_length": 33,
                "contig": "dummyII",
                "reference_start": 3,
                "reference_end": 36,
                "alignment_length": 33,
                "alignment_type": "primary_reverse",
                "mod_count": "T+T:1;(probabilities >= 0.5, PHRED base qual >= 0)"
            }
        ]);
        assert_eq!(parsed, expected);

        Ok(())
    }
}
