//! # Displays BAM file metadata including contigs and modification types
//!
//! This function displays contigs, contig lengths and modification types after looking at the header and the given records

use crate::{AllowedAGCTN, CurrRead, Error, ModChar};
use rust_htslib::bam;
use std::collections::HashSet;
use std::io;
use std::rc::Rc;

/// Run the peek command to display BAM file metadata
///
/// # Examples
///
/// ```
/// use nanalogue_core::peek;
/// use rust_htslib::bam::Reader;
/// use rust_htslib::bam::Read as _;
///
/// let mut bam = Reader::from_path("examples/example_1.bam").unwrap();
/// let header = bam.header().clone();
/// let mut output = Vec::new();
///
/// peek::run(
///     &mut output,
///     &header,
///     bam.rc_records().take(100),
/// )
/// .unwrap();
///
/// let expected = "contigs_and_lengths:\ndummyI\t22\ndummyII\t48\ndummyIII\t76\n\nmodifications:\nG-7200\nT+T\n";
/// assert_eq!(String::from_utf8(output).unwrap(), expected);
/// ```
///
/// # Errors
///
/// Returns an error if:
/// - Writing to the output handle fails
/// - BAM header parsing fails
/// - Reading or parsing BAM records fails
/// - Converting modification data fails
pub fn run<W, D>(handle: &mut W, header: &bam::HeaderView, records: D) -> Result<(), Error>
where
    W: io::Write,
    D: Iterator<Item = Result<Rc<bam::Record>, rust_htslib::errors::Error>>,
{
    // 1. Display contigs
    writeln!(handle, "contigs_and_lengths:")?;
    for tid in 0..header.target_count() {
        let name = std::str::from_utf8(
            header
                .target_names()
                .get(tid as usize)
                .ok_or_else(|| Error::InvalidSeqLength(format!("tid {tid} out of bounds")))?,
        )?;
        let length = header.target_len(tid).ok_or_else(|| {
            Error::InvalidSeqLength(format!("target_len returned None for tid {tid}"))
        })?;
        writeln!(handle, "{name}\t{length}")?;
    }
    writeln!(handle)?;

    // 2. Collect modifications from records
    let mut modifications = HashSet::new();

    for record_result in records {
        let record = record_result?;

        // Convert to CurrRead to extract modification data
        let curr_read = CurrRead::try_from(record)?;
        let base_mods = &curr_read.mod_data().0.base_mods;

        for base_mod in base_mods {
            // Convert modified_base from u8 to AllowedAGCTN
            let base = AllowedAGCTN::try_from(base_mod.modified_base)?;

            // Convert modification_type char to ModChar
            let mod_char = ModChar::from(base_mod.modification_type);

            // Create the modification string: base+strand+modification_type
            let mod_string = format!("{}{}{}", base, base_mod.strand, mod_char);
            let _: bool = modifications.insert(mod_string);
        }
    }

    // 3. Display modifications
    writeln!(handle, "modifications:")?;
    if modifications.is_empty() {
        writeln!(handle, "None")?;
    } else {
        // Sort modifications for consistent output
        let mut sorted_mods: Vec<_> = modifications.into_iter().collect();
        sorted_mods.sort();
        for mod_string in sorted_mods {
            writeln!(handle, "{mod_string}")?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BamRcRecords, InputBamBuilder, InputMods, OptionalTag, PathOrURLOrStdin};
    use std::fs;

    #[test]
    fn peek_shows_example_1() {
        // Read the BAM file
        let mut input_bam = InputBamBuilder::default()
            .bam_path(PathOrURLOrStdin::Path("examples/example_1.bam".into()))
            .build()
            .expect("should build InputBam");

        let mut reader =
            bam::Reader::from_path("examples/example_1.bam").expect("should open BAM file");

        let bam_rc_records = BamRcRecords::new(
            &mut reader,
            &mut input_bam,
            &mut InputMods::<OptionalTag>::default(),
        )
        .expect("should create BamRcRecords");

        // Run peek
        let mut output = Vec::new();
        run(
            &mut output,
            &bam_rc_records.header,
            bam_rc_records.rc_records.take(100),
        )
        .expect("peek should succeed");

        // Read expected output
        let expected = fs::read_to_string("examples/example_1_peek")
            .expect("should read expected output file");

        // Compare
        assert_eq!(
            String::from_utf8(output).expect("output should be valid UTF-8"),
            expected
        );
    }

    #[test]
    fn peek_shows_example_3() {
        // Read the BAM file
        let mut input_bam = InputBamBuilder::default()
            .bam_path(PathOrURLOrStdin::Path("examples/example_3.bam".into()))
            .build()
            .expect("should build InputBam");

        let mut reader =
            bam::Reader::from_path("examples/example_3.bam").expect("should open BAM file");

        let bam_rc_records = BamRcRecords::new(
            &mut reader,
            &mut input_bam,
            &mut InputMods::<OptionalTag>::default(),
        )
        .expect("should create BamRcRecords");

        // Run peek
        let mut output = Vec::new();
        run(
            &mut output,
            &bam_rc_records.header,
            bam_rc_records.rc_records.take(100),
        )
        .expect("peek should succeed");

        // Read expected output
        let expected = fs::read_to_string("examples/example_3_peek")
            .expect("should read expected output file");

        // Compare
        assert_eq!(
            String::from_utf8(output).expect("output should be valid UTF-8"),
            expected
        );
    }

    #[test]
    fn peek_shows_simulated_example() {
        use crate::simulate_mod_bam::{SimulationConfig, TempBamSimulation};

        let config_json = r#"{
            "contigs": {
                "number": 1,
                "len_range": [100, 100],
                "repeated_seq": "ACGTACGTACGTACGT"
            },
            "reads": [{
                "number": 20,
                "mapq_range": [20, 30],
                "base_qual_range": [20, 30],
                "len_range": [0.8, 1.0],
                "mods": [
                    {
                        "base": "C",
                        "is_strand_plus": true,
                        "mod_code": "m",
                        "win": [2],
                        "mod_range": [[0.7, 0.9]]
                    },
                    {
                        "base": "A",
                        "is_strand_plus": true,
                        "mod_code": "a",
                        "win": [3],
                        "mod_range": [[0.3, 0.5]]
                    },
                    {
                        "base": "T",
                        "is_strand_plus": false,
                        "mod_code": "t",
                        "win": [1],
                        "mod_range": [[0.5, 0.6]]
                    }
                ]
            }]
        }"#;

        let config: SimulationConfig = serde_json::from_str(config_json).unwrap();
        let sim = TempBamSimulation::new(config).unwrap();
        let mut reader = bam::Reader::from_path(sim.bam_path()).unwrap();

        let mut input_bam = InputBamBuilder::default()
            .bam_path(PathOrURLOrStdin::Path(sim.bam_path().into()))
            .build()
            .expect("should build InputBam");

        let bam_rc_records = BamRcRecords::new(
            &mut reader,
            &mut input_bam,
            &mut InputMods::<OptionalTag>::default(),
        )
        .expect("should create BamRcRecords");

        // Run peek
        let mut output = Vec::new();
        run(
            &mut output,
            &bam_rc_records.header,
            bam_rc_records.rc_records.take(100),
        )
        .expect("peek should succeed");

        let output_str = String::from_utf8(output).expect("output should be valid UTF-8");

        let expected = "contigs_and_lengths:\ncontig_00000\t100\n\nmodifications:\nA+a\nC+m\nT-t\n";

        assert_eq!(output_str, expected);
    }

    #[test]
    fn peek_empty_file() {
        use crate::simulate_mod_bam::{SimulationConfig, TempBamSimulation};

        let config_json = r#"{
            "contigs": {
                "number": 1,
                "len_range": [100, 100],
                "repeated_seq": "ACGTACGTACGTACGT"
            },
            "reads": []
        }"#;

        let config: SimulationConfig = serde_json::from_str(config_json).unwrap();
        let sim = TempBamSimulation::new(config).unwrap();
        let mut reader = bam::Reader::from_path(sim.bam_path()).unwrap();

        let mut input_bam = InputBamBuilder::default()
            .bam_path(PathOrURLOrStdin::Path(sim.bam_path().into()))
            .build()
            .expect("should build InputBam");

        let bam_rc_records = BamRcRecords::new(
            &mut reader,
            &mut input_bam,
            &mut InputMods::<OptionalTag>::default(),
        )
        .expect("should create BamRcRecords");

        // Run peek
        let mut output = Vec::new();
        run(
            &mut output,
            &bam_rc_records.header,
            bam_rc_records.rc_records.take(100),
        )
        .expect("peek should succeed");

        let output_str = String::from_utf8(output).expect("output should be valid UTF-8");

        let expected = "contigs_and_lengths:\ncontig_00000\t100\n\nmodifications:\nNone\n";

        assert_eq!(output_str, expected);
    }
}
