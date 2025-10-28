//! Tests for `read_utils.rs` extracted from doctests

use bedrs::prelude::StrandedBed3;
use bedrs::{Bed3, Coordinates, Strand};
use nanalogue_core::simulate_mod_bam::TempBamSimulation;
use nanalogue_core::{
    CurrRead, Error, Intersects, ModChar, ReadState, ThresholdState, nanalogue_bam_reader,
    read_utils::OnlyAlignData,
};
use rust_htslib::bam::Read;
use std::collections::{HashMap, hash_map::Entry};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn set_read_state_example_1() -> Result<(), Error> {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        for (count, record) in reader.records().enumerate() {
            let r = record?;
            let curr_read = CurrRead::default().set_read_state(&r)?;
            match count {
                0 | 1 => assert_eq!(curr_read.read_state(), ReadState::PrimaryFwd),
                2 => assert_eq!(curr_read.read_state(), ReadState::PrimaryRev),
                3 => assert_eq!(curr_read.read_state(), ReadState::Unmapped),
                _ => unreachable!(),
            }
        }
        Ok(())
    }

    #[test]
    fn set_read_state_example_3() -> Result<(), Error> {
        let mut reader = nanalogue_bam_reader("examples/example_3.bam")?;
        let mut count = 1; // NOTE that we start the counter from 1 here
        // as reads are called 001, 002, ..., 010 here.
        // so it is easier for us to read code when counter starts from 1.
        for record in reader.records() {
            let r = record?;
            let curr_read = CurrRead::default().set_read_state(&r)?;
            match count {
                1 | 4 | 5 | 8 | 10 => assert_eq!(curr_read.read_state(), ReadState::PrimaryFwd),
                2 => assert_eq!(curr_read.read_state(), ReadState::SecondaryFwd),
                3 => assert_eq!(curr_read.read_state(), ReadState::PrimaryRev),
                6 => assert_eq!(curr_read.read_state(), ReadState::SecondaryRev),
                7 => assert_eq!(curr_read.read_state(), ReadState::SupplementaryFwd),
                9 => assert_eq!(curr_read.read_state(), ReadState::SupplementaryRev),
                _ => unreachable!(),
            }
            count += 1;
        }
        Ok(())
    }

    #[test]
    fn set_seq_len() -> Result<(), Error> {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        for (count, record) in reader.records().enumerate() {
            let r = record?;
            let curr_read = CurrRead::default().set_read_state(&r)?.set_seq_len(&r)?;
            let Ok(len) = curr_read.seq_len() else {
                unreachable!()
            };
            match count {
                0 => assert_eq!(len, 8),
                1 | 3 => assert_eq!(len, 48),
                2 => assert_eq!(len, 33),
                _ => unreachable!(),
            }
        }
        Ok(())
    }

    #[test]
    fn set_seq_len_random() -> Result<(), Error> {
        // creates 2 contigs of 1000 bp each and reads of random
        // mapping, position etc. with lengths b/w 10-20% of contig size.
        let config_json = r#"{
        "contigs": {
            "number": 2,
            "len_range": [1000, 1000]
        },
        "reads": [{
            "number": 10000,
            "len_range": [0.1, 0.2]
        }]
    }"#;
        let simulated_bam = TempBamSimulation::new(config_json).unwrap();
        let mut reader = nanalogue_bam_reader(simulated_bam.bam_path())?;

        let (sum, deviation_sq) = {
            let mut sum: u64 = 0;
            let mut deviation_sq: u64 = 0;
            for record in reader.records() {
                let r = record?;
                let curr_read = CurrRead::default().set_read_state(&r)?.set_seq_len(&r)?;
                let len = curr_read.seq_len().unwrap();
                sum = sum.checked_add(len).unwrap();
                deviation_sq = deviation_sq.checked_add(len.abs_diff(150).pow(2)).unwrap();
            }
            (sum, deviation_sq)
        };

        // mean should be 150, and variance (200-100)^2/12.
        // we check these to 15% tolerance
        assert!((sum / 10000).abs_diff(150) < 22);
        assert!((deviation_sq / 10000).abs_diff(10000 / 12) < 1500 / 12);
        Ok(())
    }

    #[test]
    #[should_panic(expected = "InvalidDuplicates")]
    fn set_seq_len_duplicate_should_panic() {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam").unwrap();
        if let Some(record) = reader.records().next() {
            let r = record.unwrap();
            let _curr_read = CurrRead::default()
                .set_read_state(&r)
                .unwrap()
                .set_seq_len(&r)
                .unwrap()
                .set_seq_len(&r)
                .unwrap();
        }
    }

    #[test]
    fn set_align_len() -> Result<(), Error> {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let mut count = 0;
        for record in reader.records() {
            let r = record?;
            let curr_read = CurrRead::default().set_read_state(&r)?.set_align_len(&r)?;
            let Ok(len) = curr_read.align_len() else {
                unreachable!()
            };
            match count {
                0 => assert_eq!(len, 8),
                1 => assert_eq!(len, 48),
                2 => assert_eq!(len, 33),
                _ => unreachable!(),
            }
            count += 1;
            if count == 3 {
                break;
            }
        }
        Ok(())
    }

    #[test]
    fn set_align_len_random() -> Result<(), Error> {
        // creates 2 contigs of 1000 bp each and reads of random
        // mapping, position etc. with lengths b/w 10-20% of contig size.
        // then, with a barcode, the alignment length is the same but the
        // seq length will increase by 2 times the barcode length.
        let config_json = r#"{
        "contigs": {
            "number": 2,
            "len_range": [1000, 1000]
        },
        "reads": [{
            "number": 10000,
            "len_range": [0.1, 0.2],
            "barcode": "AAGTAA"
        }]
    }"#;
        let sim = TempBamSimulation::new(config_json).unwrap();
        let mut reader = nanalogue_bam_reader(sim.bam_path())?;

        let (sum_seq_len, deviation_sequence_len_sq, sum_align_len, deviation_align_len_sq, count) = {
            let mut sum_seq_len: u64 = 0;
            let mut deviation_sequence_len_sq: u64 = 0;
            let mut sum_align_len: u64 = 0;
            let mut deviation_align_len_sq: u64 = 0;
            let mut count: u64 = 0;

            for record in reader.records() {
                let r = record?;
                let curr_read = {
                    let curr_read_result = CurrRead::default()
                        .set_read_state(&r)?
                        .set_seq_len(&r)?
                        .set_align_len(&r);
                    match curr_read_result {
                        Err(Error::Unmapped) => continue,
                        v => v,
                    }
                }?;

                let len: u64 = curr_read.seq_len().unwrap();
                sum_seq_len = sum_seq_len.checked_add(len).unwrap();
                deviation_sequence_len_sq = deviation_sequence_len_sq
                    .checked_add(len.abs_diff(162).pow(2))
                    .unwrap();

                let len: u64 = curr_read.align_len().unwrap();
                sum_align_len = sum_align_len.checked_add(len).unwrap();
                deviation_align_len_sq = deviation_align_len_sq
                    .checked_add(len.abs_diff(150).pow(2))
                    .unwrap();

                count = count.checked_add(1).unwrap();
            }
            (
                sum_seq_len,
                deviation_sequence_len_sq,
                sum_align_len,
                deviation_align_len_sq,
                count,
            )
        };

        // alignment length:
        // mean should be 150, and variance (200-100)^2/12.
        // we check these to 15% tolerance
        assert!(sum_align_len.checked_div(count).unwrap().abs_diff(150) < 22);
        assert!(
            deviation_align_len_sq
                .checked_div(count)
                .unwrap()
                .abs_diff(10000 / 12)
                < 1500 / 12
        );

        // sequence length:
        // mean should be 162 i.e. 150 + 2*6 (AAGTAA), and variance the same (200-100)^2/12.
        // we check these to 15% tolerance
        assert!(sum_seq_len.checked_div(count).unwrap().abs_diff(162) < 24);
        assert!(
            deviation_sequence_len_sq
                .checked_div(count)
                .unwrap()
                .abs_diff(10000 / 12)
                < 1500 / 12
        );
        Ok(())
    }

    #[test]
    #[should_panic(expected = "Unmapped")]
    fn set_align_len_unmapped_should_panic() {
        // Fourth read in the following file is unmapped, so we check if
        // we hit an error on that read
        let mut reader = nanalogue_bam_reader("examples/example_1.bam").unwrap();
        let mut count = 0;
        for record in reader.records() {
            let r = record.unwrap();
            if count < 3 {
                count += 1;
            } else {
                let _: CurrRead<OnlyAlignData> = CurrRead::default()
                    .set_read_state(&r)
                    .unwrap()
                    .set_align_len(&r)
                    .unwrap();
            }
        }
    }

    #[test]
    #[should_panic(expected = "InvalidDuplicates")]
    fn set_align_len_duplicate_should_panic() {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam").unwrap();
        if let Some(record) = reader.records().next() {
            let r = record.unwrap();
            let _curr_read = CurrRead::default()
                .set_read_state(&r)
                .unwrap()
                .set_align_len(&r)
                .unwrap()
                .set_align_len(&r)
                .unwrap();
        }
    }

    #[test]
    fn set_contig_id_and_start() -> Result<(), Error> {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let mut count = 0;
        for record in reader.records() {
            let r = record?;
            let curr_read = CurrRead::default()
                .set_read_state(&r)?
                .set_contig_id_and_start(&r)?;
            match (count, curr_read.contig_id_and_start()) {
                (0, Ok((0, 9))) | (1, Ok((2, 23))) | (2, Ok((1, 3))) => {}
                _ => unreachable!(),
            }
            count += 1;
            if count == 3 {
                break;
            } // the fourth entry is unmapped, and will lead to an error.
        }
        Ok(())
    }

    #[test]
    #[should_panic(expected = "Unmapped")]
    fn set_contig_id_and_start_unmapped_should_panic() {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam").unwrap();
        let mut count = 0;
        for record in reader.records() {
            let r = record.unwrap();
            if count < 3 {
                count += 1;
                continue;
            }
            // the fourth read is unmapped
            let _curr_read = CurrRead::default()
                .set_read_state(&r)
                .unwrap()
                .set_contig_id_and_start(&r)
                .unwrap();
        }
    }

    #[test]
    #[should_panic(expected = "InvalidDuplicates")]
    fn set_contig_id_and_start_duplicate_should_panic() {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam").unwrap();
        if let Some(record) = reader.records().next() {
            let r = record.unwrap();
            let _curr_read = CurrRead::default()
                .set_read_state(&r)
                .unwrap()
                .set_contig_id_and_start(&r)
                .unwrap()
                .set_contig_id_and_start(&r)
                .unwrap();
        }
    }

    #[test]
    fn set_contig_name() -> Result<(), Error> {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let mut count = 0;
        for record in reader.records() {
            let r = record?;
            let curr_read = CurrRead::default()
                .set_read_state(&r)?
                .set_contig_name(&r)?;
            let Ok(contig_name) = curr_read.contig_name() else {
                unreachable!()
            };
            match (count, contig_name) {
                (0, "dummyI") | (1, "dummyIII") | (2, "dummyII") => {}
                _ => unreachable!(),
            }
            count += 1;
            if count == 3 {
                break;
            } // the fourth entry is unmapped, and will lead to an error.
        }
        Ok(())
    }

    #[test]
    #[should_panic(expected = "Unmapped")]
    fn set_contig_name_unmapped_should_panic() {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam").unwrap();
        let mut count = 0;
        for record in reader.records() {
            if count < 3 {
                count += 1;
                continue;
            }
            let r = record.unwrap();
            let curr_read = CurrRead::default().set_read_state(&r).unwrap();
            let _: CurrRead<OnlyAlignData> = curr_read.set_contig_name(&r).unwrap();
        }
    }

    #[test]
    #[should_panic(expected = "Unmapped")]
    fn get_contig_name_unmapped_should_panic() {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam").unwrap();
        let mut count = 0;
        for record in reader.records() {
            if count < 3 {
                count += 1;
                continue;
            }
            let r = record.unwrap();
            let curr_read = CurrRead::default().set_read_state(&r).unwrap();
            let _: &str = curr_read.contig_name().unwrap();
        }
    }

    #[test]
    #[should_panic(expected = "UnavailableData")]
    fn get_contig_name_without_setting_should_panic() {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam").unwrap();
        if let Some(record) = reader.records().next() {
            let r = record.unwrap();
            let curr_read = CurrRead::default().set_read_state(&r).unwrap();
            let _: &str = curr_read.contig_name().unwrap();
        }
    }

    #[test]
    #[should_panic(expected = "InvalidDuplicates")]
    fn set_contig_name_duplicate_should_panic() {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam").unwrap();
        if let Some(record) = reader.records().next() {
            let r = record.unwrap();
            let _curr_read = CurrRead::default()
                .set_read_state(&r)
                .unwrap()
                .set_contig_name(&r)
                .unwrap()
                .set_contig_name(&r)
                .unwrap();
        }
    }

    #[test]
    fn set_read_id() -> Result<(), Error> {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        for (count, record) in reader.records().enumerate() {
            let r = record?;
            let curr_read = CurrRead::default().set_read_state(&r)?.set_read_id(&r)?;
            let Ok(read_id) = curr_read.read_id() else {
                unreachable!()
            };
            match (count, read_id) {
                (0, "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575")
                | (1 | 3, "a4f36092-b4d5-47a9-813e-c22c3b477a0c")
                | (2, "fffffff1-10d2-49cb-8ca3-e8d48979001b") => {}
                _ => unreachable!(),
            }
        }
        Ok(())
    }

    #[test]
    #[should_panic(expected = "UnavailableData")]
    fn get_read_id_without_setting_should_panic() {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam").unwrap();
        if let Some(record) = reader.records().next() {
            let r = record.unwrap();
            let curr_read = CurrRead::default().set_read_state(&r).unwrap();
            let _: &str = curr_read.read_id().unwrap();
        }
    }

    #[test]
    #[should_panic(expected = "InvalidDuplicates")]
    fn set_read_id_duplicate_should_panic() {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam").unwrap();
        if let Some(record) = reader.records().next() {
            let r = record.unwrap();
            let _curr_read = CurrRead::default()
                .set_read_state(&r)
                .unwrap()
                .set_read_id(&r)
                .unwrap()
                .set_read_id(&r)
                .unwrap();
        }
    }

    #[test]
    fn strand() -> Result<(), Error> {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        for (count, record) in reader.records().enumerate() {
            let r = record?;
            let curr_read = CurrRead::default().set_read_state(&r)?;
            let strand = curr_read.strand();
            match (count, strand) {
                (0 | 1, '+') | (2, '-') | (3, '.') => {}
                _ => unreachable!(),
            }
        }
        Ok(())
    }

    #[test]
    fn seq_on_ref_coords() -> Result<(), Error> {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        for record in reader.records() {
            let r = record?;
            let curr_read = CurrRead::default().try_from_only_alignment(&r)?;

            // Skip unmapped reads
            if curr_read.read_state().to_string() != "unmapped" {
                let (contig_id, start) = curr_read.contig_id_and_start()?;
                let align_len = curr_read.align_len()?;

                // Create a region that overlaps with the read but is short of one bp.
                // Note that this BAM file has reads with all bases matching perfectly
                // with the reference.
                let region = Bed3::new(contig_id, start, start + align_len - 1);
                let seq_subset = curr_read.seq_on_ref_coords(&r, &region)?;

                // Check for sequence length match
                assert_eq!(curr_read.seq_len()? - 1, u64::try_from(seq_subset.len())?);

                // Create a region with no overlap at all and check we get no data
                let region = Bed3::new(contig_id, start + align_len, start + align_len + 2);
                match curr_read.seq_on_ref_coords(&r, &region) {
                    Err(Error::UnavailableData) => Ok(()),
                    _ => Err(Error::UnknownError),
                }?;
            }
        }
        Ok(())
    }

    #[test]
    fn seq_on_ref_coords_2() -> Result<(), Error> {
        // make a random BAM file but contigs are all just a 10 bp
        // sequence repeated to fill the required length.
        // so, we can easily check sequence retrieval.
        let config_json = r#"{
        "contigs": {
            "number": 2,
            "len_range": [1000, 1000],
            "repeated_seq": "AAGCTAGCTG"
        },
        "reads": [{
            "number": 10000,
            "len_range": [0.1, 0.2]
        }]
    }"#;
        let sim = TempBamSimulation::new(config_json).unwrap();
        let mut reader = nanalogue_bam_reader(sim.bam_path())?;
        let mut cnt = 0;

        for record in reader.records() {
            let r = record?;
            let curr_read = CurrRead::default().try_from_only_alignment(&r)?;

            // Skip unmapped reads, for others, check sequence match.
            // We probe first contig 225-236, so if we have AAGCTAGCTG repeated
            // over and over, we expect the following result.
            let expected_seq = "AGCTGAAGCTA";
            if curr_read.read_state().to_string() != "unmapped" {
                let region = Bed3::new(0, 225, 236);
                match curr_read.seq_on_ref_coords(&r, &region) {
                    Err(Error::UnavailableData) => {}
                    Ok(v) => {
                        assert!(
                            expected_seq.contains(str::from_utf8(&v).expect("no error")),
                            "unknown sequence"
                        );
                        cnt += 1;
                    }
                    _ => panic!("erroneous outcome"),
                }
            }
        }
        // 10000 reads with 2 contigs but 1/7 are unmapped => ~4300 reads per contig (10000/2 * 6/7).
        // if their length varies from 100-200 bp, we expect ~1/10-1/5 of reads
        // to pass through any given base i.e. ~430-860. So, 200 is actually quite lax, we've left it
        // this lax to account for (extreme) statistical outliers.
        assert!(cnt > 200);
        Ok(())
    }

    #[test]
    fn seq_on_ref_coords_2_but_barcode() {
        // make a random BAM file but contigs are all just a 10 bp
        // sequence repeated to fill the required length, but with a
        // barcode. Occasionally we will get a barcode in the region,
        // so we check that the barcode is removed.
        let config_json = r#"{
        "contigs": {
            "number": 1,
            "len_range": [300, 300],
            "repeated_seq": "AAGCTAGCTG"
        },
        "reads": [{
            "number": 10000,
            "len_range": [0.03, 0.03],
            "barcode": "CAG"
        }]
    }"#;
        let sim = TempBamSimulation::new(config_json).unwrap();
        let mut reader = nanalogue_bam_reader(sim.bam_path()).unwrap();

        // We probe first contig 225-229, so if we have AAGCTAGCTG repeated
        // over and over, we expect the sequence here to be "AGCT".
        // A read can start or stop here which gives
        // the following possibilities:
        let expected_seqs = ["AGCT", "T", "CT", "GCT", "AGC", "AG", "A"];

        // create a table of individual counts
        let mut cnt_individual: HashMap<String, u32> = HashMap::new();
        for k in expected_seqs {
            let _: Option<u32> = cnt_individual.insert(k.to_string(), 0);
        }

        for record in reader.records() {
            let r = record.unwrap();
            let curr_read = CurrRead::default().try_from_only_alignment(&r).unwrap();

            // Skip unmapped reads, for others, check sequence match.
            if curr_read.read_state().to_string() != "unmapped" {
                let region = Bed3::new(0, 225, 229);
                match curr_read.seq_on_ref_coords(&r, &region) {
                    Err(Error::UnavailableData) => {}
                    Ok(v) => {
                        if expected_seqs
                            .iter()
                            .any(|k| *k == str::from_utf8(&v).expect("no error"))
                        {
                            let _: Entry<String, u32> = cnt_individual
                                .entry(String::from_utf8(v).expect("string conversion error"))
                                .and_modify(|m| *m += 1);
                        } else {
                            panic!("unknown sequence {v:?}")
                        }
                    }
                    _ => panic!("erroneous outcome"),
                }
            }
        }

        // check that every entry in the list of possible sequences is visited at least once.
        // (we are not going to bother calculating the expected statistics of this count!)
        for k in cnt_individual.into_values() {
            assert!(k >= 1);
        }
    }

    #[test]
    fn basecount_per_mod() -> Result<(), Error> {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let first_count = Some(HashMap::from([(ModChar::new('T'), 0)]));
        let second_count = Some(HashMap::from([(ModChar::new('T'), 3)]));
        let third_count = Some(HashMap::from([(ModChar::new('T'), 1)]));
        let fourth_count = Some(HashMap::from([
            (ModChar::new('T'), 3),
            (ModChar::new('ᰠ'), 0),
        ]));
        for (count, record) in reader.records().enumerate() {
            let r = record?;
            let curr_read = CurrRead::default().set_read_state(&r)?.set_mod_data(
                &r,
                ThresholdState::GtEq(180),
                0,
            )?;
            let modcount = curr_read.base_count_per_mod();
            match (count, modcount) {
                (0, v) => assert_eq!(v, first_count),
                (1, v) => assert_eq!(v, second_count),
                (2, v) => assert_eq!(v, third_count),
                (3, v) => assert_eq!(v, fourth_count),
                _ => unreachable!(),
            }
        }
        Ok(())
    }

    #[test]
    fn try_from_stranded_bed3() -> Result<(), Error> {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        for (count, record) in reader.records().enumerate() {
            let r = record?;
            let curr_read = CurrRead::default().try_from_only_alignment(&r)?;
            let bed3_stranded = StrandedBed3::try_from(&curr_read).unwrap();
            let exp_bed3_stranded = match count {
                0 => StrandedBed3::new(0, 9, 17, Strand::Forward),
                1 => StrandedBed3::new(2, 23, 71, Strand::Forward),
                2 => StrandedBed3::new(1, 3, 36, Strand::Reverse),
                3 => StrandedBed3::empty(),
                _ => unreachable!(),
            };
            assert_eq!(*bed3_stranded.chr(), *exp_bed3_stranded.chr());
            assert_eq!(bed3_stranded.start(), exp_bed3_stranded.start());
            assert_eq!(bed3_stranded.end(), exp_bed3_stranded.end());
            assert_eq!(bed3_stranded.strand(), exp_bed3_stranded.strand());
        }
        Ok(())
    }

    #[test]
    fn range_intersects() {
        assert!((0..3).intersects(&(0..1)));
        assert!(!(0..3).intersects(&(5..7)));
        assert!(!(0..3).intersects(&(1..1)));
        assert!((1..3).intersects(&(0..2)));
        assert!((1..3).intersects(&(0..4)));
        assert!((0..4).intersects(&(1..3)));
    }
}

#[cfg(test)]
mod test_curr_read_align_and_mod_data {
    use super::*;
    use indoc::indoc;
    use nanalogue_core::read_utils::AlignAndModData;
    use nanalogue_core::{F32Bw0and1, InputWindowing};

    #[test]
    fn windowed_mod_data_restricted() -> Result<(), Error> {
        let input_json = indoc! {r#"
            {
              "alignment_type": "primary_forward",
              "alignment": {
                "start": 5,
                "end": 35,
                "contig": "chr1",
                "contig_id": 0
              },
              "mod_table": [
                {
                  "base": "T",
                  "is_strand_plus": true,
                  "mod_code": "T",
                  "implicit": false,
                  "data": [
                    [0, 10, 200],
                    [1, 15, 180],
                    [2, 20, 220],
                    [3, 25, 190],
                    [4, 30, 210]
                  ]
                },
                {
                  "base": "C",
                  "is_strand_plus": true,
                  "mod_code": "m",
                  "implicit": false,
                  "data": [
                    [5, 12, 150],
                    [6, 18, 160],
                    [7, 28, 170],
                    [8, 29, 180],
                    [9, 30, 190]
                  ]
                }
              ],
              "read_id": "test_read_123",
              "seq_len": 10
            }"#};

        let curr_read: CurrRead<AlignAndModData> = serde_json::from_str(input_json)?;
        let win_options: InputWindowing = serde_json::from_str(r#"{"win": 3, "step": 2}"#)?;
        let tag = ModChar::new('m');

        let window_function = |mod_data: &[u8]| -> Result<F32Bw0and1, Error> {
            let sum: f32 = mod_data.iter().map(|&x| f32::from(x)).sum();
            let mean = sum / f32::from(u8::try_from(mod_data.len())?);
            F32Bw0and1::new(mean / 255.0)
        };

        let result = curr_read.windowed_mod_data_restricted(&window_function, win_options, tag)?;

        let expected = vec![
            F32Bw0and1::new(160.0 / 255.0)?,
            F32Bw0and1::new(180.0 / 255.0)?,
        ];

        assert_eq!(result, expected);

        Ok(())
    }
}
