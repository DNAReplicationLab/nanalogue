//! Tests for read_utils.rs extracted from doctests

use bedrs::prelude::StrandedBed3;
use bedrs::{Bed3, Coordinates, Strand};
use nanalogue_core::{
    CurrRead, Error, Intersects, ModChar, ReadState, ThresholdState, nanalogue_bam_reader,
};
use rust_htslib::bam::Read;
use std::collections::HashMap;

#[test]
fn test_set_read_state() -> Result<(), Error> {
    let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    let mut count = 0;
    for record in reader.records() {
        let r = record?;
        let curr_read = CurrRead::default().set_read_state(&r)?;
        match count {
            0 => assert_eq!(curr_read.read_state(), ReadState::PrimaryFwd),
            1 => assert_eq!(curr_read.read_state(), ReadState::PrimaryFwd),
            2 => assert_eq!(curr_read.read_state(), ReadState::PrimaryRev),
            3 => assert_eq!(curr_read.read_state(), ReadState::Unmapped),
            _ => unreachable!(),
        }
        count = count + 1;
    }
    Ok(())
}

#[test]
fn test_set_seq_len() -> Result<(), Error> {
    let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    let mut count = 0;
    for record in reader.records() {
        let r = record?;
        let curr_read = CurrRead::default().set_read_state(&r)?.set_seq_len(&r)?;
        let Ok(len) = curr_read.seq_len() else {
            unreachable!()
        };
        match count {
            0 => assert_eq!(len, 8),
            1 => assert_eq!(len, 48),
            2 => assert_eq!(len, 33),
            3 => assert_eq!(len, 48),
            _ => unreachable!(),
        }
        count = count + 1;
    }
    Ok(())
}

#[test]
#[should_panic]
fn test_set_seq_len_duplicate_should_panic() {
    let mut reader = nanalogue_bam_reader(&"examples/example_1.bam").unwrap();
    for record in reader.records() {
        let r = record.unwrap();
        let _curr_read = CurrRead::default()
            .set_read_state(&r)
            .unwrap()
            .set_seq_len(&r)
            .unwrap()
            .set_seq_len(&r)
            .unwrap();
        break;
    }
}

#[test]
fn test_set_contig_id_and_start() -> Result<(), Error> {
    let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
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
        count = count + 1;
        if count == 3 {
            break;
        } // the fourth entry is unmapped, and will lead to an error.
    }
    Ok(())
}

#[test]
#[should_panic]
fn test_set_contig_id_and_start_unmapped_should_panic() {
    let mut reader = nanalogue_bam_reader(&"examples/example_1.bam").unwrap();
    let mut count = 0;
    for record in reader.records() {
        let r = record.unwrap();
        if count < 3 {
            count = count + 1;
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
#[should_panic]
fn test_set_contig_id_and_start_duplicate_should_panic() {
    let mut reader = nanalogue_bam_reader(&"examples/example_1.bam").unwrap();
    for record in reader.records() {
        let r = record.unwrap();
        let _curr_read = CurrRead::default()
            .set_read_state(&r)
            .unwrap()
            .set_contig_id_and_start(&r)
            .unwrap()
            .set_contig_id_and_start(&r)
            .unwrap();
        break;
    }
}

#[test]
fn test_set_contig_name() -> Result<(), Error> {
    let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
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
        count = count + 1;
        if count == 3 {
            break;
        } // the fourth entry is unmapped, and will lead to an error.
    }
    Ok(())
}

#[test]
#[should_panic]
fn test_set_contig_name_unmapped_should_panic() {
    let mut reader = nanalogue_bam_reader(&"examples/example_1.bam").unwrap();
    let mut count = 0;
    for record in reader.records() {
        if count < 3 {
            count = count + 1;
            continue;
        }
        let r = record.unwrap();
        let curr_read = CurrRead::default().set_read_state(&r).unwrap();
        curr_read.set_contig_name(&r).unwrap();
    }
}

#[test]
#[should_panic]
fn test_set_contig_name_duplicate_should_panic() {
    let mut reader = nanalogue_bam_reader(&"examples/example_1.bam").unwrap();
    for record in reader.records() {
        let r = record.unwrap();
        let _curr_read = CurrRead::default()
            .set_read_state(&r)
            .unwrap()
            .set_contig_name(&r)
            .unwrap()
            .set_contig_name(&r)
            .unwrap();
        break;
    }
}

#[test]
fn test_set_read_id() -> Result<(), Error> {
    let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    let mut count = 0;
    for record in reader.records() {
        let r = record?;
        let curr_read = CurrRead::default().set_read_state(&r)?.set_read_id(&r)?;
        let Ok(read_id) = curr_read.read_id() else {
            unreachable!()
        };
        match (count, read_id) {
            (0, "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575")
            | (1, "a4f36092-b4d5-47a9-813e-c22c3b477a0c")
            | (2, "fffffff1-10d2-49cb-8ca3-e8d48979001b")
            | (3, "a4f36092-b4d5-47a9-813e-c22c3b477a0c") => {}
            _ => unreachable!(),
        }
        count = count + 1;
    }
    Ok(())
}

#[test]
#[should_panic]
fn test_set_read_id_duplicate_should_panic() {
    let mut reader = nanalogue_bam_reader(&"examples/example_1.bam").unwrap();
    for record in reader.records() {
        let r = record.unwrap();
        let _curr_read = CurrRead::default()
            .set_read_state(&r)
            .unwrap()
            .set_read_id(&r)
            .unwrap()
            .set_read_id(&r)
            .unwrap();
        break;
    }
}

#[test]
fn test_strand() -> Result<(), Error> {
    let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    let mut count = 0;
    for record in reader.records() {
        let r = record?;
        let curr_read = CurrRead::default().set_read_state(&r)?;
        let strand = curr_read.strand();
        match (count, strand) {
            (0, '+') | (1, '+') | (2, '-') | (3, '.') => {}
            _ => unreachable!(),
        }
        count = count + 1;
    }
    Ok(())
}

#[test]
fn test_seq_on_ref_coords() -> Result<(), Error> {
    let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
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
fn test_basecount_per_mod() -> Result<(), Error> {
    let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    let mut count = 0;
    for record in reader.records() {
        let r = record?;
        let curr_read = CurrRead::default().set_read_state(&r)?.set_mod_data(
            &r,
            ThresholdState::GtEq(180),
            0,
        )?;
        let modcount = curr_read.base_count_per_mod();
        let zerocount = Some(HashMap::from([(ModChar::new('T'), 0)]));
        let a = Some(HashMap::from([(ModChar::new('T'), 3)]));
        let b = Some(HashMap::from([(ModChar::new('T'), 1)]));
        let c = Some(HashMap::from([
            (ModChar::new('T'), 3),
            (ModChar::new('ᰠ'), 0),
        ]));
        match (count, modcount) {
            (0, v) => assert_eq!(v, zerocount),
            (1, v) => assert_eq!(v, a),
            (2, v) => assert_eq!(v, b),
            (3, v) => assert_eq!(v, c),
            _ => unreachable!(),
        }
        count = count + 1;
    }
    Ok(())
}

#[test]
fn test_try_from_stranded_bed3() -> Result<(), Error> {
    let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    let mut count = 0;
    for record in reader.records() {
        let r = record?;
        let curr_read = CurrRead::default().try_from_only_alignment(&r)?;
        let Ok(bed3_stranded) = StrandedBed3::try_from(&curr_read) else {
            unreachable!()
        };
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
        count = count + 1;
    }
    Ok(())
}

#[test]
fn test_range_intersects() {
    assert!((0..3).intersects(&(0..1)));
    assert!(!(0..3).intersects(&(5..7)));
    assert!(!(0..3).intersects(&(1..1)));
    assert!((1..3).intersects(&(0..2)));
}

#[cfg(test)]
mod test_curr_read_align_and_mod_data {
    use super::*;
    use nanalogue_core::{InputWindowing, ModChar, F32Bw0and1};
    use nanalogue_core::read_utils::AlignAndModData;
    use indoc::indoc;

    #[test]
    fn test_windowed_mod_data_restricted() -> Result<(), Error> {
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
            let sum: f32 = mod_data.iter().map(|&x| x as f32).sum();
            let mean = sum / (mod_data.len() as f32);
            F32Bw0and1::new(mean / 256.0)
        };

        let result = curr_read.windowed_mod_data_restricted(
            &window_function,
            win_options,
            tag
        )?;

        let expected = vec![
            F32Bw0and1::new(160.0 / 256.0)?,
            F32Bw0and1::new(180.0 / 256.0)?,
        ];

        assert_eq!(result, expected);

        Ok(())
    }
}
