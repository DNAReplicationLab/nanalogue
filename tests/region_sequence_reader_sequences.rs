#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Public sequence-projection contracts for `RegionSequenceReader`.

#[cfg(feature = "bam-viewer")]
#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::region_sequences::{RegionSequence, RegionSequenceReader};
    use nanalogue_core::{Error, ModChar, uuid, write_bam_denovo};
    use rust_htslib::bam::Record;
    use rust_htslib::bam::record::{Aux, Cigar, CigarString};
    use std::fs;
    use std::path::{Path, PathBuf};

    /// A uniquely named directory removed together with all BAM sidecars.
    struct TempDir(PathBuf);

    impl TempDir {
        fn new(label: &str) -> Self {
            let path = std::env::temp_dir().join(format!(
                "nanalogue_region_sequences_{label}_{}",
                uuid::v4_random()
            ));
            fs::create_dir_all(&path).expect("temporary directory should be creatable");
            Self(path)
        }

        fn join(&self, name: &str) -> PathBuf {
            self.0.join(name)
        }
    }

    impl Drop for TempDir {
        fn drop(&mut self) {
            fs::remove_dir_all(&self.0).expect("temporary directory should be removable");
        }
    }

    /// Builds one mapped record at reference position 10.
    fn mapped_record(read_id: &str, reverse: bool, cigar: Vec<Cigar>, sequence: &[u8]) -> Record {
        let mut record = Record::new();
        record.set(
            read_id.as_bytes(),
            Some(&CigarString(cigar)),
            sequence,
            &vec![35; sequence.len()],
        );
        record.set_flags(if reverse { 16 } else { 0 });
        record.set_tid(0);
        record.set_pos(10);
        record.set_mapq(42);
        record
    }

    /// Writes records to a small indexed BAM with one reference.
    fn write_fixture(path: &Path, records: Vec<Record>) -> Result<(), Error> {
        write_bam_denovo(
            records,
            [(String::from("chr1"), 40)],
            [String::from("sequence-tests")],
            Vec::<String>::new(),
            path,
        )
    }

    /// Checks every public field on a projected sequence row.
    fn assert_row(
        row: &RegionSequence,
        read_id: &str,
        sequence: &str,
        sequence_with_insertions: &str,
        modifications: &[bool],
        modifications_with_insertions: &[bool],
        reverse: bool,
    ) {
        assert_eq!(row.read_id(), read_id);
        assert_eq!(row.sequence(), sequence);
        assert_eq!(row.sequence_with_insertions(), sequence_with_insertions);
        assert_eq!(row.modifications(), modifications);
        assert_eq!(
            row.modifications_with_insertions(),
            modifications_with_insertions
        );
        assert_eq!(row.is_reverse(), reverse);
    }

    /// A single public fetch combines sorting, strand state, CIGAR projection,
    /// insertion casing, omitted-sequence handling, and empty modification data.
    #[test]
    fn sequences_project_and_sort_complete_rows() -> Result<(), Error> {
        let temp = TempDir::new("projection");
        let path = temp.join("projection.bam");
        let complex_cigar = vec![
            Cigar::Match(2),
            Cigar::Ins(1),
            Cigar::Match(2),
            Cigar::Del(1),
            Cigar::RefSkip(1),
            Cigar::Match(2),
        ];
        let late = mapped_record("zeta", false, complex_cigar, b"ACATACA");
        let omitted = mapped_record("middle", false, vec![Cigar::Match(8)], b"");
        let early = mapped_record("alpha", true, vec![Cigar::Match(8)], b"TGCATGCC");
        write_fixture(&path, vec![late, omitted, early])?;
        let mut reader = RegionSequenceReader::from_path(&path)?;

        let rows = reader.sequences(0, 10, 18, Some(ModChar::new('a')))?;
        assert_eq!(rows.len(), 3);
        assert_row(
            rows.first().expect("alpha sorts first"),
            "alpha",
            "TGCATGCC",
            "TGCATGCC",
            &[false; 8],
            &[false; 8],
            true,
        );
        assert_row(
            rows.get(1).expect("middle sorts second"),
            "middle",
            "*",
            "*",
            &[false],
            &[false],
            false,
        );
        assert_row(
            rows.get(2).expect("zeta sorts last"),
            "zeta",
            "ACTA..CA",
            "ACaTA..CA",
            &[false; 8],
            &[false; 9],
            false,
        );
        Ok(())
    }

    /// The ML cutoff is inclusive at 128, insertion calls stay aligned with
    /// lowercase insertion bases, and selecting one type excludes another.
    #[test]
    fn sequences_filter_type_at_inclusive_threshold() -> Result<(), Error> {
        let temp = TempDir::new("modifications");
        let path = temp.join("modifications.bam");
        let mut record = mapped_record(
            "modified",
            false,
            vec![
                Cigar::Match(2),
                Cigar::Ins(1),
                Cigar::Match(2),
                Cigar::Del(1),
                Cigar::RefSkip(1),
                Cigar::Match(2),
            ],
            b"ACATACA",
        );
        record.push_aux(b"MM", Aux::String("A+a?,0,0,0,0;C+m?,0,0;"))?;
        record.push_aux(
            b"ML",
            Aux::ArrayU8((&[127u8, 255, 128, 200, 255, 0][..]).into()),
        )?;
        write_fixture(&path, vec![record])?;
        let mut reader = RegionSequenceReader::from_path(&path)?;

        let a_rows = reader.sequences(0, 10, 18, Some(ModChar::new('a')))?;
        assert_eq!(a_rows.len(), 1);
        assert_row(
            a_rows.first().expect("one a-modification row"),
            "modified",
            "ACTA..CA",
            "ACaTA..CA",
            &[false, false, false, true, false, false, false, true],
            &[false, false, true, false, true, false, false, false, true],
            false,
        );

        let m_rows = reader.sequences(0, 10, 18, Some(ModChar::new('m')))?;
        assert_eq!(m_rows.len(), 1);
        assert_row(
            m_rows.first().expect("one m-modification row"),
            "modified",
            "ACTA..CA",
            "ACaTA..CA",
            &[false, true, false, false, false, false, false, false],
            &[false, true, false, false, false, false, false, false, false],
            false,
        );
        Ok(())
    }

    /// Missing tags yield an unmodified row when a type is requested, while
    /// malformed paired tags are ignored only when modification parsing is omitted.
    #[test]
    fn sequences_distinguish_absent_and_malformed_tags() -> Result<(), Error> {
        let temp = TempDir::new("tag_states");
        let absent_path = temp.join("absent.bam");
        write_fixture(
            &absent_path,
            vec![mapped_record(
                "absent",
                false,
                vec![Cigar::Match(8)],
                b"ACGTACGT",
            )],
        )?;
        let mut absent_reader = RegionSequenceReader::from_path(&absent_path)?;
        let absent_rows = absent_reader.sequences(0, 10, 18, Some(ModChar::new('m')))?;
        assert_eq!(absent_rows.len(), 1);
        assert_row(
            absent_rows.first().expect("one absent-tag row"),
            "absent",
            "ACGTACGT",
            "ACGTACGT",
            &[false; 8],
            &[false; 8],
            false,
        );

        let malformed_path = temp.join("malformed.bam");
        let mut malformed = mapped_record("malformed", true, vec![Cigar::Match(8)], b"ACGTACGT");
        malformed.push_aux(b"MM", Aux::String("C+m?,0,0;"))?;
        malformed.push_aux(b"ML", Aux::ArrayU8((&[200u8][..]).into()))?;
        write_fixture(&malformed_path, vec![malformed])?;
        let mut malformed_reader = RegionSequenceReader::from_path(&malformed_path)?;

        let unparsed_rows = malformed_reader.sequences(0, 10, 18, None)?;
        assert_eq!(unparsed_rows.len(), 1);
        assert_row(
            unparsed_rows.first().expect("one unparsed malformed row"),
            "malformed",
            "ACGTACGT",
            "ACGTACGT",
            &[false; 8],
            &[false; 8],
            true,
        );
        let error = malformed_reader
            .sequences(0, 10, 18, Some(ModChar::new('m')))
            .expect_err("requesting modifications must reject unequal MM/ML counts");
        assert!(matches!(error, Error::InvalidModProbs(_)));
        Ok(())
    }
}
