#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Error propagation through the public indexed region reader APIs.

#[cfg(feature = "bam-viewer")]
#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::region_sequences::{ReadModProfile, RegionSequenceReader};
    use nanalogue_core::{Error, ModChar, uuid, write_bam_denovo};
    use rust_htslib::bam::Record;
    use rust_htslib::bam::record::{Aux, Cigar, CigarString};
    use rust_htslib::errors::Error as HtslibError;
    use std::fs;
    use std::num::NonZeroU32;
    use std::path::{Path, PathBuf};

    /// A uniquely named directory removed together with all BAM sidecars.
    struct TempDir(PathBuf);

    impl TempDir {
        fn new(label: &str) -> Self {
            let path = std::env::temp_dir().join(format!(
                "nanalogue_region_errors_{label}_{}",
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

    /// Builds a mapped record spanning the reference interval 10..18.
    fn mapped_record(read_id: &[u8], flags: u16, cigar: Vec<Cigar>, sequence: &[u8]) -> Record {
        let mut record = Record::new();
        record.set(
            read_id,
            Some(&CigarString(cigar)),
            sequence,
            &vec![35; sequence.len()],
        );
        record.set_flags(flags);
        record.set_tid(0);
        record.set_pos(10);
        record.set_mapq(42);
        record
    }

    /// Writes one record to a compact indexed BAM with one reference.
    fn write_fixture(path: &Path, record: Record) -> Result<(), Error> {
        write_bam_denovo(
            [record],
            [(String::from("chr1"), 40)],
            [String::from("region-error-tests")],
            Vec::<String>::new(),
            path,
        )
    }

    /// Calls the profiles API with the common fixture interval and options.
    fn profiles(reader: &mut RegionSequenceReader) -> Result<Vec<ReadModProfile>, Error> {
        reader.profiles(
            0,
            10,
            18,
            ModChar::new('m'),
            NonZeroU32::new(2).expect("two is nonzero"),
        )
    }

    /// Unsupported and contradictory alignment flags are rejected by both
    /// public retrieval methods after the indexed record has been fetched.
    #[test]
    fn retrieval_propagates_invalid_alignment_states() -> Result<(), Error> {
        let temp = TempDir::new("alignment_states");
        for (name, flags) in [("paired", 0x1u16), ("contradictory", 0x900)] {
            let path = temp.join(&format!("{name}.bam"));
            write_fixture(
                &path,
                mapped_record(name.as_bytes(), flags, vec![Cigar::Match(8)], b"ACGTACGT"),
            )?;

            let mut sequence_reader = RegionSequenceReader::from_path(&path)?;
            let sequence_error = sequence_reader
                .sequences(0, 10, 18, None)
                .expect_err("invalid flags must fail sequence retrieval");
            let mut profile_reader = RegionSequenceReader::from_path(&path)?;
            let profile_error = profiles(&mut profile_reader)
                .expect_err("invalid flags must fail profile retrieval");

            if flags == 0x1 {
                for error in [sequence_error, profile_error] {
                    assert!(matches!(
                        error,
                        Error::NotImplemented(message)
                            if message == concat!(
                                "paired-read/mate-read/duplicate/qual-check-failed flags ",
                                "not supported! read_id: paired"
                            )
                    ));
                }
            } else {
                for error in [sequence_error, profile_error] {
                    assert!(matches!(
                        error,
                        Error::UnknownAlignState(message)
                            if message == "invalid flag combination! read_id: contradictory"
                    ));
                }
            }
        }
        Ok(())
    }

    /// Read-ID validation is not bypassed by an otherwise valid indexed
    /// alignment, and each public retrieval method preserves the exact error.
    #[test]
    fn retrieval_propagates_invalid_read_ids() -> Result<(), Error> {
        let temp = TempDir::new("read_id");
        let path = temp.join("invalid-read-id.bam");
        write_fixture(
            &path,
            mapped_record(b"invalid,id", 0, vec![Cigar::Match(8)], b"ACGTACGT"),
        )?;

        let mut sequence_reader = RegionSequenceReader::from_path(&path)?;
        let sequence_error = sequence_reader
            .sequences(0, 10, 18, None)
            .expect_err("an unsafe read ID must fail sequence retrieval");
        let mut profile_reader = RegionSequenceReader::from_path(&path)?;
        let profile_error = profiles(&mut profile_reader)
            .expect_err("an unsafe read ID must fail profile retrieval");

        for error in [sequence_error, profile_error] {
            assert!(matches!(
                error,
                Error::InvalidReadID(message)
                    if message == "read_id contains forbidden characters"
            ));
        }
        Ok(())
    }

    /// Record decoding errors that occur after a successful indexed open are
    /// returned by both public iterators rather than mistaken for end-of-file.
    #[test]
    fn retrieval_propagates_record_decode_errors() -> Result<(), Error> {
        let temp = TempDir::new("decode_error");
        for api in ["sequences", "profiles"] {
            let path = temp.join(&format!("{api}.bam"));
            write_fixture(
                &path,
                mapped_record(api.as_bytes(), 0, vec![Cigar::Match(8)], b"ACGTACGT"),
            )?;
            let file_len = fs::metadata(&path)?.len();
            fs::OpenOptions::new().write(true).open(&path)?.set_len(
                file_len
                    .checked_sub(29)
                    .expect("fixture contains a data block and 28-byte EOF marker"),
            )?;
            let mut reader = RegionSequenceReader::from_path(&path)?;

            let error = if api == "sequences" {
                reader
                    .sequences(0, 10, 18, None)
                    .expect_err("truncated record data must fail sequence retrieval")
            } else {
                profiles(&mut reader)
                    .expect_err("truncated record data must fail profile retrieval")
            };
            assert!(matches!(
                error,
                Error::RustHtslibError(source)
                    if matches!(source.as_ref(), HtslibError::BamTruncatedRecord)
            ));
        }
        Ok(())
    }

    /// Both retrieval APIs reject a mapped record without a CIGAR rather than
    /// accepting `HTSlib`'s synthetic one-base reference span.
    #[test]
    fn missing_cigar_is_rejected_before_projection() -> Result<(), Error> {
        let temp = TempDir::new("missing_cigar");
        let path = temp.join("missing-cigar.bam");
        let mut record = Record::new();
        record.set(b"missing-cigar", None, b"A", &[35]);
        record.set_flags(0);
        record.set_tid(0);
        record.set_pos(10);
        record.set_mapq(42);
        write_fixture(&path, record)?;

        let mut sequence_reader = RegionSequenceReader::from_path(&path)?;
        let sequence_error = sequence_reader
            .sequences(0, 10, 11, None)
            .expect_err("sequence retrieval must reject a mapped record without a CIGAR");

        let mut profile_reader = RegionSequenceReader::from_path(&path)?;
        let profile_error = profile_reader
            .profiles(
                0,
                10,
                11,
                ModChar::new('m'),
                NonZeroU32::new(2).expect("two is nonzero"),
            )
            .expect_err("profile retrieval must reject a mapped record without a CIGAR");

        for error in [sequence_error, profile_error] {
            assert!(matches!(
                error,
                Error::InvalidAlignLength(message)
                    if message == "mapped read has no CIGAR, read_id: missing-cigar"
            ));
        }
        Ok(())
    }

    /// Omitting modification parsing allows a wrong-typed ML tag through
    /// sequence projection; requesting calls through either API rejects it.
    #[test]
    fn wrong_typed_ml_is_parsed_only_when_requested() -> Result<(), Error> {
        let temp = TempDir::new("wrong_ml_type");
        let path = temp.join("wrong-ml-type.bam");
        let mut record = mapped_record(b"wrong-ml", 0, vec![Cigar::Match(8)], b"ACGTACGT");
        record.push_aux(b"MM", Aux::String("C+m?,0;"))?;
        record.push_aux(b"ML", Aux::String("200"))?;
        write_fixture(&path, record)?;

        let mut reader = RegionSequenceReader::from_path(&path)?;
        let rows = reader.sequences(0, 10, 18, None)?;
        assert_eq!(rows.len(), 1);
        assert_eq!(
            rows.first().expect("one sequence row").sequence(),
            "ACGTACGT"
        );

        for error in [
            reader
                .sequences(0, 10, 18, Some(ModChar::new('m')))
                .expect_err("requesting sequence modifications must parse ML"),
            profiles(&mut RegionSequenceReader::from_path(&path)?)
                .expect_err("profiles must parse ML"),
        ] {
            assert!(matches!(
                error,
                Error::InvalidState(message)
                    if message
                        == "rust-htslib ML/Ml tag parsing failure: unexpected auxiliary tag type"
            ));
        }
        Ok(())
    }

    /// Profile retrieval preserves distinct parser failures for a wrong-typed
    /// MM tag, too few ML values, and probabilities without position calls.
    #[test]
    fn profiles_preserve_distinct_malformed_tag_errors() -> Result<(), Error> {
        let temp = TempDir::new("malformed_tags");

        let wrong_mm_path = temp.join("wrong-mm-type.bam");
        let mut wrong_mm = mapped_record(b"wrong-mm", 0, vec![Cigar::Match(8)], b"ACGTACGT");
        wrong_mm.push_aux(b"MM", Aux::I32(1))?;
        wrong_mm.push_aux(b"ML", Aux::ArrayU8((&[200u8][..]).into()))?;
        write_fixture(&wrong_mm_path, wrong_mm)?;
        let wrong_mm_error = profiles(&mut RegionSequenceReader::from_path(&wrong_mm_path)?)
            .expect_err("MM must be a string");
        assert!(matches!(
            wrong_mm_error,
            Error::InvalidState(message)
                if message
                    == "rust-htslib MM/Mm tag parsing failure: unexpected auxiliary tag type"
        ));

        let short_ml_path = temp.join("short-ml.bam");
        let mut short_ml = mapped_record(b"short-ml", 0, vec![Cigar::Match(8)], b"ACGTACGT");
        short_ml.push_aux(b"MM", Aux::String("C+m?,0,0;"))?;
        short_ml.push_aux(b"ML", Aux::ArrayU8((&[200u8][..]).into()))?;
        write_fixture(&short_ml_path, short_ml)?;
        let short_ml_error = profiles(&mut RegionSequenceReader::from_path(&short_ml_path)?)
            .expect_err("every MM call needs an ML probability");
        assert!(matches!(
            short_ml_error,
            Error::InvalidModProbs(message)
                if message == "ML tag appears to be insufficiently long!"
        ));

        let orphan_ml_path = temp.join("orphan-ml.bam");
        let mut orphan_ml = mapped_record(b"orphan-ml", 0, vec![Cigar::Match(8)], b"ACGTACGT");
        orphan_ml.push_aux(b"ML", Aux::ArrayU8((&[200u8][..]).into()))?;
        write_fixture(&orphan_ml_path, orphan_ml)?;
        let orphan_ml_error = profiles(&mut RegionSequenceReader::from_path(&orphan_ml_path)?)
            .expect_err("ML probabilities without MM calls must fail");
        assert!(matches!(
            orphan_ml_error,
            Error::InvalidModProbs(message)
                if message == "MM and ML tag lengths do not match!"
        ));
        Ok(())
    }
}
