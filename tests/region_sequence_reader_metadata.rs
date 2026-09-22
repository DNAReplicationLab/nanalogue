#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Public metadata and region-validation contracts for `RegionSequenceReader`.

#[cfg(feature = "bam-viewer")]
#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::region_sequences::RegionSequenceReader;
    use nanalogue_core::{Error, ModChar, uuid, write_bam_denovo};
    use rust_htslib::bam::Record;
    use rust_htslib::bam::record::{Cigar, CigarString};
    use rust_htslib::errors::Error as HtslibError;
    use std::fs;
    use std::num::NonZeroU32;
    use std::path::{Path, PathBuf};

    /// A uniquely named directory removed together with all BAM sidecars.
    struct TempDir(PathBuf);

    impl TempDir {
        fn new(label: &str) -> Self {
            let path = std::env::temp_dir().join(format!(
                "nanalogue_region_reader_{label}_{}",
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

    /// A mapped record that overlaps, but does not span, the interval 4..16.
    fn mapped_record() -> Record {
        let mut record = Record::new();
        record.set(
            b"partial-span",
            Some(&CigarString(vec![Cigar::Match(10)])),
            b"ACGTACGTAC",
            &[35; 10],
        );
        record.set_flags(0);
        record.set_tid(0);
        record.set_pos(5);
        record.set_mapq(42);
        record
    }

    /// Writes the shared two-reference, indexed BAM fixture.
    fn write_fixture(path: &Path) -> Result<(), Error> {
        write_bam_denovo(
            [mapped_record()],
            [(String::from("zeta"), 40), (String::from("alpha"), 12)],
            [String::from("reader-tests")],
            Vec::<String>::new(),
            path,
        )
    }

    /// Opening distinguishes an absent BAM, an absent index, and a valid
    /// indexed BAM whose header has no references.
    #[test]
    fn open_reports_input_and_header_failures() -> Result<(), Error> {
        let temp = TempDir::new("open_failures");
        let missing = temp.join("missing.bam");
        let missing_error =
            RegionSequenceReader::from_path(&missing).expect_err("an absent BAM should not open");
        assert!(matches!(
            missing_error,
            Error::RustHtslibError(source)
                if matches!(source.as_ref(), HtslibError::FileNotFound { path } if path == &missing)
        ));

        let unindexed = temp.join("unindexed.bam");
        write_fixture(&unindexed)?;
        fs::remove_file(temp.join("unindexed.bam.bai"))?;
        let index_error = RegionSequenceReader::from_path(&unindexed)
            .expect_err("a BAM without its BAI should not open");
        assert!(matches!(
            index_error,
            Error::RustHtslibError(source)
                if matches!(source.as_ref(), HtslibError::BamInvalidIndex { target }
                    if target == &unindexed.display().to_string())
        ));

        let reference_free = temp.join("reference-free.bam");
        write_bam_denovo(
            Vec::<Record>::new(),
            Vec::<(String, usize)>::new(),
            Vec::<String>::new(),
            Vec::<String>::new(),
            &reference_free,
        )?;
        let header_error = RegionSequenceReader::from_path(&reference_free)
            .expect_err("a reference-free BAM should not open");
        assert!(matches!(
            header_error,
            Error::RustHtslibError(source) if matches!(source.as_ref(), HtslibError::Fetch)
        ));
        Ok(())
    }

    /// Name, identifier, and length lookups preserve header order and return
    /// `None` immediately beyond either end of the public metadata tables.
    #[test]
    fn target_metadata_lookups_cover_boundaries() -> Result<(), Error> {
        let temp = TempDir::new("metadata");
        let path = temp.join("metadata.bam");
        write_fixture(&path)?;
        let reader = RegionSequenceReader::from_path(&path)?;

        assert_eq!(reader.target_id("zeta"), Some(0));
        assert_eq!(reader.target_id("alpha"), Some(1));
        assert_eq!(reader.target_id("ZETA"), None);
        assert_eq!(reader.target_id("missing"), None);

        assert_eq!(reader.target_name(0), Some("zeta"));
        assert_eq!(reader.target_name(1), Some("alpha"));
        assert_eq!(reader.target_name(2), None);
        assert_eq!(reader.target_name(u32::MAX), None);

        assert_eq!(reader.target_len(0), Some(40));
        assert_eq!(reader.target_len(1), Some(12));
        assert_eq!(reader.target_len(2), None);
        assert_eq!(reader.target_len(u32::MAX), None);
        Ok(())
    }

    /// Both retrieval APIs reject empty and reversed intervals before asking
    /// `HTSlib` to fetch, retaining the exact rejected coordinates in the error.
    #[test]
    fn retrieval_rejects_non_increasing_ranges() -> Result<(), Error> {
        let temp = TempDir::new("invalid_ranges");
        let path = temp.join("ranges.bam");
        write_fixture(&path)?;
        let mut reader = RegionSequenceReader::from_path(&path)?;
        let win = NonZeroU32::new(3).expect("three is nonzero");

        for (start, end) in [(7, 7), (9, 3)] {
            let expected = format!("0:{start}-{end}");
            let sequence_error = reader
                .sequences(0, start, end, None)
                .expect_err("sequences should reject a non-increasing range");
            assert!(matches!(
                sequence_error,
                Error::InvalidAlignCoords(coordinates) if coordinates == expected
            ));

            let profile_error = reader
                .profiles(0, start, end, ModChar::new('m'), win)
                .expect_err("profiles should reject a non-increasing range");
            assert!(matches!(
                profile_error,
                Error::InvalidAlignCoords(coordinates) if coordinates == expected
            ));
        }
        Ok(())
    }

    /// `HTSlib` accepts unknown nonnegative targets and coordinates beyond the
    /// BAI address space as successful fetches with no indexed records.
    #[test]
    fn out_of_range_fetches_have_no_records() -> Result<(), Error> {
        let temp = TempDir::new("out_of_range_fetches");
        let path = temp.join("fetches.bam");
        write_fixture(&path)?;
        let mut reader = RegionSequenceReader::from_path(&path)?;
        let win = NonZeroU32::new(2).expect("two is nonzero");

        assert_eq!(
            reader.sequences(2, 0, 1, None)?,
            [],
            "an unknown nonnegative target has no indexed records"
        );

        assert_eq!(
            reader.profiles(0, u32::MAX - 1, u32::MAX, ModChar::new('m'), win,)?,
            [],
            "coordinates beyond the BAI range have no indexed records"
        );
        Ok(())
    }

    /// Overlapping records must cover both interval boundaries to be returned.
    #[test]
    fn valid_interval_without_spanning_records_is_empty() -> Result<(), Error> {
        let temp = TempDir::new("no_spanning_records");
        let path = temp.join("partial.bam");
        write_fixture(&path)?;
        let mut reader = RegionSequenceReader::from_path(&path)?;

        assert_eq!(
            reader.sequences(0, 4, 14, None)?,
            [],
            "a record starting after the interval does not span it"
        );
        assert_eq!(
            reader.profiles(
                0,
                6,
                16,
                ModChar::new('m'),
                NonZeroU32::new(2).expect("two is nonzero"),
            )?,
            [],
            "a record ending before the interval does not span it"
        );

        let rows = reader.sequences(0, 5, 15, None)?;
        assert_eq!(rows.len(), 1, "the exact alignment interval is spanned");
        let row = rows.first().expect("the row count was checked");
        assert_eq!(row.read_id(), "partial-span");
        assert_eq!(row.sequence(), "ACGTACGTAC");

        let profiles = reader.profiles(
            0,
            5,
            15,
            ModChar::new('m'),
            NonZeroU32::new(2).expect("two is nonzero"),
        )?;
        assert_eq!(profiles.len(), 1, "the exact interval has one profile");
        let profile = profiles.first().expect("the profile count was checked");
        assert_eq!(profile.read_id(), "partial-span");
        assert_eq!((profile.align_start(), profile.align_end()), (5, 15));
        assert_eq!(profile.calls(), []);
        Ok(())
    }
}
