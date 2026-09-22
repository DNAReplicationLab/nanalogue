#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Public raw-call and window contracts for `RegionSequenceReader::profiles`.

#[cfg(feature = "bam-viewer")]
#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::region_sequences::{ReadModProfile, RegionSequenceReader};
    use nanalogue_core::{Error, ModChar, uuid, write_bam_denovo};
    use rust_htslib::bam::Record;
    use rust_htslib::bam::record::{Aux, Cigar, CigarString};
    use std::fs;
    use std::num::NonZeroU32;
    use std::path::{Path, PathBuf};

    /// A uniquely named directory removed together with all BAM sidecars.
    struct TempDir(PathBuf);

    impl TempDir {
        fn new(label: &str) -> Self {
            let path = std::env::temp_dir().join(format!(
                "nanalogue_region_profiles_{label}_{}",
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

    /// Builds one primary forward alignment at reference position 10.
    fn mapped_record(read_id: &str, cigar: Vec<Cigar>, sequence: &[u8]) -> Record {
        let mut record = Record::new();
        record.set(
            read_id.as_bytes(),
            Some(&CigarString(cigar)),
            sequence,
            &vec![35; sequence.len()],
        );
        record.set_flags(0);
        record.set_tid(0);
        record.set_pos(10);
        record.set_mapq(42);
        record
    }

    /// Attaches paired MM/ML tags to a record.
    fn with_modifications(mut record: Record, mm: &str, ml: &[u8]) -> Result<Record, Error> {
        record.push_aux(b"MM", Aux::String(mm))?;
        record.push_aux(b"ML", Aux::ArrayU8(ml.into()))?;
        Ok(record)
    }

    /// Writes records to a compact indexed BAM with one reference.
    fn write_fixture(path: &Path, records: Vec<Record>) -> Result<(), Error> {
        write_bam_denovo(
            records,
            [(String::from("chr1"), 40)],
            [String::from("profile-tests")],
            Vec::<String>::new(),
            path,
        )
    }

    /// Returns the sole profile, making fixture selection failures explicit.
    fn only_profile(profiles: &[ReadModProfile]) -> &ReadModProfile {
        assert_eq!(profiles.len(), 1, "the fixture has one spanning alignment");
        profiles.first().expect("the profile count was checked")
    }

    /// Checks one window's exact bounds and thresholded density.
    fn assert_window(profile: &ReadModProfile, index: usize, bounds: (u32, u32), value: f32) {
        let window = profile
            .windows()
            .get(index)
            .expect("the expected window should exist");
        assert_eq!((window.0, window.1), bounds);
        assert!((window.2.val() - value).abs() < f32::EPSILON);
    }

    /// Calls crossing an insertion, deletion, and reference skip retain their
    /// distinct raw and window semantics, including separate strand series.
    #[test]
    fn profiles_map_calls_and_complete_windows_across_cigar_gaps() -> Result<(), Error> {
        let temp = TempDir::new("cigar_windows");
        let path = temp.join("cigar-windows.bam");
        let record = mapped_record(
            "complex",
            vec![
                Cigar::Match(2),
                Cigar::Ins(1),
                Cigar::Match(2),
                Cigar::Del(1),
                Cigar::RefSkip(1),
                Cigar::Match(4),
            ],
            b"AAAATTTTC",
        );
        let tagged_record = with_modifications(
            record,
            "A+a?,0,0,0,0;T-a?,0,0,0,0;C+m?,0;",
            &[255, 0, 255, 64, 0, 128, 255, 127, 200],
        )?;
        write_fixture(&path, vec![tagged_record])?;
        let mut reader = RegionSequenceReader::from_path(&path)?;

        let a_profiles = reader.profiles(
            0,
            12,
            18,
            ModChar::new('a'),
            NonZeroU32::new(3).expect("three is nonzero"),
        )?;
        let a_profile = only_profile(&a_profiles);
        assert_eq!(a_profile.read_id(), "complex");
        assert_eq!((a_profile.align_start(), a_profile.align_end()), (10, 20));
        assert_eq!(
            a_profile.calls(),
            [
                (10, 255),
                (11, 0),
                (12, 64),
                (13, 0),
                (16, 128),
                (17, 255),
                (18, 127),
            ],
            "the insertion call is omitted while deletion and skip positions form a gap"
        );
        assert_eq!(a_profile.windows().len(), 2);
        assert_window(a_profile, 0, (10, 12), 2.0 / 3.0);
        assert_window(a_profile, 1, (13, 18), 2.0 / 3.0);
        assert_eq!(
            a_profile.window_series_starts(),
            [0, 1],
            "A+a and T-a windows remain independent"
        );

        let m_profiles = reader.profiles(
            0,
            10,
            20,
            ModChar::new('m'),
            NonZeroU32::new(1).expect("one is nonzero"),
        )?;
        let m_profile = only_profile(&m_profiles);
        assert_eq!(m_profile.calls(), [(19, 200)]);
        assert_eq!(m_profile.windows().len(), 1);
        assert_window(m_profile, 0, (19, 20), 1.0);
        assert_eq!(m_profile.window_series_starts(), [0]);
        Ok(())
    }

    /// Empty call states remain represented as profiles, while a complete
    /// insertion-only chunk is omitted without shifting later mapped windows.
    #[test]
    fn profiles_distinguish_empty_states_and_unmapped_windows() -> Result<(), Error> {
        let temp = TempDir::new("empty_states");
        let path = temp.join("empty-states.bam");
        let absent = mapped_record("absent-tags", vec![Cigar::Match(6)], b"ACGTAC");
        let no_match = with_modifications(
            mapped_record("different-type", vec![Cigar::Match(6)], b"CCCCCC"),
            "C+m?,0,0;",
            &[255, 128],
        )?;
        // The mismatched MM/ML counts are intentionally never parsed when SEQ is omitted.
        let omitted = with_modifications(
            mapped_record("omitted-sequence", vec![Cigar::Match(6)], b""),
            "A+a?,0,0;",
            &[255],
        )?;
        let insertion_first = with_modifications(
            mapped_record(
                "windowed",
                vec![Cigar::Ins(3), Cigar::Match(6)],
                b"AAAAAAAAA",
            ),
            "A+a?,0,0,0,0,0,0,0,0,0;",
            &[255, 255, 255, 0, 128, 255, 127, 128, 0],
        )?;
        write_fixture(&path, vec![absent, no_match, omitted, insertion_first])?;
        let mut reader = RegionSequenceReader::from_path(&path)?;

        let profiles = reader.profiles(
            0,
            10,
            16,
            ModChar::new('a'),
            NonZeroU32::new(3).expect("three is nonzero"),
        )?;
        assert_eq!(
            profiles
                .iter()
                .map(ReadModProfile::read_id)
                .collect::<Vec<_>>(),
            [
                "absent-tags",
                "different-type",
                "omitted-sequence",
                "windowed"
            ]
        );
        for profile in profiles.iter().take(3) {
            assert_eq!((profile.align_start(), profile.align_end()), (10, 16));
            assert!(profile.calls().is_empty());
            assert!(profile.windows().is_empty());
            assert!(profile.window_series_starts().is_empty());
        }

        let windowed = profiles.get(3).expect("windowed sorts last");
        assert_eq!(
            windowed.calls(),
            [(10, 0), (11, 128), (12, 255), (13, 127), (14, 128), (15, 0)]
        );
        assert_eq!(windowed.windows().len(), 2);
        assert_window(windowed, 0, (10, 13), 2.0 / 3.0);
        assert_window(windowed, 1, (13, 16), 1.0 / 3.0);
        assert_eq!(windowed.window_series_starts(), [0]);
        Ok(())
    }
}
