#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Public identity and ordering contracts for modification profiles.

#[cfg(feature = "bam-viewer")]
#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::region_sequences::{ReadModProfile, RegionSequenceReader};
    use nanalogue_core::{Error, ModChar, uuid, write_bam_denovo};
    use rust_htslib::bam::Record;
    use rust_htslib::bam::record::{Cigar, CigarString};
    use std::fs;
    use std::num::NonZeroU32;
    use std::path::{Path, PathBuf};

    /// A uniquely named directory removed together with its BAM sidecars.
    struct TempDir(PathBuf);

    impl TempDir {
        fn new() -> Self {
            let path = std::env::temp_dir()
                .join(format!("nanalogue_profile_identity_{}", uuid::v4_random()));
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

    /// Builds one mapped record with deliberately configurable identity fields.
    fn mapped_record(
        read_id: &str,
        start: i64,
        reference_length: u32,
        reverse: bool,
    ) -> Result<Record, Error> {
        let sequence_length = usize::try_from(reference_length)?;
        let mut record = Record::new();
        record.set(
            read_id.as_bytes(),
            Some(&CigarString(vec![Cigar::Match(reference_length)])),
            &vec![b'A'; sequence_length],
            &vec![35; sequence_length],
        );
        record.set_flags(0);
        record.set_tid(0);
        record.set_pos(start);
        record.set_mapq(42);
        if reverse {
            record.set_reverse();
        }
        Ok(record)
    }

    /// Writes coordinate-sorted records whose read IDs are not lexically sorted.
    fn write_fixture(path: &Path) -> Result<(), Error> {
        let records = [
            mapped_record("zeta", 2, 30, false)?,
            mapped_record("shared", 3, 30, false)?,
            mapped_record("shared", 3, 30, false)?,
            mapped_record("shared", 3, 31, false)?,
            mapped_record("shared", 3, 30, true)?,
            mapped_record("shared", 4, 29, false)?,
            mapped_record("alpha", 5, 30, false)?,
        ];
        write_bam_denovo(
            records,
            [(String::from("chr1"), 40)],
            [String::from("identity-tests")],
            Vec::<String>::new(),
            path,
        )
    }

    /// Returns the sole profile matching public alignment properties.
    fn unique_profile(
        profiles: &[ReadModProfile],
        start: u32,
        end: u32,
        reverse: bool,
    ) -> &ReadModProfile {
        let matches = profiles
            .iter()
            .filter(|profile| {
                profile.read_id() == "shared"
                    && profile.align_start() == start
                    && profile.align_end() == end
                    && profile.is_reverse() == reverse
            })
            .collect::<Vec<_>>();
        assert_eq!(matches.len(), 1, "the public properties select one profile");
        matches
            .first()
            .expect("the matching profile count was checked")
    }

    /// Profiles expose alignment properties, sort by read ID, distinguish
    /// same-name alignments and duplicates, and retain identity after refetch.
    #[test]
    fn profiles_preserve_public_alignment_identity_across_fetches() -> Result<(), Error> {
        let temp = TempDir::new();
        let path = temp.join("identity.bam");
        write_fixture(&path)?;
        let mut reader = RegionSequenceReader::from_path(&path)?;
        let win = NonZeroU32::new(3).expect("three is nonzero");

        let first_fetch = reader.profiles(0, 10, 20, ModChar::new('a'), win)?;
        assert_eq!(
            first_fetch
                .iter()
                .map(ReadModProfile::read_id)
                .collect::<Vec<_>>(),
            [
                "alpha", "shared", "shared", "shared", "shared", "shared", "zeta"
            ],
            "profiles sort by read ID rather than coordinate-sorted BAM input order"
        );

        let alpha = first_fetch.first().expect("alpha sorts first");
        assert_eq!(alpha.read_id(), "alpha");
        assert!(!alpha.is_reverse());
        assert_eq!((alpha.align_start(), alpha.align_end()), (5, 35));
        let zeta = first_fetch.last().expect("zeta sorts last");
        assert_eq!(zeta.read_id(), "zeta");
        assert!(!zeta.is_reverse());
        assert_eq!((zeta.align_start(), zeta.align_end()), (2, 32));

        let exact_duplicates = first_fetch
            .iter()
            .filter(|profile| {
                profile.read_id() == "shared"
                    && profile.align_start() == 3
                    && profile.align_end() == 33
                    && !profile.is_reverse()
            })
            .collect::<Vec<_>>();
        assert_eq!(exact_duplicates.len(), 2);
        let first_duplicate = exact_duplicates
            .first()
            .expect("the duplicate count was checked");
        let second_duplicate = exact_duplicates
            .get(1)
            .expect("the duplicate count was checked");
        assert!(
            !first_duplicate.is_same_alignment(second_duplicate),
            "otherwise identical records have distinct occurrence identities"
        );

        let longer = unique_profile(&first_fetch, 3, 34, false);
        let reverse = unique_profile(&first_fetch, 3, 33, true);
        let moved = unique_profile(&first_fetch, 4, 33, false);
        for distinct in [longer, reverse, moved] {
            assert!(
                !first_duplicate.is_same_alignment(distinct),
                "strand, start, and reference end distinguish same-name alignments"
            );
        }
        assert!(reverse.is_reverse());
        assert_eq!((reverse.align_start(), reverse.align_end()), (3, 33));
        assert_eq!((longer.align_start(), longer.align_end()), (3, 34));
        assert_eq!((moved.align_start(), moved.align_end()), (4, 33));

        for (left_index, left) in first_fetch.iter().enumerate() {
            for right in first_fetch.iter().skip(left_index.saturating_add(1)) {
                assert!(
                    !left.is_same_alignment(right),
                    "every asymmetric fixture record has a distinct alignment identity"
                );
            }
        }

        let second_fetch = reader.profiles(0, 11, 19, ModChar::new('a'), win)?;
        assert_eq!(second_fetch.len(), first_fetch.len());
        for original in &first_fetch {
            let matching_refetches = second_fetch
                .iter()
                .filter(|refetched| original.is_same_alignment(refetched))
                .count();
            assert_eq!(
                matching_refetches, 1,
                "each alignment, including each duplicate occurrence, matches once after refetch"
            );
        }
        for refetched in &second_fetch {
            assert_eq!(
                first_fetch
                    .iter()
                    .filter(|original| refetched.is_same_alignment(original))
                    .count(),
                1,
                "refetch matching is a one-to-one correspondence"
            );
        }
        Ok(())
    }
}
