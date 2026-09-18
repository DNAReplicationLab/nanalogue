//! Boundary coverage for BAM record pre-filters without file-backed fixtures.

#[cfg(test)]
mod tests {
    use nanalogue_core::{BamPreFilt as _, GenomicBed3, InputBam};
    use rust_htslib::bam::record::{Cigar, CigarString, Record};
    use std::collections::HashSet;

    /// A mapped record with a nonzero reference ID and a ten-base reference span.
    fn mapped_record() -> Record {
        let mut record = Record::new();
        record.set(
            b"read-7",
            Some(&CigarString(vec![Cigar::Match(10)])),
            b"ACGTACGTAC",
            &[30; 10],
        );
        record.set_flags(0);
        record.set_tid(2);
        record.set_pos(17);
        record.set_mapq(40);
        record
    }

    #[test]
    fn distinguishes_unavailable_mapq_from_a_high_quality_alignment() {
        let mut record = mapped_record();
        // A numeric threshold alone admits 255; the explicit exclusion must
        // reject it even when the minimum is zero or exactly 255.
        for (mapq, minimum, exclude_unavailable, expected) in [
            (39, 40, false, false),
            (40, 40, false, true),
            (41, 40, false, true),
            (39, 40, true, false),
            (40, 40, true, true),
            (254, 254, true, true),
            (254, 255, false, false),
            (255, 0, false, true),
            (255, 255, false, true),
            (255, 0, true, false),
            (255, 255, true, false),
        ] {
            record.set_mapq(mapq);
            let mut options = InputBam::default();
            options.mapq_filter = minimum;
            options.exclude_mapq_unavail = exclude_unavailable;
            assert_eq!(
                record.filt_by_mapq(minimum, exclude_unavailable),
                expected,
                "mapq={mapq}, minimum={minimum}, exclude={exclude_unavailable}"
            );
            assert_eq!(
                record.pre_filt(&options),
                expected,
                "the combined pre-filter must apply mapping quality"
            );
        }
    }

    #[test]
    fn alignment_length_requires_mapping_and_a_nonnegative_start() {
        let mut record = mapped_record();
        // Reference span, not query length: the deletion adds seven bases,
        // while the insertion consumes two query bases and no reference bases.
        record.set(
            b"read-7",
            Some(&CigarString(vec![
                Cigar::Match(3),
                Cigar::Ins(2),
                Cigar::Del(7),
                Cigar::Match(5),
            ])),
            b"ACGTACGTAC",
            &[30; 10],
        );
        assert!(record.filt_by_align_len(15), "exact span is included");
        assert!(!record.filt_by_align_len(16), "one above span is excluded");

        record.set_unmapped();
        assert!(!record.filt_by_align_len(0), "unmapped is not an alignment");
        record.unset_unmapped();
        record.set_pos(-1);
        assert!(
            !record.filt_by_align_len(0),
            "a positive span cannot legitimize a negative start"
        );
        record.set_pos(0);
        assert!(record.filt_by_align_len(15), "reference origin is valid");
    }

    #[test]
    fn read_id_set_rejects_non_utf8_without_lossy_matching() {
        let mut record = mapped_record();
        let ids = HashSet::from(["read-7".to_owned(), "read-\u{fffd}".to_owned()]);
        assert!(record.filt_by_read_id_set(&ids), "exact ID is present");
        record.set_qname(b"read-8");
        assert!(!record.filt_by_read_id_set(&ids), "different ID is absent");

        record.set_qname(b"read-\xff");
        assert_eq!(record.qname(), b"read-\xff");
        assert!(
            !record.filt_by_read_id_set(&ids),
            "invalid UTF-8 must not match a Unicode replacement character"
        );
        record.set_qname("read-\u{fffd}".as_bytes());
        assert!(
            record.filt_by_read_id_set(&ids),
            "an actual replacement character remains a valid exact ID"
        );
    }

    #[test]
    fn region_filter_rejects_reference_end_overflow_without_wrapping() {
        let mut record = mapped_record();
        let region = GenomicBed3::new(2, u32::MAX - 5, u32::MAX);
        // Ten matched bases end exactly at u32::MAX, which is representable.
        record.set_pos(i64::from(u32::MAX) - 10);
        for full_overlap in [false, true] {
            assert!(
                record.filt_by_region(&region, full_overlap),
                "the maximum representable end is accepted"
            );
        }
        // The start still fits u32, but the end does not. This distinguishes
        // the end conversion guard from the already-covered start guard.
        record.set_pos(i64::from(u32::MAX) - 9);
        for full_overlap in [false, true] {
            assert!(
                !record.filt_by_region(&region, full_overlap),
                "an unrepresentable end is rejected rather than truncated"
            );
        }
    }
}
