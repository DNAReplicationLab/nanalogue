#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Coordinate retrieval checks at the public BAM record boundary.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::{CurrRead, Error, GenomicBed3};
    use rust_htslib::bam::{
        Header, HeaderView, Record,
        header::HeaderRecord,
        record::{Cigar, CigarString},
    };
    use std::sync::Arc;

    /// Construct small in-memory records without depending on fixture files.
    fn record(cigar: Vec<Cigar>, seq: &[u8]) -> Record {
        let mut header = Header::new();
        let _header = header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 100),
        );
        let mut result = Record::new();
        result.set_header(Arc::new(HeaderView::from_header(&header)));
        result.set(
            b"coordinate-read",
            Some(&CigarString(cigar)),
            seq,
            &vec![30; seq.len()],
        );
        result.set_flags(0);
        result.set_tid(0);
        result.set_pos(10);
        result
    }

    #[test]
    fn rejects_more_query_coordinates_than_stored_bases() {
        // Exactly three bases succeed; dropping one stored base must not let
        // the third CIGAR coordinate escape as a usable sequence index.
        let complete = record(vec![Cigar::Match(3)], b"ACG");
        let short = record(vec![Cigar::Match(3)], b"AC");
        let region = GenomicBed3::new(0, 10, 13);
        let read = CurrRead::default()
            .try_from_only_alignment(&complete)
            .unwrap();
        assert_eq!(
            read.seq_coords_from_ref_coords(&complete, &region).unwrap(),
            vec![Some((true, 0)), Some((true, 1)), Some((true, 2))]
        );
        let short_read = CurrRead::default().try_from_only_alignment(&short).unwrap();
        let error = short_read
            .seq_coords_from_ref_coords(&short, &region)
            .unwrap_err();
        assert!(matches!(error, Error::InvalidState(message)
            if message == "incorrect number of coordinates received!"));
    }

    #[test]
    fn detects_a_record_that_does_not_cover_the_stored_interval() {
        // The state describes [10,15), but the second record ends at 13.
        // Both records are individually valid. Pairing the wrong record with
        // cached state must report missing coordinates, not a partial success.
        let full = record(vec![Cigar::Match(5)], b"ACGTA");
        let shorter = record(vec![Cigar::Match(3)], b"ACG");
        let read = CurrRead::default().try_from_only_alignment(&full).unwrap();
        let error = read
            .seq_coords_from_ref_coords(&shorter, &GenomicBed3::new(0, 10, 15))
            .unwrap_err();
        assert!(matches!(error, Error::InvalidState(message)
            if message == "failure from upstream libraries: missing sequence coordinates"));

        // A subinterval which is fully present still works with that record.
        let shorter_read = CurrRead::default()
            .try_from_only_alignment(&shorter)
            .unwrap();
        assert_eq!(
            shorter_read
                .seq_coords_from_ref_coords(&shorter, &GenomicBed3::new(0, 11, 13))
                .unwrap(),
            vec![Some((true, 1)), Some((true, 2))]
        );
    }

    #[test]
    fn keeps_internal_insertions_but_trims_boundary_insertions() {
        let aligned = record(
            vec![
                Cigar::SoftClip(1),
                Cigar::Match(2),
                Cigar::Ins(2),
                Cigar::Del(1),
                Cigar::Match(1),
                Cigar::SoftClip(1),
            ],
            b"TACGGAT",
        );
        let read = CurrRead::default()
            .try_from_only_alignment(&aligned)
            .unwrap();
        assert_eq!(
            read.seq_coords_from_ref_coords(&aligned, &GenomicBed3::new(0, 10, 14))
                .unwrap(),
            vec![
                Some((true, 1)),
                Some((true, 2)),
                Some((false, 3)),
                Some((false, 4)),
                None,
                Some((true, 5)),
            ]
        );
        // The insertion is now at the end of the requested interval; it is
        // excluded rather than attached to the preceding matched bases.
        assert_eq!(
            read.seq_coords_from_ref_coords(&aligned, &GenomicBed3::new(0, 10, 12))
                .unwrap(),
            vec![Some((true, 1)), Some((true, 2))]
        );
        assert_eq!(
            read.seq_coords_from_ref_coords(&aligned, &GenomicBed3::new(0, 12, 13))
                .unwrap(),
            vec![None]
        );
        // Ref coordinates are generated monotonically by aligned_pairs_full:
        // each ref-consuming operation advances exactly once. Consequently
        // ref_coord_count cannot exceed this half-open interval's length.
        // A deletion-only interval exercises that invariant without query bases.
    }

    #[test]
    fn partial_alignment_state_omits_unset_lengths_from_json() {
        let aligned = record(vec![Cigar::Match(3)], b"ACG");
        let state = CurrRead::default()
            .set_read_state_and_id(&aligned)
            .unwrap()
            .set_contig_id_and_start(&aligned)
            .unwrap();
        let initial: serde_json::Value = serde_json::from_str(&state.to_string()).unwrap();
        assert_eq!(initial.get("read_id").unwrap(), "coordinate-read");
        assert_eq!(initial.get("reference_start").unwrap(), 10);
        assert_eq!(initial.get("contig").unwrap(), "0");
        assert!(initial.get("sequence_length").is_none());
        assert!(initial.get("reference_end").is_none());
        assert!(initial.get("alignment_length").is_none());

        let complete = state
            .set_seq_len(&aligned)
            .unwrap()
            .set_align_len(&aligned)
            .unwrap();
        let final_value: serde_json::Value = serde_json::from_str(&complete.to_string()).unwrap();
        assert_eq!(final_value.get("sequence_length").unwrap(), 3);
        assert_eq!(final_value.get("reference_end").unwrap(), 13);
        assert_eq!(final_value.get("alignment_length").unwrap(), 3);
    }
}
