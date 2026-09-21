#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Interactions between public MM/ML filters, reverse coordinates, and parser errors.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use std::sync::Arc;

    use nanalogue_core::{
        BaseMods, CurrRead, Error, FiberAnnotation, ModChar, ThresholdState, nanalogue_mm_ml_parser,
    };
    use rust_htslib::bam::{
        Header, HeaderView, Record,
        header::HeaderRecord,
        record::{Aux, Cigar, CigarString},
    };

    /// Builds a mapped record whose CIGAR has soft-clipped and inserted query bases.
    fn complex_record(flags: u16, sequence: &[u8], qualities: &[u8]) -> Record {
        let mut raw_header = Header::new();
        let header = raw_header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 1_000),
        );
        let mut record = Record::new();
        let cigar = CigarString::from(vec![
            Cigar::SoftClip(1),
            Cigar::Match(2),
            Cigar::Ins(1),
            Cigar::Match(2),
        ]);
        record.set(b"filter_coordinate_read", Some(&cigar), sequence, qualities);
        record.set_header(Arc::new(HeaderView::from_header(header)));
        record.set_flags(flags);
        record.set_tid(0);
        record.set_pos(100);
        record
    }

    /// Adds one MM group and its ML values through rust-htslib's public record API.
    fn add_mod_tags(record: &mut Record, mm: &str, ml: &[u8]) {
        record
            .push_aux(b"MM", Aux::String(mm))
            .expect("MM text fits in the record");
        record
            .push_aux(b"ML", Aux::ArrayU8(ml.into()))
            .expect("ML values fit in the record");
    }

    /// Builds the parser's optimized single-operation alignment shape.
    fn contiguous_record(position: i64, sequence: &[u8]) -> Record {
        let mut record = Record::new();
        let cigar = CigarString::from(vec![Cigar::Match(
            u32::try_from(sequence.len()).expect("tiny fixture length fits u32"),
        )]);
        record.set(
            b"contiguous_coordinate_read",
            Some(&cigar),
            sequence,
            &vec![30; sequence.len()],
        );
        record.set_flags(0);
        record.set_tid(0);
        record.set_pos(position);
        record
    }

    #[test]
    fn reverse_fallback_mapping_applies_all_filters_in_forward_space() -> Result<(), Error> {
        // For a reverse record the qualities are traversed as
        // [255, 50, 30, 40, 20, 50]. Across the six N calls, index 0 tests
        // unavailable quality, 5 probability rejection, 4 low quality, 3
        // position rejection, 2 inclusive quality acceptance, and 1 a mapped
        // reference coordinate. Incorrect quality orientation changes the result.
        let mut record = complex_record(16, b"ACGTAC", &[50, 20, 40, 30, 50, 255]);
        add_mod_tags(
            &mut record,
            "N+n?,0,0,0,0,0,0;",
            &[107, 104, 103, 105, 106, 102],
        );

        let parsed = nanalogue_mm_ml_parser(
            &record,
            |probability| *probability >= 103,
            |forward_position| *forward_position != 3,
            |base, strand, tag| *base == b'N' && *strand == '+' && *tag == ModChar::new('n'),
            30,
        )?;

        let group = parsed.base_mods.first().expect("one accepted N+n group");
        assert_eq!(parsed.base_mods.len(), 1);
        assert_eq!(group.modified_base, b'N');
        assert_eq!(group.strand, '+');
        assert_eq!(group.modification_type, 'n');
        assert!(group.record_is_reverse);
        assert!(group.ranges.reverse);
        assert_eq!(group.ranges.seq_len, 6);
        assert_eq!(
            group.ranges.annotations,
            [
                FiberAnnotation {
                    // Forward index 2 becomes stored reverse-read coordinate 3.
                    pos: 3,
                    qual: 103,
                    // The insertion is third in MM's original forward orientation.
                    ref_pos: None,
                },
                FiberAnnotation {
                    pos: 4,
                    qual: 104,
                    ref_pos: Some(102),
                },
            ]
        );

        // Rejecting a valid group still consumes its calls, so it must not
        // cause an ML-length mismatch.
        let rejected = nanalogue_mm_ml_parser(&record, |_| true, |_| true, |_, _, _| false, 0)?;
        assert!(rejected.base_mods.is_empty());
        Ok(())
    }

    #[test]
    fn implicit_calls_respect_the_zero_probability_filter() -> Result<(), Error> {
        let mut record = complex_record(0, b"NNNNNN", &[30; 6]);
        add_mod_tags(&mut record, "N+n.,1,1;", &[80, 200]);

        let BaseMods { base_mods } = nanalogue_mm_ml_parser(
            &record,
            |probability| *probability >= 100,
            |_| true,
            |_, _, _| true,
            0,
        )?;
        assert_eq!(
            base_mods.first().expect("one N+n group").ranges.annotations,
            [FiberAnnotation {
                pos: 3,
                qual: 200,
                ref_pos: None,
            }],
            "zero-valued implicit calls and the low-probability explicit call are filtered"
        );

        let with_zeros = nanalogue_mm_ml_parser(&record, |_| true, |_| true, |_, _, _| true, 0)?;
        let annotations = &with_zeros
            .base_mods
            .first()
            .expect("one N+n group")
            .ranges
            .annotations;
        assert_eq!(
            annotations
                .iter()
                .map(|annotation| (annotation.pos, annotation.qual, annotation.ref_pos))
                .collect::<Vec<_>>(),
            [
                (0, 0, None),
                (1, 80, Some(100)),
                (2, 0, Some(101)),
                (3, 200, None),
                (4, 0, Some(102)),
                (5, 0, Some(103)),
            ]
        );
        Ok(())
    }

    #[test]
    fn curr_read_setters_preserve_insufficient_ml_errors() {
        let mut record = complex_record(0, b"NNNNNN", &[30; 6]);
        record
            .push_aux(b"MM", Aux::String("N+n?,0;"))
            .expect("MM text fits in the record");

        let error = CurrRead::default()
            .try_from_only_alignment(&record)
            .expect("alignment fields are valid")
            .set_mod_data(&record, ThresholdState::default(), 0)
            .expect_err("MM calls without ML probabilities must fail");
        assert!(matches!(error, Error::InvalidModProbs(message)
                if message == "ML tag appears to be insufficiently long!"));

        let restricted_error = CurrRead::default()
            .try_from_only_alignment(&record)
            .expect("alignment fields are valid")
            .set_mod_data_restricted(
                &record,
                ThresholdState::default(),
                |_| true,
                |_, _, _| true,
                0,
            )
            .expect_err("restricted parsing must preserve malformed-tag errors");
        assert!(matches!(restricted_error, Error::InvalidModProbs(message)
                if message == "ML tag appears to be insufficiently long!"));
    }

    #[test]
    fn contiguous_mapping_rejects_invalid_reference_boundaries() {
        let mut negative = contiguous_record(-1, b"N");
        add_mod_tags(&mut negative, "N+n?,0;", &[200]);
        let negative_error =
            nanalogue_mm_ml_parser(&negative, |_| true, |_| true, |_, _, _| true, 0)
                .expect_err("a mapped alignment cannot start before reference zero");
        assert!(matches!(negative_error, Error::InvalidModCoords(message)
                if message.starts_with("reference start coordinate is invalid:")));

        let mut overflowing = contiguous_record(i64::from(u32::MAX), b"NN");
        add_mod_tags(&mut overflowing, "N+n?,0;", &[200]);
        let overflow_error =
            nanalogue_mm_ml_parser(&overflowing, |_| true, |_| true, |_, _, _| true, 0)
                .expect_err("the final reference coordinate cannot exceed u32");
        assert!(matches!(overflow_error, Error::InvalidModCoords(message)
                if message == "reference coordinate exceeds u32 capacity"));
    }

    #[test]
    fn parser_distinguishes_coordinate_exhaustion_from_orphan_probabilities() {
        let mut exhausted = contiguous_record(10, b"A");
        add_mod_tags(&mut exhausted, "A+a?,1;", &[200]);
        let coordinate_error =
            nanalogue_mm_ml_parser(&exhausted, |_| true, |_| true, |_, _, _| true, 0)
                .expect_err("the second A requested by the distance does not exist");
        assert!(matches!(coordinate_error, Error::InvalidModCoords(message)
                if message == "Problem with parsing MM/ML data, counts do not match 0 != 1"));

        let mut orphan_probability = contiguous_record(10, b"A");
        orphan_probability
            .push_aux(b"ML", Aux::ArrayU8((&[200u8][..]).into()))
            .expect("one orphan probability fits in the record");
        let probability_error =
            nanalogue_mm_ml_parser(&orphan_probability, |_| true, |_| true, |_, _, _| true, 0)
                .expect_err("ML data without MM calls must not be ignored");
        assert!(matches!(probability_error, Error::InvalidModProbs(message)
                if message == "MM and ML tag lengths do not match!"));
    }
}
