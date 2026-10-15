#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Record-shape validation at the public MM/ML parsing boundary.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::{BaseMods, Error, FiberAnnotation, nanalogue_mm_ml_parser};
    use rust_htslib::bam::{
        Record,
        record::{Aux, Cigar, CigarString},
    };

    /// Parses without filtering so malformed shapes cannot hide behind a filter.
    fn parse(record: &Record) -> Result<BaseMods, Error> {
        nanalogue_mm_ml_parser(record, |_| true, |_| true, |_, _, _| true, 0)
    }

    /// Builds a mapped record through the public API, including malformed CIGARs.
    fn record(sequence: &[u8], cigar: &[Cigar]) -> Record {
        let mut result = Record::new();
        result.set(
            b"shape_read",
            Some(&CigarString(cigar.to_vec())),
            sequence,
            &vec![30; sequence.len()],
        );
        result.set_flags(0);
        result.set_tid(0);
        result.set_pos(17);
        result
    }

    #[test]
    fn empty_sequence_needs_no_shape_validation_without_mm_tag() {
        let mut empty = record(b"", &[]);
        empty.set_flags(4);
        empty.set_tid(-1);
        empty.set_pos(-1);
        assert_eq!(empty.seq_len(), 0);
        assert!(
            parse(&empty)
                .expect("no modifications")
                .base_mods
                .is_empty(),
            "a sequence-less record with no MM data has no calls"
        );
    }

    #[test]
    fn mn_is_ignored_without_mm_tag() {
        for mn_tag in [Aux::I32(14), Aux::String("stale")] {
            let mut unmodified = record(b"ACGT", &[Cigar::Match(4)]);
            unmodified.push_aux(b"MN", mn_tag).expect("MN tag");
            assert!(
                parse(&unmodified)
                    .expect("MN is irrelevant without an MM tag")
                    .base_mods
                    .is_empty(),
                "absent MM must yield no modifications"
            );
        }
    }

    #[test]
    fn empty_mm_tag_checks_mn() {
        for tag in [b"MM", b"Mm"] {
            for (mn_tag, expected_error) in [
                (None, None),
                (Some(Aux::I32(4)), None),
                (
                    Some(Aux::I32(14)),
                    Some("MN tag sequence length 14 does not match current sequence length 4"),
                ),
                (
                    Some(Aux::String("stale")),
                    Some("MN tag must contain an integer sequence length"),
                ),
            ] {
                let mut unmodified = record(b"ACGT", &[Cigar::Match(4)]);
                unmodified
                    .push_aux(tag, Aux::String(""))
                    .expect("empty MM tag");
                if let Some(value) = mn_tag {
                    unmodified.push_aux(b"MN", value).expect("MN tag");
                }
                let parsed = parse(&unmodified);
                if let Some(expected) = expected_error {
                    assert!(
                        matches!(parsed, Err(Error::InvalidModCoords(message)) if message == expected),
                        "empty MM must reject invalid MN"
                    );
                } else {
                    assert!(
                        parsed.expect("valid empty MM").base_mods.is_empty(),
                        "valid empty MM must yield no modifications"
                    );
                }
            }
        }
    }

    #[test]
    fn nonempty_mm_group_without_calls_still_checks_mn() {
        for tag in [b"MM", b"Mm"] {
            let mut modified = record(b"ACGT", &[Cigar::Match(4)]);
            modified
                .push_aux(tag, Aux::String("C+m?;"))
                .expect("MM group without calls");
            modified.push_aux(b"MN", Aux::I32(14)).expect("MN tag");
            let error = parse(&modified).expect_err("nonempty MM data requires a matching MN");
            assert!(
                matches!(error, Error::InvalidModCoords(message)
                if message == "MN tag sequence length 14 does not match current sequence length 4"),
                "MN validation must not depend on ML probabilities being present"
            );
        }
    }

    #[test]
    fn invalid_mn_precedes_mm_calls_exceeding_sequence() {
        for (mm_tag, ml_tag) in [(b"MM", b"ML"), (b"Mm", b"Ml")] {
            let mut malformed = record(b"CCGT", &[Cigar::Match(4)]);
            // Three calls cannot fit on the sequence's two cytosines.
            malformed
                .push_aux(mm_tag, Aux::String("C+m?,0,0,0;"))
                .expect("MM tag");
            malformed
                .push_aux(ml_tag, Aux::ArrayU8((&[51u8, 102, 204][..]).into()))
                .expect("ML tag");
            malformed.push_aux(b"MN", Aux::I32(14)).expect("stale MN");

            let length_error = parse(&malformed).expect_err("both MN and MM data are invalid");
            assert!(
                matches!(length_error, Error::InvalidModCoords(message)
                if message == "MN tag sequence length 14 does not match current sequence length 4"),
                "MN mismatch must be reported before processing MM calls"
            );

            malformed.remove_aux(b"MN").expect("remove stale MN");
            malformed.push_aux(b"MN", Aux::I32(4)).expect("matching MN");
            let calls_error = parse(&malformed).expect_err("MM calls still exceed the sequence");
            assert!(
                matches!(calls_error, Error::InvalidModCoords(message)
                if message == "Problem with parsing MM/ML data, counts do not match 2 != 3"),
                "correcting only MN must expose the independently invalid MM data"
            );
        }
    }

    #[test]
    fn present_mm_tag_rejects_missing_sequence_even_without_calls() {
        // BAM records can omit SEQ, but a present MM tag requires a sequence,
        // even if it is empty or has no calls and therefore no ML probabilities.
        for tag in [b"MM", b"Mm"] {
            for mm_text in ["", "C+m?;", "C+m.;", "C+m;"] {
                for mn in [None, Some(0), Some(14)] {
                    let mut empty = record(b"", &[]);
                    empty.set_flags(4);
                    empty.set_tid(-1);
                    empty.set_pos(-1);
                    empty
                        .push_aux(tag, Aux::String(mm_text))
                        .expect("syntactically valid MM tag");
                    if let Some(value) = mn {
                        empty.push_aux(b"MN", Aux::I32(value)).expect("MN tag");
                    }
                    let error = parse(&empty).expect_err("MM data requires a sequence");
                    assert!(
                        matches!(error, Error::ZeroSeqLen(message)
                        if message == "zero length sequences cannot be used for mod tag parsing"),
                        "missing sequence must be rejected before MN validation: MM: {mm_text:?}, MN: {mn:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn mapped_cigar_must_cover_exactly_the_query_sequence() {
        // Record::set permits inconsistent query consumption. The parser must
        // reject both a truncated map and an overlong map before indexing it.
        for cigar in [
            vec![Cigar::Match(2)],
            vec![Cigar::Match(4)],
            vec![Cigar::SoftClip(1), Cigar::Match(1)],
            vec![Cigar::Match(2), Cigar::Ins(2)],
        ] {
            for (mm_text, probabilities) in [("", &[][..]), ("C+m?,0,0,0;", &[51u8, 102, 204][..])]
            {
                for (mm_tag, ml_tag) in [(b"MM", b"ML"), (b"Mm", b"Ml")] {
                    for mn in [None, Some(3), Some(14)] {
                        let mut malformed = record(b"CCC", &cigar);
                        malformed
                            .push_aux(mm_tag, Aux::String(mm_text))
                            .expect("MM tag");
                        malformed
                            .push_aux(ml_tag, Aux::ArrayU8(probabilities.into()))
                            .expect("matching probabilities");
                        if let Some(value) = mn {
                            malformed.push_aux(b"MN", Aux::I32(value)).expect("MN tag");
                        }
                        let error = parse(&malformed).expect_err("invalid record");
                        if mn == Some(14) {
                            assert!(
                                matches!(error, Error::InvalidModCoords(message)
                                if message == "MN tag sequence length 14 does not match current sequence length 3"),
                                "invalid MN must be rejected before CIGAR validation"
                            );
                        } else {
                            assert!(
                                matches!(error, Error::InvalidState(message)
                                if message == "rust_htslib failure! seq coordinates malformed"),
                                "malformed CIGAR must reach the shape check: {cigar:?}, MM: {mm_text:?}, MN: {mn:?}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn exact_query_consumption_preserves_clips_insertions_and_reference_gaps() {
        for (cigar, expected_refs) in [
            (vec![Cigar::Match(3)], [Some(17), Some(18), Some(19)]),
            (
                vec![Cigar::SoftClip(1), Cigar::Match(2)],
                [None, Some(17), Some(18)],
            ),
            (
                vec![Cigar::Match(1), Cigar::Ins(1), Cigar::Match(1)],
                [Some(17), None, Some(18)],
            ),
            (
                vec![Cigar::Match(1), Cigar::Del(5), Cigar::Match(2)],
                [Some(17), Some(23), Some(24)],
            ),
        ] {
            let mut valid = record(b"CCC", &cigar);
            valid
                .push_aux(b"MM", Aux::String("C+m?,0,0,0;"))
                .expect("three calls");
            valid
                .push_aux(b"ML", Aux::ArrayU8((&[51u8, 102, 204][..]).into()))
                .expect("three probabilities");
            let parsed = parse(&valid).expect("exact query consumption is valid");
            assert_eq!(parsed.base_mods.len(), 1);
            let group = parsed.base_mods.first().expect("one modification group");
            let expected: Vec<_> = [(0, 51), (1, 102), (2, 204)]
                .into_iter()
                .zip(expected_refs)
                .map(|((pos, qual), ref_pos)| FiberAnnotation::try_new(pos, qual, ref_pos).unwrap())
                .collect();
            assert_eq!(group.ranges.seq_len(), 3);
            assert_eq!(
                group.ranges.annotations(),
                expected,
                "incorrect map for {cigar:?}"
            );
        }
        // Quality/sequence length disagreement cannot be constructed with
        // Record::set: rust-htslib enforces equal lengths. Do not mutate raw
        // record storage to reach the parser's defensive quality-length guard.
    }
}
