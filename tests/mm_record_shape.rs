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
    fn empty_sequence_needs_no_shape_validation_without_modifications() {
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
        empty
            .push_aux(b"MM", Aux::String(""))
            .expect("empty MM tag");
        assert!(
            parse(&empty).expect("empty MM data").base_mods.is_empty(),
            "an empty tag must behave like an absent tag"
        );
    }

    #[test]
    fn nonempty_mm_tag_rejects_missing_sequence_even_without_calls() {
        // BAM records can omit SEQ. A nonempty MM group still needs a sequence,
        // even if it has no distances and therefore no ML probabilities.
        for tag in [b"MM", b"Mm"] {
            for suffix in ["?", ".", ""] {
                let mut empty = record(b"", &[]);
                empty.set_flags(4);
                empty.set_tid(-1);
                empty.set_pos(-1);
                empty
                    .push_aux(tag, Aux::String(&format!("C+m{suffix};")))
                    .expect("syntactically valid empty group");
                let error = parse(&empty).expect_err("MM data requires a sequence");
                assert!(
                    matches!(error, Error::ZeroSeqLen(message)
                    if message == "zero length sequences cannot be used for mod tag parsing"),
                    "missing sequence must not produce a silently empty result"
                );
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
            let mut malformed = record(b"CCC", &cigar);
            malformed
                .push_aux(b"MM", Aux::String("C+m?,0,0,0;"))
                .expect("three calls");
            malformed
                .push_aux(b"ML", Aux::ArrayU8((&[51u8, 102, 204][..]).into()))
                .expect("three probabilities");
            let error = parse(&malformed).expect_err("query consumption mismatch");
            assert!(
                matches!(error, Error::InvalidState(message)
                if message == "rust_htslib failure! seq coordinates malformed"),
                "malformed CIGAR must reach the shape check: {cigar:?}"
            );
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
                .map(|((pos, qual), ref_pos)| FiberAnnotation { pos, qual, ref_pos })
                .collect();
            assert_eq!(group.ranges.seq_len, 3);
            assert_eq!(
                group.ranges.annotations, expected,
                "incorrect map for {cigar:?}"
            );
        }
        // Quality/sequence length disagreement cannot be constructed with
        // Record::set: rust-htslib enforces equal lengths. Do not mutate raw
        // record storage to reach the parser's defensive quality-length guard.
    }
}
