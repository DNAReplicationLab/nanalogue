#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Rejection tests for unsupported BAM flags and malformed MM/ML tags.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::read_utils::AlignAndModData;
    use nanalogue_core::{
        BaseMods, CurrRead, Error, FiberAnnotation, ReadState, nanalogue_mm_ml_parser,
    };
    use rust_htslib::bam::Record;
    use rust_htslib::bam::record::{Aux, Cigar, CigarString};
    use std::rc::Rc;

    /// Read id shared by all in-memory records; rejection messages must name it.
    const READ_ID: &str = "flag_and_tag_read";

    /// Builds a forward, mapped record with a four-base sequence and the given flags.
    fn mapped_record(flags: u16) -> Record {
        let mut record = Record::new();
        let cigar = CigarString::from(vec![Cigar::Match(4)]);
        record.set(READ_ID.as_bytes(), Some(&cigar), b"ACGT", &[30u8; 4]);
        record.set_flags(flags);
        record.set_tid(0);
        record.set_pos(17);
        record
    }

    /// Builds an unmapped record, which skips the alignment-field lookups.
    fn unmapped_record() -> Record {
        let mut record = mapped_record(4);
        record.set_tid(-1);
        record.set_pos(-1);
        record
    }

    /// Builds a mapped record carrying an MM tag and ML probabilities.
    fn tagged_record(mm: &str, ml: &[u8]) -> Record {
        let mut record = mapped_record(0);
        record
            .push_aux(b"MM", Aux::String(mm))
            .expect("MM tag fits in the record");
        record
            .push_aux(b"ML", Aux::ArrayU8(ml.into()))
            .expect("ML tag fits in the record");
        record
    }

    /// Parses MM/ML data without filters so malformed tags cannot hide behind a filter.
    fn parse(record: &Record) -> Result<BaseMods, Error> {
        nanalogue_mm_ml_parser(record, |_| true, |_| true, |_, _, _| true, 0)
    }

    /// Consumes an error and returns its `InvalidState` payload, if that is the variant.
    fn invalid_state_message(error: Error) -> Option<String> {
        if let Error::InvalidState(message) = error {
            Some(message)
        } else {
            None
        }
    }

    /// Asserts that `error` is the unsupported-flag rejection naming [`READ_ID`].
    fn assert_unsupported_flag(error: Error, context: &str) {
        let rendered = format!("{error:?}");
        assert!(
            matches!(error, Error::NotImplemented(message) if message.contains(READ_ID)),
            "{context} must be rejected as NotImplemented naming {READ_ID}, got {rendered}"
        );
    }

    #[test]
    fn unsupported_flags_are_rejected_with_read_id() {
        // Each bit marks a read shape the library explicitly does not support.
        // Testing them one by one also checks that no flag is forgotten in the
        // short-circuit condition.
        for flag in [0x1u16, 0x2, 0x40, 0x80, 0x20, 0x8, 0x400, 0x200] {
            let record = mapped_record(flag);
            let error = CurrRead::default()
                .set_read_state_and_id(&record)
                .expect_err("unsupported flags must be rejected");
            assert_unsupported_flag(error, &format!("flag {flag:#x}"));
        }
    }

    #[test]
    fn supported_flags_still_parse() {
        // Controls: plain, unmapped, and reverse reads are all supported.
        for (flag, expected) in [
            (0u16, ReadState::PrimaryFwd),
            (4, ReadState::Unmapped),
            (16, ReadState::PrimaryRev),
        ] {
            let record = mapped_record(flag);
            let curr_read = CurrRead::default()
                .set_read_state_and_id(&record)
                .expect("flags 0, 4, and 16 are supported");
            assert_eq!(curr_read.read_state(), expected, "flag {flag:#x}");
        }
    }

    #[test]
    fn rc_record_conversion_rejects_unsupported_flags() {
        // The Rc<Record> conversion must short-circuit on the same flags.
        let result: Result<CurrRead<AlignAndModData>, Error> =
            CurrRead::try_from(Rc::new(mapped_record(0x1)));
        assert_unsupported_flag(
            result.expect_err("paired reads must be rejected through Rc<Record>"),
            "Rc<Record> with flag 0x1",
        );

        // Control: an unmapped record without mod tags is accepted through the same path.
        let curr_read: CurrRead<AlignAndModData> = CurrRead::try_from(Rc::new(unmapped_record()))
            .expect("unmapped records without MM tags are supported");
        assert_eq!(
            curr_read.read_state(),
            ReadState::Unmapped,
            "unmapped control"
        );
        assert_eq!(curr_read.read_id(), READ_ID, "read id must be preserved");
        assert!(
            curr_read.base_count_per_mod().is_empty(),
            "a record without MM tags must have no modification counts"
        );
    }

    #[test]
    fn duplicate_modern_and_legacy_tags_are_rejected_distinctly() {
        let mut both_probability_tags = mapped_record(0);
        both_probability_tags
            .push_aux(b"ML", Aux::ArrayU8((&[101u8][..]).into()))
            .expect("modern ML tag");
        both_probability_tags
            .push_aux(b"Ml", Aux::ArrayU8((&[101u8][..]).into()))
            .expect("legacy Ml tag");
        let probability_message = invalid_state_message(
            parse(&both_probability_tags).expect_err("both ML and Ml must be rejected"),
        )
        .expect("duplicate probability tags must be InvalidState");

        let mut both_position_tags = tagged_record("C+m,0;", &[101]);
        both_position_tags
            .push_aux(b"Mm", Aux::String("C+m,0;"))
            .expect("legacy Mm tag");
        let position_message = invalid_state_message(
            parse(&both_position_tags).expect_err("both MM and Mm must be rejected"),
        )
        .expect("duplicate position tags must be InvalidState");

        assert!(
            probability_message.contains("`ML`") && probability_message.contains("`Ml`"),
            "duplicate probabilities must name both ML spellings: {probability_message}"
        );
        assert!(
            position_message.contains("`MM`") && position_message.contains("`Mm`"),
            "duplicate positions must name both MM spellings: {position_message}"
        );
        assert_ne!(
            probability_message, position_message,
            "the two ambiguous-tag errors must be distinguishable"
        );
    }

    #[test]
    fn unexpected_ml_and_mm_tag_types_are_rejected_distinctly() {
        let mut string_ml = mapped_record(0);
        string_ml
            .push_aux(b"ML", Aux::String("101"))
            .expect("string-valued ML tag");
        let string_ml_message =
            invalid_state_message(parse(&string_ml).expect_err("ML must be an integer array"))
                .expect("a string ML tag must be InvalidState");

        let mut integer_mm = mapped_record(0);
        integer_mm
            .push_aux(b"MM", Aux::I32(1))
            .expect("integer-valued MM tag");
        let integer_mm_message =
            invalid_state_message(parse(&integer_mm).expect_err("MM must be a string"))
                .expect("an integer MM tag must be InvalidState");

        assert_eq!(
            string_ml_message,
            "rust-htslib ML/Ml tag parsing failure: unexpected auxiliary tag type",
            "a wrong-typed ML tag must report the ML/Ml aux-tag failure"
        );
        assert_eq!(
            integer_mm_message,
            "rust-htslib MM/Mm tag parsing failure: unexpected auxiliary tag type",
            "a wrong-typed MM tag must report the MM/Mm aux-tag failure"
        );
    }

    #[test]
    fn legacy_tag_spellings_are_accepted() {
        // All four valid spellings parse, but mixing spellings of one tag kind
        // (checked above) does not.
        for (mm_tag, ml_tag) in [("MM", "ML"), ("MM", "Ml"), ("Mm", "ML"), ("Mm", "Ml")] {
            let mut record = mapped_record(0);
            record
                .push_aux(mm_tag.as_bytes(), Aux::String("C+m,0;"))
                .expect("valid MM tag");
            record
                .push_aux(ml_tag.as_bytes(), Aux::ArrayU8((&[203u8][..]).into()))
                .expect("valid ML tag");
            let mods = parse(&record).expect("legacy spellings are accepted");
            let group = mods.base_mods.first().expect("one C+m group");
            assert_eq!(group.modified_base, b'C', "{mm_tag}/{ml_tag}");
            assert_eq!(
                group.ranges.annotations,
                vec![FiberAnnotation {
                    pos: 1,
                    qual: 203,
                    ref_pos: Some(18)
                }],
                "{mm_tag}/{ml_tag}"
            );
        }
    }

    #[test]
    fn invalid_mm_base_and_strand_are_rejected() {
        let error = parse(&tagged_record("X+m,0;", &[101])).expect_err("X is not a valid MM base");
        let rendered = format!("{error:?}");
        assert!(
            matches!(error, Error::InvalidBase(message) if message == "invalid MM base `X`"),
            "an invalid base letter must be InvalidBase, got {rendered}"
        );

        let strand_error =
            parse(&tagged_record("C*m,0;", &[101])).expect_err("* is not a valid MM strand");
        let strand_rendered = format!("{strand_error:?}");
        assert!(
            matches!(strand_error, Error::InvalidModType(message)
                if message.starts_with("invalid MM strand")),
            "an invalid strand byte must be InvalidModType, got {strand_rendered}"
        );
    }

    #[test]
    fn n_base_annotates_all_positions_but_c_base_only_cytosines() {
        // N matches every base, so two zero-distance calls annotate positions 0
        // and 1 even though position 0 holds an A.
        let n_mods = parse(&tagged_record("N+m?,0,0;", &[101, 202])).expect("N group parses");
        let n_group = n_mods.base_mods.first().expect("one N group");
        assert_eq!(
            n_group.modified_base, b'N',
            "the group must record the N base"
        );
        assert_eq!(
            n_group.ranges.annotations,
            vec![
                FiberAnnotation {
                    pos: 0,
                    qual: 101,
                    ref_pos: Some(17)
                },
                FiberAnnotation {
                    pos: 1,
                    qual: 202,
                    ref_pos: Some(18)
                },
            ],
            "N must annotate positions 0 and 1 regardless of the underlying bases"
        );

        // C matches only the single cytosine at position 1, so one call lands there.
        let c_mods = parse(&tagged_record("C+m?,0;", &[203])).expect("C group parses");
        let c_group = c_mods.base_mods.first().expect("one C group");
        assert_eq!(
            c_group.modified_base, b'C',
            "the group must record the C base"
        );
        assert_eq!(
            c_group.ranges.annotations,
            vec![FiberAnnotation {
                pos: 1,
                qual: 203,
                ref_pos: Some(18)
            }],
            "C must skip the A at position 0 and annotate only the cytosine"
        );
    }
}
