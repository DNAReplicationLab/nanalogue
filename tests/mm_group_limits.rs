//! Boundary contracts for the public MM-group parser, without constructing BAM files.

#[cfg(test)]
mod tests {
    use std::fmt::Write as _;

    use nanalogue_core::{Error, constants::shared::MAX_MM_TAG_LENGTH, mm_groups};

    #[test]
    fn text_limit_is_inclusive_and_checked_before_group_syntax() {
        let limit = usize::try_from(MAX_MM_TAG_LENGTH).expect("supported platform");
        // No terminator: at the limit this must reach the syntax check, whereas
        // one byte beyond it must fail the length guard before scanning groups.
        // Reuse the allocation to keep the peak memory near the 100 MB limit.
        let mut text = "x".repeat(limit);
        text.push('x');
        let too_long = mm_groups(&text).expect_err("over-limit text must fail");
        assert!(
            matches!(too_long, Error::InvalidState(message)
                if message == "MM tag is too long to process"),
            "length rejection must precede malformed-group rejection"
        );
        assert_eq!(text.pop(), Some('x'));
        let at_limit = mm_groups(&text).expect_err("unterminated text must fail");
        assert!(
            matches!(at_limit, Error::InvalidModType(message)
                if message == "MM tag must end with `;` terminators for each group"),
            "exact-limit input must get past the length guard"
        );
        // A distance needs at least a comma plus a digit; the prefix and final
        // semicolon also consume bytes. Thus mod_dists.len() < text.len() for
        // every valid group, making the internal distance-count assertion's
        // failure unreachable once text.len() <= MAX_MM_TAG_LENGTH is checked.
    }

    #[test]
    fn exactly_one_hundred_distinct_types_preserve_order_and_payload() {
        let mut text = String::new();
        for code in 0..100 {
            write!(text, "C+{code}?,2,0;").expect("writing to a String cannot fail");
        }
        let groups = mm_groups(&text).expect("100 distinct types are accepted");
        assert_eq!(groups.len(), 100);
        for (expected_code, group) in (0u32..100).zip(&groups) {
            assert_eq!(group.mod_base, b'C');
            assert_eq!(group.mod_strand, '+');
            assert_eq!(u32::from(group.modification_type.val()), expected_code);
            assert!(!group.is_implicit, "question mark means explicit calls");
            assert_eq!(group.mod_dists, [2, 0]);
        }
        let over_limit = format!("{text}A+100,1;");
        let error = mm_groups(&over_limit).expect_err("101 distinct types must fail");
        assert!(
            matches!(error, Error::InvalidState(message)
                if message == "max types of mods exceeded 100"),
            "the first excess type must reach the mod-count guard"
        );
    }

    #[test]
    fn duplicate_identity_ignores_base_and_suffix() {
        // Numeric 109 and literal m encode the same modification. The base and
        // implicit/explicit suffix do not make that strand/type pair unique.
        for text in ["C+m?,2;A+109.,0;", "A+109,0;C+m?,2;"] {
            let error = mm_groups(text).expect_err("same strand/type is a duplicate");
            assert!(
                matches!(error, Error::InvalidDuplicates(_)),
                "base or suffix must not hide a duplicate: {text}"
            );
        }
        let groups = mm_groups("C+m?,2;A-109.,0;").expect("opposite strands are distinct");
        assert_eq!(groups.len(), 2);
        let first = groups.first().expect("first group");
        let second = groups.get(1).expect("second group");
        assert_eq!(first.mod_strand, '+');
        assert_eq!(second.mod_strand, '-');
        assert_eq!(first.modification_type, second.modification_type);
        assert!(!first.is_implicit, "first group is explicit");
        assert!(second.is_implicit, "second group is implicit");
        assert_eq!(first.mod_dists, [2]);
        assert_eq!(second.mod_dists, [0]);
    }

    #[test]
    fn distance_boundaries_are_not_confused_with_sequence_coordinates() {
        // The text parser has no sequence. It accepts the full u32 distance
        // range; the record parser separately checks that positions exist.
        let groups = mm_groups("N+n?,0,4294967295;").expect("u32 distances fit");
        let group = groups.first().expect("one group");
        assert_eq!(group.mod_dists, [0, u32::MAX]);
        for distance in ["4294967296", "-1", "1.0", " 1", "1 ", ""] {
            let text = format!("N+n?,{distance};");
            let error = mm_groups(&text).expect_err("invalid distance must fail");
            assert!(
                matches!(error, Error::InvalidModCoords(_)),
                "distance syntax must not be accepted: {distance}"
            );
        }
        let empty_groups = mm_groups("N+n?;").expect("no comma means no distances");
        assert!(
            empty_groups
                .first()
                .expect("one group")
                .mod_dists
                .is_empty(),
            "an absent distance list differs from an empty distance field"
        );
    }
}
