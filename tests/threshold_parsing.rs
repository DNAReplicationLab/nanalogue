//! Boundary coverage for fraction-band parsing and its byte-length guard.

#[cfg(test)]
mod tests {
    use nanalogue_core::{Contains as _, Error, OrdPair, ThresholdState};

    #[test]
    fn accepts_fifty_bytes_but_rejects_fifty_one() {
        // Leading whitespace is valid inside a component. Keep the numeric
        // values identical so only the length guard distinguishes these cases.
        let at_limit = format!("{}0.2,0.8", " ".repeat(43));
        let over_limit = format!(" {at_limit}");
        assert_eq!(at_limit.len(), 50);
        assert_eq!(over_limit.len(), 51);
        assert_eq!(
            ThresholdState::from_str_ordpair_fraction(&at_limit).expect("50 bytes are allowed"),
            ThresholdState::InvertGtEqLtEq(OrdPair::new(51, 204).expect("ordered bytes"))
        );
        let error = ThresholdState::from_str_ordpair_fraction(&over_limit)
            .expect_err("51 bytes exceed the guard");
        assert!(
            matches!(&error, Error::InvalidState(message)
                if message == "OrdPair conversion from very long string attempted!"),
            "length rejection must precede fraction parsing: {error}"
        );
    }

    #[test]
    fn counts_bytes_rather_than_unicode_characters() {
        // EM SPACE is trimmed by the component parser, but occupies three
        // UTF-8 bytes. Both strings have far fewer than 50 characters.
        let at_limit = format!("{} 0.2,0.8", "\u{2003}".repeat(14));
        let over_limit = format!("{}0.2,0.8", "\u{2003}".repeat(15));
        assert_eq!(at_limit.len(), 50);
        assert_eq!(over_limit.len(), 52);
        assert!(
            over_limit.chars().count() < 50,
            "byte and character limits differ"
        );
        assert_eq!(
            ThresholdState::from_str_ordpair_fraction(&at_limit).expect("trimmed Unicode spaces"),
            ThresholdState::InvertGtEqLtEq(OrdPair::new(51, 204).expect("ordered bytes"))
        );
        assert!(
            matches!(
                ThresholdState::from_str_ordpair_fraction(&over_limit).unwrap_err(),
                Error::InvalidState(_)
            ),
            "the guard must reject oversized UTF-8 input"
        );
    }

    #[test]
    fn distinguishes_empty_input_from_whitespace_and_malformed_pairs() {
        let unfiltered =
            ThresholdState::from_str_ordpair_fraction("").expect("empty disables filter");
        assert_eq!(unfiltered, ThresholdState::GtEq(0));
        for value in 0..=255 {
            assert!(unfiltered.contains(&value), "empty filter rejected {value}");
        }
        for input in [" ", "\t\n", ",", "0.2,", ",0.8", "0.2,0.8,1", "0.2;0.8"] {
            let error = ThresholdState::from_str_ordpair_fraction(input).unwrap_err();
            assert!(
                matches!(error, Error::OrdPairConversion(_)),
                "malformed pair {input:?} should report conversion failure: {error}"
            );
        }
        let reversed = ThresholdState::from_str_ordpair_fraction("0.8,0.2").unwrap_err();
        assert!(
            matches!(reversed, Error::WrongOrder(_)),
            "valid fractions in reverse order must retain the ordering error"
        );
    }

    #[test]
    fn rejects_invalid_fractions_in_either_component() {
        for invalid in ["NaN", "inf", "-inf", "-0.01", "1.01", "invalid"] {
            for input in [format!("{invalid},1"), format!("0,{invalid}")] {
                let error = ThresholdState::from_str_ordpair_fraction(&input).unwrap_err();
                assert!(
                    matches!(error, Error::OrdPairConversion(_)),
                    "invalid fraction in {input:?} escaped component validation: {error}"
                );
            }
        }
    }

    #[test]
    fn preserves_inclusive_band_edges_after_quantization() {
        // Expected byte bounds are independently calculated: round(255*p).
        // Equal endpoints and distinct fractions rounding to the same byte
        // must remain valid one-byte exclusion bands, not become empty bands.
        for (input, low, high) in [
            ("0,1", 0u8, 255u8),
            ("0,0", 0, 0),
            ("1,1", 255, 255),
            ("0.5,0.5", 128, 128),
            ("0.5001,0.5002", 128, 128),
            ("0.2,0.8", 51, 204),
            (" 2e-1 , 8e-1 ", 51, 204),
        ] {
            let threshold = ThresholdState::from_str_ordpair_fraction(input).expect("valid band");
            assert_eq!(
                threshold,
                ThresholdState::InvertGtEqLtEq(OrdPair::new(low, high).expect("ordered bounds")),
                "incorrect quantization for {input}"
            );
            for value in 0..=255 {
                assert_eq!(
                    threshold.contains(&value),
                    value < low || value > high,
                    "incorrect membership for {value} in {input}"
                );
            }
        }
    }
}
