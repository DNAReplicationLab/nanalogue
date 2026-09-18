#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Modification-code parsing boundaries and validation through the JSON wire format.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::{Error, ModChar};

    #[test]
    fn rejects_trailing_data_after_either_letter_range() {
        // The uppercase and lowercase alternatives each have a length guard.
        // Check all letters: accepting just the first character would silently
        // turn malformed multi-code input into a different modification.
        for letter in ('A'..='Z').chain('a'..='z') {
            let valid = letter
                .to_string()
                .parse::<ModChar>()
                .expect("single letter");
            assert_eq!(valid.val(), letter);
            assert_eq!(valid.to_string(), letter.to_string());
            for suffix in ["m", "7", " ", "\n", "\0", "\u{e9}"] {
                let input = format!("{letter}{suffix}");
                let error = input.parse::<ModChar>().unwrap_err();
                assert!(
                    matches!(&error, Error::InvalidModType(value) if value == &input),
                    "must reject the whole multi-character code {input:?}: {error}"
                );
            }
        }
    }

    #[test]
    fn accepts_scalar_boundaries_and_canonicalizes_numeric_letters() {
        // Expected scalars and display strings are literal, independent of
        // the parser. Adjacent values straddle both ASCII letter ranges,
        // the surrogate gap, and the final valid Unicode scalar.
        for (input, scalar, display) in [
            ("0", '\0', "0"),
            ("64", '@', "64"),
            ("65", 'A', "A"),
            ("90", 'Z', "Z"),
            ("91", '[', "91"),
            ("96", '`', "96"),
            ("97", 'a', "a"),
            ("122", 'z', "z"),
            ("123", '{', "123"),
            ("55295", '\u{d7ff}', "55295"),
            ("57344", '\u{e000}', "57344"),
            ("1114111", '\u{10ffff}', "1114111"),
            ("00097", 'a', "a"),
        ] {
            let code = input.parse::<ModChar>().expect("valid scalar code");
            assert_eq!(code.val(), scalar, "wrong scalar for {input}");
            assert_eq!(code.to_string(), display, "wrong display for {input}");
            let json = serde_json::to_string(&code).expect("serialize code");
            assert_eq!(json, format!("\"{display}\""));
            let decoded = serde_json::from_str::<ModChar>(&json).expect("decode code");
            assert_eq!(decoded, code);
        }
    }

    #[test]
    fn rejects_non_scalars_without_confusing_them_with_integer_overflow() {
        // Surrogates and values above Unicode's maximum fit in u32, but
        // cannot inhabit Rust's char type. Do not allocate or forge a char.
        for input in ["55296", "56319", "56320", "57343", "1114112", "4294967295"] {
            let error = input.parse::<ModChar>().unwrap_err();
            assert!(
                matches!(&error, Error::InvalidModType(value) if value == input),
                "non-scalar {input} must retain its original code: {error}"
            );
        }
        for input in [
            "4294967296",
            "18446744073709551616",
            "1m",
            "1.0",
            "1 ",
            "0x61",
        ] {
            let error = input.parse::<ModChar>().unwrap_err();
            assert!(
                matches!(error, Error::IntParseError(_)),
                "digit-prefixed invalid integer {input:?} must fail integer parsing: {error}"
            );
        }
    }

    #[test]
    fn rejects_non_ascii_letters_and_unrecognized_prefixes() {
        for input in ["\u{e9}", "\u{ff4d}", "\u{661}", "-1", "+97", " 97", "\0"] {
            let error = input.parse::<ModChar>().unwrap_err();
            assert!(
                matches!(&error, Error::InvalidModType(value) if value == input),
                "only ASCII letters or digit-prefixed numbers are accepted: {input:?}"
            );
        }
        assert!(
            matches!("".parse::<ModChar>().unwrap_err(), Error::EmptyModType(value) if value.is_empty()),
            "empty input must have its distinct error"
        );
    }

    #[test]
    fn json_deserialization_cannot_bypass_parser_validation() {
        for input in [
            "Mm",
            "mm",
            "55296",
            "1114112",
            "4294967296",
            "1m",
            "",
            "\u{e9}",
        ] {
            let parser_error = input.parse::<ModChar>().unwrap_err().to_string();
            let json = serde_json::to_string(input).expect("encode invalid input as a string");
            let wire_error = serde_json::from_str::<ModChar>(&json).unwrap_err();
            assert!(wire_error.is_data(), "must reject invalid code as data");
            assert!(
                wire_error.to_string().contains(&parser_error),
                "JSON must preserve the parser diagnosis for {input:?}: {wire_error}"
            );
        }
        // Numeric codes still use a JSON string; a bare number or array must
        // not accidentally be accepted as an alternative representation.
        for json in ["97", "null", "true", "[]", "[97]", "{}"] {
            let error = serde_json::from_str::<ModChar>(json).unwrap_err();
            assert!(error.is_data(), "unexpected JSON representation: {json}");
        }
        let numeric = serde_json::from_str::<ModChar>("\"97\"").expect("string code");
        let letter = serde_json::from_str::<ModChar>("\"a\"").expect("letter code");
        assert_eq!(numeric, letter);
        assert_eq!(numeric.val(), 'a');
    }
}
