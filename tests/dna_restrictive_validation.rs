#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Constructor and JSON boundaries for the validated, nonempty DNA wrapper.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::{DNARestrictive, Error, GetDNARestrictive as _};

    #[test]
    fn rejects_empty_input_through_every_constructor() {
        // This rejection is distinct from InvalidBase: there is no offending
        // byte, and callers must not receive a valid zero-length sequence.
        for result in [
            "".parse::<DNARestrictive>(),
            DNARestrictive::try_from(Vec::<u8>::new()),
        ] {
            assert!(
                matches!(result, Err(Error::InvalidSeq(message))
                    if message == "empty sequence supplied!"),
                "empty input must fail sequence validation"
            );
        }
        let error = serde_json::from_str::<DNARestrictive>(r#""""#).unwrap_err();
        assert!(
            error.to_string().contains("empty sequence supplied!"),
            "JSON must use the validating constructor: {error}"
        );
    }

    #[test]
    fn accepts_only_the_eight_canonical_ascii_base_bytes() {
        // The full byte domain includes invalid UTF-8, which FromStr alone
        // cannot supply. Expected bases are literal rather than produced by
        // the implementation's normalization routine.
        for byte in u8::MIN..=u8::MAX {
            let expected = match byte {
                b'A' | b'a' => Some(b'A'),
                b'C' | b'c' => Some(b'C'),
                b'G' | b'g' => Some(b'G'),
                b'T' | b't' => Some(b'T'),
                _ => None,
            };
            // Prefix/suffix placement catches implementations that validate
            // only the first byte, or silently discard an invalid last byte.
            for (input, normalized) in [
                (vec![byte, b'c'], expected.map(|base| vec![base, b'C'])),
                (vec![b'g', byte], expected.map(|base| vec![b'G', base])),
            ] {
                let result = DNARestrictive::try_from(input);
                if let Some(bases) = normalized {
                    let sequence = result.expect("canonical bases are valid");
                    assert_eq!(sequence.get(), bases, "byte {byte}");
                } else {
                    assert!(
                        matches!(result, Err(Error::InvalidBase(message))
                            if message == char::from(byte).to_string()),
                        "must reject byte {byte} rather than repair or discard it"
                    );
                }
            }
        }
    }

    #[test]
    fn reports_the_first_invalid_byte_without_repairing_input() {
        for input in ["NACGT", "aNgX", "acgtN", "a c", "a\0c", "a\nc"] {
            let expected = if input.contains('N') {
                "N"
            } else if input.contains(' ') {
                " "
            } else if input.contains('\0') {
                "\0"
            } else {
                "\n"
            };
            assert!(
                matches!(input.parse::<DNARestrictive>(),
                    Err(Error::InvalidBase(message)) if message == expected),
                "input {input:?} must fail on its first invalid byte"
            );
            let json = serde_json::to_string(input).expect("encode input string");
            assert!(
                serde_json::from_str::<DNARestrictive>(&json).is_err(),
                "JSON must not bypass validation for {input:?}"
            );
        }
    }

    #[test]
    fn requires_a_json_string_and_emits_canonical_sequence() {
        for json in ["null", "true", "7", "[65,67,71,84]", r#"{"0":[65]}"#] {
            assert!(
                serde_json::from_str::<DNARestrictive>(json).is_err(),
                "the wire format is a string, not {json}"
            );
        }
        let sequence: DNARestrictive =
            serde_json::from_str(r#""tGaaCcT""#).expect("mixed-case DNA string");
        assert_eq!(sequence.get(), b"TGAACCT");
        assert_eq!(sequence.get_dna_restrictive(), &sequence);
        assert_eq!(sequence.to_string(), "TGAACCT");
        assert_eq!(serde_json::to_string(&sequence).unwrap(), r#""TGAACCT""#);
        assert_eq!(sequence, "TGAACCT".parse::<DNARestrictive>().unwrap());
        // Serialization's invalid-UTF-8 error is unreachable through public
        // constructors: every stored byte is one of A/C/G/T. Do not forge an
        // invalid private representation merely to cover that error path.
        // The length-overflow guard is reachable on 64-bit systems only with
        // an allocation above 4 GiB; that is not a practical unit-test input.
    }
}
