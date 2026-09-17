//! Deterministic boundary tests for the UUIDs used as read identifiers.

#[cfg(test)]
mod tests {
    use nanalogue_core::uuid::{is_valid_v4, v4_from_bytes};

    /// A valid UUID with distinct groups makes misplaced separators visible.
    const VALID: &str = "01234567-89ab-4cde-8f01-23456789abcd";

    #[test]
    fn rejects_non_hex_in_every_payload_position() {
        // Each match arm covers a different group of hexadecimal positions.
        // A malformed last group alone cannot exercise the earlier guards.
        let positions = (0..8)
            .chain(9..13)
            .chain(15..18)
            .chain(20..23)
            .chain(24..36);
        for position in positions {
            for replacement in [b'g', b'G', b'/', b':', b' ', b'-', 0, 0x7f] {
                let mut bytes = VALID.as_bytes().to_vec();
                *bytes.get_mut(position).expect("payload position") = replacement;
                let candidate = String::from_utf8(bytes).expect("ASCII input");
                assert!(
                    !is_valid_v4(&candidate),
                    "accepted non-hex byte {replacement} at {position}"
                );
            }
        }
    }

    #[test]
    fn accepts_every_ascii_hex_digit_in_payload_positions() {
        // Only payload positions accept arbitrary hex digits: version and
        // variant positions are intentionally absent from this list.
        for position in [0, 7, 9, 12, 15, 17, 20, 22, 24, 35] {
            for replacement in b"0123456789abcdefABCDEF" {
                let mut bytes = VALID.as_bytes().to_vec();
                *bytes.get_mut(position).expect("payload position") = *replacement;
                let candidate = String::from_utf8(bytes).expect("ASCII input");
                assert!(
                    is_valid_v4(&candidate),
                    "rejected hex byte {replacement} at {position}"
                );
            }
        }
    }

    #[test]
    fn requires_exact_length_and_separator_positions() {
        for length in 0..36 {
            let prefix = VALID.get(..length).expect("ASCII boundary");
            assert!(!is_valid_v4(prefix), "accepted prefix of length {length}");
        }
        for suffix in ["0", "-", "\n", "\0"] {
            assert!(
                !is_valid_v4(&format!("{VALID}{suffix}")),
                "accepted trailing bytes"
            );
        }
        for position in [8, 13, 18, 23] {
            let mut bytes = VALID.as_bytes().to_vec();
            *bytes.get_mut(position).expect("separator position") = b'0';
            let candidate = String::from_utf8(bytes).expect("ASCII input");
            assert!(
                !is_valid_v4(&candidate),
                "accepted missing separator at {position}"
            );
        }
        assert!(is_valid_v4(VALID), "control must be valid");
    }

    #[test]
    fn constrains_version_and_variant_independently() {
        for replacement in b"0123456789abcdefABCDEF-g" {
            let mut version_bytes = VALID.as_bytes().to_vec();
            *version_bytes.get_mut(14).expect("version position") = *replacement;
            let version = String::from_utf8(version_bytes).expect("ASCII input");
            assert_eq!(
                is_valid_v4(&version),
                *replacement == b'4',
                "incorrect version acceptance for {replacement}"
            );

            let mut variant_bytes = VALID.as_bytes().to_vec();
            *variant_bytes.get_mut(19).expect("variant position") = *replacement;
            let variant = String::from_utf8(variant_bytes).expect("ASCII input");
            assert_eq!(
                is_valid_v4(&variant),
                b"89abAB".contains(replacement),
                "incorrect variant acceptance for {replacement}"
            );
        }
    }

    #[test]
    fn formatting_sets_reserved_bits_and_preserves_payload() {
        // Asymmetric bytes catch ordering errors; high bits in the version
        // and variant bytes must be cleared, not merely ORed with the mask.
        let bytes = [
            0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xfe, 0xdc, 0x7f, 0x10, 0x32, 0x54, 0x76, 0x98,
            0xba, 0xcd,
        ];
        let formatted = v4_from_bytes(bytes);
        assert_eq!(formatted, "01234567-89ab-4edc-bf10-32547698bacd");
        assert!(is_valid_v4(&formatted), "formatter must produce valid v4");
        assert_eq!(
            v4_from_bytes([0xff; 16]),
            "ffffffff-ffff-4fff-bfff-ffffffffffff"
        );
    }
}
