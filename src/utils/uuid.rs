//! Minimal UUID v4 helpers used within nanalogue.
//!
//! This module intentionally implements only the small subset of UUID
//! functionality that nanalogue needs:
//! - generate random UUID v4 strings
//! - format UUID v4 strings from caller-supplied random bytes
//! - validate lowercase/uppercase hyphenated UUID strings
//!
//! It does not provide a general-purpose UUID type.

use rand::random;
use std::fmt::Write as _;

/// Returns a random UUID v4 formatted as a hyphenated lowercase string.
#[must_use]
pub fn v4_random() -> String {
    v4_from_bytes(random())
}

/// Formats caller-supplied random bytes as a UUID v4 string.
///
/// The version and variant bits are set according to RFC 9562 before the
/// hyphenated lowercase string is produced.
#[must_use]
pub fn v4_from_bytes(mut bytes: [u8; 16]) -> String {
    // Set UUID version to 4.
    bytes[6] = (bytes[6] & 0x0F) | 0x40;
    // Set RFC 4122 / RFC 9562 variant bits.
    bytes[8] = (bytes[8] & 0x3F) | 0x80;

    let mut out = String::with_capacity(36);
    for (idx, byte) in bytes.into_iter().enumerate() {
        if matches!(idx, 4 | 6 | 8 | 10) {
            out.push('-');
        }
        write!(out, "{byte:02x}").expect("writing to a String should not fail");
    }
    out
}

/// Returns whether the provided string is a valid hyphenated UUID.
///
/// This accepts uppercase and lowercase ASCII hexadecimal digits in the
/// canonical 8-4-4-4-12 hyphenated form.
#[must_use]
pub fn is_valid_v4(s: &str) -> bool {
    for (idx, byte) in s.bytes().enumerate() {
        match (idx, byte) {
            (8 | 13 | 18 | 23, b'-')
            | (14, b'4')
            | (19, b'8' | b'9' | b'a' | b'b' | b'A' | b'B') => {}
            (0..=7 | 9..=12 | 15..=17 | 20..=22 | 24..=35, value) if value.is_ascii_hexdigit() => {}
            _ => return false,
        }
    }

    s.len() == 36
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formats_random_bytes_as_v4_uuid() {
        let uuid = v4_from_bytes([0; 16]);
        assert_eq!(uuid, "00000000-0000-4000-8000-000000000000");
    }

    #[test]
    fn validates_hyphenated_uuid() {
        assert!(is_valid_v4("550e8400-e29b-41d4-a716-446655440000"));
        assert!(is_valid_v4("550E8400-E29B-41D4-A716-446655440000"));
        assert!(!is_valid_v4("not-a-uuid"));
        assert!(!is_valid_v4("550e8400e29b41d4a716446655440000"));
        assert!(!is_valid_v4("550e8400-e29b-11d4-a716-446655440000"));
        assert!(!is_valid_v4("550e8400-e29b-41d4-c716-446655440000"));
    }
}
