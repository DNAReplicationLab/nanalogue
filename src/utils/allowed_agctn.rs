//! `AllowedAGCTN` enum for DNA bases A, G, C, T, and N
//! Provides type-safe representation of allowed DNA modification bases
//! in the modBAM format

use crate::Error;
use rand::Rng;
use rand::distr::{Distribution, StandardUniform};
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

/// Represents the allowed DNA bases for modification: A, G, C, T, or N
/// in the modBAM format
#[expect(
    clippy::exhaustive_enums,
    reason = "A, G, C, T, and N are the only modBAM DNA bases; this set is fixed and will never change"
)]
#[derive(
    Debug, Clone, Default, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize,
)]
pub enum AllowedAGCTN {
    /// Adenine
    A,
    /// Guanine
    G,
    /// Cytosine
    C,
    /// Thymine
    T,
    /// Any base (N)
    #[default]
    N,
}

// Implements random pick of a variant
impl Distribution<AllowedAGCTN> for StandardUniform {
    /// Allows us to randomly pick a variant
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> AllowedAGCTN {
        match rng.random_range(0..5) {
            0 => AllowedAGCTN::A,
            1 => AllowedAGCTN::G,
            2 => AllowedAGCTN::C,
            3 => AllowedAGCTN::T,
            4 => AllowedAGCTN::N,
            _ => unreachable!(),
        }
    }
}

/// Implements conversion from `AllowedAGCTN` to `char`
impl From<AllowedAGCTN> for char {
    fn from(base: AllowedAGCTN) -> Self {
        match base {
            AllowedAGCTN::A => 'A',
            AllowedAGCTN::G => 'G',
            AllowedAGCTN::C => 'C',
            AllowedAGCTN::T => 'T',
            AllowedAGCTN::N => 'N',
        }
    }
}

/// Implements conversion from `AllowedAGCTN` to `u8`
impl From<AllowedAGCTN> for u8 {
    fn from(base: AllowedAGCTN) -> Self {
        match base {
            AllowedAGCTN::A => b'A',
            AllowedAGCTN::G => b'G',
            AllowedAGCTN::C => b'C',
            AllowedAGCTN::T => b'T',
            AllowedAGCTN::N => b'N',
        }
    }
}

/// Implements parsing from string
///
/// ```
/// use nanalogue_core::AllowedAGCTN;
/// use std::str::FromStr;
///
/// assert_eq!(AllowedAGCTN::from_str("A")?, AllowedAGCTN::A);
/// assert_eq!(AllowedAGCTN::from_str("G")?, AllowedAGCTN::G);
/// assert_eq!(AllowedAGCTN::from_str("C")?, AllowedAGCTN::C);
/// assert_eq!(AllowedAGCTN::from_str("T")?, AllowedAGCTN::T);
/// assert_eq!(AllowedAGCTN::from_str("N")?, AllowedAGCTN::N);
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
///
/// ```should_panic
/// # use nanalogue_core::AllowedAGCTN;
/// # use std::str::FromStr;
/// // Invalid base should error
/// let base = AllowedAGCTN::from_str("X")?;
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
impl FromStr for AllowedAGCTN {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "A" => Ok(AllowedAGCTN::A),
            "G" => Ok(AllowedAGCTN::G),
            "C" => Ok(AllowedAGCTN::C),
            "T" => Ok(AllowedAGCTN::T),
            "N" => Ok(AllowedAGCTN::N),
            v => Err(Error::InvalidBase(v.to_owned())),
        }
    }
}

/// Implements conversion from `char`
///
/// ```
/// use nanalogue_core::AllowedAGCTN;
///
/// assert_eq!(AllowedAGCTN::try_from('A')?, AllowedAGCTN::A);
/// assert_eq!(AllowedAGCTN::try_from('G')?, AllowedAGCTN::G);
/// assert_eq!(AllowedAGCTN::try_from('C')?, AllowedAGCTN::C);
/// assert_eq!(AllowedAGCTN::try_from('T')?, AllowedAGCTN::T);
/// assert_eq!(AllowedAGCTN::try_from('N')?, AllowedAGCTN::N);
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
///
/// ```should_panic
/// # use nanalogue_core::AllowedAGCTN;
/// // Invalid base should error
/// let base = AllowedAGCTN::try_from('X')?;
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
impl TryFrom<char> for AllowedAGCTN {
    type Error = Error;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            'A' => Ok(AllowedAGCTN::A),
            'G' => Ok(AllowedAGCTN::G),
            'C' => Ok(AllowedAGCTN::C),
            'T' => Ok(AllowedAGCTN::T),
            'N' => Ok(AllowedAGCTN::N),
            v => Err(Error::InvalidBase(v.to_string())),
        }
    }
}

/// Implements conversion from `u8`
///
/// ```
/// use nanalogue_core::AllowedAGCTN;
///
/// assert_eq!(AllowedAGCTN::try_from(b'A')?, AllowedAGCTN::A);
/// assert_eq!(AllowedAGCTN::try_from(b'G')?, AllowedAGCTN::G);
/// assert_eq!(AllowedAGCTN::try_from(b'C')?, AllowedAGCTN::C);
/// assert_eq!(AllowedAGCTN::try_from(b'T')?, AllowedAGCTN::T);
/// assert_eq!(AllowedAGCTN::try_from(b'N')?, AllowedAGCTN::N);
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
///
/// ```should_panic
/// # use nanalogue_core::AllowedAGCTN;
/// // Invalid base should error
/// let base = AllowedAGCTN::try_from(b'X')?;
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
impl TryFrom<u8> for AllowedAGCTN {
    type Error = Error;

    fn try_from(c: u8) -> Result<Self, Self::Error> {
        match c {
            b'A' => Ok(AllowedAGCTN::A),
            b'G' => Ok(AllowedAGCTN::G),
            b'C' => Ok(AllowedAGCTN::C),
            b'T' => Ok(AllowedAGCTN::T),
            b'N' => Ok(AllowedAGCTN::N),
            v => Err(Error::InvalidBase(v.to_string())),
        }
    }
}

/// Implements printing of base
impl fmt::Display for AllowedAGCTN {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        char::from(*self).to_string().fmt(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tests default value is N
    #[test]
    fn default_is_n() {
        assert_eq!(AllowedAGCTN::default(), AllowedAGCTN::N);
    }

    /// Tests `From<AllowedAGCTN>` for `char`
    #[test]
    fn from_allowed_agctn_to_char() {
        assert_eq!(char::from(AllowedAGCTN::A), 'A');
        assert_eq!(char::from(AllowedAGCTN::G), 'G');
        assert_eq!(char::from(AllowedAGCTN::C), 'C');
        assert_eq!(char::from(AllowedAGCTN::T), 'T');
        assert_eq!(char::from(AllowedAGCTN::N), 'N');
    }

    /// Tests `From<AllowedAGCTN>` for `u8`
    #[test]
    fn from_allowed_agctn_to_u8() {
        assert_eq!(u8::from(AllowedAGCTN::A), b'A');
        assert_eq!(u8::from(AllowedAGCTN::G), b'G');
        assert_eq!(u8::from(AllowedAGCTN::C), b'C');
        assert_eq!(u8::from(AllowedAGCTN::T), b'T');
        assert_eq!(u8::from(AllowedAGCTN::N), b'N');
    }

    /// Tests `FromStr` for valid bases
    #[test]
    fn from_str_valid_bases() {
        assert_eq!(AllowedAGCTN::from_str("A").unwrap(), AllowedAGCTN::A);
        assert_eq!(AllowedAGCTN::from_str("G").unwrap(), AllowedAGCTN::G);
        assert_eq!(AllowedAGCTN::from_str("C").unwrap(), AllowedAGCTN::C);
        assert_eq!(AllowedAGCTN::from_str("T").unwrap(), AllowedAGCTN::T);
        assert_eq!(AllowedAGCTN::from_str("N").unwrap(), AllowedAGCTN::N);
    }

    /// Tests `FromStr` for invalid bases
    #[test]
    #[should_panic(expected = "InvalidBase")]
    fn from_str_invalid_base_x() {
        let _: AllowedAGCTN = AllowedAGCTN::from_str("X").unwrap();
    }

    /// Tests `FromStr` for lowercase (should fail)
    #[test]
    #[should_panic(expected = "InvalidBase")]
    fn from_str_lowercase_fails() {
        let _: AllowedAGCTN = AllowedAGCTN::from_str("a").unwrap();
    }

    /// Tests `FromStr` for empty string
    #[test]
    #[should_panic(expected = "InvalidBase")]
    fn from_str_empty_string() {
        let _: AllowedAGCTN = AllowedAGCTN::from_str("").unwrap();
    }

    /// Tests `TryFrom<char>` for valid bases
    #[test]
    fn try_from_char_valid_bases() {
        assert_eq!(AllowedAGCTN::try_from('A').unwrap(), AllowedAGCTN::A);
        assert_eq!(AllowedAGCTN::try_from('G').unwrap(), AllowedAGCTN::G);
        assert_eq!(AllowedAGCTN::try_from('C').unwrap(), AllowedAGCTN::C);
        assert_eq!(AllowedAGCTN::try_from('T').unwrap(), AllowedAGCTN::T);
        assert_eq!(AllowedAGCTN::try_from('N').unwrap(), AllowedAGCTN::N);
    }

    /// Tests `TryFrom<char>` for invalid bases
    #[test]
    #[should_panic(expected = "InvalidBase")]
    fn try_from_char_invalid_base() {
        let _: AllowedAGCTN = AllowedAGCTN::try_from('X').unwrap();
    }

    /// Tests `TryFrom<char>` for lowercase (should fail)
    #[test]
    #[should_panic(expected = "InvalidBase")]
    fn try_from_char_lowercase_fails() {
        let _: AllowedAGCTN = AllowedAGCTN::try_from('a').unwrap();
    }

    /// Tests `Display` implementation
    #[test]
    fn display_works() {
        assert_eq!(format!("{}", AllowedAGCTN::A), "A");
        assert_eq!(format!("{}", AllowedAGCTN::G), "G");
        assert_eq!(format!("{}", AllowedAGCTN::C), "C");
        assert_eq!(format!("{}", AllowedAGCTN::T), "T");
        assert_eq!(format!("{}", AllowedAGCTN::N), "N");
    }

    /// Tests deserialization failure for invalid base
    #[test]
    #[should_panic(expected = "unknown variant")]
    fn deserialize_invalid_base() {
        let _: AllowedAGCTN = serde_json::from_str(r#""X""#).unwrap();
    }

    /// Tests deserialization failure for lowercase base
    #[test]
    #[should_panic(expected = "unknown variant")]
    fn deserialize_lowercase_base_fails() {
        let _: AllowedAGCTN = serde_json::from_str(r#""a""#).unwrap();
    }

    #[test]
    fn try_from_u8_valid_bases() {
        // Test all valid uppercase bases
        assert_eq!(
            AllowedAGCTN::try_from(b'A').expect("should convert"),
            AllowedAGCTN::A
        );
        assert_eq!(
            AllowedAGCTN::try_from(b'G').expect("should convert"),
            AllowedAGCTN::G
        );
        assert_eq!(
            AllowedAGCTN::try_from(b'C').expect("should convert"),
            AllowedAGCTN::C
        );
        assert_eq!(
            AllowedAGCTN::try_from(b'T').expect("should convert"),
            AllowedAGCTN::T
        );
        assert_eq!(
            AllowedAGCTN::try_from(b'N').expect("should convert"),
            AllowedAGCTN::N
        );
    }

    #[test]
    fn try_from_u8_invalid_bases() {
        // Test invalid bases
        let result_x: Result<AllowedAGCTN, _> = AllowedAGCTN::try_from(b'X');
        let _: Error = result_x.unwrap_err();

        let result_lower_a: Result<AllowedAGCTN, _> = AllowedAGCTN::try_from(b'a');
        let _: Error = result_lower_a.unwrap_err();

        let result_lower_g: Result<AllowedAGCTN, _> = AllowedAGCTN::try_from(b'g');
        let _: Error = result_lower_g.unwrap_err();

        let result_digit: Result<AllowedAGCTN, _> = AllowedAGCTN::try_from(b'1');
        let _: Error = result_digit.unwrap_err();

        let result_space: Result<AllowedAGCTN, _> = AllowedAGCTN::try_from(b' ');
        let _: Error = result_space.unwrap_err();
    }

    #[test]
    #[expect(
        clippy::panic,
        reason = "panic is appropriate in tests for wrong error type"
    )]
    fn try_from_u8_error_type() {
        // Verify error contains the invalid base
        let result: Result<AllowedAGCTN, _> = AllowedAGCTN::try_from(b'Z');
        let err = result.unwrap_err();
        if let Error::InvalidBase(s) = err {
            // The error contains the numeric representation of the byte
            assert_eq!(s, (b'Z').to_string());
        } else {
            panic!("Expected InvalidBase error, got {err:?}");
        }
    }

    #[test]
    fn try_from_u8_roundtrip() {
        // Test converting to u8 and back
        for base in [
            AllowedAGCTN::A,
            AllowedAGCTN::G,
            AllowedAGCTN::C,
            AllowedAGCTN::T,
            AllowedAGCTN::N,
        ] {
            let as_u8: u8 = base.into();
            let converted_back = AllowedAGCTN::try_from(as_u8).expect("should convert back");
            assert_eq!(converted_back, base);
        }
    }

    /// Tests random `AllowedAGCTN` generation from `StandardUniform` produces all variants
    #[test]
    fn allowed_agctn_random_generation_all_variants() {
        let mut rng = rand::rng();

        // Generate many random bases to ensure all variants appear
        let mut generated_bases = std::collections::HashSet::new();
        for _ in 0..1000 {
            let base: AllowedAGCTN = rng.random();
            let _: bool = generated_bases.insert(base);
        }

        // Verify all 5 variants can be generated
        assert_eq!(generated_bases.len(), 5);
        assert!(generated_bases.contains(&AllowedAGCTN::A));
        assert!(generated_bases.contains(&AllowedAGCTN::G));
        assert!(generated_bases.contains(&AllowedAGCTN::C));
        assert!(generated_bases.contains(&AllowedAGCTN::T));
        assert!(generated_bases.contains(&AllowedAGCTN::N));
    }
}
