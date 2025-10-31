//! `AllowedAGCTN` enum for DNA bases A, G, C, T, and N
//! Provides type-safe representation of allowed DNA modification bases
//! in the modBAM format

use crate::Error;
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
            _ => Err(Error::InvalidBase),
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
            _ => Err(Error::InvalidBase),
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
}
