//! `ModChar` struct for handling DNA modification tags
//! Handles both letter and numeric modification codes from BAM files

use crate::Error;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

/// Our struct to hold a modification tag.
/// The BAM file format uses the syntax `base+mod_code` in its ML tag
/// to show which modification is represented e.g. C+m, A+a, T+T, ...
/// This can be a letter or a number e.g. T+472232 represents `BrdU` as that is its `CheBI` code.
/// As we rely on a fibertools-rs data structure to store mod information (`BaseMod`),
/// which uses a char datatype to represent the mod code, representing a one letter
/// mod code is easy, but numbers need to be converted to char first.
/// Fortunately, rust's char datatype is almost equivalent to u32, so we just
/// convert numbers to char before storing.
/// NOTE: the above conversion has some problems e.g. A+a is equivalent to A+97 etc.
/// (as 97 is the ascii code of a),
/// and the rust char datatype does not allow a set of u32s somewhere between 55000
/// and 59000. We have chosen to live with this problem. I think the probability of
/// having a DNA modification with a `CheBI` code overlapping with ASCII values or
/// within this narrow range of values near 59000 is very small.
#[derive(Debug, Clone, Copy, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct ModChar(char);

/// Defaults to mod tag 'N'
impl Default for ModChar {
    fn default() -> Self {
        ModChar::new('N')
    }
}

impl ModChar {
    /// We initialize with a character
    #[must_use]
    pub fn new(val: char) -> Self {
        ModChar(val)
    }
    /// Returns the character
    ///
    /// ```
    /// use nanalogue_core::ModChar;
    /// for y in vec!['a','b','c','\u{D000}']{
    ///     let x = ModChar::new(y.clone());
    ///     assert_eq!(x.val(), y);
    /// }
    /// ```
    #[must_use]
    pub fn val(&self) -> char {
        self.0
    }
}

impl From<char> for ModChar {
    fn from(value: char) -> Self {
        ModChar::new(value)
    }
}

impl From<u8> for ModChar {
    fn from(value: u8) -> Self {
        ModChar::new(char::from(value))
    }
}

impl FromStr for ModChar {
    type Err = Error;

    /// process the modification type from a string,
    /// returning the first character if it is a letter,
    /// or converting it to a character if the first character is a number
    ///
    /// ```
    /// use nanalogue_core::ModChar;
    /// use std::str::FromStr;
    ///
    /// // Single letter modification codes
    /// let mod_char = ModChar::from_str("m")?;
    /// assert_eq!(mod_char.val(), 'm');
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```
    /// # use nanalogue_core::ModChar;
    /// # use std::str::FromStr;
    /// #
    /// // CheBI code for BrdU (5-bromo-2'-deoxyuridine)
    /// let mod_char = ModChar::from_str("472232")?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```
    /// # use nanalogue_core::ModChar;
    /// # use std::str::FromStr;
    /// #
    /// // Small numeric codes
    /// let mod_char = ModChar::from_str("123")?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```should_panic
    /// # use nanalogue_core::ModChar;
    /// # use std::str::FromStr;
    /// #
    /// // Invalid: starts with special character
    /// let mod_char = ModChar::from_str("@123")?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    fn from_str(mod_type: &str) -> Result<Self, Self::Err> {
        let first_char = mod_type
            .chars()
            .next()
            .ok_or(Error::EmptyModType(String::new()))?;
        match first_char {
            'A'..='Z' | 'a'..='z' if mod_type.len() == 1 => Ok(ModChar(first_char)),
            '0'..='9' => {
                let val = char::from_u32(mod_type.parse()?)
                    .ok_or(Error::InvalidModType(mod_type.to_owned()))?;
                Ok(ModChar(val))
            }
            _ => Err(Error::InvalidModType(mod_type.to_owned())),
        }
    }
}

impl fmt::Display for ModChar {
    /// converts to string for display. If the value is in the alphabet,
    /// display it. Otherwise, display the equivalent number.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.val() {
            w @ ('A'..='Z' | 'a'..='z') => w.to_string(),
            w => (w as u32).to_string(),
        }
        .fmt(f)
    }
}

impl Serialize for ModChar {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        // Use the Display implementation to serialize as a string
        serializer.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for ModChar {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        ModChar::from_str(&s).map_err(serde::de::Error::custom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tests if `ModChar` is displayed correctly
    #[test]
    fn display_mod_char() {
        assert_eq!(
            format!("{}", ModChar::from_str("a").expect("no failure")),
            "a"
        );
        assert_eq!(
            format!("{}", ModChar::from_str("T").expect("no failure")),
            "T"
        );
        assert_eq!(
            format!("{}", ModChar::from_str("77000").expect("no failure")),
            "77000"
        );
    }

    /// Tests `ModChar` numeric conversion and edge cases
    #[expect(
        clippy::shadow_unrelated,
        reason = "repetition is fine; each block is clearly separated"
    )]
    #[test]
    fn modchar_numeric_conversion() {
        // Test letter codes
        let mod_char = ModChar::from_str("m").expect("should parse");
        assert_eq!(mod_char.val(), 'm');
        assert_eq!(format!("{mod_char}"), "m");

        let mod_char = ModChar::from_str("T").expect("should parse");
        assert_eq!(mod_char.val(), 'T');
        assert_eq!(format!("{mod_char}"), "T");

        // Test small numeric codes
        let mod_char = ModChar::from_str("123").expect("should parse");
        assert_eq!(format!("{mod_char}"), "123");

        // Test CheBI code for BrdU
        let mod_char = ModChar::from_str("472232").expect("should parse");
        assert_eq!(format!("{mod_char}"), "472232");

        // Test ASCII boundary - 97 is 'a'
        let mod_char = ModChar::from_str("97").expect("should parse");
        assert_eq!(mod_char.val(), 'a');
        // When the char value is in alphabet range, it displays as the letter, not the number
        assert_eq!(format!("{mod_char}"), "a");

        // Test very large numbers that are valid unicode
        let mod_char = ModChar::from_str("65536").expect("should parse");
        assert_eq!(format!("{mod_char}"), "65536");
    }

    #[test]
    #[should_panic(expected = "EmptyModType")]
    fn modchar_empty_string_panics() {
        let _: ModChar = ModChar::from_str("").unwrap();
    }

    #[test]
    #[should_panic(expected = "InvalidModType")]
    fn modchar_special_char_at_panics() {
        let _: ModChar = ModChar::from_str("@123").unwrap();
    }

    #[test]
    #[should_panic(expected = "InvalidModType")]
    fn modchar_special_char_hash_panics() {
        let _: ModChar = ModChar::from_str("#abc").unwrap();
    }

    /// Tests `ModChar` display format consistency
    #[test]
    fn modchar_display_consistency() {
        // Letters should display as letters
        for letter in ['a', 'b', 'z', 'A', 'B', 'Z'] {
            let mod_char = ModChar::new(letter);
            assert_eq!(format!("{mod_char}"), letter.to_string());
        }

        // Numbers converted to char should display as their numeric value
        let test_numbers = vec![123, 456, 789, 472_232];
        for num in test_numbers {
            let mod_char = ModChar::from_str(&num.to_string()).expect("should parse");
            assert_eq!(format!("{mod_char}"), num.to_string());
        }
    }

    /// Tests `From<char>` implementation for `ModChar`
    #[expect(
        clippy::shadow_unrelated,
        reason = "repetition is fine; each block is clearly separated"
    )]
    #[test]
    fn from_char() {
        // Test lowercase letters
        let mod_char = ModChar::from('a');
        assert_eq!(mod_char.val(), 'a');

        let mod_char = ModChar::from('z');
        assert_eq!(mod_char.val(), 'z');

        // Test uppercase letters
        let mod_char = ModChar::from('A');
        assert_eq!(mod_char.val(), 'A');

        let mod_char = ModChar::from('Z');
        assert_eq!(mod_char.val(), 'Z');

        // Test numeric characters
        let mod_char = ModChar::from('0');
        assert_eq!(mod_char.val(), '0');

        let mod_char = ModChar::from('9');
        assert_eq!(mod_char.val(), '9');

        // Test special characters
        let mod_char = ModChar::from('@');
        assert_eq!(mod_char.val(), '@');

        let mod_char = ModChar::from('#');
        assert_eq!(mod_char.val(), '#');

        // Test unicode characters
        let mod_char = ModChar::from('\u{D000}');
        assert_eq!(mod_char.val(), '\u{D000}');

        let mod_char = ModChar::from('\u{1F600}'); // emoji
        assert_eq!(mod_char.val(), '\u{1F600}');
    }

    /// Tests `From<u8>` implementation for `ModChar`
    #[test]
    fn from_u8() {
        // Test that u8 conversion is equivalent to char::from for all possible u8 values
        for byte_val in 0u8..=255u8 {
            let mod_char = ModChar::from(byte_val);
            assert_eq!(mod_char.val(), char::from(byte_val));
        }
    }
}
