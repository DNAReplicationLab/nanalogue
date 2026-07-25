//! `MmSuffix` enum for controlling the trailing mark of generated MM tags.
//!
//! The SAM MM tag format allows three trailing forms on each group:
//! - `?` — explicit: omitted positions are *not* assumed unmodified.
//! - `.` — implicit: omitted positions are assumed unmodified.
//! - no mark — implicit (same semantics as `.`), emitted with no trailing
//!   character.
//!
//! This module provides a validated, serde-compatible type so the simulation
//! feature can offer these as a per-mod configuration option.

use crate::Error;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

/// Trailing mark style for a generated MM tag group.
///
/// See the [module docs](self) for the semantics of each variant.
#[derive(Debug, Clone, Copy, Eq, Hash, PartialEq, PartialOrd)]
#[non_exhaustive]
pub enum MmSuffix {
    /// `?` — explicit (omitted positions are not assumed unmodified).
    QuestionMark,
    /// `.` — implicit (omitted positions are assumed unmodified).
    Dot,
    /// No trailing mark — implicit, emitted with no suffix character.
    None,
}

/// Defaults to `?` (explicit), preserving the simulator's historical behavior.
impl Default for MmSuffix {
    fn default() -> Self {
        Self::QuestionMark
    }
}

impl MmSuffix {
    /// Returns the literal suffix string to append when formatting an MM group.
    ///
    /// ```rust
    /// use nanalogue_core::MmSuffix;
    /// use std::str::FromStr;
    /// assert_eq!(MmSuffix::from_str("?").unwrap().suffix_str(), "?");
    /// assert_eq!(MmSuffix::from_str(".").unwrap().suffix_str(), ".");
    /// assert_eq!(MmSuffix::from_str("none").unwrap().suffix_str(), "");
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    #[must_use]
    pub const fn suffix_str(&self) -> &'static str {
        match *self {
            Self::QuestionMark => "?",
            Self::Dot => ".",
            Self::None => "",
        }
    }
}

impl FromStr for MmSuffix {
    type Err = Error;

    /// Parse an MM suffix style from a string.
    ///
    /// Accepts `"?"`, `"."`, or `"none"` (case-sensitive). The input length is
    /// validated to be greater than zero and at most four bytes before any
    /// further processing, to guard against unexpectedly large inputs.
    ///
    /// ```rust
    /// use nanalogue_core::MmSuffix;
    /// use std::str::FromStr;
    /// assert_eq!(MmSuffix::from_str("?")?, MmSuffix::QuestionMark);
    /// assert_eq!(MmSuffix::from_str(".")?, MmSuffix::Dot);
    /// assert_eq!(MmSuffix::from_str("none")?, MmSuffix::None);
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```should_panic(expected = "InvalidMmSuffix")
    /// # use nanalogue_core::MmSuffix;
    /// # use std::str::FromStr;
    /// let _: MmSuffix = MmSuffix::from_str("!").unwrap();
    /// ```
    ///
    /// ```should_panic(expected = "InvalidMmSuffix")
    /// # use nanalogue_core::MmSuffix;
    /// # use std::str::FromStr;
    /// let _: MmSuffix = MmSuffix::from_str("").unwrap();
    /// ```
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let len = s.len();
        if len == 0 || len > 4 {
            return Err(Error::InvalidMmSuffix(s.to_owned()));
        }
        match s {
            "?" => Ok(Self::QuestionMark),
            "." => Ok(Self::Dot),
            "none" => Ok(Self::None),
            _ => Err(Error::InvalidMmSuffix(s.to_owned())),
        }
    }
}

impl fmt::Display for MmSuffix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.suffix_str())
    }
}

impl Serialize for MmSuffix {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        match *self {
            Self::QuestionMark => serializer.serialize_str("?"),
            Self::Dot => serializer.serialize_str("."),
            Self::None => serializer.serialize_str("none"),
        }
    }
}

impl<'de> Deserialize<'de> for MmSuffix {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        Self::from_str(&s).map_err(serde::de::Error::custom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_str_valid_values() {
        assert_eq!(MmSuffix::from_str("?").unwrap(), MmSuffix::QuestionMark);
        assert_eq!(MmSuffix::from_str(".").unwrap(), MmSuffix::Dot);
        assert_eq!(MmSuffix::from_str("none").unwrap(), MmSuffix::None);
    }

    #[test]
    fn from_str_invalid_values() {
        assert!(matches!(
            MmSuffix::from_str("!"),
            Err(Error::InvalidMmSuffix(_))
        ));
        assert!(matches!(
            MmSuffix::from_str("question"),
            Err(Error::InvalidMmSuffix(_))
        ));
        assert!(matches!(
            MmSuffix::from_str("NONE"),
            Err(Error::InvalidMmSuffix(_))
        ));
    }

    #[test]
    fn from_str_rejects_empty_and_overlong() {
        assert!(matches!(
            MmSuffix::from_str(""),
            Err(Error::InvalidMmSuffix(_))
        ));
        assert!(matches!(
            MmSuffix::from_str("nonee"),
            Err(Error::InvalidMmSuffix(_))
        ));
    }

    #[test]
    fn default_is_question_mark() {
        assert_eq!(MmSuffix::default(), MmSuffix::QuestionMark);
    }

    #[test]
    fn display_matches_suffix_str() {
        assert_eq!(format!("{}", MmSuffix::QuestionMark), "?");
        assert_eq!(format!("{}", MmSuffix::Dot), ".");
        assert_eq!(format!("{}", MmSuffix::None), "");
    }

    #[test]
    fn suffix_str_values() {
        assert_eq!(MmSuffix::QuestionMark.suffix_str(), "?");
        assert_eq!(MmSuffix::Dot.suffix_str(), ".");
        assert_eq!(MmSuffix::None.suffix_str(), "");
    }

    #[test]
    fn serde_round_trip() {
        for (json, expected) in [
            ("\"?\"", MmSuffix::QuestionMark),
            ("\".\"", MmSuffix::Dot),
            ("\"none\"", MmSuffix::None),
        ] {
            let parsed: MmSuffix = serde_json::from_str(json).unwrap();
            assert_eq!(parsed, expected);
            let reserialized = serde_json::to_string(&expected).unwrap();
            assert_eq!(reserialized, json);
        }
    }

    #[test]
    fn serde_rejects_invalid() {
        let _: serde_json::Error = serde_json::from_str::<MmSuffix>("\"!\"").unwrap_err();
        let _: serde_json::Error = serde_json::from_str::<MmSuffix>("\"\"").unwrap_err();
        let _: serde_json::Error = serde_json::from_str::<MmSuffix>("\"NONE\"").unwrap_err();
    }
}
