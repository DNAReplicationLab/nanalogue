//! ReadState enum for representing BAM alignment states
//! Handles conversion between internal representation and BAM flags

use crate::Error;
use rand::Rng;
use rand::distr::{Distribution, StandardUniform};
use serde::{Deserialize, Serialize};
use std::convert::TryFrom;
use std::fmt;
use std::str::FromStr;

/// Alignment state of a read; seven possibilities + one unknown state
#[derive(Debug, Clone, Default, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReadState {
    #[default]
    /// Primary alignment to the reference strand
    #[serde(rename = "primary_forward")]
    PrimaryFwd,
    /// Primary alignment opposite the reference strand
    #[serde(rename = "primary_reverse")]
    PrimaryRev,
    /// Secondary alignment to the reference strand
    #[serde(rename = "secondary_forward")]
    SecondaryFwd,
    /// Secondary alignment opposite the reference strand
    #[serde(rename = "secondary_reverse")]
    SecondaryRev,
    /// Supplementary alignment to the reference strand
    #[serde(rename = "supplementary_forward")]
    SupplementaryFwd,
    /// Supplementary alignment opposite the reference strand
    #[serde(rename = "supplementary_reverse")]
    SupplementaryRev,
    /// Marked as unmapped in the BAM file. We are assuming
    /// that unmapped sequences will not be stored as reversed
    /// complements, as what would be the point of that?
    #[serde(rename = "unmapped")]
    Unmapped,
}

// Implements random pick of a variant
impl Distribution<ReadState> for StandardUniform {
    /// Allows us to randomly pick a variant
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> ReadState {
        match rng.random_range(0..7) {
            0 => ReadState::PrimaryFwd,
            1 => ReadState::PrimaryRev,
            2 => ReadState::SecondaryFwd,
            3 => ReadState::SecondaryRev,
            4 => ReadState::SupplementaryFwd,
            5 => ReadState::SupplementaryRev,
            6 => ReadState::Unmapped,
            _ => unreachable!(),
        }
    }
}

// Implements conversion of ReadState into the standard BAM flag format
impl TryFrom<ReadState> for u16 {
    type Error = Error;
    /// converts our internal representation to the BAM flag format
    fn try_from(value: ReadState) -> Result<u16, Error> {
        match value {
            ReadState::PrimaryFwd => Ok(0),
            ReadState::Unmapped => Ok(4),
            ReadState::PrimaryRev => Ok(16),
            ReadState::SecondaryFwd => Ok(256),
            ReadState::SecondaryRev => Ok(272),
            ReadState::SupplementaryFwd => Ok(2048),
            ReadState::SupplementaryRev => Ok(2064),
        }
    }
}

// Implements conversion of the standard BAM flag format into ReadState
impl TryFrom<u16> for ReadState {
    type Error = Error;
    /// converts our internal representation to the BAM flag format
    fn try_from(value: u16) -> Result<ReadState, Error> {
        match value {
            0 => Ok(ReadState::PrimaryFwd),
            4 => Ok(ReadState::Unmapped),
            16 => Ok(ReadState::PrimaryRev),
            256 => Ok(ReadState::SecondaryFwd),
            272 => Ok(ReadState::SecondaryRev),
            2048 => Ok(ReadState::SupplementaryFwd),
            2064 => Ok(ReadState::SupplementaryRev),
            _ => Err(Error::UnknownAlignState),
        }
    }
}

/// Implements from string for ReadState
///
/// ```
/// use nanalogue_core::ReadState;
/// use std::str::FromStr;
///
/// // Primary alignments
/// let state = ReadState::from_str("primary_forward")?;
/// assert_eq!(state, ReadState::PrimaryFwd);
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
///
/// ```
/// # use nanalogue_core::ReadState;
/// # use std::str::FromStr;
/// #
/// // Secondary alignments
/// let state = ReadState::from_str("secondary_reverse")?;
/// assert_eq!(state, ReadState::SecondaryRev);
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
///
/// ```
/// # use nanalogue_core::ReadState;
/// # use std::str::FromStr;
/// #
/// // Supplementary alignments
/// let state = ReadState::from_str("supplementary_forward")?;
/// assert_eq!(state, ReadState::SupplementaryFwd);
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
///
/// ```
/// # use nanalogue_core::ReadState;
/// # use std::str::FromStr;
/// #
/// // Unmapped reads
/// let state = ReadState::from_str("unmapped")?;
/// assert_eq!(state, ReadState::Unmapped);
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
///
/// ```should_panic
/// # use nanalogue_core::ReadState;
/// # use std::str::FromStr;
/// #
/// // Invalid string should error
/// let state = ReadState::from_str("invalid_state")?;
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
impl FromStr for ReadState {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "primary_forward" => Ok(ReadState::PrimaryFwd),
            "primary_reverse" => Ok(ReadState::PrimaryRev),
            "secondary_forward" => Ok(ReadState::SecondaryFwd),
            "secondary_reverse" => Ok(ReadState::SecondaryRev),
            "supplementary_forward" => Ok(ReadState::SupplementaryFwd),
            "supplementary_reverse" => Ok(ReadState::SupplementaryRev),
            "unmapped" => Ok(ReadState::Unmapped),
            _ => Err(Error::UnknownAlignState),
        }
    }
}

/// Implements printing of read state
impl fmt::Display for ReadState {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let printable = match *self {
            ReadState::PrimaryFwd => "primary_forward",
            ReadState::SecondaryFwd => "secondary_forward",
            ReadState::SupplementaryFwd => "supplementary_forward",
            ReadState::PrimaryRev => "primary_reverse",
            ReadState::SecondaryRev => "secondary_reverse",
            ReadState::SupplementaryRev => "supplementary_reverse",
            ReadState::Unmapped => "unmapped",
        };
        write!(f, "{printable}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tests ReadState u16 conversion round-trip
    #[test]
    fn test_readstate_u16_conversion_roundtrip() {
        let states = vec![
            ReadState::PrimaryFwd,
            ReadState::PrimaryRev,
            ReadState::SecondaryFwd,
            ReadState::SecondaryRev,
            ReadState::SupplementaryFwd,
            ReadState::SupplementaryRev,
            ReadState::Unmapped,
        ];

        for state in states {
            // Convert to u16 and back
            let flag: u16 = state.try_into().expect("conversion to u16 should work");
            let recovered_state: ReadState =
                flag.try_into().expect("conversion from u16 should work");
            assert_eq!(state, recovered_state);
        }
    }

    /// Tests specific ReadState u16 flag values
    #[test]
    fn test_readstate_specific_flag_values() {
        assert_eq!(u16::try_from(ReadState::PrimaryFwd).unwrap(), 0);
        assert_eq!(u16::try_from(ReadState::Unmapped).unwrap(), 4);
        assert_eq!(u16::try_from(ReadState::PrimaryRev).unwrap(), 16);
        assert_eq!(u16::try_from(ReadState::SecondaryFwd).unwrap(), 256);
        assert_eq!(u16::try_from(ReadState::SecondaryRev).unwrap(), 272);
        assert_eq!(u16::try_from(ReadState::SupplementaryFwd).unwrap(), 2048);
        assert_eq!(u16::try_from(ReadState::SupplementaryRev).unwrap(), 2064);
    }

    /// Tests ReadState from invalid u16 values
    #[test]
    fn test_readstate_invalid_u16_values() {
        // Test various invalid flag combinations
        let invalid_flags = vec![1, 2, 8, 32, 64, 128, 512, 1024, 4096, 8192];
        for flag in invalid_flags {
            assert!(matches!(
                ReadState::try_from(flag),
                Err(Error::UnknownAlignState)
            ));
        }
    }

    /// Tests ReadState string parsing and display consistency
    #[test]
    fn test_readstate_string_consistency() {
        let states = vec![
            ReadState::PrimaryFwd,
            ReadState::PrimaryRev,
            ReadState::SecondaryFwd,
            ReadState::SecondaryRev,
            ReadState::SupplementaryFwd,
            ReadState::SupplementaryRev,
            ReadState::Unmapped,
        ];

        for state in states {
            let string_repr = format!("{}", state);
            let parsed_state = ReadState::from_str(&string_repr).expect("should parse");
            assert_eq!(state, parsed_state);
        }
    }

    /// Tests ReadState from_str error cases
    #[test]
    fn test_readstate_from_str_errors() {
        assert!(matches!(
            ReadState::from_str("invalid_state"),
            Err(Error::UnknownAlignState)
        ));
        assert!(matches!(
            ReadState::from_str(""),
            Err(Error::UnknownAlignState)
        ));
        assert!(matches!(
            ReadState::from_str("primary"), // Incomplete
            Err(Error::UnknownAlignState)
        ));
    }
}
