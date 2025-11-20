//! `ReadStates` struct for representing a collection of BAM alignment states
//! Handles conversion between string representation and BAM flags

use crate::Error;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::str::FromStr;

use super::read_state::ReadState;

/// Implements a collection-of-states of `ReadState`
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReadStates(Vec<u16>);

impl FromStr for ReadStates {
    type Err = Error;

    /// converts a comma-separated list of read states to `ReadStates`
    ///
    /// # Examples
    ///
    /// ```
    /// use nanalogue_core::{Error, ReadStates};
    /// use std::str::FromStr;
    /// let mut op = ReadStates::from_str("primary_forward")?;
    /// assert_eq!(op.bam_flags(), &[0]);
    /// let mut op = ReadStates::from_str("unmapped,secondary_reverse,supplementary_reverse")?;
    /// assert_eq!(op.bam_flags(), &[4, 256 + 16, 2048 + 16]);
    /// op = ReadStates::from_str("primary_reverse")?;
    /// assert_eq!(op.bam_flags(), &[16]);
    /// op = ReadStates::from_str("primary_reverse,primary_forward")?;
    /// assert_eq!(op.bam_flags(), &[0, 16]);
    /// op = ReadStates::from_str("primary_reverse,primary_forward,primary_reverse")?;
    /// assert_eq!(op.bam_flags(), &[0, 16]);
    /// # Ok::<(), Error>(())
    /// ```
    ///
    /// Bad strings cause an error
    ///
    /// ```should_panic
    /// use nanalogue_core::{Error, ReadStates};
    /// use std::str::FromStr;
    /// let mut op = ReadStates::from_str("random")?;
    /// # Ok::<(), Error>(())
    /// ```
    ///
    /// Empty strings cause an error
    ///
    /// ```should_panic
    /// use nanalogue_core::{Error, ReadStates};
    /// use std::str::FromStr;
    /// let mut op = ReadStates::from_str("")?;
    /// # Ok::<(), Error>(())
    /// ```
    fn from_str(s: &str) -> Result<ReadStates, Self::Err> {
        let mut states = {
            let mut temp_states = HashSet::<u16>::new();
            for part in s.split(',') {
                let _: bool = temp_states.insert(u16::from(ReadState::from_str(part.trim())?));
            }
            temp_states.into_iter().collect::<Vec<u16>>()
        };
        if states.is_empty() {
            Err(Error::UnknownAlignState(
                "Set of allowed read states cannot be empty!".to_owned(),
            ))
        } else {
            states.sort_unstable();
            Ok(ReadStates(states))
        }
    }
}

impl TryFrom<Vec<u16>> for ReadStates {
    type Error = Error;

    /// Converts from a vector of `u16`, only numbers corresponding to valid [`crate::ReadState`]
    /// are permitted
    ///
    /// # Example
    /// ```
    /// use nanalogue_core::{Error, ReadStates};
    /// let val_1 = ReadStates::try_from(vec![0, 16, 256]).unwrap();
    /// let val_2: Error = ReadStates::try_from(vec![700]).unwrap_err();
    /// let val_3: Error = ReadStates::try_from(vec![16, 400]).unwrap_err();
    /// let val_4: Error = ReadStates::try_from(vec![]).unwrap_err();
    /// let val_5 = ReadStates::try_from(vec![16, 16]).unwrap();
    /// assert_eq!(val_5.bam_flags(), &[16]);
    /// ```
    fn try_from(s: Vec<u16>) -> Result<Self, Self::Error> {
        if s.is_empty() {
            Err(Error::UnknownAlignState(
                "Set of allowed read states cannot be empty!".to_owned(),
            ))
        } else {
            let temp_output = s
                .into_iter()
                .map(|x| match x {
                    v @ (0 | 4 | 16 | 256 | 272 | 2048 | 2064) => Ok(v),
                    v => Err(Error::UnknownAlignState(format!(
                        "{v} is disallowed either by the BAM format or by us"
                    ))),
                })
                .collect::<Result<HashSet<u16>, _>>()?;

            let mut output = temp_output.into_iter().collect::<Vec<u16>>();
            output.sort_unstable();
            Ok(ReadStates(output))
        }
    }
}

impl ReadStates {
    /// Returns the flags contained within
    #[must_use]
    pub fn bam_flags(&self) -> &Vec<u16> {
        &self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[expect(
        clippy::shadow_unrelated,
        reason = "repetition is fine; each block is clearly separated"
    )]
    #[test]
    fn read_states() -> Result<(), Error> {
        let op = ReadStates::from_str("primary_forward")?;
        assert_eq!(op.bam_flags(), &[0]);

        let op = ReadStates::from_str("unmapped,secondary_reverse,supplementary_reverse")?;
        assert_eq!(op.bam_flags(), &[4, 256 + 16, 2048 + 16]);

        let op = ReadStates::from_str("primary_reverse")?;
        assert_eq!(op.bam_flags(), &[16]);

        let op = ReadStates::from_str("primary_reverse,primary_forward")?;
        assert_eq!(op.bam_flags(), &[0, 16]);

        let op = ReadStates::from_str("primary_reverse,primary_forward,primary_reverse")?;
        assert_eq!(op.bam_flags(), &[0, 16]);

        Ok(())
    }

    #[test]
    #[expect(
        clippy::shadow_unrelated,
        reason = "repetition is fine; each block is clearly separated"
    )]
    fn read_states_from_vec() {
        let val = ReadStates::try_from(vec![16, 256, 0]).unwrap();
        assert_eq!(val.bam_flags(), &[0, 16, 256]);

        let _: Error = ReadStates::try_from(vec![700]).unwrap_err();
        let _: Error = ReadStates::try_from(vec![16, 400]).unwrap_err();
        let _: Error = ReadStates::try_from(vec![]).unwrap_err();

        let val = ReadStates::try_from(vec![16, 16]).unwrap();
        assert_eq!(val.bam_flags(), &[16]);
    }

    #[test]
    #[should_panic(expected = "UnknownAlignState")]
    fn invalid_string() {
        let _op = ReadStates::from_str("random").unwrap();
    }

    #[test]
    fn empty_string() {
        let _op: Error = ReadStates::from_str("").unwrap_err();
    }
}
