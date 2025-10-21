//! ReadStates struct for representing a collection of BAM alignment states
//! Handles conversion between string representation and BAM flags

use crate::Error;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::str::FromStr;

use super::read_state::ReadState;

/// Implements a collection-of-states of ReadState
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReadStates(Vec<u16>);

impl FromStr for ReadStates {
    type Err = Error;

    /// converts a comma-separated list of read states to ReadStates
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
    /// ```should_panic
    /// use nanalogue_core::{Error, ReadStates};
    /// use std::str::FromStr;
    /// let mut op = ReadStates::from_str("random")?;
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
        states.sort();
        Ok(ReadStates(states))
    }
}

impl ReadStates {
    /// Returns the flags contained within
    pub fn bam_flags(&self) -> &Vec<u16> {
        &self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_states() -> Result<(), Error> {
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
    #[should_panic]
    fn test_invalid_state() {
        let _op = ReadStates::from_str("random").unwrap();
    }
}
