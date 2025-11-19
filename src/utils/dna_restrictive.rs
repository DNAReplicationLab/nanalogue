//! # DNA Restrictive
//!
//! Validated DNA sequence wrapper that guarantees only valid bases (A, C, G, T).
//! Provides type-safe handling of DNA sequences.

use crate::Error;
use serde::{Deserialize, Serialize};
use std::{fmt, str::FromStr};

/// Validated DNA sequence wrapper that guarantees only valid bases (A, C, G, T).
/// Stores sequences in uppercase.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct DNARestrictive(Vec<u8>);

impl DNARestrictive {
    /// Returns a reference to the underlying DNA sequence bytes
    #[must_use]
    pub fn get(&self) -> &[u8] {
        &self.0
    }
}

/// Trait that returns a `DNARestrictive` object
pub trait GetDNARestrictive {
    /// Returns a `DNARestrictive` object
    fn get_dna_restrictive(&self) -> &DNARestrictive;
}

impl GetDNARestrictive for DNARestrictive {
    fn get_dna_restrictive(&self) -> &DNARestrictive {
        self
    }
}

impl FromStr for DNARestrictive {
    type Err = Error;

    /// Convert from string with only valid bases (A, C, G, T, or lowercase).
    /// Does not accept ambiguous bases like 'N'.
    ///
    /// # Examples
    /// ```
    /// use nanalogue_core::{DNARestrictive, Error};
    /// use std::str::FromStr;
    ///
    /// let val_1 = DNARestrictive::from_str("ACGT").unwrap();
    /// let val_2 = DNARestrictive::from_str("acgt").unwrap();
    /// let val_3: Error = DNARestrictive::from_str("ACGTN").unwrap_err();
    /// let val_4: Error = DNARestrictive::from_str("").unwrap_err();
    /// ```
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        DNARestrictive::try_from(s.as_bytes().to_vec())
    }
}

impl TryFrom<Vec<u8>> for DNARestrictive {
    type Error = Error;

    /// Converts from a vector of `u8`, only `ACGT` upper or lowercases allowed.
    ///
    /// # Example
    /// ```
    /// use nanalogue_core::{DNARestrictive, Error};
    /// let val_1 = DNARestrictive::try_from(vec![b'A', b'C', b'G', b'T', b'a', b'c', b'g', b't']).unwrap();
    /// let val_2: Error = DNARestrictive::try_from(vec![b'h']).unwrap_err();
    /// let val_3: Error = DNARestrictive::try_from(vec![b'N']).unwrap_err();
    /// let val_4: Error = DNARestrictive::try_from(vec![]).unwrap_err();
    /// ```
    fn try_from(s: Vec<u8>) -> Result<Self, Self::Error> {
        if s.is_empty() {
            Err(Error::InvalidSeq("empty sequence supplied!".to_owned()))
        } else {
            Ok(DNARestrictive(
                s.into_iter()
                    .map(|x| match x {
                        b'A' | b'a' => Ok(b'A'),
                        b'C' | b'c' => Ok(b'C'),
                        b'G' | b'g' => Ok(b'G'),
                        b'T' | b't' => Ok(b'T'),
                        v => Err(Error::InvalidBase(char::from(v).to_string())),
                    })
                    .collect::<Result<Vec<u8>, _>>()?,
            ))
        }
    }
}

impl fmt::Display for DNARestrictive {
    /// Standard display function.
    ///
    /// We are fine with unsafe as we guarantee only AGCT are allowed.
    ///
    /// # Example
    /// ```
    /// use nanalogue_core::DNARestrictive;
    /// let val_1 = DNARestrictive::try_from(vec![b'A', b'C', b'G', b'T', b'a', b'c', b'g', b't']).unwrap();
    /// assert_eq!(val_1.to_string(), String::from("ACGTACGT"));
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        unsafe { String::from_utf8_unchecked(self.0.clone()).fmt(f) }
    }
}

impl<'de> Deserialize<'de> for DNARestrictive {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        DNARestrictive::from_str(&s).map_err(serde::de::Error::custom)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tests `DNARestrictive` parsing with invalid barcode
    #[test]
    #[should_panic(expected = "InvalidBase")]
    fn dna_restrictive_invalid() {
        let invalid_barcode = "ACGTN";
        let _: DNARestrictive = DNARestrictive::from_str(invalid_barcode).unwrap();
    }
}
