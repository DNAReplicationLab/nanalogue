//! # DNA Restrictive
//!
//! Validated DNA sequence wrapper that guarantees only valid bases (A, C, G, T).
//! Provides type-safe handling of DNA sequences with compile-time validation.

use crate::Error;
use serde::{Deserialize, Serialize};
use std::str::FromStr;

/// Validated DNA sequence wrapper that guarantees only valid bases (A, C, G, T).
/// Stores sequences in uppercase.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct DNARestrictive(Vec<u8>);

impl DNARestrictive {
    /// Returns a reference to the underlying DNA sequence bytes
    pub fn get(&self) -> &[u8] {
        &self.0
    }
}

impl FromStr for DNARestrictive {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if !is_valid_dna_restrictive(s) {
            return Err(Error::InvalidSeq);
        }
        let bytes: Vec<u8> = s.bytes().map(|b| b.to_ascii_uppercase()).collect();
        Ok(DNARestrictive(bytes))
    }
}

impl<'de> Deserialize<'de> for DNARestrictive {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        DNARestrictive::from_str(&s).map_err(|e| serde::de::Error::custom(e.to_string()))
    }
}

/// Validates that a DNA sequence contains only valid bases (A, C, G, T).
/// Does not accept ambiguous bases like 'N'.
///
/// # Examples
/// ```
/// use nanalogue_core::utils::is_valid_dna_restrictive;
///
/// assert!(is_valid_dna_restrictive("ACGT"));
/// assert!(is_valid_dna_restrictive("acgt"));
/// assert!(!is_valid_dna_restrictive("ACGTN"));
/// assert!(!is_valid_dna_restrictive(""));
/// ```
pub fn is_valid_dna_restrictive(seq: &str) -> bool {
    (!seq.is_empty())
        && seq
            .bytes()
            .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tests DNARestrictive parsing with invalid barcode
    #[test]
    #[should_panic(expected = "InvalidSeq")]
    fn test_dna_restrictive_invalid() {
        let invalid_barcode = "ACGTN";
        let _ = DNARestrictive::from_str(invalid_barcode).unwrap();
    }
}
