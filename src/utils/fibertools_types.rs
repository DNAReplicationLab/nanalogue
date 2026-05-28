//! # Vendored types from fibertools-rs
//!
//! This module contains data structures and helper functions adapted from the
//! published crate [`fibertools-rs` v0.8.2](https://crates.io/crates/fibertools-rs).
//! The published crate metadata declares the MIT license. We vendor only the
//! small subset of types
//! we use (`FiberAnnotation`, `FiberAnnotations`/`Ranges`, `BaseMod`,
//! `BaseMods`, `convert_seq_uppercase`, `get_u8_tag`) so that our crate
//! avoids compiling the full fibertools-rs dependency and its heavy build
//! script.
//!
//! See `THIRD_PARTY_NOTICES.md` for the full upstream license text.

use rust_htslib::bam::record::Aux;

// ---------------------------------------------------------------------------
// bamannotations types
// ---------------------------------------------------------------------------

/// A single genomic annotation with query and reference coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, PartialOrd)]
#[expect(
    clippy::exhaustive_structs,
    reason = "vendored type constructed directly in user code and doctests"
)]
pub struct FiberAnnotation {
    /// Position on the query sequence (0-based, inclusive)
    pub pos: i64,
    /// Quality / probability value (0–255)
    pub qual: u8,
    /// Position on the reference (0-based, inclusive), if mapped
    pub ref_pos: Option<i64>,
}

/// A collection of [`FiberAnnotation`] items along a single read.
#[derive(Debug, Clone, PartialEq, Eq, Ord, PartialOrd)]
#[expect(
    clippy::exhaustive_structs,
    reason = "vendored type constructed directly in user code and doctests"
)]
pub struct FiberAnnotations {
    /// Sorted annotations along the read
    pub annotations: Vec<FiberAnnotation>,
    /// Length of the query sequence
    pub seq_len: i64,
    /// Whether the read is on the reverse strand
    pub reverse: bool,
}

/// Backward-compatible alias used throughout the codebase.
pub type Ranges = FiberAnnotations;

impl FiberAnnotations {
    /// Creates `FiberAnnotations` from a vector of annotations, sorting by start position.
    #[must_use]
    pub fn from_annotations(
        mut annotations: Vec<FiberAnnotation>,
        seq_len: i64,
        reverse: bool,
    ) -> Self {
        annotations.sort_by_key(|a| a.pos);
        Self {
            annotations,
            seq_len,
            reverse,
        }
    }

    /// Query start positions.
    #[must_use]
    pub fn starts(&self) -> Vec<i64> {
        self.annotations.iter().map(|a| a.pos).collect()
    }

    /// Query end positions.
    #[must_use]
    pub fn ends(&self) -> Vec<i64> {
        self.annotations
            .iter()
            .map(|a| a.pos.saturating_add(1))
            .collect()
    }

    /// Query lengths.
    #[must_use]
    pub fn lengths(&self) -> Vec<i64> {
        self.annotations.iter().map(|_| 1).collect()
    }

    /// Quality values.
    #[must_use]
    pub fn qual(&self) -> Vec<u8> {
        self.annotations.iter().map(|a| a.qual).collect()
    }

    /// Reference start positions.
    #[must_use]
    pub fn reference_starts(&self) -> Vec<Option<i64>> {
        self.annotations.iter().map(|a| a.ref_pos).collect()
    }

    /// Reference end positions.
    #[must_use]
    pub fn reference_ends(&self) -> Vec<Option<i64>> {
        self.annotations.iter().map(|a| a.ref_pos).collect()
    }

    /// Reference lengths.
    #[must_use]
    pub fn reference_lengths(&self) -> Vec<Option<i64>> {
        self.annotations
            .iter()
            .map(|a| a.ref_pos.map(|_| 0))
            .collect()
    }
}

// ---------------------------------------------------------------------------
// basemods types
// ---------------------------------------------------------------------------

/// A single base-modification type on a read (e.g. C+m on the + strand).
#[derive(Eq, PartialEq, Debug, PartialOrd, Ord, Clone)]
#[expect(
    clippy::exhaustive_structs,
    reason = "vendored type constructed directly in user code and doctests"
)]
pub struct BaseMod {
    /// The canonical base that is modified (e.g. `b'C'`)
    pub modified_base: u8,
    /// Strand indicator (`+` or `-`)
    pub strand: char,
    /// Single-character modification code (e.g. `m` for 5mC)
    pub modification_type: char,
    /// Per-position annotations (coordinates + probabilities)
    pub ranges: Ranges,
    /// Whether the originating BAM record is reverse-complemented
    pub record_is_reverse: bool,
}

/// Collection of all base-modification types found on a single read.
#[derive(Eq, PartialEq, Debug, Clone)]
#[expect(
    clippy::exhaustive_structs,
    reason = "vendored type constructed directly in user code and doctests"
)]
pub struct BaseMods {
    /// One entry per distinct modification type
    pub base_mods: Vec<BaseMod>,
}

// ---------------------------------------------------------------------------
// bio_io helpers
// ---------------------------------------------------------------------------

/// Converts sequence bases to uppercase, leaving non-ACGTN characters unchanged.
///
/// # Examples
///
/// ```
/// use nanalogue_core::convert_seq_uppercase;
/// let input = vec![b'a', b'c', b'g', b't', b'n', b'A', b'='];
/// let output = convert_seq_uppercase(input);
/// assert_eq!(output, vec![b'A', b'C', b'G', b'T', b'N', b'A', b'=']);
/// ```
#[must_use]
pub fn convert_seq_uppercase(mut seq: Vec<u8>) -> Vec<u8> {
    for base in &mut seq {
        match *base {
            b'a' => *base = b'A',
            b'c' => *base = b'C',
            b'g' => *base = b'G',
            b't' => *base = b'T',
            b'n' => *base = b'N',
            _ => {}
        }
    }
    seq
}

/// Extracts a `u8` array auxiliary tag from a BAM record.
///
/// Returns an empty `Vec` if the tag is absent or not of the expected type.
///
/// # Examples
///
/// ```
/// use nanalogue_core::get_u8_tag;
/// use rust_htslib::bam::Record;
/// let record = Record::new();
/// assert_eq!(get_u8_tag(&record, b"ML"), Vec::<u8>::new());
/// ```
#[must_use]
pub fn get_u8_tag(record: &rust_htslib::bam::Record, tag: &[u8; 2]) -> Vec<u8> {
    if let Ok(Aux::ArrayU8(array)) = record.aux(tag) {
        array.iter().collect()
    } else {
        vec![]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn convert_seq_uppercase_mixed() {
        let input = vec![
            b'A', b'C', b'G', b'T', b'N', b'a', b'c', b'g', b't', b'n', b'=',
        ];
        let expected = vec![
            b'A', b'C', b'G', b'T', b'N', b'A', b'C', b'G', b'T', b'N', b'=',
        ];
        assert_eq!(convert_seq_uppercase(input), expected);
    }

    #[test]
    fn convert_seq_uppercase_already_upper() {
        let input = vec![b'A', b'C', b'G', b'T'];
        assert_eq!(convert_seq_uppercase(input.clone()), input);
    }

    #[test]
    fn convert_seq_uppercase_empty() {
        let input: Vec<u8> = vec![];
        assert_eq!(convert_seq_uppercase(input), Vec::<u8>::new());
    }

    #[test]
    fn get_u8_tag_missing() {
        let record = rust_htslib::bam::Record::new();
        assert_eq!(get_u8_tag(&record, b"ML"), Vec::<u8>::new());
    }

    #[test]
    fn fiber_annotations_accessors() {
        let annotations = FiberAnnotations {
            annotations: vec![
                FiberAnnotation {
                    pos: 5,
                    qual: 100,
                    ref_pos: Some(50),
                },
                FiberAnnotation {
                    pos: 20,
                    qual: 150,
                    ref_pos: None,
                },
            ],
            seq_len: 50,
            reverse: false,
        };
        assert_eq!(annotations.starts(), vec![5, 20]);
        assert_eq!(annotations.ends(), vec![6, 21]);
        assert_eq!(annotations.lengths(), vec![1, 1]);
        assert_eq!(annotations.qual(), vec![100, 150]);
        assert_eq!(annotations.reference_starts(), vec![Some(50), None]);
        assert_eq!(annotations.reference_ends(), vec![Some(50), None]);
        assert_eq!(annotations.reference_lengths(), vec![Some(0), None]);
    }

    #[test]
    fn fiber_annotations_empty() {
        let annotations = FiberAnnotations {
            annotations: vec![],
            seq_len: 0,
            reverse: false,
        };
        assert!(annotations.starts().is_empty());
        assert!(annotations.ends().is_empty());
        assert!(annotations.lengths().is_empty());
        assert!(annotations.qual().is_empty());
        assert!(annotations.reference_starts().is_empty());
        assert!(annotations.reference_ends().is_empty());
        assert!(annotations.reference_lengths().is_empty());
    }
}
