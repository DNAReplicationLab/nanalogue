//! # Vendored types from fibertools-rs
//!
//! This module contains data structures and helper functions adapted from the
//! published crate [`fibertools-rs` v0.8.2](https://crates.io/crates/fibertools-rs).
//! The published crate metadata declares the MIT license. We vendor only the
//! small subset of types
//! we use (`FiberAnnotation`, `FiberAnnotations`/`Ranges`, `BaseMod`,
//! `BaseMods`, `convert_seq_uppercase`) so that our crate
//! avoids compiling the full fibertools-rs dependency and its heavy build
//! script.
//!
//! Unlike upstream, [`FiberAnnotation`] is bit-packed into 8 bytes; its fields
//! are private and accessed through methods.
//!
//! See `THIRD_PARTY_NOTICES.md` for the full upstream license text.

// ---------------------------------------------------------------------------
// bamannotations types
// ---------------------------------------------------------------------------

use crate::Error;
use crate::constants::shared::{MAX_CONTIG_LEN, MAX_SEQ_LEN};
use core::cmp::Ordering;
use core::num::NonZeroU32;

/// Number of low bits of [`FiberAnnotation`]'s packed word that hold the quality.
const QUAL_BITS: u32 = 8;

/// Mask selecting the quality bits of [`FiberAnnotation`]'s packed word.
const QUAL_MASK: u32 = (1 << QUAL_BITS) - 1;

/// A single genomic annotation with query and reference coordinates.
///
/// The annotation is packed into 8 bytes (two `u32` words) rather than the
/// 16 bytes the three fields would occupy side by side:
///
/// - the query position and the quality share one `u32`: the position is in
///   the upper 24 bits and the quality in the lower 8 bits, which restricts
///   query positions to [`FiberAnnotation::MAX_POS`] (see also
///   [`MAX_SEQ_LEN`]);
/// - the reference position is stored as `Option<NonZeroU32>` holding
///   `ref_pos + 1`, so that reference position 0 is stored as 1, 1 as 2 and so
///   on, and the compiler's zero niche encodes an unmapped position (`None`)
///   without any extra space. Reference positions are restricted to be below
///   [`MAX_CONTIG_LEN`], which is below `u32::MAX`, so the increment never
///   overflows.
///
/// Values are validated on construction via [`FiberAnnotation::try_new`] and
/// read back through [`FiberAnnotation::pos`], [`FiberAnnotation::qual`] and
/// [`FiberAnnotation::ref_pos`].
///
/// ```
/// use nanalogue_core::FiberAnnotation;
/// let annotation = FiberAnnotation::try_new(12, 200, Some(0))?;
/// assert_eq!(annotation.pos(), 12);
/// assert_eq!(annotation.qual(), 200);
/// assert_eq!(annotation.ref_pos(), Some(0));
/// assert_eq!(size_of::<FiberAnnotation>(), 8);
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct FiberAnnotation {
    /// Query position (upper 24 bits) and quality (lower 8 bits).
    pos_qual: u32,
    /// Reference position plus one; `None` if the query position is unmapped.
    ref_pos_plus_one: Option<NonZeroU32>,
}

const _: () = assert!(
    size_of::<FiberAnnotation>() == 8,
    "FiberAnnotation must pack into two u32 words"
);
const _: () = assert!(
    MAX_SEQ_LEN <= FiberAnnotation::MAX_POS + 1,
    "every query position below MAX_SEQ_LEN must fit in the packed position bits"
);
const _: () = assert!(
    MAX_CONTIG_LEN < u32::MAX,
    "reference positions below MAX_CONTIG_LEN must be storable as ref_pos + 1"
);

impl FiberAnnotation {
    /// Largest query position that can be stored (`2^24 - 1`).
    pub const MAX_POS: u32 = u32::MAX >> QUAL_BITS;

    /// Creates an annotation from a query position, a quality and an optional
    /// reference position.
    ///
    /// # Errors
    /// Returns [`Error::InvalidModCoords`] if `pos` exceeds
    /// [`FiberAnnotation::MAX_POS`], and [`Error::InvalidAlignCoords`] if
    /// `ref_pos` is at least [`MAX_CONTIG_LEN`].
    ///
    /// ```
    /// use nanalogue_core::{Error, FiberAnnotation};
    /// assert!(FiberAnnotation::try_new(FiberAnnotation::MAX_POS, 0, None).is_ok());
    /// assert!(matches!(
    ///     FiberAnnotation::try_new(FiberAnnotation::MAX_POS + 1, 0, None),
    ///     Err(Error::InvalidModCoords(_))
    /// ));
    /// assert!(matches!(
    ///     FiberAnnotation::try_new(0, 0, Some(u32::MAX)),
    ///     Err(Error::InvalidAlignCoords(_))
    /// ));
    /// ```
    #[inline]
    pub fn try_new(pos: u32, qual: u8, ref_pos: Option<u32>) -> Result<Self, Error> {
        if pos > Self::MAX_POS {
            return Err(Error::InvalidModCoords(format!(
                "query position {pos} exceeds maximum {}",
                Self::MAX_POS
            )));
        }
        let ref_pos_plus_one = match ref_pos {
            None => None,
            Some(v) if v < MAX_CONTIG_LEN => {
                #[expect(
                    clippy::arithmetic_side_effects,
                    reason = "v < MAX_CONTIG_LEN < u32::MAX so v + 1 cannot overflow"
                )]
                let stored = v + 1;
                NonZeroU32::new(stored)
            }
            Some(v) => {
                return Err(Error::InvalidAlignCoords(format!(
                    "reference position {v} must be below maximum contig length {MAX_CONTIG_LEN}"
                )));
            }
        };
        Ok(Self {
            pos_qual: (pos << QUAL_BITS) | u32::from(qual),
            ref_pos_plus_one,
        })
    }

    /// Position on the query sequence (0-based).
    #[inline]
    #[must_use]
    pub const fn pos(&self) -> u32 {
        self.pos_qual >> QUAL_BITS
    }

    /// Quality / probability value (0–255).
    #[inline]
    #[must_use]
    pub const fn qual(&self) -> u8 {
        // Masked to the lower 8 bits, so the cast cannot truncate.
        (self.pos_qual & QUAL_MASK) as u8
    }

    /// Position on the reference (0-based), if mapped.
    #[inline]
    #[must_use]
    pub const fn ref_pos(&self) -> Option<u32> {
        match self.ref_pos_plus_one {
            // The stored value is non-zero, so subtracting one cannot underflow.
            Some(v) => Some(v.get() - 1),
            None => None,
        }
    }
}

/// Formats the unpacked fields so that output matches the original struct.
impl core::fmt::Debug for FiberAnnotation {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("FiberAnnotation")
            .field("pos", &self.pos())
            .field("ref_pos", &self.ref_pos())
            .field("qual", &self.qual())
            .finish()
    }
}

/// Orders by query position, then reference position (unmapped first), then
/// quality, matching the field order of the original unpacked struct.
impl Ord for FiberAnnotation {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        self.pos()
            .cmp(&other.pos())
            .then_with(|| self.ref_pos_plus_one.cmp(&other.ref_pos_plus_one))
            .then_with(|| self.qual().cmp(&other.qual()))
    }
}

impl PartialOrd for FiberAnnotation {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
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
    pub seq_len: u32,
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
        seq_len: u32,
        reverse: bool,
    ) -> Self {
        annotations.sort_by_key(FiberAnnotation::pos);
        Self {
            annotations,
            seq_len,
            reverse,
        }
    }

    /// Query positions.
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn pos(&self) -> impl Iterator<Item = u32> + '_ {
        self.annotations.iter().map(FiberAnnotation::pos)
    }

    /// Reference positions.
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn ref_pos(&self) -> impl Iterator<Item = Option<u32>> + '_ {
        self.annotations.iter().map(FiberAnnotation::ref_pos)
    }

    /// Quality values.
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn qual(&self) -> impl Iterator<Item = u8> + '_ {
        self.annotations.iter().map(FiberAnnotation::qual)
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

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
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
    fn fiber_annotations_accessors() {
        let annotations = FiberAnnotations {
            annotations: vec![
                FiberAnnotation::try_new(5, 100, Some(50)).unwrap(),
                FiberAnnotation::try_new(20, 150, None).unwrap(),
            ],
            seq_len: 50,
            reverse: false,
        };
        assert_eq!(annotations.pos().collect::<Vec<_>>(), vec![5, 20]);
        assert_eq!(annotations.qual().collect::<Vec<_>>(), vec![100, 150]);
        assert_eq!(
            annotations.ref_pos().collect::<Vec<_>>(),
            vec![Some(50), None]
        );
    }

    #[test]
    fn fiber_annotation_is_eight_bytes() {
        assert_eq!(size_of::<FiberAnnotation>(), 8);
        assert_eq!(align_of::<FiberAnnotation>(), 4);
        assert_eq!(size_of::<Option<NonZeroU32>>(), 4);
    }

    #[test]
    fn fiber_annotation_round_trips_boundaries() {
        for pos in [0, 1, 255, 256, MAX_SEQ_LEN - 1, FiberAnnotation::MAX_POS] {
            for qual in [0, 1, 127, 128, 254, 255] {
                for ref_pos in [None, Some(0), Some(1), Some(MAX_CONTIG_LEN - 1)] {
                    let annotation = FiberAnnotation::try_new(pos, qual, ref_pos).unwrap();
                    assert_eq!(annotation.pos(), pos, "pos for {pos}, {qual}, {ref_pos:?}");
                    assert_eq!(
                        annotation.qual(),
                        qual,
                        "qual for {pos}, {qual}, {ref_pos:?}"
                    );
                    assert_eq!(
                        annotation.ref_pos(),
                        ref_pos,
                        "ref_pos for {pos}, {qual}, {ref_pos:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn fiber_annotation_ref_pos_zero_is_distinct_from_unmapped() {
        let mapped = FiberAnnotation::try_new(3, 9, Some(0)).unwrap();
        let unmapped = FiberAnnotation::try_new(3, 9, None).unwrap();
        assert_ne!(mapped, unmapped);
        assert_eq!(mapped.ref_pos(), Some(0));
        assert_eq!(unmapped.ref_pos(), None);
    }

    #[test]
    fn fiber_annotation_rejects_out_of_range_values() {
        assert!(matches!(
            FiberAnnotation::try_new(FiberAnnotation::MAX_POS + 1, 0, None),
            Err(Error::InvalidModCoords(_))
        ));
        assert!(matches!(
            FiberAnnotation::try_new(u32::MAX, 0, None),
            Err(Error::InvalidModCoords(_))
        ));
        assert!(matches!(
            FiberAnnotation::try_new(0, 0, Some(MAX_CONTIG_LEN)),
            Err(Error::InvalidAlignCoords(_))
        ));
        assert!(matches!(
            FiberAnnotation::try_new(0, 0, Some(u32::MAX)),
            Err(Error::InvalidAlignCoords(_))
        ));
    }

    #[test]
    fn fiber_annotation_orders_like_the_unpacked_struct() {
        // The unpacked struct derived `Ord` over (pos, ref_pos, qual); a naive
        // derive over the packed word would order by (pos, qual, ref_pos).
        let tuples = [
            (5, 200, Some(1)),
            (5, 10, Some(2)),
            (5, 255, None),
            (4, 0, Some(9)),
            (5, 10, None),
            (6, 0, Some(0)),
            (5, 200, Some(2)),
        ];
        let mut expected = tuples.map(|(pos, qual, ref_pos)| (pos, ref_pos, qual));
        expected.sort_unstable();
        let mut annotations = tuples
            .map(|(pos, qual, ref_pos)| FiberAnnotation::try_new(pos, qual, ref_pos).unwrap());
        annotations.sort_unstable();
        assert_eq!(
            annotations.map(|a| (a.pos(), a.ref_pos(), a.qual())),
            expected
        );
    }

    #[test]
    fn fiber_annotation_debug_shows_unpacked_fields() {
        let annotation = FiberAnnotation::try_new(7, 42, Some(99)).unwrap();
        assert_eq!(
            format!("{annotation:?}"),
            "FiberAnnotation { pos: 7, ref_pos: Some(99), qual: 42 }"
        );
    }

    #[test]
    fn fiber_annotations_empty() {
        let annotations = FiberAnnotations {
            annotations: vec![],
            seq_len: 0,
            reverse: false,
        };
        assert!(annotations.pos().next().is_none());
        assert!(annotations.qual().next().is_none());
        assert!(annotations.ref_pos().next().is_none());
    }
}
