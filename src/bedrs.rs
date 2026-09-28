//! Minimal compatibility subset adapted from `bedrs` 0.2.26.

/// A genomic strand.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, serde::Serialize, serde::Deserialize)]
#[expect(
    clippy::exhaustive_enums,
    reason = "public compatibility with bedrs 0.2.26"
)]
pub enum Strand {
    /// Forward (`+`) strand.
    #[serde(rename = "+")]
    Forward,
    /// Reverse (`-`) strand.
    #[serde(rename = "-")]
    Reverse,
    /// Unknown (`.`) strand.
    #[default]
    #[serde(rename = ".")]
    Unknown,
}

/// A three-column BED interval.
#[derive(Clone, Copy, Debug, Eq, PartialEq, serde::Serialize)]
pub struct Bed3<C, T> {
    /// Chromosome identifier.
    chr: C,
    /// Zero-based start coordinate.
    start: T,
    /// Exclusive end coordinate.
    end: T,
}

impl<C, T: PartialOrd> Bed3<C, T> {
    /// Constructs an interval.
    ///
    /// # Panics
    /// Panics if `start` is not less than or equal to `end`.
    pub fn new(chr: C, start: T, end: T) -> Self {
        assert!(start <= end, "interval start cannot exceed end");
        Self { chr, start, end }
    }
}

impl<C: Default, T: Default + PartialOrd> Default for Bed3<C, T> {
    fn default() -> Self {
        Self::new(C::default(), T::default(), T::default())
    }
}

impl<'de, C, T> serde::Deserialize<'de> for Bed3<C, T>
where
    C: serde::Deserialize<'de>,
    T: serde::Deserialize<'de> + PartialOrd,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        #[derive(serde::Deserialize)]
        #[serde(rename = "Bed3")]
        struct SerializedBed3<C, T> {
            chr: C,
            start: T,
            end: T,
        }

        let serialized = SerializedBed3::deserialize(deserializer)?;
        if serialized.start <= serialized.end {
            Ok(Self {
                chr: serialized.chr,
                start: serialized.start,
                end: serialized.end,
            })
        } else {
            Err(serde::de::Error::custom("interval start cannot exceed end"))
        }
    }
}

impl<C: Default, T: Default + PartialOrd> Bed3<C, T> {
    /// Constructs an empty/default interval.
    #[must_use]
    pub fn empty() -> Self {
        Self::default()
    }
}

/// A stranded three-column BED interval.
#[derive(Clone, Copy, Debug, Eq, PartialEq, serde::Serialize)]
pub struct StrandedBed3<C, T> {
    /// Chromosome identifier.
    chr: C,
    /// Zero-based start coordinate.
    start: T,
    /// Exclusive end coordinate.
    end: T,
    /// Genomic strand.
    strand: Strand,
}

impl<C, T: PartialOrd> StrandedBed3<C, T> {
    /// Constructs a stranded interval.
    ///
    /// # Panics
    /// Panics if `start` is not less than or equal to `end`.
    pub fn new(chr: C, start: T, end: T, strand: Strand) -> Self {
        assert!(start <= end, "interval start cannot exceed end");
        Self {
            chr,
            start,
            end,
            strand,
        }
    }
}

impl<C: Default, T: Default + PartialOrd> Default for StrandedBed3<C, T> {
    fn default() -> Self {
        Self::new(C::default(), T::default(), T::default(), Strand::default())
    }
}

impl<'de, C, T> serde::Deserialize<'de> for StrandedBed3<C, T>
where
    C: serde::Deserialize<'de>,
    T: serde::Deserialize<'de> + PartialOrd,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        #[derive(serde::Deserialize)]
        #[serde(rename = "StrandedBed3")]
        struct SerializedStrandedBed3<C, T> {
            chr: C,
            start: T,
            end: T,
            strand: Strand,
        }

        let serialized = SerializedStrandedBed3::deserialize(deserializer)?;
        if serialized.start <= serialized.end {
            Ok(Self {
                chr: serialized.chr,
                start: serialized.start,
                end: serialized.end,
                strand: serialized.strand,
            })
        } else {
            Err(serde::de::Error::custom("interval start cannot exceed end"))
        }
    }
}

impl<C: Default, T: Default + PartialOrd> StrandedBed3<C, T> {
    /// Constructs an empty/default interval.
    #[must_use]
    pub fn empty() -> Self {
        Self::default()
    }
}

/// Shared access to interval coordinates.
pub trait Coordinates {
    /// Chromosome identifier type.
    type Chrom;
    /// Coordinate type.
    type Coord: Copy;

    /// Returns the chromosome identifier.
    fn chr(&self) -> &Self::Chrom;
    /// Returns the zero-based start coordinate.
    fn start(&self) -> Self::Coord;
    /// Returns the exclusive end coordinate.
    fn end(&self) -> Self::Coord;
    /// Returns the strand, if represented by this record.
    fn strand(&self) -> Option<Strand>;
}

impl<C, T: Copy> Coordinates for Bed3<C, T> {
    type Chrom = C;
    type Coord = T;

    fn chr(&self) -> &C {
        &self.chr
    }
    fn start(&self) -> T {
        self.start
    }
    fn end(&self) -> T {
        self.end
    }
    fn strand(&self) -> Option<Strand> {
        None
    }
}

impl<C, T: Copy> Coordinates for StrandedBed3<C, T> {
    type Chrom = C;
    type Coord = T;

    fn chr(&self) -> &C {
        &self.chr
    }
    fn start(&self) -> T {
        self.start
    }
    fn end(&self) -> T {
        self.end
    }
    fn strand(&self) -> Option<Strand> {
        Some(self.strand)
    }
}

/// Computes the half-open intersection of two genomic intervals.
pub trait Intersect<Rhs = Self> {
    /// Intersection record type.
    type Output;
    /// Returns the overlap, or `None` for different chromosomes or touching/disjoint intervals.
    fn intersect(&self, other: &Rhs) -> Option<Self::Output>;
}

impl<C: Clone + PartialEq, T: Copy + Ord> Intersect<StrandedBed3<C, T>> for Bed3<C, T> {
    type Output = StrandedBed3<C, T>;

    fn intersect(&self, other: &StrandedBed3<C, T>) -> Option<Self::Output> {
        let start = self.start.max(other.start);
        let end = self.end.min(other.end);
        (self.chr == other.chr && self.start < other.end && self.end > other.start)
            .then(|| StrandedBed3::new(self.chr.clone(), start, end, other.strand))
    }
}

/// Commonly used record and operation imports.
pub mod prelude {
    pub use super::{Bed3, Coordinates, Intersect, Strand, StrandedBed3};
}

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use super::{Bed3, Coordinates as _, Intersect as _, Strand, StrandedBed3};
    use std::sync::atomic::{AtomicI32, Ordering};

    // Supplies successive coordinate defaults of 1, then 0. This models a valid,
    // though unusual, `Default` implementation whose return value changes between calls.
    static NEXT_DEFAULT_COORD: AtomicI32 = AtomicI32::new(1);

    /// A test coordinate whose default counts down on each call.
    ///
    /// Constructing a default interval calls `T::default()` separately for its start
    /// and end, so this type attempts to produce the reversed interval `1..0`.
    #[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
    struct ChangingDefaultCoord(i32);

    impl Default for ChangingDefaultCoord {
        fn default() -> Self {
            Self(NEXT_DEFAULT_COORD.fetch_sub(1, Ordering::Relaxed))
        }
    }

    #[test]
    fn bed3_serde_matches_upstream() {
        let bed = Bed3::new(0, 0, 10);
        assert_eq!(
            serde_json::to_string(&bed).unwrap(),
            r#"{"chr":0,"start":0,"end":10}"#
        );
        assert_eq!(
            serde_json::from_str::<Bed3<i32, u32>>("[0,0,10]").unwrap(),
            bed
        );
    }

    #[test]
    fn intersection_is_half_open() {
        let bed = Bed3::new(0, 5, 10);
        let overlapping = StrandedBed3::new(0, 8, 12, Strand::Forward);
        let touching = StrandedBed3::new(0, 10, 12, Strand::Reverse);
        let other_chromosome = StrandedBed3::new(1, 8, 12, Strand::Unknown);

        let overlap = bed.intersect(&overlapping).unwrap();
        assert_eq!((overlap.start(), overlap.end()), (8, 10));
        assert_eq!(overlap.strand(), Some(Strand::Forward));
        assert_eq!(bed.intersect(&touching), None);
        assert_eq!(bed.intersect(&other_chromosome), None);
    }

    #[test]
    fn intersection_preserves_empty_interval() {
        let interior_empty = StrandedBed3::new(0, 3, 3, Strand::Reverse);
        let containing_region = Bed3::new(0, 0, 20);

        let empty_overlap = containing_region.intersect(&interior_empty).unwrap();
        assert_eq!((empty_overlap.start(), empty_overlap.end()), (3, 3));
    }

    #[test]
    #[should_panic(expected = "interval start cannot exceed end")]
    fn bed3_rejects_reversed_coordinates() {
        let _: Bed3<i32, u32> = Bed3::new(0, 10, 5);
    }

    #[test]
    fn bed3_deserialization_rejects_reversed_coordinates() {
        let err = serde_json::from_str::<Bed3<i32, u32>>("[0,10,5]").unwrap_err();
        assert!(err.to_string().contains("interval start cannot exceed end"));
    }

    #[test]
    #[should_panic(expected = "interval start cannot exceed end")]
    fn stranded_bed3_rejects_reversed_coordinates() {
        let _: StrandedBed3<i32, u32> = StrandedBed3::new(0, 10, 5, Strand::Forward);
    }

    #[test]
    fn stranded_bed3_deserialization_rejects_reversed_coordinates() {
        let err = serde_json::from_str::<StrandedBed3<i32, u32>>(r#"[0,10,5,"+"]"#).unwrap_err();
        assert!(err.to_string().contains("interval start cannot exceed end"));
    }

    #[test]
    fn defaults_validate_coordinates() {
        // Reset the counter so the start defaults to 1 and the end defaults to 0.
        // The resulting panic proves that `Bed3::default()` routes through `new`
        // rather than bypassing its start <= end invariant.
        NEXT_DEFAULT_COORD.store(1, Ordering::Relaxed);
        let _bed_panic =
            std::panic::catch_unwind(Bed3::<(), ChangingDefaultCoord>::default).unwrap_err();

        // Verify the same invariant for the stranded interval type.
        NEXT_DEFAULT_COORD.store(1, Ordering::Relaxed);
        let _stranded_bed_panic =
            std::panic::catch_unwind(StrandedBed3::<(), ChangingDefaultCoord>::default)
                .unwrap_err();
    }

    #[test]
    fn strand_serde_uses_bed_symbols() {
        assert_eq!(serde_json::to_string(&Strand::Forward).unwrap(), "\"+\"");
        assert_eq!(serde_json::to_string(&Strand::Reverse).unwrap(), "\"-\"");
        assert_eq!(serde_json::to_string(&Strand::Unknown).unwrap(), "\".\"");
    }
}
