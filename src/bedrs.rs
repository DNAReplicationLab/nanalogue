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
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Bed3<C, T> {
    /// Chromosome identifier.
    chr: C,
    /// Zero-based start coordinate.
    start: T,
    /// Exclusive end coordinate.
    end: T,
}

impl<C, T> Bed3<C, T> {
    /// Constructs an interval.
    pub const fn new(chr: C, start: T, end: T) -> Self {
        Self { chr, start, end }
    }
}

impl<C: Default, T: Default> Bed3<C, T> {
    /// Constructs an empty/default interval.
    #[must_use]
    pub fn empty() -> Self {
        Self::default()
    }
}

/// A stranded three-column BED interval.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, serde::Serialize, serde::Deserialize)]
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

impl<C, T> StrandedBed3<C, T> {
    /// Constructs a stranded interval.
    pub const fn new(chr: C, start: T, end: T, strand: Strand) -> Self {
        Self {
            chr,
            start,
            end,
            strand,
        }
    }
}

impl<C: Default, T: Default> StrandedBed3<C, T> {
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
mod tests {
    use super::{Bed3, Coordinates as _, Intersect as _, Strand, StrandedBed3};

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
    fn intersection_preserves_upstream_malformed_interval_behavior() {
        let reversed = Bed3::new(0, 10, 5);
        let containing = StrandedBed3::new(0, 0, 20, Strand::Forward);
        let interior_empty = StrandedBed3::new(0, 3, 3, Strand::Reverse);
        let containing_region = Bed3::new(0, 0, 20);

        let malformed_overlap = reversed.intersect(&containing).unwrap();
        assert_eq!(
            (malformed_overlap.start(), malformed_overlap.end()),
            (10, 5)
        );

        let empty_overlap = containing_region.intersect(&interior_empty).unwrap();
        assert_eq!((empty_overlap.start(), empty_overlap.end()), (3, 3));
    }

    #[test]
    fn strand_serde_uses_bed_symbols() {
        assert_eq!(serde_json::to_string(&Strand::Forward).unwrap(), "\"+\"");
        assert_eq!(serde_json::to_string(&Strand::Reverse).unwrap(), "\"-\"");
        assert_eq!(serde_json::to_string(&Strand::Unknown).unwrap(), "\".\"");
    }
}
