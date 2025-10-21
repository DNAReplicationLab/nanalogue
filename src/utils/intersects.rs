//! Intersects trait for testing if an interval intersects with self
//! Used for geometric and interval intersection testing

use std::ops::Range;

/// Implements test if an interval intersects with self
pub trait Intersects<T> {
    /// see if interval intersects with self
    fn intersects(&self, val: &T) -> bool;
}

impl Intersects<Range<u64>> for Range<u64> {
    /// Check if two ranges intersect (overlap)
    ///
    /// ```
    /// use nanalogue_core::Intersects;
    /// assert!((0..3).intersects(&(0..1)));
    /// assert!(!(0..3).intersects(&(5..7)));
    /// assert!(!(0..3).intersects(&(1..1)));
    /// assert!((1..3).intersects(&(0..2)));
    /// ```
    fn intersects(&self, val: &Range<u64>) -> bool {
        (self.start < val.end && self.end > val.start) && !(self.is_empty() || val.is_empty())
    }
}
