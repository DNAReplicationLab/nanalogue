//! Contains trait for testing if a value is within some interval
//! Used by various types to implement containment checking

/// Trait for testing whether a value is contained within an interval.
pub trait Contains<T> {
    /// Returns true if the value is contained within this interval.
    fn contains(&self, val: &T) -> bool;
}
