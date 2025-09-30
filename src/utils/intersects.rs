//! Intersects trait for testing if an interval intersects with self
//! Used for geometric and interval intersection testing

/// Implements test if an interval intersects with self
pub trait Intersects<T> {
    /// see if interval intersects with self
    fn intersects(&self, val: &T) -> bool;
}
