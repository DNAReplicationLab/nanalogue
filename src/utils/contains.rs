//! Contains trait for testing if a value is within some interval
//! Used by various types to implement containment checking

/// Implements test if a value is within some interval
pub trait Contains<T> {
    /// see if value is contained within
    fn contains(&self, val: &T) -> bool;
}
