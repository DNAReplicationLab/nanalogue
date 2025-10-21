//! F32AbsValBelow1 struct for constrained float between -1 and 1
//! Ensures floating-point values are within valid range at construction

use crate::Error;
use std::fmt;
use std::str::FromStr;

/// Datatype holding a float (f32) between -1 and 1 (both inclusive) guaranteed at creation.
#[derive(Debug, Clone, Default, Copy, PartialOrd, PartialEq)]
pub struct F32AbsValBelow1(f32);

impl F32AbsValBelow1 {
    /// Constructor, will fail if float is not between -1 and 1
    ///
    /// ```should_panic
    /// use nanalogue_core::Error;
    /// use nanalogue_core::F32AbsValBelow1;
    /// let x = F32AbsValBelow1::new(-1.1).unwrap();
    /// ```
    /// ```should_panic
    /// # use nanalogue_core::Error;
    /// # use nanalogue_core::F32AbsValBelow1;
    /// let x = F32AbsValBelow1::new(1.1).unwrap();
    /// ```
    /// ```
    /// # use nanalogue_core::Error;
    /// # use nanalogue_core::F32AbsValBelow1;
    /// let x = F32AbsValBelow1::new(0.1)?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    /// ```
    /// # use nanalogue_core::Error;
    /// # use nanalogue_core::F32AbsValBelow1;
    /// let x = F32AbsValBelow1::new(-0.5)?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    pub fn new(val: f32) -> Result<Self, Error> {
        if (-1.0..=1.0).contains(&val) {
            Ok(F32AbsValBelow1(val))
        } else {
            Err(Error::InvalidState("Num not b/w -1 and 1!".to_string()))
        }
    }
    /// Returns the value of the float.
    ///
    /// ```
    /// use nanalogue_core::F32AbsValBelow1;
    /// for y in vec![0.0,0.1,-0.7,1.0,-1.0]{
    ///     let x = F32AbsValBelow1::new(y.clone())?;
    ///     assert_eq!(x.val(), y);
    /// }
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    pub fn val(&self) -> f32 {
        self.0
    }
}

impl FromStr for F32AbsValBelow1 {
    type Err = Error;

    /// Parse a string to obtain float and then convert if b/w -1 and 1
    ///
    /// ```
    /// use nanalogue_core::F32AbsValBelow1;
    /// use std::str::FromStr;
    ///
    /// // Boundary values - exactly -1.0, 0.0, and 1.0 should work
    /// let neg_one = F32AbsValBelow1::from_str("-1.0")?;
    /// assert_eq!(neg_one.val(), -1.0);
    /// let zero = F32AbsValBelow1::from_str("0.0")?;
    /// assert_eq!(zero.val(), 0.0);
    /// let one = F32AbsValBelow1::from_str("1.0")?;
    /// assert_eq!(one.val(), 1.0);
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```
    /// # use nanalogue_core::F32AbsValBelow1;
    /// # use std::str::FromStr;
    /// #
    /// // Near-boundary values
    /// let near_neg_one = F32AbsValBelow1::from_str("-0.999999")?;
    /// let near_one = F32AbsValBelow1::from_str("0.999999")?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```should_panic
    /// # use nanalogue_core::F32AbsValBelow1;
    /// # use std::str::FromStr;
    /// #
    /// // Just outside boundaries should fail
    /// let outside = F32AbsValBelow1::from_str("-1.000001")?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    fn from_str(val_str: &str) -> Result<Self, Self::Err> {
        Self::new(f32::from_str(val_str)?)
    }
}

impl fmt::Display for F32AbsValBelow1 {
    /// converts to string for display.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.val())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_f32_abs_val_below1_basic() {
        // Test boundary values
        assert!(F32AbsValBelow1::new(-1.0).is_ok());
        assert!(F32AbsValBelow1::new(1.0).is_ok());
        assert!(F32AbsValBelow1::new(0.0).is_ok());

        // Test near-boundary values
        assert!(F32AbsValBelow1::new(-0.999999).is_ok());
        assert!(F32AbsValBelow1::new(0.999999).is_ok());

        // Test outside boundaries
        assert!(F32AbsValBelow1::new(-1.000001).is_err());
        assert!(F32AbsValBelow1::new(1.000001).is_err());
    }

    #[test]
    fn test_f32_abs_val_below1_from_str() {
        // Valid strings
        assert!(F32AbsValBelow1::from_str("-1.0").is_ok());
        assert!(F32AbsValBelow1::from_str("1.0").is_ok());
        assert!(F32AbsValBelow1::from_str("0.0").is_ok());
        assert!(F32AbsValBelow1::from_str("-0.5").is_ok());
        assert!(F32AbsValBelow1::from_str("0.5").is_ok());

        // Invalid strings
        assert!(F32AbsValBelow1::from_str("-1.1").is_err());
        assert!(F32AbsValBelow1::from_str("1.1").is_err());
        assert!(F32AbsValBelow1::from_str("abc").is_err());
        assert!(F32AbsValBelow1::from_str("").is_err());
    }

    #[test]
    fn test_f32_abs_val_below1_display() {
        let val = F32AbsValBelow1::new(-0.5).expect("should create");
        assert_eq!(format!("{}", val), "-0.5");

        let val = F32AbsValBelow1::new(0.75).expect("should create");
        assert_eq!(format!("{}", val), "0.75");
    }

    #[test]
    fn test_f32_abs_val_below1_val() {
        for test_val in vec![0.0, 0.1, -0.7, 1.0, -1.0] {
            let val = F32AbsValBelow1::new(test_val).expect("should create");
            assert_eq!(val.val(), test_val);
        }
    }
}
