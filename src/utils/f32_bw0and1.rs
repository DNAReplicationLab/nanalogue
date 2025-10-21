//! F32Bw0and1 struct for constrained float between 0 and 1
//! Ensures floating-point values are within valid range at construction

use super::f32_abs_val_below1::F32AbsValBelow1;
use crate::Error;
use serde::{Deserialize, Serialize};
use std::convert::From;
use std::fmt;
use std::str::FromStr;

/// Datatype holding a float (f32) between 0 and 1 (both inclusive) guaranteed at creation.
#[derive(Debug, Clone, Default, Copy, PartialOrd, PartialEq, Serialize, Deserialize)]
pub struct F32Bw0and1(f32);

impl F32Bw0and1 {
    /// Constructor, will fail if float is not between 0 and 1
    ///
    /// ```should_panic
    /// use nanalogue_core::Error;
    /// use nanalogue_core::F32Bw0and1;
    /// let x = F32Bw0and1::new(-0.1).unwrap();
    /// ```
    /// ```
    /// # use nanalogue_core::Error;
    /// # use nanalogue_core::F32Bw0and1;
    /// let x = F32Bw0and1::new(0.1)?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    pub fn new(val: f32) -> Result<Self, Error> {
        if (0.0..=1.0).contains(&val) {
            Ok(F32Bw0and1(val))
        } else {
            Err(Error::InvalidState("Num not b/w 0 and 1!".to_string()))
        }
    }
    /// Returns the value of the float.
    ///
    /// ```
    /// use nanalogue_core::F32Bw0and1;
    /// for y in vec![0.0,0.1,0.7,1.0]{
    ///     let x = F32Bw0and1::new(y.clone())?;
    ///     assert_eq!(x.val(), y);
    /// }
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    pub fn val(&self) -> f32 {
        self.0
    }
    /// Shortcut for 1.0
    pub fn one() -> Self {
        F32Bw0and1::new(1.0).expect("no error")
    }
    /// Shortcut for 0.0
    pub fn zero() -> Self {
        F32Bw0and1::new(0.0).expect("no error")
    }
    /// Converts from F32AbsValBelow1 using the absolute value
    pub fn abs_f32_abs_val_below_1(val: F32AbsValBelow1) -> Self {
        F32Bw0and1::new(f32::abs(val.val())).expect("no error")
    }
}

impl FromStr for F32Bw0and1 {
    type Err = Error;

    /// Parse a string to obtain float and then convert if b/w 0 and 1
    ///
    /// ```
    /// use nanalogue_core::F32Bw0and1;
    /// use std::str::FromStr;
    ///
    /// // Boundary values - exactly 0.0 and 1.0 should work
    /// let zero = F32Bw0and1::from_str("0.0")?;
    /// assert_eq!(zero.val(), 0.0);
    /// let one = F32Bw0and1::from_str("1.0")?;
    /// assert_eq!(one.val(), 1.0);
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```
    /// # use nanalogue_core::F32Bw0and1;
    /// # use std::str::FromStr;
    /// #
    /// // Near-boundary values
    /// let near_zero = F32Bw0and1::from_str("0.000001")?;
    /// let near_one = F32Bw0and1::from_str("0.999999")?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```should_panic
    /// # use nanalogue_core::F32Bw0and1;
    /// # use std::str::FromStr;
    /// #
    /// // Just outside boundaries should fail
    /// let outside = F32Bw0and1::from_str("1.000001").unwrap();
    /// ```
    fn from_str(val_str: &str) -> Result<Self, Self::Err> {
        Self::new(f32::from_str(val_str)?)
    }
}

impl From<u8> for F32Bw0and1 {
    /// Convert from a u8 i.e. a number >= 0 and <= 255
    fn from(value: u8) -> Self {
        F32Bw0and1::new((value as f32) / (u8::MAX as f32 + 1.0)).expect("no F32 conversion error")
    }
}

impl fmt::Display for F32Bw0and1 {
    /// converts to string for display.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.val())
    }
}

impl From<F32Bw0and1> for F32AbsValBelow1 {
    /// Convert between the two types of floats
    fn from(value: F32Bw0and1) -> Self {
        F32AbsValBelow1::new(value.val()).expect("no F32 conversion error")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_f32_bw0and1_basic() {
        // Test boundary values
        assert!(F32Bw0and1::new(0.0).is_ok());
        assert!(F32Bw0and1::new(1.0).is_ok());

        // Test near-boundary values
        assert!(F32Bw0and1::new(0.000001).is_ok());
        assert!(F32Bw0and1::new(0.999999).is_ok());

        // Test outside boundaries
        assert!(F32Bw0and1::new(-0.000001).is_err());
        assert!(F32Bw0and1::new(1.000001).is_err());
    }

    #[test]
    fn test_f32_bw0and1_from_str() {
        // Valid strings
        assert!(F32Bw0and1::from_str("0.0").is_ok());
        assert!(F32Bw0and1::from_str("1.0").is_ok());
        assert!(F32Bw0and1::from_str("0.5").is_ok());

        // Invalid strings
        assert!(F32Bw0and1::from_str("-0.1").is_err());
        assert!(F32Bw0and1::from_str("1.1").is_err());
        assert!(F32Bw0and1::from_str("abc").is_err());
        assert!(F32Bw0and1::from_str("").is_err());
    }

    #[test]
    fn test_f32_bw0and1_shortcuts() {
        let zero = F32Bw0and1::zero();
        assert_eq!(zero.val(), 0.0);

        let one = F32Bw0and1::one();
        assert_eq!(one.val(), 1.0);
    }

    #[test]
    fn test_f32_bw0and1_from_u8() {
        // Test boundary values
        let zero = F32Bw0and1::from(0u8);
        assert_eq!(zero.val(), 0.0);

        let max = F32Bw0and1::from(255u8);
        assert!(max.val() < 1.0); // Should be 255/256 = 0.99609375

        // Test some intermediate values
        let half = F32Bw0and1::from(128u8);
        assert!(half.val() > 0.49 && half.val() < 0.51);

        // Test exact calculation
        let test_val = 100u8;
        let converted = F32Bw0and1::from(test_val);
        let expected = (test_val as f32) / (u8::MAX as f32 + 1.0);
        assert_eq!(converted.val(), expected);
    }

    #[test]
    fn test_f32_types_integration() {
        // Test conversion from F32Bw0and1 to F32AbsValBelow1
        let pos_values = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        for val in pos_values {
            let bw_val = F32Bw0and1::new(val).expect("should create");
            let abs_val: F32AbsValBelow1 = bw_val.into();
            assert_eq!(abs_val.val(), val);
        }

        // Test conversion via absolute value function
        let neg_val = F32AbsValBelow1::new(-0.5).expect("should create");
        let abs_converted = F32Bw0and1::abs_f32_abs_val_below_1(neg_val);
        assert_eq!(abs_converted.val(), 0.5);

        let pos_val = F32AbsValBelow1::new(0.7).expect("should create");
        let abs_converted = F32Bw0and1::abs_f32_abs_val_below_1(pos_val);
        assert_eq!(abs_converted.val(), 0.7);
    }

    #[test]
    fn test_f32_bw0and1_display() {
        let val = F32Bw0and1::new(0.5).expect("should create");
        assert_eq!(format!("{}", val), "0.5");

        let val = F32Bw0and1::new(0.75).expect("should create");
        assert_eq!(format!("{}", val), "0.75");
    }
}
