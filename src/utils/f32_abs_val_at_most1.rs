//! `F32AbsValAtMost1` struct for constrained float between -1 and 1
//! Ensures floating-point values are within valid range at construction

use crate::Error;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

/// Datatype holding a float (f32) between -1 and 1 (both inclusive) guaranteed at creation.
#[derive(Debug, Clone, Default, Copy, PartialOrd, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "f32")]
pub struct F32AbsValAtMost1(f32);

impl F32AbsValAtMost1 {
    /// Constructor, will fail if float is not between -1 and 1
    ///
    /// ```should_panic
    /// use nanalogue_core::Error;
    /// use nanalogue_core::F32AbsValAtMost1;
    /// let x = F32AbsValAtMost1::new(-1.1).unwrap();
    /// ```
    /// ```should_panic
    /// # use nanalogue_core::Error;
    /// # use nanalogue_core::F32AbsValAtMost1;
    /// let x = F32AbsValAtMost1::new(1.1).unwrap();
    /// ```
    /// ```
    /// # use nanalogue_core::Error;
    /// # use nanalogue_core::F32AbsValAtMost1;
    /// let x = F32AbsValAtMost1::new(0.1)?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    /// ```
    /// # use nanalogue_core::Error;
    /// # use nanalogue_core::F32AbsValAtMost1;
    /// let x = F32AbsValAtMost1::new(-0.5)?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// # Errors
    /// Returns an error if the value is not between -1.0 and 1.0 (inclusive).
    pub fn new(val: f32) -> Result<Self, Error> {
        if (-1.0..=1.0).contains(&val) {
            Ok(F32AbsValAtMost1(val))
        } else {
            Err(Error::InvalidState("Num not b/w -1 and 1!".to_string()))
        }
    }
    /// Returns the value of the float.
    ///
    /// ```
    /// use nanalogue_core::F32AbsValAtMost1;
    /// for y in vec![0.0,0.1,-0.7,1.0,-1.0]{
    ///     let x = F32AbsValAtMost1::new(y.clone())?;
    ///     assert_eq!(x.val(), y);
    /// }
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    #[must_use]
    pub fn val(&self) -> f32 {
        self.0
    }
}

impl FromStr for F32AbsValAtMost1 {
    type Err = Error;

    /// Parse a string to obtain float and then convert if b/w -1 and 1
    ///
    /// ```
    /// use nanalogue_core::F32AbsValAtMost1;
    /// use std::str::FromStr;
    ///
    /// // Boundary values - exactly -1.0, 0.0, and 1.0 should work
    /// let neg_one = F32AbsValAtMost1::from_str("-1.0")?;
    /// assert_eq!(neg_one.val(), -1.0);
    /// let zero = F32AbsValAtMost1::from_str("0.0")?;
    /// assert_eq!(zero.val(), 0.0);
    /// let one = F32AbsValAtMost1::from_str("1.0")?;
    /// assert_eq!(one.val(), 1.0);
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```
    /// # use nanalogue_core::F32AbsValAtMost1;
    /// # use std::str::FromStr;
    /// #
    /// // Near-boundary values
    /// let near_neg_one = F32AbsValAtMost1::from_str("-0.999999")?;
    /// let near_one = F32AbsValAtMost1::from_str("0.999999")?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```should_panic
    /// # use nanalogue_core::F32AbsValAtMost1;
    /// # use std::str::FromStr;
    /// #
    /// // Just outside boundaries should fail
    /// let outside = F32AbsValAtMost1::from_str("-1.000001")?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    fn from_str(val_str: &str) -> Result<Self, Self::Err> {
        Self::new(f32::from_str(val_str)?)
    }
}

impl fmt::Display for F32AbsValAtMost1 {
    /// converts to string for display.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.val().fmt(f)
    }
}

impl TryFrom<f32> for F32AbsValAtMost1 {
    type Error = Error;

    /// attempts conversion from `f32`, will succeed if -1 <= value <= 1
    ///
    /// # Errors
    /// If conversion doesn't work.
    ///
    /// # Examples
    /// ```
    /// use nanalogue_core::{Error, F32AbsValAtMost1};
    ///
    /// let val1: F32AbsValAtMost1 = 0.5f32.try_into().unwrap();
    /// let val2: F32AbsValAtMost1 = (-0.7_f32).try_into().unwrap();
    /// let val3: Error = F32AbsValAtMost1::try_from(1.7).unwrap_err();
    /// ```
    fn try_from(value: f32) -> Result<Self, Self::Error> {
        F32AbsValAtMost1::new(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f32_abs_val_at_most1_basic() {
        // Test boundary values
        let _: F32AbsValAtMost1 = F32AbsValAtMost1::new(-1.0).unwrap();
        let _: F32AbsValAtMost1 = F32AbsValAtMost1::new(1.0).unwrap();
        let _: F32AbsValAtMost1 = F32AbsValAtMost1::new(0.0).unwrap();

        // Test near-boundary values
        let _: F32AbsValAtMost1 = F32AbsValAtMost1::new(-0.999_999).unwrap();
        let _: F32AbsValAtMost1 = F32AbsValAtMost1::new(0.999_999).unwrap();

        // Test outside boundaries
        let _: Error = F32AbsValAtMost1::new(-1.000_001).unwrap_err();
        let _: Error = F32AbsValAtMost1::new(1.000_001).unwrap_err();
    }

    #[test]
    fn f32_abs_val_at_most1_from_str() {
        // Valid strings
        let _: F32AbsValAtMost1 = F32AbsValAtMost1::from_str("-1.0").unwrap();
        let _: F32AbsValAtMost1 = F32AbsValAtMost1::from_str("1.0").unwrap();
        let _: F32AbsValAtMost1 = F32AbsValAtMost1::from_str("0.0").unwrap();
        let _: F32AbsValAtMost1 = F32AbsValAtMost1::from_str("-0.5").unwrap();
        let _: F32AbsValAtMost1 = F32AbsValAtMost1::from_str("0.5").unwrap();

        // Invalid strings
        let _: Error = F32AbsValAtMost1::from_str("-1.1").unwrap_err();
        let _: Error = F32AbsValAtMost1::from_str("1.1").unwrap_err();
        let _: Error = F32AbsValAtMost1::from_str("abc").unwrap_err();
        let _: Error = F32AbsValAtMost1::from_str("").unwrap_err();
    }

    #[expect(
        clippy::shadow_unrelated,
        reason = "repetition is fine; each block is clearly separated"
    )]
    #[test]
    fn f32_abs_val_at_most1_display() {
        let val = F32AbsValAtMost1::new(-0.5).expect("should create");
        assert_eq!(format!("{val}"), "-0.5");

        let val = F32AbsValAtMost1::new(0.75).expect("should create");
        assert_eq!(format!("{val}"), "0.75");
    }

    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "exact compare ok as (1) very few significant digits, and (2) no arithmetic"
    )]
    fn f32_abs_val_at_most1_val() {
        for test_val in [0.0, 0.1, -0.7, 1.0, -1.0] {
            let val = F32AbsValAtMost1::new(test_val).expect("should create");
            assert_eq!(val.val(), test_val);
        }
    }

    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "exact comparison ok for boundary values and simple fractions"
    )]
    fn f32_abs_val_at_most1_try_from_f32_valid() {
        // Test boundary values
        let neg_one: F32AbsValAtMost1 = (-1.0f32).try_into().expect("should convert");
        assert_eq!(neg_one.val(), -1.0);

        let zero: F32AbsValAtMost1 = 0.0f32.try_into().expect("should convert");
        assert_eq!(zero.val(), 0.0);

        let one: F32AbsValAtMost1 = 1.0f32.try_into().expect("should convert");
        assert_eq!(one.val(), 1.0);

        // Test intermediate positive values
        let half: F32AbsValAtMost1 = 0.5f32.try_into().expect("should convert");
        assert_eq!(half.val(), 0.5);

        // Test intermediate negative values
        let neg_half: F32AbsValAtMost1 = (-0.5f32).try_into().expect("should convert");
        assert_eq!(neg_half.val(), -0.5);

        let neg_three_quarters: F32AbsValAtMost1 = (-0.75f32).try_into().expect("should convert");
        assert_eq!(neg_three_quarters.val(), -0.75);

        // Test near-boundary values
        let near_neg_one: F32AbsValAtMost1 = (-0.999_999f32).try_into().expect("should convert");
        assert_eq!(near_neg_one.val(), -0.999_999);

        let near_one: F32AbsValAtMost1 = 0.999_999f32.try_into().expect("should convert");
        assert_eq!(near_one.val(), 0.999_999);
    }

    #[test]
    fn f32_abs_val_at_most1_try_from_f32_invalid() {
        // Test below lower boundary
        let below_neg_one: Result<F32AbsValAtMost1, _> = (-1.000_001f32).try_into();
        let _: Error = below_neg_one.unwrap_err();

        let too_negative: Result<F32AbsValAtMost1, _> = (-1.5f32).try_into();
        let _: Error = too_negative.unwrap_err();

        // Test above upper boundary
        let above_one: Result<F32AbsValAtMost1, _> = 1.000_001f32.try_into();
        let _: Error = above_one.unwrap_err();

        let too_large: Result<F32AbsValAtMost1, _> = 1.5f32.try_into();
        let _: Error = too_large.unwrap_err();

        // Test special float values
        let infinity: Result<F32AbsValAtMost1, _> = f32::INFINITY.try_into();
        let _: Error = infinity.unwrap_err();

        let neg_infinity: Result<F32AbsValAtMost1, _> = f32::NEG_INFINITY.try_into();
        let _: Error = neg_infinity.unwrap_err();
    }
}
