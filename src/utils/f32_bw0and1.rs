//! `F32Bw0and1` struct for constrained float between 0 and 1
//! Ensures floating-point values are within valid range at construction

use super::f32_abs_val_at_most1::F32AbsValAtMost1;
use crate::Error;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

/// Datatype holding a float (f32) between 0 and 1 (both inclusive) guaranteed at creation.
#[derive(Debug, Clone, Default, Copy, PartialOrd, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "f32")]
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
    ///
    /// # Errors
    /// Returns an error if the value is not between 0.0 and 1.0 (inclusive).
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
    #[must_use]
    pub fn val(&self) -> f32 {
        self.0
    }
    /// Shortcut for 1.0
    ///
    #[must_use]
    #[expect(clippy::missing_panics_doc, reason = "no error possible here")]
    pub fn one() -> Self {
        F32Bw0and1::new(1.0).expect("no error")
    }
    /// Shortcut for 0.0
    ///
    #[expect(clippy::missing_panics_doc, reason = "no error possible here")]
    #[must_use]
    pub fn zero() -> Self {
        F32Bw0and1::new(0.0).expect("no error")
    }
    /// Converts from `F32AbsValAtMost1` using the absolute value
    ///
    #[expect(
        clippy::missing_panics_doc,
        reason = "absolute value of a number between -1 and 1 is always between 0 and 1"
    )]
    #[must_use]
    pub fn abs_f32_abs_val_at_most_1(val: F32AbsValAtMost1) -> Self {
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
        F32Bw0and1::new(f32::from(value) / f32::from(u8::MAX)).expect("no F32 conversion error")
    }
}

impl From<F32Bw0and1> for u8 {
    /// Convert into a u8 i.e. a number >= 0 and <= 255
    #[expect(
        clippy::cast_possible_truncation,
        reason = "float to non-negative int involves loss, we limit this with round()"
    )]
    #[expect(clippy::cast_sign_loss, reason = "these are positive numbers")]
    fn from(value: F32Bw0and1) -> Self {
        (value.val() * 255.0).round() as u8
    }
}

impl fmt::Display for F32Bw0and1 {
    /// converts to string for display.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.val().fmt(f)
    }
}

impl TryFrom<f32> for F32Bw0and1 {
    type Error = Error;

    /// attempts conversion from `f32`, will succeed if 0 <= value <= 1
    ///
    /// # Errors
    /// If conversion doesn't work.
    ///
    /// # Examples
    /// ```
    /// use nanalogue_core::{Error, F32Bw0and1};
    ///
    /// let val1: F32Bw0and1 = 0.5f32.try_into().unwrap();
    /// let val2: F32Bw0and1 = 0.7_f32.try_into().unwrap();
    /// let val3: Error = F32Bw0and1::try_from(1.7).unwrap_err();
    /// ```
    fn try_from(value: f32) -> Result<Self, Self::Error> {
        F32Bw0and1::new(value)
    }
}

impl From<F32Bw0and1> for F32AbsValAtMost1 {
    /// Convert between the two types of floats
    fn from(value: F32Bw0and1) -> Self {
        F32AbsValAtMost1::new(value.val()).expect("no F32 conversion error")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f32_bw0and1_basic() {
        // Test boundary values
        let _: F32Bw0and1 = F32Bw0and1::new(0.0).unwrap();
        let _: F32Bw0and1 = F32Bw0and1::new(1.0).unwrap();

        // Test near-boundary values
        let _: F32Bw0and1 = F32Bw0and1::new(0.000_001).unwrap();
        let _: F32Bw0and1 = F32Bw0and1::new(0.999_999).unwrap();

        // Test outside boundaries
        let _: Error = F32Bw0and1::new(-0.000_001).unwrap_err();
        let _: Error = F32Bw0and1::new(1.000_001).unwrap_err();
    }

    #[test]
    fn f32_bw0and1_from_str() {
        // Valid strings
        let _: F32Bw0and1 = F32Bw0and1::from_str("0.0").unwrap();
        let _: F32Bw0and1 = F32Bw0and1::from_str("1.0").unwrap();
        let _: F32Bw0and1 = F32Bw0and1::from_str("0.5").unwrap();

        // Invalid strings
        let _: Error = F32Bw0and1::from_str("-0.1").unwrap_err();
        let _: Error = F32Bw0and1::from_str("1.1").unwrap_err();
        let _: Error = F32Bw0and1::from_str("abc").unwrap_err();
        let _: Error = F32Bw0and1::from_str("").unwrap_err();
    }

    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "0.0, 1.0 generated without computation, so can compare"
    )]
    fn f32_bw0and1_shortcuts() {
        let zero = F32Bw0and1::zero();
        assert_eq!(zero.val(), 0.0);

        let one = F32Bw0and1::one();
        assert_eq!(one.val(), 1.0);
    }

    #[test]
    #[expect(clippy::float_cmp, reason = "comparing exactly is ok for 0 and 1")]
    fn f32_bw0and1_from_u8() {
        // Test boundary values
        let zero = F32Bw0and1::from(0u8);
        assert_eq!(zero.val(), 0.0);

        let max = F32Bw0and1::from(255u8);
        assert_eq!(max.val(), 1.0);

        // Test some intermediate values
        let half = F32Bw0and1::from(128u8);
        // 128/255 is approx 0.5020
        assert_eq!(format!("{:.4}", half.val()), "0.5020");

        // Test exact calculation
        let test_val = 100u8;
        let converted = F32Bw0and1::from(test_val);
        // 100/255 is approx 0.3922
        assert_eq!(format!("{:.4}", converted.val()), "0.3922");
    }

    #[test]
    fn f32_bw0and1_into_u8() {
        // Test boundary values
        let zero = F32Bw0and1::new(0.0).expect("should create");
        let zero_u8: u8 = zero.into();
        assert_eq!(zero_u8, 0u8);

        let one = F32Bw0and1::new(1.0).expect("should create");
        let one_u8: u8 = one.into();
        assert_eq!(one_u8, 255u8);

        // Test some intermediate values
        let half = F32Bw0and1::new(0.5).expect("should create");
        let half_u8: u8 = half.into();
        assert!((127..=128).contains(&half_u8)); // Should be ~127.5

        // Test exact calculation
        let test_val = F32Bw0and1::new(0.392_156_87).expect("should create");
        let converted_u8: u8 = test_val.into();
        assert_eq!(converted_u8, 100u8);

        // Test near-boundary values
        let near_zero = F32Bw0and1::new(0.001).expect("should create");
        let near_zero_u8: u8 = near_zero.into();
        assert_eq!(near_zero_u8, 0u8);

        let near_one = F32Bw0and1::new(0.999).expect("should create");
        let near_one_u8: u8 = near_one.into();
        assert_eq!(near_one_u8, 255u8);
    }

    #[test]
    fn f32_bw0and1_u8_roundtrip() {
        // Test that converting u8 -> F32Bw0and1 -> u8 gives reasonable results
        for val in [0u8, 1, 50, 100, 128, 200, 254, 255] {
            let f32_val = F32Bw0and1::from(val);
            let converted_back: u8 = f32_val.into();
            assert_eq!(converted_back, val);
        }
    }

    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "conversion to abs values shouldn't result in floating point problems"
    )]
    fn f32_types_integration() {
        // Test conversion from F32Bw0and1 to F32AbsValAtMost1
        let pos_values = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        for val in pos_values {
            let bw_val = F32Bw0and1::new(val).expect("should create");
            let abs_val: F32AbsValAtMost1 = bw_val.into();
            assert_eq!(abs_val.val(), val);
        }

        // Test conversion via absolute value function
        let neg_val = F32AbsValAtMost1::new(-0.5).expect("should create");
        let abs_converted_neg = F32Bw0and1::abs_f32_abs_val_at_most_1(neg_val);
        assert_eq!(abs_converted_neg.val(), 0.5);

        let pos_val = F32AbsValAtMost1::new(0.7).expect("should create");
        let abs_converted_pos = F32Bw0and1::abs_f32_abs_val_at_most_1(pos_val);
        assert_eq!(abs_converted_pos.val(), 0.7);
    }

    #[expect(
        clippy::shadow_unrelated,
        reason = "repetition is fine; each block is clearly separated"
    )]
    #[test]
    fn f32_bw0and1_display() {
        let val = F32Bw0and1::new(0.5).expect("should create");
        assert_eq!(format!("{val}"), "0.5");

        let val = F32Bw0and1::new(0.75).expect("should create");
        assert_eq!(format!("{val}"), "0.75");
    }

    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "exact comparison ok for boundary values and simple fractions"
    )]
    fn f32_bw0and1_try_from_f32_valid() {
        // Test boundary values
        let zero: F32Bw0and1 = 0.0f32.try_into().expect("should convert");
        assert_eq!(zero.val(), 0.0);

        let one: F32Bw0and1 = 1.0f32.try_into().expect("should convert");
        assert_eq!(one.val(), 1.0);

        // Test intermediate values
        let half: F32Bw0and1 = 0.5f32.try_into().expect("should convert");
        assert_eq!(half.val(), 0.5);

        let three_quarters: F32Bw0and1 = 0.75f32.try_into().expect("should convert");
        assert_eq!(three_quarters.val(), 0.75);

        // Test near-boundary values
        let near_zero: F32Bw0and1 = 0.000_001f32.try_into().expect("should convert");
        assert_eq!(near_zero.val(), 0.000_001);

        let near_one: F32Bw0and1 = 0.999_999f32.try_into().expect("should convert");
        assert_eq!(near_one.val(), 0.999_999);
    }

    #[test]
    fn f32_bw0and1_try_from_f32_invalid() {
        // Test below lower boundary
        let below_zero: Result<F32Bw0and1, _> = (-0.000_001f32).try_into();
        let _: Error = below_zero.unwrap_err();

        let negative: Result<F32Bw0and1, _> = (-0.5f32).try_into();
        let _: Error = negative.unwrap_err();

        // Test above upper boundary
        let above_one: Result<F32Bw0and1, _> = 1.000_001f32.try_into();
        let _: Error = above_one.unwrap_err();

        let too_large: Result<F32Bw0and1, _> = 1.5f32.try_into();
        let _: Error = too_large.unwrap_err();

        // Test special float values
        let infinity: Result<F32Bw0and1, _> = f32::INFINITY.try_into();
        let _: Error = infinity.unwrap_err();

        let neg_infinity: Result<F32Bw0and1, _> = f32::NEG_INFINITY.try_into();
        let _: Error = neg_infinity.unwrap_err();
    }
}
