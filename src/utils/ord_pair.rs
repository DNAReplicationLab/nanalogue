//! `OrdPair` struct for ordered pairs with low <= high guarantee
//! Provides ordered pair datatype with validation and interval operations

use super::contains::Contains;
use crate::Error;
use serde::{Deserialize, Serialize};
use std::convert::From;
use std::fmt;
use std::fmt::Debug;
use std::ops::RangeInclusive;
use std::str::FromStr;

/// Datatype holding two values low, high such that low <= high is guaranteed at creation.
#[derive(Debug, Clone, Copy, PartialOrd, PartialEq, Serialize, Deserialize)]
pub struct OrdPair<T: Clone + Copy + Debug> {
    low: T,
    high: T,
}

impl<T: Clone + Copy + Debug + Default + Ord> Default for OrdPair<T> {
    fn default() -> Self {
        OrdPair {
            low: T::default(),
            high: T::default(),
        }
    }
}

impl<T: Clone + Copy + Debug + PartialEq + PartialOrd> OrdPair<T> {
    /// Constructor with two values, will fail if ordering in input is not respected.
    ///
    /// ```should_panic
    /// use nanalogue_core::OrdPair;
    /// let x = OrdPair::<f32>::new(1.0,0.0).unwrap();
    /// ```
    /// ```
    /// # use nanalogue_core::OrdPair;
    /// let x = OrdPair::<f32>::new(0.0, 1.0)?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    pub fn new(low: T, high: T) -> Result<Self, Error> {
        if low <= high {
            Ok(OrdPair { low, high })
        } else {
            Err(Error::WrongOrder)
        }
    }
    /// Gets the low value
    ///
    /// ```
    /// use nanalogue_core::OrdPair;
    /// let x = OrdPair::<u8>::new(10,11).expect("no failure");
    /// assert_eq!(x.get_low(),10);
    /// ```
    pub fn get_low(&self) -> T {
        self.low
    }
    /// Gets the high value
    ///
    /// ```
    /// use nanalogue_core::OrdPair;
    /// let x = OrdPair::<u8>::new(10,11).expect("no failure");
    /// assert_eq!(x.get_high(),11);
    /// ```
    pub fn get_high(&self) -> T {
        self.high
    }
}

impl OrdPair<u64> {
    /// Parse an interval string specifically for genomic regions.
    /// Supports formats like "1000-2000" and "1000-" (where end defaults to `u64::MAX`).
    /// Enforces strict inequality (start < end).
    ///
    /// ```
    /// use nanalogue_core::OrdPair;
    ///
    /// // Standard interval
    /// let interval = OrdPair::<u64>::from_interval("1000-2000")?;
    /// assert_eq!(interval.get_low(), 1000);
    /// assert_eq!(interval.get_high(), 2000);
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```
    /// # use nanalogue_core::OrdPair;
    /// #
    /// // Open-ended interval (end defaults to u64::MAX)
    /// let interval = OrdPair::<u64>::from_interval("1000-")?;
    /// assert_eq!(interval.get_low(), 1000);
    /// assert_eq!(interval.get_high(), u64::MAX);
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```should_panic
    /// # use nanalogue_core::OrdPair;
    /// #
    /// // Equal start and end should fail
    /// let interval = OrdPair::<u64>::from_interval("1000-1000")?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    pub fn from_interval(interval_str: &str) -> Result<Self, Error> {
        let parts: Vec<&str> = interval_str.split('-').collect();

        match parts.len() {
            2 => {
                let start = parts[0].trim().parse::<u64>().map_err(|_| {
                    Error::OrdPairConversionError(
                        "Invalid start coordinate in interval!".to_string(),
                    )
                })?;

                let end = if parts[1].trim().is_empty() {
                    // Open-ended interval: "1000-"
                    u64::MAX
                } else {
                    // Closed interval: "1000-2000"
                    parts[1].trim().parse::<u64>().map_err(|_| {
                        Error::OrdPairConversionError(
                            "Invalid end coordinate in interval!".to_string(),
                        )
                    })?
                };

                // Enforce strict inequality (start < end)
                if start < end {
                    Ok(OrdPair {
                        low: start,
                        high: end,
                    })
                } else {
                    Err(Error::OrdPairConversionError(
                        "Genomic intervals require start < end (strict inequality)".to_string(),
                    ))
                }
            }
            _ => Err(Error::OrdPairConversionError(
                "Invalid interval format! Expected 'start-end' or 'start-'".to_string(),
            )),
        }
    }
}

impl<T: Clone + Copy + Debug + PartialEq + PartialOrd + FromStr> FromStr for OrdPair<T> {
    type Err = Error;

    /// Parse a string to obtain an Ordered Pair, return Error if cannot be done.
    fn from_str(val_str: &str) -> Result<Self, Self::Err> {
        macro_rules! parse_error {
            () => {
                Err(Error::OrdPairConversionError(
                    "Bad ordered pair inputs!".to_string(),
                ))
            };
        }
        let v: Vec<&str> = val_str.split(',').map(str::trim).collect();
        match v.len() {
            2 => {
                let Ok(low) = T::from_str(v[0]) else {
                    parse_error!()?
                };
                let Ok(high) = T::from_str(v[1]) else {
                    parse_error!()?
                };
                OrdPair::<T>::new(low, high)
            }
            _ => parse_error!(),
        }
    }
}

impl<T: Clone + Copy + Debug + PartialEq + PartialOrd> From<OrdPair<T>> for RangeInclusive<T> {
    /// Convert the `OrdPair` into a `RangeInclusive` i.e. (start..=end)
    fn from(value: OrdPair<T>) -> Self {
        RangeInclusive::<T>::new(value.get_low(), value.get_high())
    }
}

impl<T: Clone + Copy + Debug + PartialEq + PartialOrd> Contains<T> for OrdPair<T> {
    /// Check if the provided value is within the Range of the `OrdPair`
    fn contains(&self, val: &T) -> bool {
        RangeInclusive::<T>::from(*self).contains(val)
    }
}

impl<T: Clone + Copy + Debug + fmt::Display + PartialEq + PartialOrd> fmt::Display for OrdPair<T> {
    /// converts to string for display i.e. "low, high"
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}, {}", self.get_low(), self.get_high())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tests if our Ordered Pair struct can be correctly obtained from strings
    #[test]
    fn test_ord_pair_from_str() {
        assert_eq!(
            OrdPair::<f32>::from_str("1.0,2.0")
                .expect("no failure")
                .get_low(),
            1.0
        );
        assert_eq!(
            OrdPair::<f32>::from_str("1.0,2.0")
                .expect("no failure")
                .get_high(),
            2.0
        );
        assert_eq!(
            OrdPair::<u8>::from_str("1, 2")
                .expect("no failure")
                .get_low(),
            1
        );
        assert_eq!(
            OrdPair::<u8>::from_str("1, 2")
                .expect("no failure")
                .get_high(),
            2
        );
    }

    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn test_ord_pair_from_str_empty_first_value_panics() {
        let _ = OrdPair::<u8>::from_str(",2").unwrap();
    }

    #[test]
    #[should_panic(expected = "WrongOrder")]
    fn test_ord_pair_from_str_wrong_order_panics() {
        let _ = OrdPair::<u8>::from_str("2,1").unwrap();
    }

    /// Tests if `OrdPair` can be converted into a range
    #[test]
    fn test_ord_pair_to_range() {
        assert_eq!(
            (3..=5),
            RangeInclusive::from(OrdPair::new(3, 5).expect("no failure"))
        );
    }

    /// Tests `OrdPair::from_interval` method for genomic intervals
    #[test]
    fn test_ord_pair_from_interval() {
        // Standard interval
        let interval = OrdPair::<u64>::from_interval("1000-2000").expect("should parse");
        assert_eq!(interval.get_low(), 1000);
        assert_eq!(interval.get_high(), 2000);

        // Open-ended interval
        let interval = OrdPair::<u64>::from_interval("1000-").expect("should parse");
        assert_eq!(interval.get_low(), 1000);
        assert_eq!(interval.get_high(), u64::MAX);
    }

    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn test_ord_pair_from_interval_equal_start_end_panics() {
        let _ = OrdPair::<u64>::from_interval("1000-1000").unwrap();
    }

    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn test_ord_pair_from_interval_start_greater_than_end_panics() {
        let _ = OrdPair::<u64>::from_interval("2000-1000").unwrap();
    }

    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn test_ord_pair_from_interval_no_dash_panics() {
        let _ = OrdPair::<u64>::from_interval("1000").unwrap();
    }

    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn test_ord_pair_from_interval_too_many_dashes_panics() {
        let _ = OrdPair::<u64>::from_interval("1000-2000-3000").unwrap();
    }

    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn test_ord_pair_from_interval_invalid_start_panics() {
        let _ = OrdPair::<u64>::from_interval("abc-2000").unwrap();
    }

    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn test_ord_pair_from_interval_invalid_end_panics() {
        let _ = OrdPair::<u64>::from_interval("1000-xyz").unwrap();
    }

    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn test_ord_pair_from_interval_empty_string_panics() {
        let _ = OrdPair::<u64>::from_interval("").unwrap();
    }

    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn test_ord_pair_from_interval_just_dash_panics() {
        let _ = OrdPair::<u64>::from_interval("-").unwrap();
    }

    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn test_ord_pair_from_interval_negative_numbers_panics() {
        let _ = OrdPair::<u64>::from_interval("-100-200").unwrap();
    }

    #[test]
    fn test_ord_pair_contains() {
        let pair = OrdPair::new(10, 20).expect("should create");
        assert!(pair.contains(&15));
        assert!(pair.contains(&10)); // inclusive lower bound
        assert!(pair.contains(&20)); // inclusive upper bound
        assert!(!pair.contains(&5));
        assert!(!pair.contains(&25));
    }

    #[test]
    fn test_ord_pair_display() {
        let pair = OrdPair::new(10, 20).expect("should create");
        assert_eq!(format!("{pair}"), "10, 20");

        let float_pair = OrdPair::new(1.5, 2.5).expect("should create");
        assert_eq!(format!("{float_pair}"), "1.5, 2.5");
    }

    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn test_ord_pair_from_str_empty_string_panics() {
        let _ = OrdPair::<i32>::from_str("").unwrap();
    }

    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn test_ord_pair_from_str_single_value_panics() {
        let _ = OrdPair::<i32>::from_str("1").unwrap();
    }

    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn test_ord_pair_from_str_too_many_values_panics() {
        let _ = OrdPair::<i32>::from_str("1,2,3").unwrap();
    }

    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn test_ord_pair_from_str_non_numeric_panics() {
        let _ = OrdPair::<i32>::from_str("abc,def").unwrap();
    }
}
