//! ThresholdState enum for modification probability thresholds
//! Handles different threshold types for modification data filtering

use super::contains::Contains;
use super::f32_bw0and1::F32Bw0and1;
use super::ord_pair::OrdPair;
use serde::{Deserialize, Serialize};
use std::fmt;

/// Types of thresholds on modification level that can be applied to modification data.
/// Two possible use cases: (1) to specify that reading mod data should be restricted
/// to bases at least this level of modified, or (2) to specify that only bases
/// in this range should be regarded as modified.
/// Values are 0 to 255 below as that's how they are stored in a modBAM file and
/// this struct is expected to be used in contexts dealing directly with this data.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum ThresholdState {
    /// modification probability >= this value, values are 0 to 255
    GtEq(u8),
    /// modification probability not within this range.
    /// We expect this to be used to filter out modification calls
    /// around 0.5 i.e. ones with the most uncertainty, although
    /// users of this crate are free to set this to an interval
    /// not including 0.5
    InvertGtEqLtEq(OrdPair<u8>),
    /// modification probability >= first value, and mod prob
    /// not within the second range i.e. the 'and' combination
    /// of the two possibilities above
    Both((u8, OrdPair<u8>)),
}

/// default threshold is >= 0 i.e. all mods are allowed
impl Default for ThresholdState {
    fn default() -> Self {
        ThresholdState::GtEq(0)
    }
}

/// Displays thresholds but using floating point numbers between 0 and 1
///
/// Example 1:
/// ```
/// use nanalogue_core::{ThresholdState, OrdPair};
/// let b = ThresholdState::GtEq(100);
/// assert_eq!("probabilities >= 0.390625", format!("{}", b));
/// ```
/// Example 2:
/// ```
/// # use nanalogue_core::{ThresholdState, OrdPair};
/// let b = ThresholdState::InvertGtEqLtEq(OrdPair::new(200, 220).expect("no error"));
/// assert_eq!("probabilities < 0.78125 or > 0.859375", format!("{}", b));
/// ```
///
/// Example 3:
/// ```
/// # use nanalogue_core::{ThresholdState, OrdPair};
/// let b = ThresholdState::Both((100, OrdPair::new(200, 220).expect("no error")));
/// assert_eq!("probabilities >= 0.390625 and (probabilities < 0.78125 or > 0.859375)", format!("{}", b));
/// ```
impl fmt::Display for ThresholdState {
    /// display the u8 thresholds as a floating point number between 0 and 1
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let printable = match &self {
            ThresholdState::GtEq(v) => format!("probabilities >= {}", F32Bw0and1::from(*v)),
            ThresholdState::InvertGtEqLtEq(v) => {
                format!(
                    "probabilities < {} or > {}",
                    F32Bw0and1::from(v.get_low()),
                    F32Bw0and1::from(v.get_high())
                )
            }
            ThresholdState::Both((a, b)) => {
                format!(
                    "{} and ({})",
                    ThresholdState::GtEq(*a),
                    ThresholdState::InvertGtEqLtEq(*b)
                )
            }
        };
        write!(f, "{printable}")
    }
}

/// Check if a given u8 is within the interval covered
///
/// Example 1:
/// ```
/// use nanalogue_core::{Error, OrdPair, ThresholdState, Contains};
/// let b = ThresholdState::GtEq(100);
/// assert!(b.contains(&101));
/// assert!(b.contains(&100));
/// assert!(!b.contains(&99));
/// ```
/// Example 2:
/// ```
/// # use nanalogue_core::{Error, OrdPair, ThresholdState, Contains};
/// let b = ThresholdState::InvertGtEqLtEq(OrdPair::new(200, 220).expect("no error"));
/// assert!(b.contains(&100));
/// assert!(!b.contains(&200));
/// assert!(!b.contains(&210));
/// assert!(!b.contains(&220));
/// assert!(b.contains(&250));
/// ```
/// Example 3:
/// ```
/// # use nanalogue_core::{Error, OrdPair, ThresholdState, Contains};
/// let b = ThresholdState::Both((100, OrdPair::new(200, 220).expect("no error")));
/// assert!(!b.contains(&99));
/// assert!(b.contains(&100));
/// assert!(b.contains(&101));
/// assert!(!b.contains(&200));
/// assert!(!b.contains(&210));
/// assert!(!b.contains(&220));
/// assert!(b.contains(&250));
/// ```
impl Contains<u8> for ThresholdState {
    /// see if value is contained within the interval
    /// specified by the threshold state
    fn contains(&self, val: &u8) -> bool {
        match &self {
            ThresholdState::GtEq(v) => val >= v,
            ThresholdState::InvertGtEqLtEq(w) => !w.contains(val),
            ThresholdState::Both((a, b)) => {
                ThresholdState::GtEq(*a).contains(val)
                    && ThresholdState::InvertGtEqLtEq(*b).contains(val)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_threshold_state_gt_eq() {
        let threshold = ThresholdState::GtEq(100);

        // Test contains functionality
        assert!(threshold.contains(&101));
        assert!(threshold.contains(&100));
        assert!(!threshold.contains(&99));

        // Test display
        let display_str = format!("{}", threshold);
        assert!(display_str.contains("probabilities >="));
        assert!(display_str.contains("0.390625")); // 100/256 = 0.390625
    }

    #[test]
    fn test_threshold_state_invert_gt_eq_lt_eq() {
        let pair = OrdPair::new(200, 220).expect("should create");
        let threshold = ThresholdState::InvertGtEqLtEq(pair);

        // Test contains functionality
        assert!(threshold.contains(&100)); // outside range (below)
        assert!(!threshold.contains(&200)); // within range (boundary)
        assert!(!threshold.contains(&210)); // within range (middle)
        assert!(!threshold.contains(&220)); // within range (boundary)
        assert!(threshold.contains(&250)); // outside range (above)

        // Test display
        let display_str = format!("{}", threshold);
        assert!(display_str.contains("probabilities <"));
        assert!(display_str.contains("or >"));
    }

    #[test]
    fn test_threshold_state_both() {
        let pair = OrdPair::new(200, 220).expect("should create");
        let threshold = ThresholdState::Both((100, pair));

        // Test contains functionality
        assert!(!threshold.contains(&99)); // fails first condition
        assert!(threshold.contains(&100)); // meets both conditions
        assert!(threshold.contains(&101)); // meets both conditions
        assert!(!threshold.contains(&200)); // fails second condition
        assert!(!threshold.contains(&210)); // fails second condition
        assert!(!threshold.contains(&220)); // fails second condition
        assert!(threshold.contains(&250)); // meets both conditions

        // Test display
        let display_str = format!("{}", threshold);
        assert!(display_str.contains("and"));
        assert!(display_str.contains("probabilities >="));
    }

    #[test]
    fn test_threshold_state_default() {
        let default_threshold = ThresholdState::default();
        assert!(matches!(default_threshold, ThresholdState::GtEq(0)));

        // Default should accept all values
        for val in 0..=255u8 {
            assert!(default_threshold.contains(&val));
        }
    }

    #[test]
    fn test_threshold_state_display_consistency() {
        // Test that display format is consistent and meaningful
        let thresholds = vec![
            ThresholdState::GtEq(128),
            ThresholdState::InvertGtEqLtEq(OrdPair::new(100, 150).expect("should create")),
            ThresholdState::Both((50, OrdPair::new(120, 140).expect("should create"))),
        ];

        for threshold in thresholds {
            let display_str = format!("{}", threshold);
            assert!(display_str.contains("probabilities"));
            assert!(!display_str.is_empty());
        }
    }

    #[test]
    fn test_threshold_state_edge_cases() {
        // Test boundary conditions
        let threshold_255 = ThresholdState::GtEq(255);
        assert!(threshold_255.contains(&255));
        assert!(!threshold_255.contains(&254));

        let threshold_0 = ThresholdState::GtEq(0);
        assert!(threshold_0.contains(&0));
        assert!(threshold_0.contains(&255));

        // Test single-value range
        let single_val_pair = OrdPair::new(128, 129).expect("should create");
        let threshold_single = ThresholdState::InvertGtEqLtEq(single_val_pair);
        assert!(threshold_single.contains(&127));
        assert!(!threshold_single.contains(&128));
        assert!(!threshold_single.contains(&129));
        assert!(threshold_single.contains(&130));
    }
}
