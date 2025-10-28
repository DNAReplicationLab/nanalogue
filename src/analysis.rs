//! # Analysis functions for modification data processing.
//!
//! Contains `threshold_and_mean` and `threshold_and_gradient` calculations.
use crate::{Contains, Error, F32AbsValAtMost1, F32Bw0and1, ThresholdState};

/// Threshold and calculate mean modification density per window
///
/// # Errors
///
/// Returns `Error::EmptyWindow` if the input slice is empty,
/// or error associated with `F32Bw0and1` creation.
///
/// # Examples
///
/// ```
/// use nanalogue_core::analysis::threshold_and_mean;
/// use nanalogue_core::F32Bw0and1;
///
/// // Test with all modified bases (values >= 128)
/// let mod_data = [128, 150, 200, 255];
/// let result = threshold_and_mean(&mod_data).unwrap();
/// assert_eq!(result.val(), 1.0);
/// ```
pub fn threshold_and_mean(mod_list: &[u8]) -> Result<F32Bw0and1, Error> {
    let win_size: usize = match mod_list.len() {
        0 => Err(Error::EmptyWindow),
        v => Ok(v),
    }?;
    let count_mod: usize = mod_list
        .iter()
        .filter(|x| ThresholdState::GtEq(128).contains(x))
        .count();
    #[expect(
        clippy::cast_precision_loss,
        reason = "if two nums below have >6 digits, we suffer precision loss. \
may fix in future. this means we are using ~Mbp windows, which is unlikely..."
    )]
    F32Bw0and1::new(count_mod as f32 / win_size as f32)
}

/// Threshold and calculate gradient of mod density per window.
///
/// NOTE: to calculate gradient, we need (x, y) data where x is the
/// coordinate along the read and y is 0 or 1 depending on whether the
/// base is modified or not. But, in the calculation below, we use only
/// the y values i.e. the modified values and assume even spacing
/// along the x direction i.e. along the read coordinate. This assumption
/// will break down if modifications occur in bursts, or if positions
/// are missing along the read, or in other scenarios. As our goal
/// here is to use a simple method to calculate gradients and select
/// interesting reads, we are o.k. with an approximate calculation.
///
/// # Errors
///
/// Returns `Error::EmptyWindow` if the input slice is empty, or
/// `Error::InsufficientDataSize` if the input slice has only one element, or
/// an error associated with `F32AbsValAtMost1` creation.
///
/// # Examples
///
/// ```
/// use nanalogue_core::analysis::threshold_and_gradient;
/// use nanalogue_core::F32AbsValAtMost1;
///
/// // Test with increasing modification pattern
/// // For [0, 0, 128, 200]: x_mean = 2.5, modified positions are 3,4
/// // numerator = (3-2.5) + (4-2.5) = 2.0
/// // denominator = (1-2.5)² + (2-2.5)² + (3-2.5)² + (4-2.5)² = 2.25 + 0.25 + 0.25 + 2.25 = 5.0
/// // gradient = 2.0 / 5.0 = 0.4
/// let mod_data = [0, 0, 128, 200];
/// let result = threshold_and_gradient(&mod_data).unwrap();
/// assert_eq!(result.val(), 0.4);
/// ```
#[expect(
    clippy::cast_precision_loss,
    reason = "if win_size or mod_list length has too many digits (>6), we suffer precision loss. \
may fix in future. this means we are using ~Mbp windows, which is unlikely..."
)]
pub fn threshold_and_gradient(mod_list: &[u8]) -> Result<F32AbsValAtMost1, Error> {
    let win_size = match mod_list.len() {
        0 => Err(Error::EmptyWindow),
        1 => Err(Error::InsufficientDataSize),
        v => Ok(v),
    }?;
    let x_mean = f32::midpoint(win_size as f32, 1.0);
    let numerator: f32 = mod_list
        .iter()
        .enumerate()
        .map(|(i, x)| {
            if ThresholdState::GtEq(128).contains(x) {
                i as f32 + 1.0 - x_mean
            } else {
                0.0
            }
        })
        .sum();
    let denominator: f32 = (1..=win_size)
        .map(|x| (x as f32 - x_mean) * (x as f32 - x_mean))
        .sum();
    F32AbsValAtMost1::new(numerator / denominator)
}

/// Threshold, calculate mean modification density per window, and check against threshold
///
/// This function calculates the mean modification density and returns an error if the
/// mean is below the provided threshold value.
///
/// # Errors
///
/// Returns `Error::EmptyWindow` if the input slice is empty, or
/// `Error::WindowDensBelowThres` if the calculated density is below the threshold, or
/// error associated with `F32Bw0and1` creation.
///
/// # Examples
///
/// ```
/// use nanalogue_core::analysis::threshold_and_mean_and_thres_win;
/// use nanalogue_core::F32Bw0and1;
///
/// // Create a threshold of 0.5 (50% of bases must be modified)
/// let threshold = F32Bw0and1::new(0.5).unwrap();
///
/// // Test with 75% modified bases (values >= 128)
/// let mod_data = [0, 128, 200, 255];
/// let result = threshold_and_mean_and_thres_win(&mod_data, threshold).unwrap();
/// assert_eq!(result.val(), 0.75);
///
/// // Test with 25% modified bases (should return WindowDensBelowThres error)
/// let mod_data = [0, 0, 0, 128];
/// let result = threshold_and_mean_and_thres_win(&mod_data, threshold);
/// assert!(result.is_err());
/// ```
pub fn threshold_and_mean_and_thres_win(
    mod_list: &[u8],
    threshold: F32Bw0and1,
) -> Result<F32Bw0and1, Error> {
    let density = threshold_and_mean(mod_list)?;

    if density.val() < threshold.val() {
        Err(Error::WindowDensBelowThres { density, threshold })
    } else {
        Ok(density)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "divide-by-four unlikely to give floating point errors"
    )]
    fn test_threshold_and_mean_no_modified() {
        // Test with no modified bases (values < 128)
        let mod_data = [0, 50, 100, 127];
        let result = threshold_and_mean(&mod_data).unwrap();
        assert_eq!(result.val(), 0.0);
    }

    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "divide-by-four unlikely to give floating point errors"
    )]
    fn test_threshold_and_mean_half_modified() {
        // Test with half modified bases
        let mod_data = [100, 128, 50, 200];
        let result = threshold_and_mean(&mod_data).unwrap();
        assert_eq!(result.val(), 0.5);
    }

    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "divide-by-five unlikely to give floating point errors"
    )]
    fn test_threshold_and_gradient_decreasing() {
        // Test with decreasing modification pattern
        // For [200, 128, 0, 0]: x_mean = 2.5, modified positions are 1,2
        // numerator = (1-2.5) + (2-2.5) = -1.5 + (-0.5) = -2.0
        // denominator = (1-2.5)^2 + (2-2.5)^2 + (3-2.5)^2 + (4-2.5)^2 = 5.0
        // gradient = -2.0 / 5.0 = -0.4
        let mod_data = [200, 128, 0, 0];
        let result = threshold_and_gradient(&mod_data).unwrap();
        assert_eq!(result.val(), -0.4);
    }

    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "divide-by-five unlikely to give floating point errors"
    )]
    fn test_threshold_and_gradient_alternating() {
        // Test with alternating pattern
        // For [128, 0, 128, 0]: modified positions are 1,3
        // numerator = (1-2.5) + (3-2.5) = -1.5 + 0.5 = -1.0
        // denominator = (1-2.5)^2 + (2-2.5)^2 + (3-2.5)^2 + (4-2.5)^2 = 5.0
        // gradient = -1.0 / 5.0 = -0.2
        let mod_data = [128, 0, 128, 0];
        let result = threshold_and_gradient(&mod_data).unwrap();
        assert_eq!(result.val(), -0.2);
    }

    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "divide-by-four unlikely to give floating point errors"
    )]
    fn test_threshold_and_mean_and_thres_win_above_threshold() {
        // Test with density above threshold
        let mod_data = [0, 128, 200, 255];
        let threshold = F32Bw0and1::new(0.5).unwrap();
        let result = threshold_and_mean_and_thres_win(&mod_data, threshold).unwrap();
        assert_eq!(result.val(), 0.75);
    }

    #[test]
    #[expect(
        clippy::float_cmp,
        reason = "divide-by-four unlikely to give floating point errors"
    )]
    fn test_threshold_and_mean_and_thres_win_equal_threshold() {
        // Test with density equal to threshold
        let mod_data = [0, 0, 128, 128];
        let threshold = F32Bw0and1::new(0.5).unwrap();
        let result = threshold_and_mean_and_thres_win(&mod_data, threshold).unwrap();
        assert_eq!(result.val(), 0.5);
    }

    #[test]
    #[should_panic(expected = "WindowDensBelowThres")]
    fn test_threshold_and_mean_and_thres_win_below_threshold() {
        // Test with density below threshold
        let mod_data = [0, 0, 0, 128];
        let threshold = F32Bw0and1::new(0.5).unwrap();

        // This should panic with a WindowDensBelowThres error
        let _ = threshold_and_mean_and_thres_win(&mod_data, threshold).unwrap();
    }

    #[test]
    #[should_panic(expected = "EmptyWindow")]
    fn test_threshold_and_mean_empty_array() {
        let mod_data: [u8; 0] = [];
        let _ = threshold_and_mean(&mod_data).unwrap();
    }

    #[test]
    #[should_panic(expected = "EmptyWindow")]
    fn test_threshold_and_gradient_empty_array() {
        let mod_data: [u8; 0] = [];
        let _ = threshold_and_gradient(&mod_data).unwrap();
    }

    #[test]
    #[should_panic(expected = "InsufficientDataSize")]
    fn test_threshold_and_gradient_single_element() {
        let mod_data = [128];
        let _ = threshold_and_gradient(&mod_data).unwrap();
    }
}
