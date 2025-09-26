//! # Analysis functions for modification data processing.
//!
//! Contains threshold_and_mean and threshold_and_gradient calculations.
use crate::{Contains, Error, F32AbsValBelow1, F32Bw0and1, ThresholdState};

/// Threshold and calculate mean modification density per window
pub fn threshold_and_mean(mod_list: &[u8]) -> Result<F32Bw0and1, Error> {
    let win_size: usize = mod_list.len();
    let count_mod: usize = mod_list
        .iter()
        .filter(|x| ThresholdState::GtEq(128).contains(x))
        .count();
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
pub fn threshold_and_gradient(mod_list: &[u8]) -> Result<F32AbsValBelow1, Error> {
    let win_size = mod_list.len();
    let x_mean: f32 = (win_size as f32 + 1.0) / 2.0;
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
    F32AbsValBelow1::new(numerator / denominator)
}
