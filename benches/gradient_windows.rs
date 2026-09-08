//! Benchmark repeated gradient calculation over overlapping windows.

use std::hint::black_box;
use std::time::{Duration, Instant};

use nanalogue_core::analysis::threshold_and_gradient;

/// Number of deterministic candidate values in the input.
const DATA_LEN: usize = 16_384;
/// Number of candidates in each overlapping window.
const WINDOW_SIZE: usize = 1_024;
/// Number of complete window scans per timed sample.
const REPETITIONS: usize = 4;
/// Odd sample count so the reported median is an observed duration.
const SAMPLES: usize = 11;

/// Calculates every overlapping window repeatedly and returns a checksum.
fn run_windows(data: &[u8]) -> f32 {
    let mut checksum = 0.0;
    for _ in 0..REPETITIONS {
        for window in data.windows(WINDOW_SIZE) {
            checksum += black_box(threshold_and_gradient(black_box(window)).unwrap().val());
        }
    }
    checksum
}

#[expect(
    clippy::cast_precision_loss,
    reason = "benchmark reporting converts bounded durations and work counts to f64"
)]
fn main() {
    let data: Vec<u8> = (0..DATA_LEN)
        .map(|index| {
            if index.wrapping_mul(17).wrapping_add(index / 23) % 11 < 5 {
                200
            } else {
                20
            }
        })
        .collect();

    let _: f32 = black_box(run_windows(&data));

    let mut timings = Vec::with_capacity(SAMPLES);
    let mut checksum = 0.0;
    for _ in 0..SAMPLES {
        let start = Instant::now();
        checksum = black_box(run_windows(&data));
        timings.push(start.elapsed());
    }
    timings.sort_unstable();

    let windows_per_sample = (DATA_LEN - WINDOW_SIZE + 1) * REPETITIONS;
    let median = timings.get(SAMPLES / 2).copied().unwrap_or(Duration::ZERO);
    let minimum = timings.first().copied().unwrap_or(Duration::ZERO);
    let maximum = timings.last().copied().unwrap_or(Duration::ZERO);
    println!(
        "gradient_windows: median={:.3} ms, {:.1} ns/window ({} overlapping windows/sample, {SAMPLES} samples, checksum={checksum})",
        median.as_secs_f64() * 1_000.0,
        median.as_secs_f64() * 1_000_000_000.0 / windows_per_sample as f64,
        windows_per_sample,
    );
    println!(
        "range: {:.3}..{:.3} ms",
        minimum.as_secs_f64() * 1_000.0,
        maximum.as_secs_f64() * 1_000.0,
    );
}
