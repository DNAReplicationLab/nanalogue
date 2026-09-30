#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Benchmark dense, heavily overlapping windows through the public window-reads path.
//!
//! # Workload
//!
//! The fixture is one in-memory, unmapped 20,000-base `ModBAM` record with a modification
//! candidate at every base. It uses a 10,000-candidate window and step size 1, producing 10,001
//! windows whose neighbours overlap by 9,999 candidates. This deliberately stresses work that
//! can be shared while a window slides forward. The benchmark verifies that all synthetic
//! candidates were parsed before timing starts.
//!
//! # Timed work
//!
//! Each of the nine samples calls the public [`window_reads::run`] path once. A sample includes
//! record validation, MM/ML parsing, window construction, reference-bound and base-quality
//! calculations, and formatting every TSV row into a reusable byte buffer. Reusing the buffer
//! excludes output-buffer growth after warm-up, while black-boxing its bytes ensures formatting
//! is performed. Filesystem I/O, BAM decompression, record construction, and CLI startup are not
//! measured. The reported duration is the observed median sample. Allocation statistics come
//! from a separate pass so their atomic counters do not affect the timing samples.
//!
//! # Interpretation
//!
//! The window callback intentionally reads only its first modification value. This isolates the
//! windowing infrastructure and exposes repeated scans of overlapping reference and base-quality
//! spans. It is not an end-to-end benchmark of analysis functions such as `threshold_and_mean`,
//! which scan every modification value in every window and can become the dominant cost. Treat
//! this benchmark as the attainable overhead of the surrounding window machinery, not as expected
//! throughput for every real analysis.

use nanalogue_core::{
    F32AbsValAtMost1, InputMods, InputWindowing, OptionalTag, nanalogue_mm_ml_parser, window_reads,
};
use rust_htslib::bam::{Record, record::Aux};
use std::alloc::{GlobalAlloc, Layout, System};
use std::hint::black_box;
use std::rc::Rc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::time::{Duration, Instant};

/// Number of modification candidates in the synthetic read.
const CANDIDATES: usize = 20_000;
/// Number of candidates in each window.
const WINDOW_SIZE: usize = 10_000;
/// Number of independently timed samples.
const SAMPLES: usize = 9;

const _: () = {
    assert!(
        CANDIDATES > WINDOW_SIZE,
        "read must contain multiple windows"
    );
    assert!(WINDOW_SIZE > 0, "window must be nonzero");
    assert!(SAMPLES >= 3, "fewer than three samples found");
    assert!(SAMPLES & 1 == 1, "sample count must be odd");
};

/// Whether allocator activity should be counted.
static COUNT_ALLOCATIONS: AtomicBool = AtomicBool::new(false);
/// Number of allocations observed during measurement.
static ALLOCATIONS: AtomicUsize = AtomicUsize::new(0);
/// Number of reallocations observed during measurement.
static REALLOCATIONS: AtomicUsize = AtomicUsize::new(0);
/// Total bytes requested by allocations and reallocations during measurement.
static ALLOCATED_BYTES: AtomicUsize = AtomicUsize::new(0);

/// System allocator wrapper that records allocation activity.
struct CountingAllocator;

// SAFETY: Every operation delegates directly to the system allocator with the
// original pointer and layout. The counters do not affect allocation behavior.
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        if COUNT_ALLOCATIONS.load(Ordering::Relaxed) {
            let _previous_allocations = ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
            let _previous_bytes = ALLOCATED_BYTES.fetch_add(layout.size(), Ordering::Relaxed);
        }
        // SAFETY: The caller supplies the layout required by `GlobalAlloc`.
        unsafe { System.alloc(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        // SAFETY: The caller supplies the pointer and layout from allocation.
        unsafe { System.dealloc(ptr, layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        if COUNT_ALLOCATIONS.load(Ordering::Relaxed) {
            let _previous_reallocations = REALLOCATIONS.fetch_add(1, Ordering::Relaxed);
            let _previous_bytes = ALLOCATED_BYTES.fetch_add(new_size, Ordering::Relaxed);
        }
        // SAFETY: The caller supplies the original allocation and new size.
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
/// Allocator used to count benchmark allocation activity.
static GLOBAL: CountingAllocator = CountingAllocator;

/// Builds one unmapped read with an annotated candidate at every base.
#[expect(
    clippy::integer_division_remainder_used,
    reason = "remainders generate deterministic benchmark qualities"
)]
fn dense_record() -> Record {
    let sequence = vec![b'A'; CANDIDATES];
    let qualities: Vec<u8> = (0..CANDIDATES)
        .map(|position| u8::try_from(position % 94).expect("quality fits u8"))
        .collect();
    let modification_qualities: Vec<u8> = (0..CANDIDATES)
        .map(|position| u8::try_from(position % 256).expect("quality fits u8"))
        .collect();
    let mut mm_tag = String::with_capacity(CANDIDATES.saturating_mul(2).saturating_add(5));
    mm_tag.push_str("A+a");
    for _ in 0..CANDIDATES {
        mm_tag.push_str(",0");
    }
    mm_tag.push(';');

    let mut record = Record::new();
    record.set(b"dense-overlap", None, &sequence, &qualities);
    record.set_flags(4);
    record
        .push_aux(b"MM", Aux::String(&mm_tag))
        .expect("MM tag should fit the record");
    record
        .push_aux(b"ML", Aux::ArrayU8((&modification_qualities).into()))
        .expect("ML tag should fit the record");
    record
}

/// Runs one complete dense-window pass and retains formatted output for observation.
fn run_once(
    output: &mut Vec<u8>,
    record: &Rc<Record>,
    options: InputWindowing,
    mods: &InputMods<OptionalTag>,
) {
    output.clear();
    window_reads::run(
        output,
        [Ok::<Rc<Record>, rust_htslib::errors::Error>(Rc::clone(
            record,
        ))],
        options,
        mods,
        |values| {
            F32AbsValAtMost1::new(f32::from(*values.first().expect("window is non-empty")) / 255.0)
        },
    )
    .expect("dense fixture should produce windows");
    let _observed_output: &[u8] = black_box(output.as_slice());
}

#[expect(
    clippy::cast_precision_loss,
    reason = "benchmark reporting converts bounded work counts to f64"
)]
#[expect(
    clippy::integer_division,
    clippy::integer_division_remainder_used,
    reason = "integer division selects the observed median sample"
)]
#[expect(
    clippy::print_stdout,
    reason = "the benchmark reports its measurements to stdout"
)]
fn main() {
    let record = Rc::new(dense_record());
    let parsed = nanalogue_mm_ml_parser(&record, |&_| true, |&_| true, |&_, &_, &_| true, 0)
        .expect("dense fixture should parse");
    let parsed_candidates = parsed
        .base_mods
        .first()
        .expect("fixture should contain one modification type")
        .ranges
        .annotations
        .len();
    assert_eq!(
        parsed_candidates, CANDIDATES,
        "benchmark must parse every synthetic candidate"
    );
    drop(parsed);
    let options: InputWindowing =
        serde_json::from_str(&format!("{{\"win\": {WINDOW_SIZE}, \"step\": 1}}"))
            .expect("window options should parse");
    let mods = InputMods::default();
    let windows_per_sample = CANDIDATES - WINDOW_SIZE + 1;
    let mut output = Vec::new();

    run_once(&mut output, &record, options, &mods);

    let mut timings = Vec::with_capacity(SAMPLES);
    for _ in 0..SAMPLES {
        let start = Instant::now();
        run_once(&mut output, black_box(&record), options, black_box(&mods));
        timings.push(start.elapsed());
    }
    timings.sort_unstable();
    let median = timings.get(SAMPLES / 2).copied().unwrap_or(Duration::ZERO);

    ALLOCATIONS.store(0, Ordering::Relaxed);
    REALLOCATIONS.store(0, Ordering::Relaxed);
    ALLOCATED_BYTES.store(0, Ordering::Relaxed);
    COUNT_ALLOCATIONS.store(true, Ordering::Relaxed);
    run_once(&mut output, &record, options, &mods);
    COUNT_ALLOCATIONS.store(false, Ordering::Relaxed);

    println!("window_reads dense-overlap benchmark");
    println!(
        "candidates={CANDIDATES} width={WINDOW_SIZE} step=1 windows_per_sample={windows_per_sample} samples={SAMPLES}"
    );
    println!(
        "median_ms={:.3} throughput={:.1} windows/s",
        median.as_secs_f64() * 1_000.0,
        windows_per_sample as f64 / median.as_secs_f64()
    );
    println!(
        "allocations={} reallocations={} requested_bytes={}",
        ALLOCATIONS.load(Ordering::Relaxed),
        REALLOCATIONS.load(Ordering::Relaxed),
        ALLOCATED_BYTES.load(Ordering::Relaxed)
    );
    let minimum = timings.first().copied().unwrap_or(Duration::ZERO);
    let maximum = timings.last().copied().unwrap_or(Duration::ZERO);
    println!(
        "range_ms={:.3}..{:.3}",
        minimum.as_secs_f64() * 1_000.0,
        maximum.as_secs_f64() * 1_000.0
    );
}
