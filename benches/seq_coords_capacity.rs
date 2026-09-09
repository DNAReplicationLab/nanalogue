//! Benchmarks successful coordinate extraction when the output retains spare capacity.

use nanalogue_core::bedrs::Bed3;
use nanalogue_core::read_utils::OnlyAlignDataComplete;
use nanalogue_core::{CurrRead, GenomicBed3, nanalogue_bam_reader};
use rust_htslib::bam::{Read as _, Record};
use std::alloc::{GlobalAlloc, Layout, System};
use std::hint::black_box;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::time::Instant;

/// Odd number of timing samples collected.
const SAMPLES: usize = 21;
/// Calls measured in each timing sample.
const CALLS_PER_SAMPLE: usize = 250_000;
/// Calls made before allocator counters and timing samples start.
const WARMUP_CALLS: usize = 100_000;
/// Calls used for the separate allocation measurement.
const ALLOCATION_CALLS: usize = 100_000;

const _: () = {
    assert!(SAMPLES >= 3, "fewer than three samples found");
    assert!(SAMPLES & 1 == 1, "sample count must be odd");
    assert!(CALLS_PER_SAMPLE > 0, "zero calls per sample found");
    assert!(WARMUP_CALLS > 0, "zero warmup calls found");
    assert!(ALLOCATION_CALLS > 0, "zero allocation calls found");
};

/// Whether allocator activity should be counted.
static COUNT_ALLOCATIONS: AtomicBool = AtomicBool::new(false);
/// Number of allocations observed during measurement.
static ALLOCATIONS: AtomicU64 = AtomicU64::new(0);
/// Number of reallocations observed during measurement.
static REALLOCATIONS: AtomicU64 = AtomicU64::new(0);
/// Sum of requested new sizes for reallocations during measurement.
static REALLOCATED_BYTES: AtomicU64 = AtomicU64::new(0);

/// System allocator wrapper that records allocation activity.
struct CountingAllocator;

// SAFETY: Every operation delegates directly to the system allocator with the
// original pointer and layout. The counters do not affect allocation behavior.
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        if COUNT_ALLOCATIONS.load(Ordering::Relaxed) {
            let _previous_allocations = ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
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
            let _previous_reallocated_bytes = REALLOCATED_BYTES.fetch_add(
                u64::try_from(new_size).expect("allocation size should fit u64"),
                Ordering::Relaxed,
            );
        }
        // SAFETY: The caller supplies the original allocation and new size.
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
/// Allocator used to count benchmark allocation activity.
static GLOBAL: CountingAllocator = CountingAllocator;

/// Repeatedly extracts and consumes coordinates from the fixture.
fn run_calls(
    curr_read: &CurrRead<OnlyAlignDataComplete>,
    record: &Record,
    region: &GenomicBed3,
    calls: usize,
) {
    for _ in 0..calls {
        let coords = curr_read
            .seq_coords_from_ref_coords(black_box(record), black_box(region))
            .expect("fixture should produce coordinates");
        drop(black_box(coords));
    }
}

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
    let mut reader = nanalogue_bam_reader("examples/example_7.sam").expect("fixture should open");
    let record = reader
        .records()
        .next()
        .expect("fixture should contain a record")
        .expect("fixture record should parse");
    let curr_read = CurrRead::default()
        .try_from_only_alignment(&record)
        .expect("fixture alignment should parse");
    let region = Bed3::new(0, 9, 13);
    let expected = [
        Some((true, 0)),
        None,
        None,
        Some((false, 1)),
        Some((true, 2)),
    ];
    let probe = curr_read
        .seq_coords_from_ref_coords(&record, &region)
        .expect("fixture should produce coordinates");
    assert_eq!(probe, expected, "benchmark fixture output changed");
    assert!(
        probe.capacity() > probe.len(),
        "benchmark fixture must produce output with excess capacity"
    );
    let output_len = probe.len();
    let output_capacity = probe.capacity();

    run_calls(&curr_read, &record, &region, WARMUP_CALLS);

    let mut elapsed = Vec::with_capacity(SAMPLES);
    for _ in 0..SAMPLES {
        let start = Instant::now();
        run_calls(&curr_read, &record, &region, CALLS_PER_SAMPLE);
        elapsed.push(start.elapsed());
    }

    elapsed.sort_unstable();
    let median = *elapsed
        .get(SAMPLES / 2)
        .expect("sample count should contain its median");
    let calls_per_second = u128::try_from(CALLS_PER_SAMPLE)
        .expect("call count should fit u128")
        .checked_mul(1_000_000_000)
        .expect("throughput numerator should not overflow")
        .checked_div(median.as_nanos())
        .expect("median duration should be nonzero");
    ALLOCATIONS.store(0, Ordering::Relaxed);
    REALLOCATIONS.store(0, Ordering::Relaxed);
    REALLOCATED_BYTES.store(0, Ordering::Relaxed);
    COUNT_ALLOCATIONS.store(true, Ordering::Relaxed);
    run_calls(&curr_read, &record, &region, ALLOCATION_CALLS);
    COUNT_ALLOCATIONS.store(false, Ordering::Relaxed);

    let allocation_calls =
        u64::try_from(ALLOCATION_CALLS).expect("allocation call count should fit u64");
    let allocations = ALLOCATIONS.load(Ordering::Relaxed);
    let reallocations = REALLOCATIONS.load(Ordering::Relaxed);
    let reallocated_bytes = REALLOCATED_BYTES.load(Ordering::Relaxed);

    println!("seq_coords_from_ref_coords capacity benchmark");
    println!("samples={SAMPLES} calls_per_sample={CALLS_PER_SAMPLE}");
    println!(
        "median_ns={} throughput={calls_per_second} calls/s",
        median.as_nanos()
    );
    println!("output_len={output_len} output_capacity={output_capacity}");
    println!(
        "allocations/call={} reallocations/call={} reallocated_bytes/call={}",
        allocations
            .checked_div(allocation_calls)
            .expect("calls are nonzero"),
        reallocations
            .checked_div(allocation_calls)
            .expect("calls are nonzero"),
        reallocated_bytes
            .checked_div(allocation_calls)
            .expect("calls are nonzero"),
    );

    let min = elapsed.first().expect("samples should not be empty");
    let max = elapsed.last().expect("samples should not be empty");
    let spread = max.saturating_sub(*min);
    println!(
        "sample_min_ns={} sample_max_ns={} spread_ns={}",
        min.as_nanos(),
        max.as_nanos(),
        spread.as_nanos()
    );
}
