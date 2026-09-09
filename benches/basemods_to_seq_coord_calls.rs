//! Long-read batch benchmark for sequence-coordinate base modification conversion.

use std::alloc::{GlobalAlloc, Layout, System};
use std::hint::black_box;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::time::Instant;

use nanalogue_core::{BaseMod, BaseMods, FiberAnnotation, Ranges, SeqCoordCalls};

/// Sequence length representative of a long nanopore read.
const SEQ_LEN: u32 = 100_000;
/// Number of modification type/strand combinations in each position row.
const MOD_TYPES: u8 = 4;
/// Distance between annotations for each modification type.
const ANNOTATION_STEP: u8 = 4;
/// Number of conversions timed together to reduce timer noise.
const BATCH_SIZE: u32 = 24;
/// Number of independently reported timed batches.
const SAMPLES: u8 = 9;

const _: () = {
    assert!(SEQ_LEN > 0, "zero sequence length found");
    assert!(MOD_TYPES > 0, "zero modification types found");
    assert!(ANNOTATION_STEP > 0, "zero annotation step found");
    assert!(
        MOD_TYPES == ANNOTATION_STEP,
        "interleaved annotations must cover every position"
    );
    assert!(BATCH_SIZE > 0, "zero batch size found");
    assert!(SAMPLES >= 3, "fewer than three samples found");
};

/// System allocator wrapper that can count allocation calls on demand.
struct CountingAllocator;

/// Whether allocator calls should currently be counted.
static COUNT_ALLOCATIONS: AtomicBool = AtomicBool::new(false);
/// Number of allocator calls made while counting was enabled.
static ALLOCATIONS: AtomicUsize = AtomicUsize::new(0);

// SAFETY: Every operation delegates directly to the system allocator with the
// original pointer and layout. The counter does not affect allocation behavior.
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        if COUNT_ALLOCATIONS.load(Ordering::Relaxed) {
            let _previous_allocation_count = ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        }
        // SAFETY: This allocator delegates the unchanged layout to `System`.
        unsafe { System.alloc(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        // SAFETY: This allocator delegates the pointer and its original layout to `System`.
        unsafe { System.dealloc(ptr, layout) }
    }
}

#[global_allocator]
/// Allocator used to measure calls made by conversion without timing the counter.
static GLOBAL: CountingAllocator = CountingAllocator;

/// Builds a long read whose four interleaved annotation sets cover every position.
fn fixture() -> BaseMods {
    let base_mods = (0..MOD_TYPES)
        .map(|mod_type| BaseMod {
            modified_base: b'C',
            strand: if mod_type & 1 == 0 { '+' } else { '-' },
            modification_type: char::from(
                b'a'.checked_add(mod_type).expect("four types fit in ASCII"),
            ),
            ranges: Ranges {
                annotations: (usize::from(mod_type)..usize::try_from(SEQ_LEN).expect("fits usize"))
                    .step_by(usize::from(ANNOTATION_STEP))
                    .map(|pos| FiberAnnotation {
                        pos: u32::try_from(pos).expect("fixture positions fit u32"),
                        ref_pos: None,
                        qual: mod_type.checked_add(1).expect("four quality values fit u8"),
                    })
                    .collect(),
                seq_len: SEQ_LEN,
                reverse: false,
            },
            record_is_reverse: false,
        })
        .collect();
    BaseMods { base_mods }
}

/// Converts one batch and returns a checksum that keeps the result observable.
fn convert_batch(input: &BaseMods) -> usize {
    let mut checksum = 0usize;
    for _ in 0..BATCH_SIZE {
        let calls = black_box(SeqCoordCalls::try_from(black_box(input)).expect("valid fixture"));
        checksum ^= usize::from(
            calls
                .mod_calls(
                    usize::try_from(SEQ_LEN.checked_sub(1).expect("nonzero")).expect("fits usize"),
                )
                .expect("last position exists")
                .get(usize::from(MOD_TYPES.checked_sub(1).expect("nonzero")))
                .copied()
                .expect("last mod type exists"),
        );
        drop(black_box(calls));
    }
    checksum
}

#[expect(
    clippy::print_stdout,
    reason = "the benchmark reports its measurements to stdout"
)]
fn main() {
    let input = fixture();
    // Sanity-check the fixture; this is not exhaustive benchmark validation.
    let probe = SeqCoordCalls::try_from(&input).expect("valid fixture");
    let last_position =
        usize::try_from(SEQ_LEN.checked_sub(1).expect("nonzero")).expect("fits usize");
    let expected_last_row = [0, 0, 0, 4];
    assert_eq!(
        probe.mod_calls(last_position),
        Some(expected_last_row.as_slice()),
        "benchmark fixture output changed"
    );
    drop(probe);
    // Warm up before measuring.
    let _: usize = black_box(convert_batch(&input));

    ALLOCATIONS.store(0, Ordering::Relaxed);
    COUNT_ALLOCATIONS.store(true, Ordering::Relaxed);
    let _: usize = black_box(convert_batch(&input));
    COUNT_ALLOCATIONS.store(false, Ordering::Relaxed);
    let allocations = ALLOCATIONS.load(Ordering::Relaxed);

    println!(
        "config,seq_len={SEQ_LEN},mod_types={MOD_TYPES},annotation_step={ANNOTATION_STEP},batch_size={BATCH_SIZE},samples={SAMPLES}"
    );
    println!(
        "allocations,batch,{allocations},per_read,{}",
        allocations
            .checked_div(usize::try_from(BATCH_SIZE).expect("batch size fits usize"))
            .expect("batch size is nonzero")
    );
    for sample in 1..=SAMPLES {
        let start = Instant::now();
        let checksum = convert_batch(&input);
        let elapsed = start.elapsed();
        let reads_per_second = f64::from(BATCH_SIZE) / elapsed.as_secs_f64();
        println!(
            "sample,{sample},elapsed_ms,{:.3},reads_per_second,{reads_per_second:.3},checksum,{checksum}",
            elapsed.as_secs_f64() * 1_000.0
        );
    }
}
