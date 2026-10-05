#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Deterministic in-memory benchmark for MM/ML parsing.
//!
//! The scenarios cover short reads, dense and sparse calls, implicit calls, multiple
//! modifications sharing one canonical base, four canonical bases, and sparse `N`
//! annotations on an ultra-long read. Record construction and fixture validation happen
//! before timing. Each reported sample measures only repeated calls to
//! [`nanalogue_mm_ml_parser`], including its output allocations.

use nanalogue_core::nanalogue_mm_ml_parser;
use rust_htslib::bam::{Record, record::Aux};
use std::hint::black_box;
use std::time::{Duration, Instant};

/// Number of independently timed samples per scenario.
const SAMPLES: u8 = 9;

/// One MM group in a generated record.
#[derive(Clone, Copy)]
struct GroupSpec {
    /// Canonical base named by the MM group.
    base: u8,
    /// Modification code named by the MM group.
    code: &'static str,
    /// MM suffix (`?`, `.`, or empty).
    suffix: &'static str,
    /// Retain one call for every this many candidate bases.
    stride: usize,
    /// ML probability assigned to every retained call.
    probability: u8,
}

/// One parser workload.
#[derive(Clone, Copy)]
struct FixtureSpec {
    /// Name printed in benchmark output.
    name: &'static str,
    /// Read length.
    seq_len: u32,
    /// Repeated sequence motif.
    pattern: &'static [u8],
    /// MM groups attached to the read.
    groups: &'static [GroupSpec],
    /// Parses performed in each timed sample.
    batch_size: u32,
}

/// Dense C+m calls.
const DENSE_C: GroupSpec = GroupSpec {
    base: b'C',
    code: "m",
    suffix: "?",
    stride: 1,
    probability: 200,
};
/// Sparse C+m calls.
const SPARSE_C: GroupSpec = GroupSpec {
    base: b'C',
    code: "m",
    suffix: "?",
    stride: 1_000,
    probability: 200,
};
/// Sparse implicit C+m calls.
const IMPLICIT_C: GroupSpec = GroupSpec {
    base: b'C',
    code: "m",
    suffix: ".",
    stride: 1_000,
    probability: 200,
};
/// Dense C+h calls used with [`DENSE_C`].
const DENSE_H: GroupSpec = GroupSpec {
    base: b'C',
    code: "h",
    suffix: "?",
    stride: 1,
    probability: 50,
};
/// Dense A+a calls.
const DENSE_A: GroupSpec = GroupSpec {
    base: b'A',
    code: "a",
    suffix: "?",
    stride: 1,
    probability: 200,
};
/// Dense G+o calls.
const DENSE_G: GroupSpec = GroupSpec {
    base: b'G',
    code: "o",
    suffix: "?",
    stride: 1,
    probability: 200,
};
/// Dense T+T calls.
const DENSE_T: GroupSpec = GroupSpec {
    base: b'T',
    code: "T",
    suffix: "?",
    stride: 1,
    probability: 200,
};
/// Sparse N+n calls.
const SPARSE_N: GroupSpec = GroupSpec {
    base: b'N',
    code: "n",
    suffix: "?",
    stride: 1_000,
    probability: 200,
};

/// Groups for single-modification dense scenarios.
static DENSE_GROUPS: [GroupSpec; 1] = [DENSE_C];
/// Groups for a sparse explicit scenario.
static SPARSE_GROUPS: [GroupSpec; 1] = [SPARSE_C];
/// Groups for a sparse implicit scenario.
static IMPLICIT_GROUPS: [GroupSpec; 1] = [IMPLICIT_C];
/// Two modifications sharing one canonical base.
static SAME_BASE_GROUPS: [GroupSpec; 2] = [DENSE_C, DENSE_H];
/// Dense modifications on all four canonical bases.
static FOUR_BASE_GROUPS: [GroupSpec; 4] = [DENSE_A, DENSE_C, DENSE_G, DENSE_T];
/// Sparse all-base annotations.
static N_GROUPS: [GroupSpec; 1] = [SPARSE_N];

/// Benchmark scenarios. Batch sizes keep timer noise low without making a run excessively long.
static FIXTURES: [FixtureSpec; 7] = [
    FixtureSpec {
        name: "short-dense-C",
        seq_len: 150,
        pattern: b"ACGT",
        groups: &DENSE_GROUPS,
        batch_size: 20_000,
    },
    FixtureSpec {
        name: "long-dense-C",
        seq_len: 100_000,
        pattern: b"ACGT",
        groups: &DENSE_GROUPS,
        batch_size: 20,
    },
    FixtureSpec {
        name: "long-sparse-C",
        seq_len: 100_000,
        pattern: b"ACGT",
        groups: &SPARSE_GROUPS,
        batch_size: 100,
    },
    FixtureSpec {
        name: "long-implicit-sparse-C",
        seq_len: 100_000,
        pattern: b"ACGT",
        groups: &IMPLICIT_GROUPS,
        batch_size: 20,
    },
    FixtureSpec {
        name: "long-two-mods-on-C",
        seq_len: 100_000,
        pattern: b"ACGT",
        groups: &SAME_BASE_GROUPS,
        batch_size: 10,
    },
    FixtureSpec {
        name: "long-four-bases",
        seq_len: 100_000,
        pattern: b"ACGT",
        groups: &FOUR_BASE_GROUPS,
        batch_size: 5,
    },
    FixtureSpec {
        name: "ultra-sparse-N",
        seq_len: 500_000,
        pattern: b"ACGT",
        groups: &N_GROUPS,
        batch_size: 100,
    },
];

/// Counts candidate bases for one MM group.
#[expect(
    clippy::naive_bytecount,
    reason = "fixture construction is untimed and does not justify another dependency"
)]
fn candidate_count(sequence: &[u8], base: u8) -> usize {
    if base == b'N' {
        sequence.len()
    } else {
        sequence
            .iter()
            .filter(|candidate| **candidate == base)
            .count()
    }
}

/// Returns the number of output annotations expected for one group.
fn expected_annotations(candidates: usize, group: GroupSpec) -> usize {
    if group.suffix == "?" {
        candidates.div_ceil(group.stride)
    } else {
        candidates
    }
}

/// Adds one generated MM group and its ML values.
fn add_group(
    position_tag: &mut String,
    probability_tag: &mut Vec<u8>,
    candidates: usize,
    group: GroupSpec,
) {
    position_tag.push(char::from(group.base));
    position_tag.push('+');
    position_tag.push_str(group.code);
    position_tag.push_str(group.suffix);

    let mut previous: Option<usize> = None;
    for occurrence in (0..candidates).step_by(group.stride) {
        let distance = previous.map_or(occurrence, |prior| {
            occurrence
                .checked_sub(prior)
                .and_then(|difference| difference.checked_sub(1))
                .expect("selected candidate positions increase by at least one")
        });
        position_tag.push(',');
        position_tag.push_str(&distance.to_string());
        probability_tag.push(group.probability);
        previous = Some(occurrence);
    }
    position_tag.push(';');
}

/// Builds and validates one deterministic in-memory `ModBAM` record.
fn fixture(spec: FixtureSpec) -> (Record, usize) {
    let seq_len = usize::try_from(spec.seq_len).expect("benchmark length fits usize");
    let sequence: Vec<u8> = spec.pattern.iter().copied().cycle().take(seq_len).collect();
    let qualities = vec![30; seq_len];
    let mut position_tag = String::new();
    let mut probability_tag = Vec::new();
    let mut expected_total = 0usize;

    for &group in spec.groups {
        let candidates = candidate_count(&sequence, group.base);
        add_group(&mut position_tag, &mut probability_tag, candidates, group);
        expected_total = expected_total
            .checked_add(expected_annotations(candidates, group))
            .expect("fixture annotation total fits usize");
    }

    let mut record = Record::new();
    record.set(spec.name.as_bytes(), None, &sequence, &qualities);
    record.set_flags(4);
    record
        .push_aux(b"MM", Aux::String(&position_tag))
        .expect("MM tag should fit the record");
    record
        .push_aux(b"ML", Aux::ArrayU8((&probability_tag).into()))
        .expect("ML tag should fit the record");

    let parsed = nanalogue_mm_ml_parser(&record, |&_| true, |&_| true, |&_, &_, &_| true, 0)
        .expect("generated fixture should parse");
    let observed_total = parsed
        .base_mods
        .iter()
        .map(|base_mod| base_mod.ranges.annotations().len())
        .sum::<usize>();
    assert_eq!(
        observed_total, expected_total,
        "generated fixture annotation count changed"
    );
    assert_eq!(
        parsed.base_mods.len(),
        spec.groups.len(),
        "generated fixture group count changed"
    );
    (record, expected_total)
}

/// Parses one fixture repeatedly and returns an observable annotation checksum.
fn run_batch(record: &Record, batch_size: u32) -> usize {
    let mut checksum = 0usize;
    for _ in 0..batch_size {
        let parsed = nanalogue_mm_ml_parser(
            black_box(record),
            |&_| true,
            |&_| true,
            |&_, &_, &_| true,
            0,
        )
        .expect("generated fixture should parse");
        checksum = checksum
            .checked_add(
                parsed
                    .base_mods
                    .iter()
                    .map(|base_mod| base_mod.ranges.annotations().len())
                    .sum::<usize>(),
            )
            .expect("benchmark checksum fits usize");
        drop(black_box(parsed));
    }
    checksum
}

#[expect(
    clippy::print_stdout,
    reason = "the benchmark reports its measurements to stdout"
)]
fn main() {
    println!("MM/ML parser benchmark: samples={SAMPLES}");
    println!(
        "scenario\tseq_len\tgroups\tannotations\tbatch\tmedian_us/read\tthroughput_Mb/s\trange_us/read"
    );

    for &spec in &FIXTURES {
        let (record, annotations) = fixture(spec);
        let expected_checksum = annotations
            .checked_mul(usize::try_from(spec.batch_size).expect("batch size fits usize"))
            .expect("expected checksum fits usize");
        assert_eq!(
            black_box(run_batch(&record, spec.batch_size)),
            expected_checksum,
            "warm-up checksum changed"
        );

        let mut timings = Vec::with_capacity(usize::from(SAMPLES));
        for _ in 0..SAMPLES {
            let start = Instant::now();
            let checksum = black_box(run_batch(&record, spec.batch_size));
            let elapsed = start.elapsed();
            assert_eq!(checksum, expected_checksum, "timed checksum changed");
            timings.push(elapsed);
        }
        timings.sort_unstable();
        let median_index = usize::from(SAMPLES.checked_div(2).expect("sample count is nonzero"));
        let median = timings
            .get(median_index)
            .copied()
            .expect("median sample exists");
        let minimum = timings.first().copied().unwrap_or(Duration::ZERO);
        let maximum = timings.last().copied().unwrap_or(Duration::ZERO);
        let batch = f64::from(spec.batch_size);
        let median_per_read = median.as_secs_f64() / batch;
        let minimum_per_read = minimum.as_secs_f64() / batch;
        let maximum_per_read = maximum.as_secs_f64() / batch;
        let throughput = f64::from(spec.seq_len) * batch / median.as_secs_f64() / 1_000_000.0;

        println!(
            "{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}..{:.3}",
            spec.name,
            spec.seq_len,
            spec.groups.len(),
            annotations,
            spec.batch_size,
            median_per_read * 1_000_000.0,
            throughput,
            minimum_per_read * 1_000_000.0,
            maximum_per_read * 1_000_000.0,
        );
    }
}
