# Benchmarks

This directory contains reproducible performance measurements rather than
pass/fail performance tests. Results depend on the machine, storage, cache
state, compiler, and enabled features, so compare builds on the same system
instead of treating any recorded runtime as a fixed threshold.

## Rust microbenchmarks

The Rust microbenchmarks generate deterministic inputs in memory and report
multiple timing samples. Run one with:

```console
cargo bench --bench mm_ml_parser
```

Available benchmarks:

| Benchmark | Timed work |
| --- | --- |
| `mm_ml_parser` | MM/ML parsing across short, long, sparse, dense, implicit, multi-modification, and ultra-long records. |
| `basemods_to_seq_coord_calls` | Conversion of long-read base modifications into sequence-coordinate calls, including separate allocation measurements. |
| `gradient_windows` | Gradient calculation over repeated, overlapping windows. |
| `seq_coords_capacity` | Coordinate extraction when output vectors retain spare capacity, including allocation and reallocation measurements. |
| `window_reads_dense` | The public `window_reads` path over dense, heavily overlapping windows. It excludes BAM decompression and filesystem I/O. |

Run all registered benchmarks with `cargo bench`, or replace
`mm_ml_parser` in the example above with another benchmark name.

## End-to-end MM/ML CLI benchmark

`mm_ml_cli` is a harness-free Cargo benchmark so it can accept workload options.
It generates deterministic ModBAM fixtures, runs real `nanalogue` commands,
and writes timings and descriptive statistics as JSON.

Show every fixture and option:

```console
cargo bench --bench mm_ml_cli -- --help
```

Run the default representative fixtures with warm caches:

```console
cargo bench --bench mm_ml_cli
```

Run one workload:

```console
cargo bench --bench mm_ml_cli -- \
  --fixture dense-single \
  --command read-info \
  --runs 5
```

Compare an existing baseline executable with the current checkout:

```console
cargo bench --bench mm_ml_cli -- \
  --baseline-bin /path/to/old/nanalogue \
  --fixture dense-single \
  --command read-info \
  --command read-table-show-mods \
  --runs 7 \
  --cache cold
```

The baseline executable must not be the current driver's release-build output;
copy it elsewhere first. Baseline and candidate runs alternate order. Reported
speedup is baseline mean divided by candidate mean.

By default, fixtures and `results.json` are stored under
`target/nanalogue-mm-ml-cli-benchmark/`. Matching fixtures are reused; pass
`--regenerate` to rebuild them. Use `--work-dir` and `--output` to choose other
locations.

### Cache modes

- `--cache warm` performs an unmeasured warm-up before timed runs.
- `--cache cold` uses `POSIX_FADV_DONTNEED` and `mincore` to evict and verify
  the BAM, index, and executables before each timed run. This mode is supported
  on 64-bit x86 and ARM Linux, requires files owned by the current user, and
  does not require root.

### Large `window-dens` workload

The one-billion-base `window-dens` workload is available through the same
benchmark driver:

```console
cargo bench --bench mm_ml_cli -- \
  --fixture window-dens-e2e \
  --command window-dens \
  --cache cold
```

This large benchmark generates one billion read bases. Before timing, it
validates the read count and lengths, contigs, `T+T` modification type, and the
expected number of output lines. Results include median throughput in Mb/s.
