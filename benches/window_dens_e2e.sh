#!/usr/bin/env bash
set -euo pipefail

# End-to-end benchmark for dense window-dens processing.
#
# Environment overrides:
#   THREADS=2            worker threads passed to nanalogue (default: 2)
#   RUNS=5               number of timed runs (default: 5)
#   BENCH_DIR=/path      new disk-backed directory for the fixture and results

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"
export LC_ALL=C

threads="${THREADS:-2}"
runs="${RUNS:-5}"

case "${1:-}" in
  "") ;;
  -h|--help)
    sed -n '3,10p' "$0"
    exit 0
    ;;
  *)
    echo "usage: $0 [--help]" >&2
    exit 2
    ;;
esac

for value_name in threads runs; do
  value="${!value_name}"
  if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
    echo "$value_name must be a positive integer, got: $value" >&2
    exit 2
  fi
done

for command in cargo fincore python3 rustc; do
  if ! command -v "$command" >/dev/null 2>&1; then
    echo "required command not found: $command" >&2
    exit 1
  fi
done

if [[ -n "${BENCH_DIR:-}" ]]; then
  bench_dir="$BENCH_DIR"
  if ! mkdir "$bench_dir"; then
    echo "BENCH_DIR must not already exist: $bench_dir" >&2
    exit 1
  fi
else
  mkdir -p "$repo_root/target"
  bench_dir="$(mktemp -d "$repo_root/target/nanalogue-window-dens-benchmark.XXXXXX")"
fi
config="$bench_dir/config.json"
bam="$bench_dir/simulated.bam"
fasta="$bench_dir/simulated.fasta"
results="$bench_dir/window-dens-times.tsv"

cache_probe="$bench_dir/cache-probe"
python3 - "$cache_probe" <<'PY'
import os
import sys

if not hasattr(os, "posix_fadvise") or not hasattr(os, "POSIX_FADV_DONTNEED"):
    raise SystemExit("Python does not expose POSIX_FADV_DONTNEED on this platform")

with open(sys.argv[1], "wb") as handle:
    handle.write(b"0" * 4096)
    handle.flush()
    os.fsync(handle.fileno())
    os.posix_fadvise(handle.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)
PY
probe_resident="$(fincore --bytes --noheadings --output RES "$cache_probe" | tr -d '[:space:]')"
rm "$cache_probe"
if [[ "$probe_resident" != 0 ]]; then
  echo "BENCH_DIR does not support cache eviction; choose a disk-backed filesystem: $bench_dir" >&2
  exit 1
fi

cat > "$config" <<'JSON'
{
  "contigs": {
    "number": 12,
    "len_range": [1000000, 1000000],
    "repeated_seq": "ACGT"
  },
  "reads": [{
    "number": 100000,
    "mapq_range": [60, 60],
    "base_qual_range": [30, 30],
    "len_range": [0.01, 0.01],
    "mods": [{
      "base": "T",
      "is_strand_plus": true,
      "mod_code": "T",
      "mm_suffix": "?",
      "win": [2500],
      "mod_range": [[0.8, 1.0]]
    }]
  }],
  "seed": 42
}
JSON

echo "Building release binaries..."
host_target="$(rustc -vV | awk '$1 == "host:" { print $2 }')"
if [[ -z "$host_target" ]]; then
  echo "could not determine the host Rust target" >&2
  exit 1
fi
cargo build --locked --quiet --release --bins --target "$host_target" --target-dir "$repo_root/target"
nanalogue="$repo_root/target/$host_target/release/nanalogue"
simulator="$repo_root/target/$host_target/release/nanalogue_sim_bam"

echo "Simulating 12 x 1 Mb contigs and 100,000 x 10 kb reads..."
"$simulator" "$config" "$bam" "$fasta"

echo "Validating fixture..."
stats="$("$nanalogue" read-stats --threads "$threads" "$bam")"
read_count="$(awk -F '\t' '$1 ~ /^n_(primary_alignments|secondary_alignments|supplementary_alignments|unmapped_reads)$/ { total += $2 } END { print total + 0 }' <<< "$stats")"
if [[ "$read_count" != 100000 ]]; then
  echo "expected 100000 reads, found: $read_count" >&2
  exit 1
fi
for key in seq_len_mean seq_len_min seq_len_max; do
  value="$(awk -F '\t' -v key="$key" '$1 == key { print $2 }' <<< "$stats")"
  if [[ "$value" != 10000 ]]; then
    echo "expected $key=10000, found: ${value:-missing}" >&2
    exit 1
  fi
done

peek="$("$nanalogue" peek "$bam")"
contig_count="$(awk -F '\t' '$1 ~ /^contig_[0-9]+$/ && $2 == 1000000 { count++ } END { print count + 0 }' <<< "$peek")"
if [[ "$contig_count" != 12 ]] || ! grep -qx 'T+T' <<< "$peek"; then
  echo "fixture does not contain twelve 1 Mb contigs and T+T modifications" >&2
  exit 1
fi

echo "Validating window output..."
line_count="$("$nanalogue" window-dens --win 300 --step 300 --threads "$threads" "$bam" | wc -l)"
line_count="${line_count//[[:space:]]/}"
if [[ "$line_count" != 800001 ]]; then
  echo "expected 800001 output lines including the header, found: $line_count" >&2
  exit 1
fi

printf 'run\tthreads\telapsed_s\n' > "$results"
echo "Running $runs cold-cache trials with $threads threads..."
for ((run = 1; run <= runs; run++)); do
  python3 - "$bam" "$bam.bai" "$nanalogue" <<'PY'
import os
import sys

if not hasattr(os, "posix_fadvise") or not hasattr(os, "POSIX_FADV_DONTNEED"):
    raise SystemExit("Python does not expose POSIX_FADV_DONTNEED on this platform")

for path in sys.argv[1:]:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.posix_fadvise(descriptor, 0, 0, os.POSIX_FADV_DONTNEED)
    finally:
        os.close(descriptor)
PY

  resident="$(fincore --bytes --noheadings --output RES "$bam" | tr -d '[:space:]')"
  if [[ "$resident" != 0 ]]; then
    echo "BAM cache eviction failed before run $run; resident size: $resident" >&2
    exit 1
  fi

  elapsed="$(python3 - "$nanalogue" "$threads" "$bam" <<'PY'
import subprocess
import sys
import time

command = [
    sys.argv[1],
    "window-dens",
    "--win", "300",
    "--step", "300",
    "--threads", sys.argv[2],
    sys.argv[3],
]
start = time.perf_counter()
subprocess.run(command, stdout=subprocess.DEVNULL, check=True)
print(f"{time.perf_counter() - start:.6f}")
PY
)"
  printf '%d\t%s\t%s\n' "$run" "$threads" "$elapsed" >> "$results"
  printf 'run %d/%d: %s s\n' "$run" "$runs" "$elapsed"
done

python3 - "$results" <<'PY'
import csv
import statistics
import sys

with open(sys.argv[1], newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

timings = [float(row["elapsed_s"]) for row in rows]
median = statistics.median(timings)
print(f"median: {median:.3f} s")
print(f"throughput: {1_000.0 / median:.2f} Mb/s")
print(f"results: {sys.argv[1]}")
PY
printf 'benchmark directory: %s\n' "$bench_dir"
