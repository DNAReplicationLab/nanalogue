# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- Contig names now support PanSN names like `HG002#1#chr1` and bare names with
  colons like `HLA-A*01:01:01:01`. Names starting with `#` or containing a
  hyphen after the last colon are unsupported; the latter clashes with
  region syntax. Read-ID rules are unchanged.

### Fixed
- `find-modified-reads dens-range-above` allows an absolute `f32::EPSILON`
  tolerance when comparing density ranges to `--min-range`, so subtraction
  rounding does not exclude ranges at the requested boundary. Positive
  thresholds always reject zero density range, including thresholds at or
  below epsilon; `--min-range 0` still accepts zero range.
- Record capacity can now reach exactly 32 MiB, so HTSlib's rounded allocations
  are no longer rejected for being one byte over the old limit.
- Applied modification-region filtering even when the requested region
  contains the entire aligned span of a read, including an exact match.
  Terminal soft-clipped modification calls are now excluded in these cases;
  insertions between retained calls are still kept.
- Invalid MM-strand errors now display the offending character and explain
  that the strand must be `+` or `-`.
- Corrected the `--include-zero-len` documentation: zero-length sequences
  bypass minimum sequence-length filtering only; minimum alignment-length
  and all other read filters still apply.
- Reject reversed `--low`/`--high` thresholds.

## [0.2.0] - 2026-10-06

### Added
- Added CRAM output to `nanalogue_sim_bam`. The simulator now selects BAM or
  CRAM from the output extension and creates the corresponding BAI or CRAI;
  CRAM 3.1 output uses a generated external FASTA reference and FAI. The public
  `write_cram_denovo` API and detailed simulator CLI configuration help are
  also new.
- Added `mm_suffix` and `drop` simulation options. Simulated MM groups can use
  `?`, `.`, or no suffix, and selected canonical bases can be omitted to
  generate non-zero MM distances and omit their corresponding ML entries.
- Added the public `mm_groups`/`ParsedMmGroup` MM-tag parser, resource-limit
  constants, and reusable input-validation helpers.
- Added reusable MM/ML parser and end-to-end CLI benchmarks with deterministic
  ModBAM fixtures, cache controls, baseline comparisons, and JSON results.

### Changed
- Sped up modification parsing by scanning BAM's packed sequence directly,
  parsing MM distances without intermediate strings, flattening stored calls,
  caching candidate positions, and using a sparse CIGAR-segment coordinate map.
  End-to-end benchmarks improved representative sparse long-read commands by
  more than 2×, with gains varying by read and modification structure.
- Sped up dense, overlapping `window-dens` and `window-grad` workloads by
  reusing rolling reference-coordinate and base-quality calculations. Window
  `basecall_qual` values can differ by one Phred unit at rounding boundaries
  because of the changed floating-point summation order.
- Made Polars an optional, default-off dependency. Library consumers using
  `curr_reads_to_dataframe`, `reads_table::run_df`, `window_reads::run_df`, or
  Polars-specific error handling must explicitly enable the `polars` Cargo
  feature. CLI functionality is unchanged. Several Polars output columns,
  including sequence lengths, alignment bounds, and query positions, now use
  `UInt32` instead of `UInt64`.
- Reduced the dependency surface by replacing `thiserror`, `regex`, `csv`, and
  `itertools` usage, and by incorporating only the required `bedrs` and
  `openssl-probe` functionality. Enabled libdeflate in `rust-htslib`.
- Hardened analysis of untrusted or pathological inputs with explicit limits
  for record capacity/count, CIGAR operations, sequence and contig lengths,
  modification types and annotations, MM/ML data, identifiers, paths, regions,
  sequencing-summary files, and `peek` records. This includes a maximum of 20
  modification groups per record and stricter ASCII/character rules for read
  IDs, contigs, and paths; displayed user-controlled error context is truncated
  after 300 Unicode characters. `find-modified-reads`, table/window commands,
  and `peek` now error on zero input records, including when applicable filters
  remove every record. `read-info` continues to return an empty JSON array, and
  `read-stats` returns zero-valued summaries. TSV/DataFrame window output also
  errors when records are present but no windows are emitted.
- Applied the HTTP, HTTPS, and FTP URL-scheme allow-list consistently to public
  URL-reader helpers and deserialization. Hardened the install instructions and
  script to require HTTPS and TLS 1.2 when using curl (and HTTPS when using
  wget).
- Sequencing-summary parsing now uses bounded plain TSV input rather than CSV
  quoting rules. Comments are accepted only before the header, and supplied
  summaries must contain a header and at least one read row.
- Simulator output paths must end case-insensitively in `.bam` or `.cram`, and
  simulated unmapped reads now use MAPQ 0 rather than 255.
- Breaking for library users: many coordinates and lengths now use bounded
  `u32` types in `InputBam`, `InputWindowing`, `GenomicRegion`, simulation
  configuration and helpers, `BamPreFilt`, `OrdPair`, alignments, modification
  annotations, and relevant error payloads. Conversion from `&GenomicRegion`
  to `bam::FetchDefinition` is now infallible (`From` rather than `TryFrom`).
  The crate now exposes its own `Bed3`/`StrandedBed3` types and `GenomicBed3`
  aliases; interval constructors are fallible and reject reversed coordinates.
  `GenomicRegion` and `AlignmentInfo` no longer implement `Default`. Only
  32-bit and 64-bit targets are supported.
- Breaking for library users: `FiberAnnotation` is now an 8-byte packed value
  representing single query/reference positions rather than intervals; its old
  end, length, and extra-column fields were removed. Construct it with
  `FiberAnnotation::try_new(pos, qual, ref_pos)` and use `pos()`, `qual()`, and
  `ref_pos()` accessors. `FiberAnnotations` fields are also private and its
  coordinate/quality accessors return iterators (`starts()` and
  `reference_starts()` become `pos()` and `ref_pos()`); `from_annotations` now
  returns `Result`. The maximum supported read length is now 2^24 - 1 bases
  (about 16.7 megabases), and the maximum contig length is `u32::MAX - 1`.
- Breaking for library users: callbacks passed to `window_reads::run`,
  `run_json`, and `run_df` now return `Result<Option<_>, Error>`: `Some` emits
  a window, `None` omits it, and errors are propagated instead of logged and
  skipped. Accordingly,
  `threshold_and_mean_and_thres_win` now returns `Ok(None)` rather than an
  error for a window below the requested density, and the obsolete public
  `Error::WindowDensBelowThres` variant was removed.
- Breaking for library users: `TempBamSimulation::new` now takes an
  `AlignmentFormat`, and `TempBamSimulation` no longer implements
  `Deserialize`. `SeqCoordCalls::mod_calls` now returns `&[u8]` rather than
  `Option<&[u8]>` and panics for an out-of-range position.
- Breaking for library users: external errors wrapped by `Error` now use boxed
  payloads, affecting direct construction and nested pattern matching;
  `Error::CsvError` was removed.
- (Project tooling, not code) Expanded locked and unlocked dependency CI to
  macOS and both default/all-feature builds, strengthened repository guardrails,
  and updated the pinned Zig release toolchain with download verification and
  a libdeflate compatibility workaround.

### Removed
- Removed the public `get_u8_tag` and
  `utils::filter_by_ref_coords::WindowState` APIs. The simulation helpers
  `generate_random_dna_modification` and `generate_reads_denovo` are now
  private.

### Fixed
- Corrected `read-stats` median and N50 boundary calculations, including exact
  halfway totals and zero-length values.
- Made MM/ML parsing reject malformed records instead of panicking or silently
  accepting them. This includes non-canonical numeric modification codes,
  oversized headers, gaps above 128,000,000 bases, zero-length sequences,
  mismatched MM/ML counts, MN tags that disagree with sequence length, and
  mapped reads without a CIGAR. Canonical `U` is handled as BAM-encoded `T`,
  and MM/ML parsing handles padding (`P`) CIGAR operations without panicking.
- Enforced sorted query coordinates and strictly monotonic mapped reference
  coordinates in modification annotations, including during deserialization,
  so reference filtering and window bounds cannot silently use invalid order.
- Applied complete validation consistently to parsed, converted, built, and
  deserialized genomic regions and intervals. Alignment intervals must be
  non-empty; BED intervals may be empty but cannot be reversed.
- Fixed probability-exclusion bands to reject exactly the representable
  byte-valued ML probabilities within the inclusive bounds; a band containing
  no representable probability now excludes nothing.
- Deferred each `read-info` JSON delimiter until its record has been parsed and
  rendered successfully, avoiding a dangling object for a failing first
  record. Moved missing-index warnings from stdout to stderr so they no longer
  corrupt structured output.
- Gradient commands now reject one-element windows, and invalid
  sequence-display flag combinations are rejected by the CLI.
- Fixed simulator validation and cleanup: empty window/modification schedules
  are rejected, failed temporary simulations remove their directories, and
  aliased or hard-linked output paths cannot overwrite one another.
- BAM/CRAM writers now reject unusable index paths and out-of-range CRAM thread
  counts before creating alignment output.
- Explicitly flushed subcommand output and propagated flush failures, preventing
  buffered-output errors from being silently lost.
- Read-ID filter files now permit comments only before the first ID and reject
  blank/whitespace-only lines and malformed IDs.
- Fixed HTTPS certificate discovery on macOS by checking the standard macOS
  certificate bundle path.
- (Project tooling, not code) Fixed `install.sh` to use the ARM release archive
  name produced by the release workflow.
- (GitHub workflow, not code) Docker publication dry runs now authenticate with
  Docker Hub so invalid release credentials fail before a real push.
- Updated the locked networking dependencies to `h2` 0.4.16, `rustls` 0.23.45,
  and `rustls-webpki` 0.103.15 to address security advisories.
- Updated locked `event-listener` from 5.4.1 to 5.4.2 to address
  `RUSTSEC-2026-0221`, and replaced the yanked `chacha20` 0.10.0 with 0.10.2.

## [0.1.11] - 2026-05-16

### Changed
- (GitHub workflow, not code) Switched release binary builds to `cargo zigbuild` and expanded the matrix to cover additional Linux targets/compatibility tiers
- (Project tooling, not code) Updated `install.sh` to match the expanded artifact names and architecture aliases, and refreshed the install-script checksum file
- (Documentation, not code) Updated `README.md` install guidance to point at the GitHub Actions artifact section and list the current artifact names/suffixes
- Lowered the reverse-complement pre-reservation cap from 4 GiB to 3 GiB
  to avoid overflowing `usize` on 32-bit platforms
- Unvendored `hts-sys` and `rust-htslib`, updated third-party notices, and
  refreshed `Cargo.lock` to use crates.io directly

### Fixed
- (Project packaging, not code) Stopped ignoring `CLAUDE.md` in git and
  excluded both `AGENTS.md` and `CLAUDE.md` from crate packaging
- (GitHub workflow, not code) Exchanged GitHub OIDC tokens for temporary
  crates.io publish tokens during crate publish dry runs and real publishes

## [0.1.10] - 2026-05-14

### Added
- `CurrRead` now exposes mapping quality (`mapq`) in library output, the
  `Display` implementation, and the `read-info` subcommand; `read-info` also
  reports when MAPQ is unavailable
- Remote BAM retrieval coverage now includes a dedicated test plus CI
  execution for remote URL handling
- (GitHub workflow, not code) CI workflow to test sister packages
  (pynanalogue, nanalogue-node) against current commit using Cargo's
  `[patch.crates-io]` mechanism
- (GitHub workflow, not code) Weekly documentation audit automation plus
  follow-up workflow hardening to better constrain automated doc updates
- (Project tooling, not code) Repository guardrails including git hooks and
  install-script integrity checks

### Changed
- Vendored selected fibertools-derived types, DNA complement helpers, and UUID
  helpers to reduce external dependency surface and improve build provenance
- Renamed reference-coordinate filtering APIs from
  `FilterByRefCoords`/`filter_by_ref_pos` to
  `FilterModsByRefCoords`/`filter_mods_by_ref_pos`
- Hardened dependency and workflow safety by pinning GitHub Actions to commit
  SHAs, tightening dependency safety checks, and applying additional workflow
  cleanups from zizmor reviews
- Applied additional hardening fixes using reports from Opus, and updated
  `ethnum` to avoid a docs.rs nightly crash problem
- Hoisted `target_names()` out of the `peek` loop for lower repeated
  overhead, and silence broken-pipe exits instead of surfacing them as errors
- Refined documentation and README links, naming, and usage details through
  repeated doc audit passes

### Fixed
- Removed multiple panic paths across BAM parsing, read processing,
  simulation code, and missing-quality handling; also tightened trait
  contracts and broader error handling
- Fixed rand compatibility issues encountered during vendoring and toolchain
  updates, including follow-up warning cleanups
- Made libcurl/OpenSSL environment variable initialization safer and
  tightened the SSL initialization contract for remote access paths
- Fixed an unlikely `i64` overflow path in `read_utils`

## [0.1.9] - 2026-02-18

### Added
- JSON output for `window-reads` via `run_json`: per-read structured output with alignment info, modification tables, and windowed data
- Stochastic tests for `run_json` covering reads with mods, without mods, non-perfectly aligned reads, two modification types, and zero-read edge cases

### Changed
- Refactored windowing logic in `window_reads` into reusable `compute_windowed_mod_data` function shared by TSV and JSON paths
- Added `Serialize` and `serde(try_from = "f32")` to `F32AbsValAtMost1` for JSON serialization and safe deserialization
- Updated packages in Cargo.lock

## [0.1.8] - 2026-02-12

### Added
- Optional `--sample-seed` for reproducible read subsampling: hash-based deterministic filtering on read name + seed ([`e8c576f`](https://github.com/DNAReplicationLab/nanalogue/commit/e8c576f3455312b7be1eb81922047319eefae5fb))
- Optional `seed` field in `SimulationConfig` for reproducible BAM simulation output ([`e8c576f`](https://github.com/DNAReplicationLab/nanalogue/commit/e8c576f3455312b7be1eb81922047319eefae5fb))
- (GitHub workflow, not code) Docker image push workflow for multi-arch images to Docker Hub ([`7c2388e`](https://github.com/DNAReplicationLab/nanalogue/commit/7c2388e22053248e0089707a26b1cd3d1c14af18))
- (GitHub workflow, not code) workflow_call trigger for publish-crates so release can invoke it ([`e6f7e2d`](https://github.com/DNAReplicationLab/nanalogue/commit/e6f7e2dab3fae7fd475a3d3b51dc9d43c5617901))
- (GitHub workflow, not code) Docker push and crate publish wired into release pipeline ([`e9901d2`](https://github.com/DNAReplicationLab/nanalogue/commit/e9901d23608d4fdc246f101e2db111724b6cd1f5))

### Fixed
- (GitHub workflow, not code) Fixed release workflow to upload artifacts as zip files instead of individual files
- (GitHub workflow, not code) Fixed publish-crates workflow to use environment variable for cargo registry token
- (GitHub workflow, not code) Release workflow now runs all CI tests before building artifacts
- (GitHub workflow, not code) Release workflow verifies Cargo.toml and Cargo.lock versions match the release tag
- (GitHub workflow, not code) Install script test now runs after release artifacts are uploaded (instead of racing with them)
- (GitHub workflow, not code) Fixed Docker push workflow zip extraction and image tag collision ([`e72e582`](https://github.com/DNAReplicationLab/nanalogue/commit/e72e582))

### Changed
- Updated packages in Cargo.lock

## [0.1.7] - 2026-02-03

### Added
- Install script (`install.sh`) for quick installation of pre-built binaries on macOS and Linux
- GitHub workflows for automated release artifact uploads and crates.io publishing

### Fixed
- Fixed `peek` to skip zero-length reads instead of returning an error

### Changed
- Improved README install, update, and Docker snippets
- Updated `clap` from 4.5.54 to 4.5.56
- Updated `openssl-probe` from 0.2.0 to 0.2.1
- Updated `uuid` from 1.18.1 to 1.20.0

## [0.1.6] - 2026-01-19

### Fixed
- Fixed `BamPreFilt` so that full region filtering works when whole contig filtering is requested
- Fixed `GenomicRegion::try_to_bed3` to return actual contig lengths instead of `u64::MAX` for open-ended regions

### Changed
- Ran `cargo update` to update dependencies

## [0.1.5] - 2026-01-18

### Changed
- Improved struct documentation to reference Builder patterns in `src/simulate_mod_bam.rs`
- Updated `openssl-probe` to 0.2.0 and switched to `try_init_openssl_env_vars()`
- Updated `bio` to 3.0.0
- Fixes `vergen` at "=9.0.6" for `fibertools-rs` building

### Fixed
- Fixed mismatch generation bug in `src/simulate_mod_bam.rs`: replaced `partial_shuffle` with `choose_multiple` to preserve position information during base mutations

### Added
- New test `mismatch_mod_check()` in `src/simulate_mod_bam.rs` to validate that simulated mismatches correctly affect modification reference positions while preserving modification quality
- Integration tests for `read-stats`, `read-info`, and `window-grad` commands with example fixtures
- Integration tests for all `find-modified-reads` subcommands using BAM simulation

## [0.1.4] - 2026-01-11

### Added
- Github actions to check build without Cargo.lock.
- Dependabot on Github.

### Changed
- Updates cargo packages: `cargo install` failing in previous version although `cargo install --locked` worked.

## [0.1.3] - 2026-01-09

### Added
- Dockerfiles based on distroless
- CI/CD workflows in github to make binaries and docker images.
- nanalogue `peek` to get (from header) contigs, contig lengths, and (from first 100 records) types of mods present.
- Support for fallback MM/ML tag variants: parser now accepts both standard MM/ML tags and Mm/Ml capitalization variants used by some sequencing technologies.

## [0.1.2] - 2026-01-02

### Added
- CI workflows to build and check code
- Tests to increase coverage to > 92%.
- Adds more documentation and link to our new [cookbook](https://www.nanalogue.com)

### Changed
- Updates cargo packages
- Use `mimalloc` if target_env is `musl`; supposed to decrease program runtimes.

### Fixed
- `https` BAM retrieval
- `read_info` was not processing mod options when `--detailed` options were absent.
