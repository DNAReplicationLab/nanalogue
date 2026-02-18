# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

## [0.1.15] - 2026-01-18

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
