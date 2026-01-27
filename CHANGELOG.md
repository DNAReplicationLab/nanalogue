# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- Fixed `peek` to skip zero-length reads instead of returning an error

### Changed
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
