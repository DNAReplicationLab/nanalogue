# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
