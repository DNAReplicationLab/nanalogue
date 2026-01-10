# `nanalogue`

Nanalogue = *N*ucleic Acid *Analogue*

Nanalogue is a tool to parse or analyse BAM/Mod BAM files with a single-molecule focus.

[![Cargo Build & Test](https://github.com/DNAReplicationLab/nanalogue/actions/workflows/ci.yml/badge.svg)](https://github.com/DNAReplicationLab/nanalogue/actions/workflows/ci.yml)
[![Code test coverage > 92\%](https://github.com/DNAReplicationLab/nanalogue/actions/workflows/cargo-llvm-cov.yml/badge.svg)](https://github.com/DNAReplicationLab/nanalogue/actions/workflows/cargo-llvm-cov.yml)
[![crates.io](https://img.shields.io/crates/v/nanalogue.svg)](https://crates.io/crates/nanalogue)
[![Documentation](https://docs.rs/nanalogue/badge.svg)](https://docs.rs/nanalogue)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A common pain point in genomics analyses is that BAM files are information-dense
which makes it difficult to gain insight from them. Nanalogue hopes to make it easy
to extract and process this information, with a particular focus on single-molecule
aspects and DNA/RNA modifications. Despite this focus, some of nanalogue's commands are
quite general and can be applied to almost any BAM file.

We can process any type of DNA/RNA modifications occurring in any pattern (single/multiple mods,
spatially-isolated/non-isolated etc.). All we require is that the data is stored
in a BAM file in the mod BAM format (i.e. using MM/ML tags as laid down in the
[specifications](https://samtools.github.io/hts-specs/SAMtags.pdf)).

## Table of Contents

- [Usage and documentation](#usage-and-documentation)
  - [Simulate BAM files](#simulate-bam-files)
- [Installation](#installation)
  - [Using Cargo](#using-cargo)
  - [Using Docker](#using-docker)
  - [Pre-built Binaries](#pre-built-binaries)
    - [GitHub Releases](#github-releases)
    - [GitHub Actions Artifacts](#github-actions-artifacts)
- [Commands](#commands)
  - [`nanalogue read-info`](#nanalogue-read-info)
  - [`nanalogue read-table-hide-mods`](#nanalogue-read-table-hide-mods)
  - [`nanalogue read-table-show-mods`](#nanalogue-read-table-show-mods)
  - [`nanalogue read-stats`](#nanalogue-read-stats)
  - [`nanalogue find-modified-reads`](#nanalogue-find-modified-reads)
  - [`nanalogue window-dens`](#nanalogue-window-dens)
  - [`nanalogue window-grad`](#nanalogue-window-grad)
  - [`nanalogue peek`](#nanalogue-peek)
- [Python wrapper](#python-wrapper)
- [Contributing](#contributing)
- [Security](#security)
- [Changelog](#changelog)
- [Acknowledgments](#acknowledgments)

# Usage and documentation

There are three ways to use nanalogue with a pre-existing BAM file:
- as a tool on the command line
- as a rust library
- as a python library

This note is largely about command line usage.
Please refer to the documentation in [docs.rs](https://docs.rs/nanalogue),
or look at this [section](#python-wrapper) for usage as a rust
and a python library respectively.

In addition to these resources, we are developing a
companion cookbook [here](https://www.nanalogue.com).

## Simulate BAM files

For developers: if you are looking to make a custom BAM file containing synthetic, simulated
DNA/RNA modification data to develop/test your tool, you may be interested in `nanalogue_sim_bam`.
This is an executable that ships with nanalogue that can create a BAM file according to your
specifications. Please run `nanalogue_sim_bam --help`. If you are a rust developer looking
to use this functionality in your library, please look at the documentation of the module
`nanalogue_core::simulate_mod_bam` in the docs.rs link [above](#usage).

# Installation

## Using Cargo

Run the following command to install or update `nanalogue` for usage on the command line.

```bash
cargo install nanalogue
```

`cargo` is the rust package manager. If you do not have `cargo`,
then follow these [instructions](https://doc.rust-lang.org/cargo/getting-started/installation.html)
to get it. On Linux and macOS systems, the install command is as simple as
`curl https://sh.rustup.rs -sSf | sh`

## Using Docker

You can also use `nanalogue` via Docker:

```bash
docker pull dockerofsat/nanalogue:latest
```

## Pre-built Binaries

Pre-built binaries for macOS and Linux are available from two sources:

### GitHub Releases

Official release binaries can be downloaded from the [Releases page](https://github.com/DNAReplicationLab/nanalogue/releases)
on the Github repository. Each release includes binaries for multiple platforms.

### GitHub Actions Artifacts

Binaries built from the latest code are available as artifacts from the 
[Build Release Binaries workflow](https://github.com/DNAReplicationLab/nanalogue/actions/workflows/build-binaries.yml)
from the github repository. To download:

1. Navigate to the workflow runs
2. Click on a successful workflow run
3. Scroll to the "Artifacts" section at the bottom
4. Download the binary artifact for your platform:
   - `binaries-macos` - macOS binaries
   - `binaries-musllinux_1_2` - Alpine/musl-based Linux (static binaries)
   - `binaries-manylinux_2_28` - Modern Linux distributions (glibc 2.28+)
   - `binaries-manylinux_2_34` - Newer Linux distributions (glibc 2.34+)
   - `binaries-manylinux2014` - Older Linux distributions (glibc 2.17+, maximum compatibility)

# Commands

All the commands below have options you can specify on the command line.
Please run `--help` with a command to learn what these are. Among other operations,
the options allow you to subsample the BAM file (`-s`),
restrict read and/or modification data to a specific genomic region (`--region` or `--mod-region`),
restrict by one or several read ids (`--read-id` or `--read-id-list`),
a specific mapping type (`--read-filter`), filter modification data suitably
(`--mod-prob-filter`) etc.

## `nanalogue read-info`
Prints information about reads in JSON. A sample output snippet follows.

```json
[
{
        "read_id": "a4f36092-b4d5-47a9-813e-c22c3b477a0c",
        "sequence_length": 48,
        "contig": "dummyIII",
        "reference_start": 23,
        "reference_end": 71,
        "alignment_length": 48,
        "alignment_type": "primary_forward",
        "mod_count": "T+T:3;(probabilities >= 0.5020, PHRED base qual >= 0)"
}
]
```

With options like `--detailed` and `--detailed-pretty`, the modification information in the BAM file
is converted to a more-usable JSON format. A sample output snippet follows.

```json
[
{
  "alignment_type": "primary_forward",
  "alignment": {
    "start": 23,
    "end": 71,
    "contig": "dummyIII",
    "contig_id": 2
  },
  "mod_table": [
    {
      "base": "T",
      "is_strand_plus": true,
      "mod_code": "T",
      "data": [
        [
          3,
          26,
          221
        ],
        [
          8,
          31,
          242
        ],
        [
          27,
          50,
          3
        ],
        [
          39,
          62,
          47
        ],
        [
          47,
          70,
          239
        ]
      ]
    }
  ],
  "read_id": "a4f36092-b4d5-47a9-813e-c22c3b477a0c",
  "seq_len": 48
}
]
```

## `nanalogue read-table-hide-mods`

Prints basecalled length, alignment length, read id,
and optionally other information such as sequence per molecule.
This command does not expect modification data at all.
A sample output snippet follows.

```text
read_id align_length    sequence_length_template        alignment_type
a4f36092-b4d5-47a9-813e-c22c3b477a0c    48, 0   46      primary_forward, unmapped
5d10eb9a-aae1-4db8-8ec6-7ebb34d32575    8       6       primary_forward
fffffff1-10d2-49cb-8ca3-e8d48979001b    33      38      primary_reverse
```

## `nanalogue read-table-show-mods`

Prints basecalled length, alignment length, read id, modification counts,
and optionally other information such as sequence per molecule.
If modification data is not available, any modification-related columns
have outputs like "NA". A sample output snippet follows.

```text
# mod-unmod threshold is 0.5
read_id align_length    sequence_length_template        alignment_type  mod_count
a4f36092-b4d5-47a9-813e-c22c3b477a0c    48, 0   48      primary_forward, unmapped       T:3, T:3;232486:0
5d10eb9a-aae1-4db8-8ec6-7ebb34d32575    8       8       primary_forward T:0
fffffff1-10d2-49cb-8ca3-e8d48979001b    33      33      primary_reverse T:1
```

## `nanalogue read-stats`
Calculates various summary statistics on all reads. A sample output follows.

```text
key     value
n_primary_alignments    3
n_secondary_alignments  0
n_supplementary_alignments      0
n_unmapped_reads        1
n_reversed_reads        1
align_len_mean  29
align_len_max   48
align_len_min   8
align_len_median        8
align_len_n50   48
seq_len_mean    34
seq_len_max     48
seq_len_min     8
seq_len_median  33
seq_len_n50     48
```

## `nanalogue find-modified-reads`
Find names of modified reads through criteria specified by sub commands
e.g.  at least one window with a modification density above
some value (`any-dens-above`). Please run
`nanalogue find-modified-reads --help` to learn more.
Output is a list of read ids that satisfy the specified criterion e.g.

```text
a4f36092-b4d5-47a9-813e-c22c3b477a0f
5d10eb9a-aae1-4db8-8ec6-7ebb34d32576
fffffff1-10d2-49cb-8ca3-e8d48979001a
```

## `nanalogue window-dens`
Output windowed densities of reads. Sample output follows.

```text
#contig ref_win_start   ref_win_end     read_id win_val strand  base    mod_strand      mod_type        win_start       win_end basecall_qual
dummyI  9       13      read1   0       +       T       +       T       0       3       2
dummyI  12      14      read1   0       +       T       +       T       2       4       3
dummyI  13      17      read1   0       +       T       +       T       3       7       32
```

## `nanalogue window-grad`
Output windowed gradients of reads.
Sample output is similar to the one above but with gradients reported instead of
mean modification densities in the `win_val` column.

## `nanalogue peek`
Read the BAM file header and the first few records, and output contig information
and modification tag information. Sample output follows.

```text
contigs_and_lengths:
dummyI  22
dummyII 48
dummyIII        76

modifications:
G-7200
T+T
```

If no modifications are found in the first few records, you will see

```text
contigs_and_lengths:
dummyI  22
dummyII 48
dummyIII        76

modifications:
None
```

# Python wrapper

A python wrapper of some of these commands `pynanalogue` is in development
[here](https://github.com/DNAReplicationLab/pynanalogue).

# Contributing

Contributions are welcome! Please see [CONTRIBUTIONS.md](CONTRIBUTIONS.md)
for guidelines on how to contribute to this project.

# Security

For security concerns and vulnerability reporting, please see [SECURITY.md](SECURITY.md).

# Changelog

Changelog of the project is at [CHANGELOG.md](CHANGELOG.md).

# Acknowledgments

This software was developed at the Earlham Institute in the UK.
This work was supported by the Biotechnology and Biological Sciences
Research Council (BBSRC), part of UK Research and Innovation,
through the Core Capability Grant BB/CCG2220/1 at the Earlham Institute
and the Earlham Institute Strategic Programme Grant Cellular Genomics
BBX011070/1 and its constituent work packages BBS/E/ER/230001B
(CellGen WP2 Consequences of somatic genome variation on traits).
The work was also supported by the following response-mode project grants:
BB/W006014/1 (Single molecule detection of DNA replication errors) and
BB/Y00549X/1 (Single molecule analysis of Human DNA replication).
This research was supported in part by NBI Research Computing
through use of the High-Performance Computing system and Isilon storage.
