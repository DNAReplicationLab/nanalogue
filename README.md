# `nanalogue`

Nanalogue = *N*ucleic Acid *Analogue*

Nanalogue is a tool to parse or analyse BAM/Mod BAM files with a single-molecule focus.

[![Cargo Build & Test](https://github.com/DNAReplicationLab/nanalogue/actions/workflows/ci.yml/badge.svg)](https://github.com/DNAReplicationLab/nanalogue/actions/workflows/ci.yml)
[![Code test coverage > 92\%](https://github.com/DNAReplicationLab/nanalogue/actions/workflows/cargo-llvm-cov.yml/badge.svg)](https://github.com/DNAReplicationLab/nanalogue/actions/workflows/cargo-llvm-cov.yml)
[![crates.io](https://img.shields.io/crates/v/nanalogue.svg)](https://crates.io/crates/nanalogue)
[![Documentation](https://docs.rs/nanalogue/badge.svg)](https://docs.rs/nanalogue)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Installation

Run the following command to install or update `nanalogue`.

```bash
cargo install nanalogue
```

`cargo` is the rust package manager. If you do not have `cargo`,
then follow these [instructions](https://doc.rust-lang.org/cargo/getting-started/installation.html)
to get it. On Linux and macOS systems, the install command is as simple as
`curl https://sh.rustup.rs -sSf | sh` 

# Commands

All the commands below have options you can specify on the command line.
Please run `--help` with a command to learn what the options are.

## `nanalogue read-table-hide-mods`

Prints basecalled length, alignment length, read id,
and optionally other information such as sequence per molecule.
This command does not expect modification data at all.

## `nanalogue read-table-show-mods`

Prints basecalled length, alignment length, read id, modification counts,
and optionally other information such as sequence per molecule.
If modification data is not available, any modification-related columns
have outputs like "NA".

## `nanalogue read-stats`
Calculates various summary statistics on all reads

## `nanalogue read-info`
Prints information about reads

## `nanalogue find-modified-reads`
Find names of modified reads through criteria specified by sub commands

## `nanalogue window-dens`
Output windowed densities of reads

## `nanalogue window-grad`
Output windowed gradients of reads

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
