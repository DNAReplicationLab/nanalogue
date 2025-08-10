//! # Cli
//!
//! This file provides some global options in the command line interface.
use clap::Args;
use std::num::NonZeroU32;

/// This struct is used to parse the input bam file and the filters that should be applied to the bam file.
/// This struct is parsed to create command line arguments and then passed to many functions.
/// We have copied and edited a similar struct from the fibertools-rs repository.
#[derive(Debug, Args)]
pub struct InputBam {
    /// Input BAM file. If no path is provided stdin is used.
    #[clap(default_value = "-")]
    pub bam_path: String,
    /// Exclude reads whose sequence length in the BAM file is
    /// below this value.
    /// NOTE: This is applied using the BAM file and not on
    /// any other files possibly used by subcommands such as the
    /// sequencing summary file.
    #[clap(long, default_value_t = 0)]
    pub min_seq_len: u64,
    /// Number of threads used during some aspects of program execution
    #[clap(long, default_value_t = NonZeroU32::new(2).expect("no error"))]
    pub threads: NonZeroU32,
    /// Exclude "zero-length" sequences e.g. sequences with "*" in the sequence
    /// field. Please use this flag if you encounter this error in our program.
    /// NOTE: due to a technical reason, we need a DNA sequence
    /// in the sequence field and cannot infer sequence length from other sources
    /// e.g. CIGAR strings.
    #[clap(long, default_value_t = false)]
    pub exclude_zero_len: bool,
}
