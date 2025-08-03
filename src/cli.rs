//! # Cli
//!
//! This file provides some global options in the command line interface.
use clap::Args;

/// This struct is used to parse the input bam file and the filters that should be applied to the bam file.
/// This struct is parsed to create command line arguments and then passed to many functions.
/// We have copied and edited a similar struct from the fibertools-rs repository.
#[derive(Debug, Args)]
pub struct InputBam {
    /// Input BAM file. If no path is provided stdin is used.
    #[clap(default_value = "-")]
    pub bam_path: String,
    /// File with one column of read ids. If provided, only these
    /// reads are considered in any operation. Defaults to unused
    /// i.e. all reads in the BAM file are considered.
    #[clap(default_value = "")]
    pub read_list: String,
}
