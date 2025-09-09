//! # Cli
//!
//! This file provides some global options in the command line interface.
use crate::{ModChar, ReadStates, RestrictModCalledStrand, ThresholdState, F32Bw0and1};
use clap::Args;
use serde::{Deserialize, Serialize};
use std::num::NonZeroU32;

/// This struct is used to parse the input bam file and the filters that should be applied to the bam file.
/// This struct is parsed to create command line arguments and then passed to many functions.
/// We have copied and edited a similar struct from the fibertools-rs repository.
#[derive(Debug, Args, Serialize, Deserialize)]
pub struct InputBam {
    /// Input BAM file. If no path is provided stdin is used.
    #[clap(default_value = "-")]
    pub bam_path: String,
    /// Exclude reads whose sequence length in the BAM file is
    /// below this value. Defaults to 0.
    #[clap(long, default_value_t = 0)]
    pub min_seq_len: u64,
    /// Exclude reads whose alignment length in the BAM file is
    /// below this value. Defaults to unused.
    #[clap(long)]
    pub min_align_len: Option<i64>,
    /// Only include this read id, defaults to unused i.e. all reads are used.
    /// NOTE: if there are multiple alignments corresponding
    /// to this read id, all of them are used.
    #[clap(long)]
    pub read_id: Option<String>,
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
    /// Only retain reads of this type. Allowed types are primary_forward,
    /// primary_reverse, secondary_forward, secondary_reverse, supplementary_forward,
    /// supplementary_reverse and unmapped. Specify more than one type if needed
    /// separated by commas, in which case reads of any type in list are retained.
    #[clap(long)]
    pub read_filter: Option<ReadStates>,
    /// Subsample BAM to retain only this fraction of total number of reads.
    /// NOTE: a new subsample is drawn every time as the seed is not fixed.
    /// If you want reproducibility, consider pipeing the output of `samtools view -s`.
    #[clap(short, long, default_value_t = F32Bw0and1::one())]
    pub sample_fraction: F32Bw0and1,
}

/// Implements a default class for InputBAM
impl Default for InputBam {
    fn default() -> Self {
        InputBam {
            bam_path: "".to_string(),
            min_seq_len: 0,
            min_align_len: None,
            read_id: None,
            threads: NonZeroU32::new(1).expect("no error"),
            exclude_zero_len: false,
            read_filter: None,
            sample_fraction: F32Bw0and1::one(),
        }
    }
}

/// This struct contains the options input to our
/// modification-data functions with restrictions on data received
#[derive(Debug, Args, Serialize, Deserialize)]
pub struct InputMods {
    /// modified tag
    #[clap(long)]
    pub tag: ModChar,
    /// modified strand, set this to bc or bc_comp, meaning
    /// on basecalled strand or its complement. Some technologies
    /// like PacBio or ONT duplex can call mod data on both a strand
    /// and its complementary DNA and store it in the record corresponding
    /// to the strand, so you can use this filter to select only for
    /// mod data on a strand or its complement. Please note that this
    /// filter is different from selecting for forward or reverse
    /// aligned reads using the BAM flags.
    #[clap(long)]
    pub mod_strand: Option<RestrictModCalledStrand>,
    /// Filter to reject mods before data operations.
    /// We allow all mods through for now and
    /// do not expose this to the user.
    #[clap(skip)]
    pub mod_prob_filter: ThresholdState,
    /// Filter this many bp at the start and
    /// end of a read before any operations.
    /// Please note that the units here are bp and
    /// not units of base being queried.
    #[clap(long, default_value_t = 0)]
    pub trim_read_ends: usize,
}

/// Implements a default for InputMods
impl Default for InputMods {
    fn default() -> Self {
        InputMods {
            tag: ModChar::new('N'),
            mod_strand: None,
            mod_prob_filter: ThresholdState::GtEq(0),
            trim_read_ends: 0,
        }
    }
}

/// This struct contains the options input to our
/// modification-data-windowing functions
#[derive(Debug, Args, Serialize, Deserialize)]
pub struct InputWindowing {
    /// size of window in units of base being queried i.e.
    /// if you are looking for cytosine modifications, then
    /// a window of a value 300 means create windows each with
    /// 300 cytosines irrespective of their modification status.
    #[clap(long)]
    pub win: NonZeroU32,
    /// step window by this size in units of base being queried.
    #[clap(long)]
    pub step: NonZeroU32,
}

/// Implements a default for InputWindowing
/// NOTE: we just choose 100 as the default value arbitrarily.
impl Default for InputWindowing {
    fn default() -> Self {
        InputWindowing {
            win: NonZeroU32::new(100).expect("no error"),
            step: NonZeroU32::new(100).expect("no error"),
        }
    }
}
