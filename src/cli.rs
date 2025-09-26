//! # Cli
//!
//! This file provides some global options in the command line interface.
use crate::{
    Error, F32Bw0and1, GenomicRegion, ModChar, ReadStates, RestrictModCalledStrand, ThresholdState,
};
use bedrs::prelude::Bed3;
use clap::{Args, FromArgMatches};
use rust_htslib::bam;
use serde::{Deserialize, Serialize};
use std::num::NonZeroU32;

/// This struct is used to parse the input bam file and the filters that should be applied to the bam file.
/// This struct is parsed to create command line arguments and then passed to many functions.
/// We have copied and edited a similar struct from the fibertools-rs repository.
#[derive(Debug, Args, Clone, Serialize, Deserialize)]
pub struct InputBam {
    /// Input BAM file. Set this to - to read from stdin.
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
    /// Include "zero-length" sequences e.g. sequences with "*" in the sequence
    /// field. By default, these sequences are excluded to avoid processing errors.
    /// If this flag is set, these reads are included irrespective of any
    /// minimum sequence or align length criteria the user may have set.
    /// WARNINGS: (1) Some functions of the codebase may break or produce incorrect
    /// results if you use this flag. (2) due to a technical reason, we need a DNA sequence
    /// in the sequence field and cannot infer sequence length from other sources
    /// e.g. CIGAR strings.
    #[clap(long, default_value_t = false)]
    pub include_zero_len: bool,
    /// Only retain reads of this type. Allowed types are primary_forward,
    /// primary_reverse, secondary_forward, secondary_reverse, supplementary_forward,
    /// supplementary_reverse and unmapped. Specify more than one type if needed
    /// separated by commas, in which case reads of any type in list are retained.
    /// Defaults to retain reads of all types.
    #[clap(long)]
    pub read_filter: Option<ReadStates>,
    /// Subsample BAM to retain only this fraction of total number of reads,
    /// defaults to 1.0.
    /// NOTE: a new subsample is drawn every time as the seed is not fixed.
    /// If you want reproducibility, consider pipeing the output of `samtools view -s`.
    #[clap(short, long, default_value_t = F32Bw0and1::one())]
    pub sample_fraction: F32Bw0and1,
    /// Exclude reads whose MAPQ (Mapping quality of position) is below this value.
    /// Defaults to zero i.e. do not exclude any read.
    #[clap(long, default_value_t = 0)]
    pub mapq_filter: u8,
    /// Exclude sequences with MAPQ unavailable.
    /// In the BAM format, a value of 255 in this column means MAPQ is unavailable.
    /// These reads are allowed by default, set this flag to exclude.
    #[clap(long, default_value_t = false)]
    pub exclude_mapq_unavail: bool,
    /// Only keep reads passing through this region.
    #[clap(long)]
    pub region: Option<GenomicRegion>,
    /// Only keep read data from this region.
    /// This is an internal option not exposed to the user, we will set it
    /// based on the other options that the user sets.
    #[clap(skip)]
    pub region_bed3: Option<Bed3<i32, u64>>,
    /// Only keep reads if they pass through the specified region in full.
    /// Related to the input `--region`; has no effect if that is not set.
    #[clap(long, default_value_t = false)]
    pub full_region: bool,
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
            include_zero_len: false,
            read_filter: None,
            sample_fraction: F32Bw0and1::one(),
            mapq_filter: 0,
            exclude_mapq_unavail: false,
            region: None,
            region_bed3: None,
            full_region: false,
        }
    }
}

/// This struct contains an optional modification tag
#[derive(Debug, Default, Args, Clone, Copy, Serialize, Deserialize)]
pub struct OptionalTag {
    /// modified tag
    #[clap(long)]
    pub tag: Option<ModChar>,
}

/// This struct contains a required modification tag
#[derive(Debug, Default, Args, Clone, Copy, Serialize, Deserialize)]
pub struct RequiredTag {
    /// modified tag
    #[clap(long)]
    pub tag: ModChar,
}

/// Trait that returns a modification tag
pub trait TagState {
    /// Returns the modification tag of the tag state in an option
    fn tag(&self) -> Option<ModChar> {
        todo!();
    }
}

impl TagState for OptionalTag {
    fn tag(&self) -> Option<ModChar> {
        self.tag
    }
}

impl TagState for RequiredTag {
    fn tag(&self) -> Option<ModChar> {
        Some(self.tag)
    }
}

/// This struct contains the options input to our
/// modification-data functions with restrictions on data received
#[derive(Debug, Args, Clone, Serialize, Deserialize)]
pub struct InputMods<S: TagState + Args + FromArgMatches> {
    /// modified tag
    #[clap(flatten)]
    pub tag: S,
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
    /// Exclude bases whose base quality is below
    /// this threshold, defaults to 0 i.e. unused.
    /// NOTE: No offsets such as +33 are needed here.
    /// NOTE: Reads with missing base quality information
    /// are rejected if this is non-zero.
    #[clap(long, default_value_t = 0)]
    pub base_qual_filter: u8,
    /// Only keep modification data from this region
    #[clap(long)]
    pub mod_region: Option<GenomicRegion>,
    /// Only keep modification data from this region.
    /// We do not expose this to the user, but infer it from
    /// the other options set by the user.
    #[clap(skip)]
    pub region_bed3: Option<Bed3<i32, u64>>,
}

/// Implements defaults for InputMods
impl Default for InputMods<OptionalTag> {
    fn default() -> Self {
        InputMods::<OptionalTag> {
            tag: OptionalTag { tag: None },
            mod_strand: None,
            mod_prob_filter: ThresholdState::GtEq(0),
            trim_read_ends: 0,
            base_qual_filter: 0,
            mod_region: None,
            region_bed3: None,
        }
    }
}

/// Retrieves options for modification input
pub trait InputModOptions {
    /// retrieves tag
    fn tag(&self) -> Option<ModChar> {
        todo!()
    }
    /// retrieves option to set basecalled strand or opposite in mod retrieval
    fn mod_strand(&self) -> Option<RestrictModCalledStrand> {
        todo!()
    }
    /// returns probability filter
    fn mod_prob_filter(&self) -> ThresholdState {
        todo!()
    }
    /// returns read end trimming
    fn trim_read_ends(&self) -> usize {
        todo!()
    }
    /// returns threshold for filtering base PHRED quality
    fn base_qual_filter(&self) -> u8 {
        todo!()
    }
}

/// Retrieves options for region
pub trait InputRegionOptions {
    /// returns region requested
    fn region_filter(&self) -> &Option<Bed3<i32, u64>> {
        todo!()
    }
    /// returns region requested but region in genomic string format
    fn region_filter_genomic_string(&self) -> Option<GenomicRegion> {
        todo!()
    }
    /// sets region requested
    fn set_region_filter(&mut self, _value: Option<Bed3<i32, u64>>) {
        todo!()
    }
    /// returns true if full overlap with region is requested as opposed to
    /// only partial overlap. defaults to false.
    fn is_full_overlap(&self) -> bool {
        false
    }
    /// converts region from genomic string representation to bed3 representation
    fn convert_region_to_bed3(&mut self, header: bam::HeaderView) -> Result<bool, Error> {
        match self.region_filter_genomic_string() {
            None => self.set_region_filter(None),
            Some(v) => self.set_region_filter(Some(v.try_to_bed3(header)?)),
        }
        Ok(true)
    }
}

impl<S: TagState + Args + FromArgMatches> InputModOptions for InputMods<S> {
    fn tag(&self) -> Option<ModChar> {
        self.tag.tag()
    }
    fn mod_strand(&self) -> Option<RestrictModCalledStrand> {
        self.mod_strand
    }
    fn mod_prob_filter(&self) -> ThresholdState {
        self.mod_prob_filter
    }
    fn trim_read_ends(&self) -> usize {
        self.trim_read_ends
    }
    fn base_qual_filter(&self) -> u8 {
        self.base_qual_filter
    }
}

impl<S: TagState + Args + FromArgMatches> InputRegionOptions for InputMods<S> {
    fn region_filter_genomic_string(&self) -> Option<GenomicRegion> {
        self.mod_region.clone()
    }
    fn region_filter(&self) -> &Option<Bed3<i32, u64>> {
        &self.region_bed3
    }
    fn set_region_filter(&mut self, value: Option<Bed3<i32, u64>>) {
        self.region_bed3 = value;
    }
}

/// can return tag without encasing in an option in the RequiredTag variant
impl InputMods<RequiredTag> {
    /// retrieves tag
    pub fn tag(&self) -> ModChar {
        self.tag.tag
    }
}

impl InputRegionOptions for InputBam {
    fn region_filter_genomic_string(&self) -> Option<GenomicRegion> {
        self.region.clone()
    }
    fn region_filter(&self) -> &Option<Bed3<i32, u64>> {
        &self.region_bed3
    }
    fn set_region_filter(&mut self, value: Option<Bed3<i32, u64>>) {
        self.region_bed3 = value;
    }
    fn is_full_overlap(&self) -> bool {
        self.full_region
    }
}

/// This struct contains the options input to our
/// modification-data-windowing functions
#[derive(Debug, Args, Clone, Copy, Serialize, Deserialize)]
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
