//! # Cli
//!
//! This file provides some global options in the command line interface.
use crate::{
    Error, F32Bw0and1, GenomicRegion, ModChar, PathOrURLOrStdin, ReadStates,
    RestrictModCalledStrand, ThresholdState,
};
use bedrs::prelude::Bed3;
use clap::{Args, FromArgMatches};
use derive_builder::Builder;
use rust_htslib::bam;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::num::{NonZeroU32, NonZeroUsize};
use std::str::FromStr;

/// Options to parse the input bam file and the filters that should be applied to the bam file.
///
/// This struct is parsed to create command line arguments and then passed to many functions.
/// We have copied and edited a similar struct from the fibertools-rs repository.
/// You can build this through [`InputBamBuilder`].
/// In CLI mode, `clap` populates this struct.
/// This and the [`InputMods`] struct are used to set almost all input options
/// to many of our functions that process BAM/modBAM files.
///
/// # Examples
///
/// Sample way to build the struct. Some of the parameters are optional
/// and can be left unset which would give them default values.
/// We do not check if the specified bam path or URL exists as there are
/// use cases where files are generated before the `InputBam` object is used.
///
/// ```
/// use nanalogue_core::{Error, InputBamBuilder, PathOrURLOrStdin};
///
/// let bam = InputBamBuilder::default()
///     .bam_path(PathOrURLOrStdin::Path("/some/path/to/bam.bam".into()))
///     .min_seq_len(30000u64)
///     .min_align_len(20000i64)
///     .read_id("some-id")
///     .read_filter("primary_forward,secondary_forward".into())
///     .sample_fraction(1.0)
///     .mapq_filter(20)
///     .exclude_mapq_unavail(true)
///     .region("chr4:1000-2000".into())
///     .full_region(true)
///     .build()?;
/// # Ok::<(), Error>(())
/// ```
///
/// This struct and the [`InputMods`] struct allow us to set input options
/// for BAM/modBAM calculations. An example is shown below where the [`crate::read_info::run`]
/// command is called to process data from a BAM file with some input options.
///
/// ```
/// use nanalogue_core::{BamRcRecords, BamPreFilt, Error, InputBamBuilder, InputModsBuilder,
///     OptionalTag, PathOrURLOrStdin, ThresholdState, nanalogue_bam_reader, read_info};
///
/// let mut bam = InputBamBuilder::default()
///     .bam_path(PathOrURLOrStdin::Path("./examples/example_1.bam".into()))
///     .region("dummyI".into())
///     .build()?;
/// let mut mods = InputModsBuilder::<OptionalTag>::default()
///     .mod_prob_filter(ThresholdState::GtEq(0))
///     .build()?;
///
/// let mut buffer = Vec::new();
/// let mut reader = nanalogue_bam_reader(&bam.bam_path.to_string())?;
/// let bam_rc_records = BamRcRecords::new(&mut reader, &mut bam, &mut mods)?;
/// read_info::run(
///     &mut buffer,
///     bam_rc_records.rc_records
///         .filter(|r| r.as_ref().map_or(true, |v| v.pre_filt(&bam))),
///     &mods,
///     None,
/// )?;
/// # Ok::<(), Error>(())
/// ```
///
/// ## Examples resulting in errors
///
/// Full region without actually setting a region
///
/// ```should_panic
/// use nanalogue_core::{Error, InputBamBuilder, PathOrURLOrStdin};
///
/// let bam = InputBamBuilder::default()
///     .bam_path(PathOrURLOrStdin::Path("/some/path/to/bam.bam".into()))
///     .read_id("some-id")
///     .full_region(true)
///     .build()?;
/// # Ok::<(), Error>(())
/// ```
///
/// Setting both `region` and `region_bed3`. `region` can be converted to
/// `region_bed3` using [`GenomicRegion::try_to_bed3`] and a BAM header.
///
/// ```should_panic
/// use bedrs::prelude::Bed3;
/// use nanalogue_core::{Error, InputBamBuilder, PathOrURLOrStdin};
///
/// let bam = InputBamBuilder::default()
///     .bam_path(PathOrURLOrStdin::Path("/some/path/to/bam.bam".into()))
///     .read_id("some-id")
///     .region("chr4:1000-2000".into())
///     .region_bed3(Bed3::<i32,u64>::new(3, 1000, 2000))
///     .build()?;
/// # Ok::<(), Error>(())
/// ```
///
/// Setting more than one of `read_id`, `read_id_list` and `read_id_set`.
/// * `read_id` means filter to retain only this read.
/// * `read_id_list` is a path to a file with a list of read ids.
/// * `read_id_set` is a set of read ids supplied directly.
///
/// ```
/// use bedrs::prelude::Bed3;
/// use nanalogue_core::{Error, InputBamBuilder, PathOrURLOrStdin};
/// use std::collections::HashSet;
///
/// let _ = InputBamBuilder::default()
///     .bam_path(PathOrURLOrStdin::Path("/some/path/to/bam.bam".into()))
///     .read_id("some-id")
///     .read_id_list("/some/file.txt")
///     .build().unwrap_err();
///
/// let mut read_id_set = HashSet::<String>::new();
/// read_id_set.insert("some-read-a".to_owned());
/// read_id_set.insert("some-read-b".to_owned());
///
/// let _ = InputBamBuilder::default()
///     .bam_path(PathOrURLOrStdin::Path("/some/path/to/bam.bam".into()))
///     .read_id_list("/some/file.txt")
///     .read_id_set(read_id_set.clone())
///     .build().unwrap_err();
///
/// let _ = InputBamBuilder::default()
///     .bam_path(PathOrURLOrStdin::Path("/some/path/to/bam.bam".into()))
///     .read_id("some-id")
///     .read_id_set(read_id_set)
///     .build().unwrap_err();
/// ```
#[derive(Builder, Debug, Args, Clone, Serialize, Deserialize)]
#[serde(default)]
#[builder(default, build_fn(error = "Error", validate = "Self::validate"))]
#[non_exhaustive]
pub struct InputBam {
    /// Input BAM file. Set to a local file path, or set to - to read from stdin,
    /// or set to a URL to read from a remote file. If using stdin and piping in
    /// from `samtools view`, always include the header with the `-h` option.
    pub bam_path: PathOrURLOrStdin,
    /// Exclude reads whose sequence length in the BAM file is
    /// below this value. Defaults to 0.
    #[clap(long, default_value_t)]
    #[builder(setter(into))]
    pub min_seq_len: u64,
    /// Exclude reads whose alignment length in the BAM file is
    /// below this value. Defaults to unused.
    #[clap(long)]
    #[builder(setter(into, strip_option))]
    pub min_align_len: Option<i64>,
    /// Only include this read id, defaults to unused i.e. all reads are used.
    /// NOTE: if there are multiple alignments corresponding
    /// to this read id, all of them are used.
    #[clap(long, conflicts_with = "read_id_list")]
    #[builder(setter(into, strip_option))]
    pub read_id: Option<String>,
    /// Path to file containing list of read IDs (one per line).
    /// Lines starting with '#' are treated as comments and ignored.
    /// Cannot be used together with --read-id.
    #[clap(long, conflicts_with = "read_id")]
    #[builder(setter(into, strip_option))]
    pub read_id_list: Option<String>,
    /// Internal `HashSet` of read IDs loaded from `read_id_list` file.
    /// This is populated automatically and not exposed to users.
    #[clap(skip)]
    #[builder(setter(into, strip_option))]
    pub read_id_set: Option<HashSet<String>>,
    /// Number of threads used during some aspects of program execution
    #[clap(long, default_value_t = NonZeroU32::new(2).expect("no error"))]
    #[builder(setter(into))]
    pub threads: NonZeroU32,
    /// Include "zero-length" sequences e.g. sequences with "*" in the sequence
    /// field. By default, these sequences are excluded to avoid processing errors.
    /// If this flag is set, these reads are included irrespective of any
    /// minimum sequence or align length criteria the user may have set.
    /// WARNINGS: (1) Some functions of the codebase may break or produce incorrect
    /// results if you use this flag. (2) due to a technical reason, we need a DNA sequence
    /// in the sequence field and cannot infer sequence length from other sources
    /// e.g. CIGAR strings.
    #[clap(long, default_value_t)]
    pub include_zero_len: bool,
    /// Only retain reads of this type. Allowed types are `primary_forward`,
    /// `primary_reverse`, `secondary_forward`, `secondary_reverse`, `supplementary_forward`,
    /// `supplementary_reverse` and unmapped. Specify more than one type if needed
    /// separated by commas, in which case reads of any type in list are retained.
    /// Defaults to retain reads of all types.
    #[clap(long)]
    #[builder(field(
        ty = "String",
        build = "(!self.read_filter.is_empty()).then(|| ReadStates::from_str(&self.read_filter)).transpose()?"
    ))]
    pub read_filter: Option<ReadStates>,
    /// Subsample BAM to retain only this fraction of total number of reads,
    /// defaults to 1.0. The sampling algorithm considers every read according
    /// to the specified probability, so due to this, you may not always get
    /// the same number of reads e.g. if you set `-s 0.05` in a file with 1000 reads,
    /// you will get 50 +- sqrt(50) reads.
    /// NOTE: a new subsample is drawn every time as the seed is not fixed.
    /// If you want reproducibility, consider piping the output of `samtools view -s`
    /// to our program.
    #[clap(short, long, default_value_t = F32Bw0and1::one())]
    #[builder(field(ty = "f32", build = "self.sample_fraction.try_into()?"))]
    pub sample_fraction: F32Bw0and1,
    /// Exclude reads whose MAPQ (Mapping quality of position) is below this value.
    /// Defaults to zero i.e. do not exclude any read.
    #[clap(long, default_value_t)]
    #[builder(setter(into))]
    pub mapq_filter: u8,
    /// Exclude sequences with MAPQ unavailable.
    /// In the BAM format, a value of 255 in this column means MAPQ is unavailable.
    /// These reads are allowed by default, set this flag to exclude.
    #[clap(long, default_value_t)]
    pub exclude_mapq_unavail: bool,
    /// Only keep reads passing through this region. If a BAM index is available
    /// with a name same as the BAM file but with the .bai suffix, the operation
    /// of selecting such reads will be faster. If you are using standard input
    /// as your input e.g. you are piping in the output from samtools, then
    /// you cannot use an index as a BAM filename is not available.
    #[clap(long)]
    #[builder(field(
        ty = "String",
        build = "(!self.region.is_empty()).then(|| GenomicRegion::from_str(&self.region)).transpose()?"
    ))]
    pub region: Option<GenomicRegion>,
    /// Only keep read data from this region.
    /// This is an internal option not exposed to the user, we will set it
    /// based on the other options that the user sets.
    #[clap(skip)]
    #[builder(setter(into, strip_option))]
    pub region_bed3: Option<Bed3<i32, u64>>,
    /// Only keep reads if they pass through the specified region in full.
    /// Related to the input `--region`; has no effect if that is not set.
    #[clap(long, default_value_t, requires = "region")]
    pub full_region: bool,
}

impl InputBamBuilder {
    /// Validate [`InputBam`] creation through the Builder method.
    ///
    /// [`InputBam`] can be created in many ways:
    /// - directly populating the fields as they are public.
    /// - through the CLI-machinery
    /// - through the builder.
    ///
    /// While we cannot check a user's usage through the first method,
    /// we can check constructions through methods 2 and 3 above.
    /// The CLI construction checking is handled by `clap`, and the builder
    /// checking is done through this function.
    /// The construction checks are slightly different as what may be feasible
    /// through a CLI and through code are different.
    /// So we have this validate method to check construction through
    /// method number 3 above.
    #[expect(
        clippy::arithmetic_side_effects,
        reason = "1 + 1 + 1 is not gonna overflow"
    )]
    fn validate(&self) -> Result<(), Error> {
        match (
            !self.region.is_empty(),
            self.region_bed3.is_some(),
            self.full_region,
        ) {
            (false, false, Some(true)) => {
                return Err(Error::BuilderValidation("InputBamBuilder: cannot set `full_region` without setting `region` or `region_bed3`".to_owned()));
            }
            (true, true, _) => {
                return Err(Error::BuilderValidation(
                    "InputBamBuilder: cannot set both `region` and `region_bed3`".to_owned(),
                ));
            }
            _ => {}
        }

        if u8::from(self.read_id.is_some())
            + u8::from(self.read_id_list.is_some())
            + u8::from(self.read_id_set.is_some())
            > 1
        {
            Err(Error::BuilderValidation(
                "InputBamBuilder: cannot set >1 of `read_id`, `read_id_list` and `read_id_set`"
                    .to_owned(),
            ))
        } else {
            Ok(())
        }
    }
}

/// Implements a default class for `InputBAM`
impl Default for InputBam {
    fn default() -> Self {
        InputBam {
            bam_path: PathOrURLOrStdin::Stdin,
            min_seq_len: 0,
            min_align_len: None,
            read_id: None,
            read_id_list: None,
            read_id_set: None,
            threads: NonZeroU32::new(2).expect("no error"),
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
#[non_exhaustive]
pub struct OptionalTag {
    /// modified tag
    #[clap(long)]
    pub tag: Option<ModChar>,
}

impl FromStr for OptionalTag {
    type Err = Error;

    /// process the modification type from a string.
    ///
    /// ```
    /// use nanalogue_core::{Error, ModChar, OptionalTag};
    /// use std::str::FromStr;
    ///
    /// let tag = OptionalTag::from_str("m")?;
    /// assert_eq!(tag.tag.unwrap().val(), 'm');
    /// # Ok::<(), Error>(())
    /// ```
    fn from_str(mod_type: &str) -> Result<Self, Self::Err> {
        Ok(OptionalTag {
            tag: Some(ModChar::from_str(mod_type)?),
        })
    }
}

/// This struct contains a required modification tag
#[derive(Debug, Default, Args, Clone, Copy, Serialize, Deserialize)]
#[non_exhaustive]
pub struct RequiredTag {
    /// modified tag
    #[clap(long)]
    pub tag: ModChar,
}

impl FromStr for RequiredTag {
    type Err = Error;

    /// process the modification type from a string.
    ///
    /// ```
    /// use nanalogue_core::{Error, ModChar, RequiredTag};
    /// use std::str::FromStr;
    ///
    /// let tag = RequiredTag::from_str("m")?;
    /// assert_eq!(tag.tag.val(), 'm');
    /// # Ok::<(), Error>(())
    /// ```
    fn from_str(mod_type: &str) -> Result<Self, Self::Err> {
        Ok(RequiredTag {
            tag: ModChar::from_str(mod_type)?,
        })
    }
}

/// Trait that returns a modification tag
pub trait TagState {
    /// Returns the modification tag of the tag state in an option
    fn tag(&self) -> Option<ModChar> {
        unimplemented!();
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
/// modification-data functions with restrictions on data received.
///
/// If you want to build this, use [`InputModsBuilder`].
/// In CLI mode, `clap` populates this struct.
/// This and the [`InputBam`] struct are used to set almost all input options
/// to many of our functions that process BAM/modBAM files.
///
/// # Examples
///
/// Following example uses many fields and is for a [`RequiredTag`] variant.
/// You can omit some of them depending on your use case. `mod_region` must
/// can be converted to `region_bed3` using [`GenomicRegion::try_to_bed3`]
/// with a suitable BAM header before this struct can be used with most of
/// our functions. We do not force this in the builder route as we may not
/// want to do this for some reason.
///
/// ```
/// use nanalogue_core::{InputModsBuilder, Error, RequiredTag, ThresholdState};
/// use std::str::FromStr;
///
/// let input_options = InputModsBuilder::default()
///     .tag(RequiredTag::from_str("m")?)
///     .mod_strand("bc".into())
///     .mod_prob_filter(ThresholdState::GtEq(0))
///     .trim_read_ends_mod(10)
///     .base_qual_filter_mod(10)
///     .mod_region("chr3:4000-5000".into())
///     .build()?;
///
/// # Ok::<(), Error>(())
/// ```
///
/// Please look at the documentation of [`InputBam`] to see a fully worked
/// example to set up inputs to perform a calculation.
///
/// A [`RequiredTag`] variant that sets `region_bed3` instead of `mod_region`.
/// If `mod_region` is set, then a downstream function must parse a BAM file,
/// and then convert `mod_region` into `region_bed3` using the header.
/// If the user is confident that they know the index number of a contig,
/// then they can directly set `region_bed3`. In the example below, the user
/// knows `chr3` corresponds to `4`, so they set the region directly.
///
/// ```
/// use bedrs::prelude::Bed3;
/// use nanalogue_core::{InputModsBuilder, Error, RequiredTag, ThresholdState};
/// use std::str::FromStr;
///
/// let input_options = InputModsBuilder::<RequiredTag>::default()
///     .tag(RequiredTag::from_str("m")?)
///     .mod_strand("bc".into())
///     .mod_prob_filter(ThresholdState::GtEq(0))
///     .trim_read_ends_mod(10)
///     .base_qual_filter_mod(10)
///     .region_bed3(Bed3::<i32, u64>::new(4, 4000, 5000))
///     .build()?;
///
/// # Ok::<(), Error>(())
/// ```
///
/// Setting both `region_bed3` and `mod_region` should result in an error,
/// even if they both refer to the same region. This is because setting both
/// can result in an undefined state, so we do not allow this.
///
/// ```should_panic
/// use bedrs::prelude::Bed3;
/// use nanalogue_core::{InputModsBuilder, Error, RequiredTag, ThresholdState};
/// use std::str::FromStr;
///
/// let input_options = InputModsBuilder::<RequiredTag>::default()
///     .mod_prob_filter(ThresholdState::GtEq(0))
///     .region_bed3(Bed3::<i32, u64>::new(4, 4000, 5000))
///     .mod_region("chr3:4000-5000".into())
///     .build()?;
///
/// # Ok::<(), Error>(())
/// ```
///
/// Following example is for an [`OptionalTag`] variant with the `tag` omitted,
/// and many other fields omitted as well, which take on a default value.
///
/// ```
/// use nanalogue_core::{InputModsBuilder, Error, OptionalTag, ThresholdState};
/// use std::str::FromStr;
///
/// let input_options = InputModsBuilder::<OptionalTag>::default()
///     .mod_prob_filter(ThresholdState::GtEq(0))
///     .build()?;
///
/// # Ok::<(), Error>(())
/// ```
#[derive(Builder, Debug, Args, Clone, Serialize, Deserialize)]
#[builder(default, build_fn(error = "Error", validate = "Self::validate"))]
#[serde(default)]
#[non_exhaustive]
pub struct InputMods<S: TagState + Args + FromArgMatches + Default> {
    /// modified tag
    #[clap(flatten)]
    pub tag: S,
    /// modified strand, set this to `bc` or `bc_comp`, meaning
    /// on basecalled strand or its complement. Some technologies
    /// like `PacBio` or `ONT` duplex can call mod data on both a strand
    /// and its complementary DNA and store it in the record corresponding
    /// to the strand, so you can use this filter to select only for
    /// mod data on a strand or its complement. Please note that this
    /// filter is different from selecting for forward or reverse
    /// aligned reads using the BAM flags.
    #[clap(long)]
    #[builder(field(
        ty = "String",
        build = "(!self.mod_strand.is_empty()).then(|| RestrictModCalledStrand::from_str(&self.mod_strand)).transpose()?"
    ))]
    pub mod_strand: Option<RestrictModCalledStrand>,
    /// Filter to reject mods before analysis. Specify as low,high where
    /// both are fractions to reject modifications where the probabilities (p)
    /// are in this range e.g. "0.4,0.6" rejects 0.4 <= p <= 0.6.
    /// You can use this to reject 'weak' modification calls before analysis
    /// i.e. those with probabilities close to 0.5.
    /// NOTE: (1) Whether this filtration is applied or not, mods < 0.5
    /// are considered unmodified and >= 0.5 are considered modified by our program.
    /// (2) mod probabilities are stored as a number from 0-255 in the modBAM format,
    /// so we internally convert 0.0-1.0 to 0-255. Default: reject nothing.
    #[clap(long, value_parser=ThresholdState::from_str_ordpair_fraction, default_value = "")]
    pub mod_prob_filter: ThresholdState,
    /// Filter this many bp at the start and
    /// end of a read before any mod operations.
    /// Please note that the units here are bp and
    /// not units of base being queried.
    #[clap(long, default_value_t)]
    pub trim_read_ends_mod: usize,
    /// Exclude bases whose base quality is below
    /// this threshold before any mod operation, defaults to 0 i.e. unused.
    /// NOTE: (1) This step is only applied before modification operations, and
    /// not before any other operations.
    /// (2) No offsets such as +33 are needed here.
    /// (3) Modifications on reads where base quality information
    /// is not available are all rejected if this is non-zero.
    #[clap(long, default_value_t)]
    pub base_qual_filter_mod: u8,
    /// Only keep modification data from this region
    #[clap(long)]
    #[builder(field(
        ty = "String",
        build = "(!self.mod_region.is_empty()).then(|| GenomicRegion::from_str(&self.mod_region)).transpose()?"
    ))]
    pub mod_region: Option<GenomicRegion>,
    /// Only keep modification data from this region.
    /// We do not expose this to the user, but infer it from
    /// the other options set by the user.
    /// We cannot populate this directly at the time of CLI parsing,
    /// as we need to look at the BAM header to convert a contig name
    /// into a numeric contig id.
    #[clap(skip)]
    #[builder(setter(strip_option))]
    pub region_bed3: Option<Bed3<i32, u64>>,
}

impl<S: TagState + Args + FromArgMatches + Default> InputModsBuilder<S> {
    /// Validate [`InputMods`] creation through the Builder method.
    ///
    /// [`InputMods`] can be created in many ways:
    /// - directly populating the fields as they are public.
    /// - through the CLI-machinery
    /// - through the builder.
    ///
    /// While we cannot check a user's usage through the first method,
    /// we can check constructions through methods 2 and 3 above.
    /// The CLI construction checking is handled by `clap`, and the builder
    /// checking is done through this function.
    /// The construction checks are slightly different as what may be feasible
    /// through a CLI and through code are different.
    /// So we have this validate method to check construction through
    /// method number 3 above.
    fn validate(&self) -> Result<(), Error> {
        // self.mod_region is a `String` for the Builder before conversion into an `Option<_>`
        if (!self.mod_region.is_empty()) && self.region_bed3.is_some() {
            Err(Error::BuilderValidation(
                "cannot set `mod_region` and `region_bed3` simultaneously!".to_owned(),
            ))
        } else {
            Ok(())
        }
    }
}

/// Implements defaults for `InputMods`
impl<S: TagState + Args + FromArgMatches + Default> Default for InputMods<S> {
    fn default() -> Self {
        InputMods::<S> {
            tag: S::default(),
            mod_strand: None,
            mod_prob_filter: ThresholdState::GtEq(0),
            trim_read_ends_mod: 0,
            base_qual_filter_mod: 0,
            mod_region: None,
            region_bed3: None,
        }
    }
}

/// Retrieves options for modification input
pub trait InputModOptions {
    /// retrieves tag
    fn tag(&self) -> Option<ModChar> {
        unimplemented!()
    }
    /// retrieves option to set basecalled strand or opposite in mod retrieval
    fn mod_strand(&self) -> Option<RestrictModCalledStrand> {
        unimplemented!()
    }
    /// returns probability filter
    fn mod_prob_filter(&self) -> ThresholdState {
        unimplemented!()
    }
    /// returns read end trimming
    fn trim_read_ends_mod(&self) -> usize {
        unimplemented!()
    }
    /// returns threshold for filtering base PHRED quality
    fn base_qual_filter_mod(&self) -> u8 {
        unimplemented!()
    }
}

/// Retrieves options for region
pub trait InputRegionOptions {
    /// returns region requested
    fn region_filter(&self) -> &Option<Bed3<i32, u64>> {
        unimplemented!()
    }
    /// returns region requested but region in genomic string format
    fn region_filter_genomic_string(&self) -> Option<GenomicRegion> {
        unimplemented!()
    }
    /// sets region requested
    fn set_region_filter(&mut self, _value: Option<Bed3<i32, u64>>) {
        unimplemented!()
    }
    /// returns true if full overlap with region is requested as opposed to
    /// only partial overlap. defaults to false.
    fn is_full_overlap(&self) -> bool {
        false
    }
    /// converts region from genomic string representation to bed3 representation
    ///
    /// # Errors
    /// Returns error if conversion fails
    fn convert_region_to_bed3(&mut self, header: bam::HeaderView) -> Result<(), Error> {
        match (
            self.region_filter_genomic_string(),
            self.region_filter().as_ref().is_some(),
        ) {
            (None, false) => self.set_region_filter(None),
            (Some(v), false) => self.set_region_filter(Some(v.try_to_bed3(&header)?)),
            (None, true) => {}
            (Some(_), true) => {
                return Err(Error::InvalidState(
                    "cannot set a region as both a `GenomicRegion` and a `Bed3`".to_owned(),
                ));
            }
        }
        Ok(())
    }
}

impl<S: TagState + Args + FromArgMatches + Default> InputModOptions for InputMods<S> {
    fn tag(&self) -> Option<ModChar> {
        self.tag.tag()
    }
    fn mod_strand(&self) -> Option<RestrictModCalledStrand> {
        self.mod_strand
    }
    fn mod_prob_filter(&self) -> ThresholdState {
        self.mod_prob_filter
    }
    fn trim_read_ends_mod(&self) -> usize {
        self.trim_read_ends_mod
    }
    fn base_qual_filter_mod(&self) -> u8 {
        self.base_qual_filter_mod
    }
}

impl<S: TagState + Args + FromArgMatches + Default> InputRegionOptions for InputMods<S> {
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

/// can return tag without encasing in an option in the `RequiredTag` variant
impl InputMods<RequiredTag> {
    /// retrieves tag
    #[must_use]
    pub fn required_tag(&self) -> ModChar {
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
///
/// # Example
/// ```
/// use nanalogue_core::{Error, InputWindowingBuilder};
///
/// let window = InputWindowingBuilder::default()
///     .win(20)
///     .step(10).build()?;
/// # Ok::<(), Error>(())
/// ```
#[derive(Builder, Debug, Args, Clone, Copy, Serialize, Deserialize)]
#[builder(default, build_fn(error = "Error"))]
#[serde(default)]
#[non_exhaustive]
pub struct InputWindowing {
    /// size of window in units of base being queried i.e.
    /// if you are looking for cytosine modifications, then
    /// a window of a value 300 means create windows each with
    /// 300 cytosines irrespective of their modification status.
    #[clap(long)]
    #[builder(field(ty = "usize", build = "NonZeroUsize::try_from(self.win)?"))]
    pub win: NonZeroUsize,
    /// step window by this size in units of base being queried.
    #[clap(long)]
    #[builder(field(ty = "usize", build = "NonZeroUsize::try_from(self.step)?"))]
    pub step: NonZeroUsize,
}

/// Implements a default for `InputWindowing`.
/// NOTE the defaults of 1 for each are just for ease of programming.
/// We do not expose these defaults to the command-line user.
impl Default for InputWindowing {
    fn default() -> Self {
        InputWindowing {
            win: NonZeroUsize::new(1).expect("no error"),
            step: NonZeroUsize::new(1).expect("no error"),
        }
    }
}

/// This enum contains inputs to display sequences with various options
#[derive(Debug, Clone, Default, Copy, Serialize, Deserialize)]
#[non_exhaustive]
pub enum SeqDisplayOptions {
    /// Do not display any sequence, default
    #[default]
    No,
    /// Display the entire sequence without any subsetting
    Full {
        /// Option to display basecalling qualities
        show_base_qual: bool,
    },
    /// Display the sequence within a region
    Region {
        /// Option to display basecalling qualities
        show_base_qual: bool,
        /// Option to show bases in lowercase if they are insertions
        show_ins_lowercase: bool,
        /// Option to show modified bases as Z (or z depending on other options)
        show_mod_z: bool,
        /// Region over which sequence must be shown
        region: Bed3<i32, u64>,
    },
}

#[cfg(test)]
mod tag_struct_tests {
    use super::*;

    #[test]
    fn optional_tag_default() {
        let optional_tag = OptionalTag::default();
        assert!(optional_tag.tag.is_none());
        assert_eq!(optional_tag.tag(), None);
    }

    #[test]
    fn optional_tag_with_value() {
        let mod_char = ModChar::new('m');
        let optional_tag = OptionalTag {
            tag: Some(mod_char),
        };
        assert_eq!(optional_tag.tag(), Some(mod_char));
    }

    #[test]
    fn required_tag_default() {
        let required_tag = RequiredTag::default();
        assert_eq!(required_tag.tag(), Some(ModChar::default()));
    }

    #[test]
    fn required_tag_with_value() {
        let mod_char = ModChar::new('C');
        let required_tag = RequiredTag { tag: mod_char };
        assert_eq!(required_tag.tag(), Some(mod_char));
    }
}

#[cfg(test)]
mod input_windowing_tests {
    use super::*;

    #[test]
    fn input_windowing_default() {
        let windowing = InputWindowing::default();
        assert_eq!(windowing.win, NonZeroUsize::new(1).unwrap());
        assert_eq!(windowing.step, NonZeroUsize::new(1).unwrap());
    }

    #[test]
    fn input_windowing_custom_values() {
        let windowing = InputWindowing {
            win: NonZeroUsize::new(300).unwrap(),
            step: NonZeroUsize::new(150).unwrap(),
        };
        assert_eq!(windowing.win, NonZeroUsize::new(300).unwrap());
        assert_eq!(windowing.step, NonZeroUsize::new(150).unwrap());
    }
}

#[cfg(test)]
mod input_mods_required_tag_tests {
    use super::*;

    #[test]
    fn input_mods_required_tag_fn_tag() {
        let mod_char = ModChar::new('C');
        let input_mods = InputMods::<RequiredTag> {
            tag: RequiredTag { tag: mod_char },
            mod_strand: None,
            mod_prob_filter: ThresholdState::GtEq(0),
            trim_read_ends_mod: 0,
            base_qual_filter_mod: 0,
            mod_region: None,
            region_bed3: None,
        };

        assert_eq!(input_mods.required_tag(), mod_char);
    }
}

#[cfg(test)]
mod input_bam_tests {
    use super::*;
    use bedrs::Coordinates as _;
    use indoc::indoc;

    #[test]
    fn input_bam_is_full_overlap() {
        // Test default (false)
        let input_bam_default = InputBam::default();
        assert!(!input_bam_default.is_full_overlap());

        // Test explicit false
        let input_bam_false = InputBam {
            full_region: false,
            ..Default::default()
        };
        assert!(!input_bam_false.is_full_overlap());

        // Test true
        let input_bam_true = InputBam {
            full_region: true,
            ..Default::default()
        };
        assert!(input_bam_true.is_full_overlap());
    }

    #[test]
    fn input_bam_convert_region_to_bed3_none() {
        let mut input_bam = InputBam::default();

        let header_view = bam::HeaderView::from_bytes(indoc! {b"@HD\tVN:1.6\tSO:coordinate
        @SQ\tSN:chr1\tLN:248956422\n"});

        input_bam.convert_region_to_bed3(header_view).unwrap();
        assert!(input_bam.region_bed3.is_none());
    }

    #[test]
    fn input_bam_convert_region_to_bed3_with_region() {
        let mut input_bam = InputBam {
            region: Some(GenomicRegion::from_str("chr2:3400-3600").unwrap()),
            ..Default::default()
        };

        let header_view = bam::HeaderView::from_bytes(indoc! {b"@HD\tVN:1.6\tSO:coordinate
                @SQ\tSN:chr1\tLN:3000
                @SQ\tSN:chr2\tLN:4000\n"});

        input_bam.convert_region_to_bed3(header_view).unwrap();
        assert!(input_bam.region_bed3.is_some());

        let bed3 = input_bam.region_bed3.unwrap();
        assert_eq!(bed3.chr(), &1);
        assert_eq!(bed3.start(), 3400);
        assert_eq!(bed3.end(), 3600);
    }

    #[test]
    #[should_panic(expected = "InvalidRegion")]
    fn input_bam_convert_region_to_bed3_invalid_region() {
        let mut input_bam = InputBam {
            region: Some(GenomicRegion::from_str("chr2:4400-4600").expect("no error")),
            ..Default::default()
        };
        let header_view = bam::HeaderView::from_bytes(indoc! {b"@HD\tVN:1.6\tSO:coordinate
                @SQ\tSN:chr1\tLN:3000
                @SQ\tSN:chr2\tLN:4000\n"});

        input_bam.convert_region_to_bed3(header_view).unwrap();
    }

    #[test]
    #[should_panic(expected = "InvalidRegion")]
    fn input_bam_convert_region_to_bed3_invalid_open_ended_region() {
        let mut input_bam = InputBam {
            region: Some(GenomicRegion::from_str("chr2:4600-").expect("no error")),
            ..Default::default()
        };
        let header_view = bam::HeaderView::from_bytes(indoc! {b"@HD\tVN:1.6\tSO:coordinate
                @SQ\tSN:chr1\tLN:3000
                @SQ\tSN:chr2\tLN:4000\n"});

        input_bam.convert_region_to_bed3(header_view).unwrap();
    }

    #[test]
    #[should_panic(expected = "InvalidAlignCoords")]
    fn input_bam_convert_region_to_bed3_invalid_contig() {
        let mut input_bam = InputBam {
            region: Some(GenomicRegion::from_str("chr3:1000-2000").expect("no error")),
            ..Default::default()
        };
        let header_view = bam::HeaderView::from_bytes(indoc! {b"@HD\tVN:1.6\tSO:coordinate
                @SQ\tSN:chr1\tLN:3000
                @SQ\tSN:chr2\tLN:4000\n"});

        input_bam.convert_region_to_bed3(header_view).unwrap();
    }
}
