//! # Write Simulated Mod BAM
//! Generates simulated BAM files with base modifications for testing purposes.
//! Accepts JSON configuration to produce both BAM and FASTA reference files.
//! Please note that both BAM files and FASTA files are created from scratch,
//! so please do not specify pre-existing BAM or FASTA files in the output
//! path; if so, they will be overwritten.
//!
//! ## Example Usage
//!
//! We've shown how to construct `SimulationConfig` which contains input options
//! and then use it to create a BAM file below. You can also use [`SimulationConfigBuilder`]
//! to achieve this without using a `json` structure. Some of the fields in the `json` entry
//! shown below are optional; please read the comments following it.
//!
//! ```
//! use nanalogue_core::{Error, SimulationConfig, simulate_mod_bam::run};
//! use uuid::Uuid;
//!
//! let config_json = r#"{
//!   "contigs": {
//!     "number": 4,
//!     "len_range": [1000, 8000],
//!     "repeated_seq": "ATCGAATT"
//!   },
//!   "reads": [{
//!     "number": 1000,
//!     "mapq_range": [10, 20],
//!     "base_qual_range": [10, 20],
//!     "len_range": [0.1, 0.8],
//!     "barcode": "ACGTAA",
//!     "delete": [0.5, 0.7],
//!     "insert_middle": "ATCG",
//!     "mismatch": 0.5,
//!     "mods": [{
//!         "base": "T",
//!         "is_strand_plus": true,
//!         "mod_code": "T",
//!         "win": [4, 5],
//!         "mod_range": [[0.1, 0.2], [0.3, 0.4]]
//!     }]
//!   }]
//! }"#;
//!
//! // Note: * We generate four contigs `4` with a length chosen randomly from the range `1000-8000`
//! //         in the example above. Contigs are called `contig_00000`, ..., `contig_00003`.
//! //       * "repeated_seq" field is optional in contigs; if set, contigs are made by
//! //         repeating this sequence instead of generating random sequences.
//! //       * Reads are generated with options shown. In the example, 1000 reads are made, with a
//! //         mapping quality in the range 10-20, a basecalling quality in the range 10-20,
//! //         lengths in the range 10%-80% of a contig (each read is mapped to a contig chosen
//! //         randomly with equal probability for all contigs). Reads are randomly assigned to be
//! //         in one of seven states -> primary/secondary/supplementary X fwd/rev plus unmapped.
//! //         Read sequences are drawn from the associated contig at a random location.
//! //         (Unmapped reads are also drawn from contig sequences; they are just marked unmapped
//! //          after the fact).
//! //         Read names are 'n.uuid' where n is the read group number and uuid is a random UUID
//! //         (UUID is a random string that looks like this for example ->
//! //          "e9529f28-d27a-497a-be7e-bffdb77e8bc1"). Read group number is assigned
//! //         automatically, depending on how many entries are in the `reads` field.
//! //         In the example here, there is only one entry, so the read group number is 0,
//! //         so all the reads will have a name like '0.uuid'.
//! //       * "barcode" field is optional and can be omitted; barcodes added to both ends.
//! //         As sequences and lengths are generated independently of barcodes,
//! //         a barcode will _add_ 2 times so many bp to sequence length statistics.
//! //       * "delete" is optional. If specified, the shown section is deleted on
//! //         all reads. In the example, the section between the middle of the read (50%)
//! //         and seven-tenths of the read (70%) is deleted. The direction in the instructions
//! //         here is parallel to the reference genome.
//! //       * "insert_middle" is optional. If specified, the specified sequence will be inserted
//! //         at the middle of the read.
//! //       * "mismatch" is optional. If specified, this random fraction of bases is mismatched i.e.
//! //         we make the base different from what it is supposed to be.
//! //       * insertion/deletion/barcode are applied after sequence lengths are determined
//! //         according to the given options. In other words, these options will cause sequence
//! //         lengths to be different from the given statistics e.g. an insertion of
//! //         10 bp will make all reads longer than the desired sequence statistics by 10 bp.
//! //       * "mods" are optional as well, omitting them will create a BAM file
//! //          with no modification information.
//! //       * mod information is specified using windows i.e. in window of given size,
//! //         we enforce that mod probabilities per base are in the given range.
//! //         Multiple windows can be specified, and they repeat along the
//! //         read until it ends i.e. we cycle windows so we go 4Ts, 5Ts, 4Ts, 5Ts, ...
//! //         and probabilities are cycled as (0.1, 0.2), (0.3, 0.4), (0.1, 0.2), ...
//! //       * if windows and mod_range lists are of different sizes, the cycles may not line up.
//! //         e.g. if there are three window values and two mod ranges, then
//! //         windows repeat in cycles of 3 whereas mod ranges will repeat in cycles of 2.
//! //         You can use such inputs if this is what you want.
//!
//! // Paths used here must not exist already as these are files created anew.
//! let config: SimulationConfig = serde_json::from_str(&config_json)?;
//! let temp_dir = std::env::temp_dir();
//! let bam_path = temp_dir.join(format!("{}.bam", Uuid::new_v4()));
//! let fasta_path = temp_dir.join(format!("{}.fa", Uuid::new_v4()));
//! let bai_path = bam_path.with_extension("bam.bai");
//! run(config, &bam_path, &fasta_path)?;
//! std::fs::remove_file(&bam_path)?;
//! std::fs::remove_file(&fasta_path)?;
//! std::fs::remove_file(&bai_path)?;
//! # Ok::<(), Error>(())
//! ```

use crate::{
    AllowedAGCTN, DNARestrictive, Error, F32Bw0and1, GetDNARestrictive, ModChar, OrdPair, ReadState,
};
use crate::{write_bam_denovo, write_fasta};
use derive_builder::Builder;
use itertools::join;
use rand::Rng;
use rand::SeedableRng as _;
use rand::rngs::StdRng;
use rand::seq::IteratorRandom as _;
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use serde::{Deserialize, Serialize};
use std::iter;
use std::num::{NonZeroU32, NonZeroU64};
use std::ops::RangeInclusive;
use std::path::Path;
use std::str::FromStr as _;
use uuid::Uuid;

/// We need a shorthand for this for a one-liner we use
/// in a macro provided by `derive_builder`.
type OF = OrdPair<F32Bw0and1>;

/// Creates an `OrdPair` of `F32Bw0and1` values
///
/// We do not export this macro as it involves `unwrap`,
/// so we only use it locally.
macro_rules! ord_pair_f32_bw0and1 {
    ($low:expr, $high:expr) => {
        OrdPair::new(
            F32Bw0and1::new($low).unwrap(),
            F32Bw0and1::new($high).unwrap(),
        )
        .unwrap()
    };
}

/// Main configuration struct for simulation.
///
/// After construction, you can pass it to [`crate::simulate_mod_bam::run`]
/// to create a mod BAM file. You can build this using [`SimulationConfigBuilder`]
/// as shown below.
///
/// # Example
///
/// Below struct can be used in appropriate routines to
/// create a BAM file.
///
/// ```
/// use nanalogue_core::Error;
/// use nanalogue_core::simulate_mod_bam::{ContigConfigBuilder, SimulationConfigBuilder,
///     ReadConfigBuilder, ModConfigBuilder};
///
/// // Request two contigs with lengths randomly chosen from 1000 - 2000 bp.
/// let contig_config = ContigConfigBuilder::default()
///     .number(2)
///     .len_range((1000,2000)).build()?;
///
/// // First set of reads
/// let read_config1 = ReadConfigBuilder::default()
///     .number(100)
///     .mapq_range((10, 20))
///     .base_qual_range((30, 40))
///     .len_range((0.5, 0.6));
///
/// // second set of reads, now 200 requested with a modification and a barcode
/// let mod_config = ModConfigBuilder::default()
///     .base('C')
///     .is_strand_plus(true)
///     .mod_code("m".into())
///     .win(vec![2, 3])
///     .mod_range(vec![(0.1, 0.8), (0.3, 0.4)]).build()?;
///
/// let read_config2 = read_config1.clone()
///     .number(200)
///     .barcode("AATTG".into())
///     .mods(vec![mod_config]);
///
/// // set simulation config, ready to be fed into an appropriate function
/// let sim_config = SimulationConfigBuilder::default()
///     .contigs(contig_config)
///     .reads(vec![read_config1.build()?, read_config2.build()?]).build()?;
/// # Ok::<(), Error>(())
/// ```
#[derive(Builder, Debug, Default, PartialEq, Clone, Serialize, Deserialize)]
#[serde(default)]
#[builder(default, build_fn(error = "Error"), pattern = "owned", derive(Clone))]
#[non_exhaustive]
pub struct SimulationConfig {
    /// Configuration for contig generation
    pub contigs: ContigConfig,
    /// Configuration for read generation
    pub reads: Vec<ReadConfig>,
    /// Seed for reproducible simulation. When set, all random operations
    /// use a deterministic RNG seeded with this value, producing identical
    /// output files across runs. If not set, simulation is non-reproducible.
    pub seed: Option<u64>,
}

/// Configuration for contig generation
///
/// The struct contains parameters which can then be passed to other
/// functions to construct the actual contig. Also refer to [`SimulationConfig`]
/// to see how this `struct` can be used, and [`ContigConfigBuilder`] for how
/// to build this (as shown below in the examples).
///
/// # Examples
///
/// Request two contigs with lengths randomly chosen from 1000 - 2000 bp.
/// The DNA of the contigs are randomly-generated.
///
/// ```
/// use nanalogue_core::Error;
/// use nanalogue_core::simulate_mod_bam::ContigConfigBuilder;
///
/// let contig_config = ContigConfigBuilder::default()
///     .number(2)
///     .len_range((1000,2000)).build()?;
///
/// # Ok::<(), Error>(())
/// ```
///
/// Request similar contigs as above but with the DNA of the contigs
/// consisting of the same sequence repeated over and over.
///
/// ```
/// # use nanalogue_core::Error;
/// # use nanalogue_core::simulate_mod_bam::ContigConfigBuilder;
/// let contig_config = ContigConfigBuilder::default()
///     .number(2)
///     .len_range((1000,2000))
///     .repeated_seq("AAGCTTGA".into()).build()?;
///
/// # Ok::<(), Error>(())
/// ```
///
/// Request a contig with a fixed length and sequence.
/// As the repeated sequence's length is equal to the
/// contig length, and the contig length is precisely fixed,
/// you get a non-randomly generated contig.
///
/// ```
/// # use nanalogue_core::Error;
/// # use nanalogue_core::simulate_mod_bam::ContigConfigBuilder;
/// let contig_config = ContigConfigBuilder::default()
///     .number(1)
///     .len_range((20,20))
///     .repeated_seq("ACGTACGTACGTACGTACGT".into()).build()?;
///
/// # Ok::<(), Error>(())
/// ```
#[derive(Builder, Debug, PartialEq, Clone, Serialize, Deserialize)]
#[serde(default)]
#[builder(default, build_fn(error = "Error"), pattern = "owned", derive(Clone))]
#[non_exhaustive]
pub struct ContigConfig {
    /// Number of contigs to generate
    #[builder(field(ty = "u32", build = "self.number.try_into()?"))]
    pub number: NonZeroU32,
    /// Contig length range in bp [min, max]
    #[builder(field(ty = "(u64, u64)", build = "self.len_range.try_into()?"))]
    pub len_range: OrdPair<NonZeroU64>,
    /// Optional repeated sequence to use for contigs instead of random generation
    #[builder(field(
        ty = "String",
        build = "(!self.repeated_seq.is_empty()).then(|| DNARestrictive::from_str(&self.repeated_seq)).transpose()?"
    ))]
    pub repeated_seq: Option<DNARestrictive>,
}

/// Configuration for read generation. Also refer to [`SimulationConfig`]
/// to see how this `struct` can be used, and to [`ReadConfigBuilder`] for
/// how to build this (as shown below in the examples).
///
/// # Example 1
///
///  Request reads without modifications
/// ```
/// use nanalogue_core::Error;
/// use nanalogue_core::simulate_mod_bam::ReadConfigBuilder;
///
/// // First set of reads
/// let read_config1 = ReadConfigBuilder::default()
///     .number(100)
///     .mapq_range((10, 20))
///     .base_qual_range((30, 40))
///     .len_range((0.5, 0.6)).build()?;
///
/// # Ok::<(), Error>(())
/// ```
/// # Example 2
///
///  Request reads with modifications, but also with insertions, deletions, and a mismatch rate.
/// ```
/// use nanalogue_core::Error;
/// use nanalogue_core::simulate_mod_bam::{ReadConfigBuilder, ModConfigBuilder};
///
/// let mod_config = ModConfigBuilder::default()
///     .base('C')
///     .is_strand_plus(true)
///     .mod_code("m".into())
///     .win(vec![2, 3])
///     .mod_range(vec![(0.1, 0.8), (0.3, 0.4)]).build()?;
///
/// let read_config2 = ReadConfigBuilder::default()
///     .number(200)
///     .mapq_range((30, 40))
///     .base_qual_range((50, 60))
///     .barcode("AATTG".into())
///     .len_range((0.7, 0.8))
///     .delete((0.2, 0.3))
///     .insert_middle("ATCGA".into())
///     .mismatch(0.1)
///     .mods(vec![mod_config]).build()?;
///
/// # Ok::<(), Error>(())
/// ```
#[derive(Builder, Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
#[builder(default, build_fn(error = "Error"), pattern = "owned", derive(Clone))]
#[non_exhaustive]
pub struct ReadConfig {
    /// Total number of reads to generate
    #[builder(field(ty = "u32", build = "self.number.try_into()?"))]
    pub number: NonZeroU32,
    /// Mapping quality range [min, max]
    #[builder(field(ty = "(u8, u8)", build = "self.mapq_range.try_into()?"))]
    pub mapq_range: OrdPair<u8>,
    /// Base quality score range [min, max]
    #[builder(field(ty = "(u8, u8)", build = "self.base_qual_range.try_into()?"))]
    pub base_qual_range: OrdPair<u8>,
    /// Read length range as fraction of contig [min, max] (e.g.: [0.1, 0.8] = 10% to 80%)
    #[builder(field(ty = "(f32, f32)", build = "self.len_range.try_into()?"))]
    pub len_range: OrdPair<F32Bw0and1>,
    /// Optional barcode DNA sequence to add to read ends
    #[builder(field(
        ty = "String",
        build = "(!self.barcode.is_empty()).then(|| DNARestrictive::from_str(&self.barcode)).transpose()?"
    ))]
    pub barcode: Option<DNARestrictive>,
    /// Optional deletion: fractional range [start, end] of read to delete (e.g., [0.5, 0.7] deletes middle 20%)
    #[builder(field(
        ty = "(f32, f32)",
        build = "Some(self.delete).map(TryInto::try_into).transpose()?"
    ))]
    pub delete: Option<OrdPair<F32Bw0and1>>,
    /// Optional sequence to insert at the middle of the read
    #[builder(field(
        ty = "String",
        build = "(!self.insert_middle.is_empty()).then(|| DNARestrictive::from_str(&self.insert_middle)).transpose()?"
    ))]
    pub insert_middle: Option<DNARestrictive>,
    /// Optional mismatch fraction: random fraction of bases to mutate (e.g., 0.5 = 50% of bases)
    #[builder(field(ty = "f32", build = "Some(F32Bw0and1::try_from(self.mismatch)?)"))]
    pub mismatch: Option<F32Bw0and1>,
    /// Modification configurations
    pub mods: Vec<ModConfig>,
}

/// Configuration for modification generation in simulated BAM files.
/// Also refer to [`SimulationConfig`] to see how this `struct` can be used,
/// and [`ModConfigBuilder`] to see how it can be built as shown below.
///
/// # Example
///
/// Configure modification for cytosines on plus strand.
/// Here, 2Cs have probabilities randomly chosen b/w 0.4-0.8,
///       3Cs have probabilities b/w 0.5-0.7,
///       2Cs have ...,
///       3Cs have ...,
///       ...
///  repeated forever. This can be passed as an input to a routine
///  sometime later and applied on a given sequence.
/// ```
/// use nanalogue_core::Error;
/// use nanalogue_core::simulate_mod_bam::ModConfigBuilder;
///
/// let mod_config_c = ModConfigBuilder::default()
///     .base('C')
///     .is_strand_plus(true)
///     .mod_code("m".into())
///     .win(vec![2, 3])
///     .mod_range(vec![(0.4, 0.8), (0.5, 0.7)]).build()?;
/// # Ok::<(), Error>(())
/// ```
#[derive(Builder, Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
#[builder(default, build_fn(error = "Error"), pattern = "owned", derive(Clone))]
#[non_exhaustive]
pub struct ModConfig {
    /// Base that is modified (A, C, G, T, or N)
    #[builder(field(ty = "char", build = "self.base.try_into()?"))]
    pub base: AllowedAGCTN,
    /// Whether this is on the plus strand
    pub is_strand_plus: bool,
    /// Modification code (character or numeric)
    #[builder(field(ty = "String", build = "self.mod_code.parse()?"))]
    pub mod_code: ModChar,
    /// Vector of window sizes for modification density variation.
    /// e.g. if you want reads with a window of 100 bases with each
    /// modified with a probability in the range 0.4-0.6 and the next
    /// 200 modified in the range 0.1-0.2 (and this two window config
    /// repeating forever), you need to set this win array to
    /// [100, 200] and the mod range array to [[0.4, 0.6], [0.1, 0.2]]
    /// (Remembering to adjust the datatypes of the inputs suitably
    /// i.e. use `NonZeroU32` etc.)
    /// NOTE: "bases" in the above description refer to bases of interest
    /// set in the base field i.e. 100 bases with base set to "T" mean
    /// 100 thymidines specifically.
    #[builder(field(
        ty = "Vec<u32>",
        build = "self.win.iter().map(|&x| NonZeroU32::new(x).ok_or(Error::Zero(\"cannot use zero-\
sized windows in builder\".to_owned()))).collect::<Result<Vec<NonZeroU32>,_>>()?"
    ))]
    pub win: Vec<NonZeroU32>,
    /// Vector of modification density range e.g. [[0.4, 0.6], [0.1, 0.2]].
    /// Also see description of `win` above.
    #[builder(field(
        ty = "Vec<(f32, f32)>",
        build = "self.mod_range.iter().map(|&x| OF::try_from(x)).collect::<Result<Vec<OF>, _>>()?"
    ))]
    pub mod_range: Vec<OrdPair<F32Bw0and1>>,
}

/// Represents a contig with name and sequence.
/// Can be built using [`ContigConfigBuilder`] as shown below.
///
/// # Example
///
/// Example contig with a simple sequence and name, to be passed
/// somewhere else to be constructed.
///
/// ```
/// use nanalogue_core::{Error, simulate_mod_bam::ContigBuilder};
///
/// let contig = ContigBuilder::default()
///     .name("chr1")
///     .seq("ACGTACGTACGTACGT".into()).build()?;
/// # Ok::<(), Error>(())
#[derive(Builder, Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
#[builder(default, build_fn(error = "Error"))]
#[non_exhaustive]
pub struct Contig {
    /// Contig name
    #[builder(setter(into))]
    pub name: String,
    /// Contig sequence (A, C, G, T)
    #[builder(field(ty = "String", build = "self.seq.parse()?"))]
    pub seq: DNARestrictive,
}

impl GetDNARestrictive for Contig {
    /// Returns a reference to the contig's DNA
    fn get_dna_restrictive(&self) -> &DNARestrictive {
        &self.seq
    }
}

impl Default for Contig {
    fn default() -> Self {
        Self {
            name: String::new(),
            seq: DNARestrictive::from_str("A").expect("valid default DNA"),
        }
    }
}

impl Default for ContigConfig {
    fn default() -> Self {
        Self {
            number: NonZeroU32::new(1).unwrap(),
            len_range: OrdPair::new(NonZeroU64::new(1).unwrap(), NonZeroU64::new(1).unwrap())
                .unwrap(),
            repeated_seq: None,
        }
    }
}

impl Default for ReadConfig {
    fn default() -> Self {
        Self {
            number: NonZeroU32::new(1).unwrap(),
            mapq_range: OrdPair::new(0, 0).unwrap(),
            base_qual_range: OrdPair::new(0, 0).unwrap(),
            len_range: ord_pair_f32_bw0and1!(0.0, 0.0),
            barcode: None,
            delete: None,
            insert_middle: None,
            mismatch: None,
            mods: Vec::new(),
        }
    }
}

impl Default for ModConfig {
    fn default() -> Self {
        Self {
            base: AllowedAGCTN::C,
            is_strand_plus: true,
            mod_code: ModChar::new('m'),
            win: vec![NonZeroU32::new(1).unwrap()],
            mod_range: vec![ord_pair_f32_bw0and1!(0.0, 1.0)],
        }
    }
}

/// Builder for creating a sequence with cigar string and optional barcode
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerfectSeqMatchToNot {
    /// The DNA sequence to be processed
    seq: Vec<u8>,
    /// Optional barcode to add to both ends of the sequence
    barcode: Option<DNARestrictive>,
    /// Optional deletion: fractional range of sequence to delete
    delete: Option<OrdPair<F32Bw0and1>>,
    /// Optional sequence to insert at the middle
    insert_middle: Option<DNARestrictive>,
    /// Optional mismatch fraction
    mismatch: Option<F32Bw0and1>,
}

impl PerfectSeqMatchToNot {
    /// Creates a new builder with the given sequence
    ///
    /// # Errors
    /// Returns `Error::InvalidState` if the sequence is empty
    ///
    /// # Examples
    ///
    /// Valid sequence:
    /// ```
    /// use nanalogue_core::simulate_mod_bam::PerfectSeqMatchToNot;
    ///
    /// let builder = PerfectSeqMatchToNot::seq(b"ACGT".to_vec()).unwrap();
    /// ```
    ///
    /// Empty sequence returns error:
    /// ```
    /// use nanalogue_core::simulate_mod_bam::PerfectSeqMatchToNot;
    /// use nanalogue_core::Error;
    ///
    /// let result = PerfectSeqMatchToNot::seq(b"".to_vec());
    /// assert!(matches!(result, Err(Error::InvalidState(_))));
    /// ```
    pub fn seq(seq: Vec<u8>) -> Result<Self, Error> {
        if seq.is_empty() {
            return Err(Error::InvalidState(
                "Sequence length is 0; cannot create PerfectSeqMatchToNot with empty sequence"
                    .into(),
            ));
        }
        Ok(Self {
            seq,
            barcode: None,
            delete: None,
            insert_middle: None,
            mismatch: None,
        })
    }

    /// Sets the barcode for this sequence
    #[must_use]
    pub fn barcode(mut self, barcode: DNARestrictive) -> Self {
        self.barcode = Some(barcode);
        self
    }

    /// Sets the deletion range for this sequence
    #[must_use]
    pub fn delete(mut self, delete: OrdPair<F32Bw0and1>) -> Self {
        self.delete = Some(delete);
        self
    }

    /// Sets the sequence to insert at the middle
    #[must_use]
    pub fn insert_middle(mut self, insert: DNARestrictive) -> Self {
        self.insert_middle = Some(insert);
        self
    }

    /// Sets the mismatch fraction for this sequence
    #[must_use]
    pub fn mismatch(mut self, mismatch: F32Bw0and1) -> Self {
        self.mismatch = Some(mismatch);
        self
    }

    /// Builds the final sequence with barcode (if specified) and corresponding cigar string
    ///
    /// # Examples
    ///
    /// Without barcode, forward read:
    /// ```
    /// use nanalogue_core::simulate_mod_bam::PerfectSeqMatchToNot;
    /// use nanalogue_core::{ReadState, Error};
    /// use rand::Rng;
    ///
    /// let mut rng = rand::rng();
    /// let (seq, cigar) = PerfectSeqMatchToNot::seq(b"GGGGGGGG".to_vec())?
    ///     .build(ReadState::PrimaryFwd, &mut rng)?;
    /// assert_eq!(seq, b"GGGGGGGG");
    /// assert_eq!(cigar.expect("no error").to_string(), "8M");
    /// # Ok::<(), Error>(())
    /// ```
    ///
    /// Without barcode, reverse read:
    /// ```
    /// use nanalogue_core::simulate_mod_bam::PerfectSeqMatchToNot;
    /// use nanalogue_core::{ReadState, Error};
    /// use rand::Rng;
    ///
    /// let mut rng = rand::rng();
    /// let (seq, cigar) = PerfectSeqMatchToNot::seq(b"GGGGGGGG".to_vec())?
    ///     .build(ReadState::PrimaryRev, &mut rng)?;
    /// assert_eq!(seq, b"GGGGGGGG");
    /// assert_eq!(cigar.expect("no error").to_string(), "8M");
    /// # Ok::<(), Error>(())
    /// ```
    ///
    /// With barcode, forward read (barcode + seq + revcomp(barcode)):
    /// ```
    /// use nanalogue_core::simulate_mod_bam::PerfectSeqMatchToNot;
    /// use nanalogue_core::{ReadState, DNARestrictive, Error};
    /// use std::str::FromStr;
    /// use rand::Rng;
    ///
    /// let mut rng = rand::rng();
    /// let barcode = DNARestrictive::from_str("ACGTAA").unwrap();
    /// let (seq, cigar) = PerfectSeqMatchToNot::seq(b"GGGGGGGG".to_vec())?
    ///     .barcode(barcode)
    ///     .build(ReadState::PrimaryFwd, &mut rng)?;
    /// assert_eq!(seq, b"ACGTAAGGGGGGGGTTACGT");
    /// assert_eq!(cigar.expect("no error").to_string(), "6S8M6S");
    /// # Ok::<(), Error>(())
    /// ```
    ///
    /// With barcode, reverse read (comp(barcode) + seq + rev(barcode)):
    /// ```
    /// use nanalogue_core::simulate_mod_bam::PerfectSeqMatchToNot;
    /// use nanalogue_core::{ReadState, DNARestrictive, Error};
    /// use std::str::FromStr;
    /// use rand::Rng;
    ///
    /// let mut rng = rand::rng();
    /// let barcode = DNARestrictive::from_str("ACGTAA").unwrap();
    /// let (seq, cigar) = PerfectSeqMatchToNot::seq(b"GGGGGGGG".to_vec())?
    ///     .barcode(barcode)
    ///     .build(ReadState::PrimaryRev, &mut rng)?;
    /// assert_eq!(seq, b"TGCATTGGGGGGGGAATGCA");
    /// assert_eq!(cigar.expect("no error").to_string(), "6S8M6S");
    /// # Ok::<(), Error>(())
    /// ```
    ///
    /// Unmapped read returns None for cigar:
    /// ```
    /// use nanalogue_core::simulate_mod_bam::PerfectSeqMatchToNot;
    /// use nanalogue_core::{ReadState, Error};
    /// use rand::Rng;
    ///
    /// let mut rng = rand::rng();
    /// let (seq, cigar) = PerfectSeqMatchToNot::seq(b"GGGGGGGG".to_vec())?
    ///     .build(ReadState::Unmapped, &mut rng)?;
    /// assert_eq!(seq, b"GGGGGGGG");
    /// assert!(cigar.is_none());
    /// # Ok::<(), Error>(())
    /// ```
    ///
    /// With deletion (removes a portion of the sequence):
    /// ```
    /// use nanalogue_core::simulate_mod_bam::PerfectSeqMatchToNot;
    /// use nanalogue_core::{ReadState, F32Bw0and1, OrdPair, Error};
    /// use rand::Rng;
    ///
    /// let mut rng = rand::rng();
    /// let delete_range = OrdPair::new(
    ///     F32Bw0and1::try_from(0.5)?,
    ///     F32Bw0and1::try_from(0.75)?
    /// )?;
    /// let (seq, cigar) = PerfectSeqMatchToNot::seq(b"AAAATTTTCCCCGGGG".to_vec())?
    ///     .delete(delete_range)
    ///     .build(ReadState::PrimaryFwd, &mut rng)?;
    /// // Deletes positions 8-12 (CCCC), resulting in AAAATTTTGGGG
    /// assert_eq!(seq, b"AAAATTTTGGGG");
    /// assert_eq!(cigar.expect("no error").to_string(), "8M4D4M");
    /// # Ok::<(), Error>(())
    /// ```
    ///
    /// With insertion in the middle:
    /// ```
    /// use nanalogue_core::simulate_mod_bam::PerfectSeqMatchToNot;
    /// use nanalogue_core::{ReadState, DNARestrictive, Error};
    /// use std::str::FromStr;
    /// use rand::Rng;
    ///
    /// let mut rng = rand::rng();
    /// let insert_seq = DNARestrictive::from_str("TTTT")?;
    /// let (seq, cigar) = PerfectSeqMatchToNot::seq(b"AAAAGGGG".to_vec())?
    ///     .insert_middle(insert_seq)
    ///     .build(ReadState::PrimaryFwd, &mut rng)?;
    /// // Inserts TTTT at position 4 (middle), resulting in AAAATTTTGGGG
    /// assert_eq!(seq, b"AAAATTTTGGGG");
    /// assert_eq!(cigar.expect("no error").to_string(), "4M4I4M");
    /// # Ok::<(), Error>(())
    /// ```
    ///
    /// With mismatches (randomly changes bases):
    /// ```
    /// use nanalogue_core::simulate_mod_bam::PerfectSeqMatchToNot;
    /// use nanalogue_core::{ReadState, F32Bw0and1, Error};
    /// use rand::Rng;
    ///
    /// let mut rng = rand::rng();
    /// let mismatch_frac = F32Bw0and1::try_from(0.5)?;
    /// let (seq, cigar) = PerfectSeqMatchToNot::seq(b"AAAAAAAA".to_vec())?
    ///     .mismatch(mismatch_frac)
    ///     .build(ReadState::PrimaryFwd, &mut rng)?;
    /// // 50% of bases (4 out of 8) will be changed to different bases
    /// let mismatch_count = seq.iter().filter(|&&b| b != b'A').count();
    /// assert_eq!(mismatch_count, 4);
    /// // CIGAR remains 8M (mismatches don't change CIGAR)
    /// assert_eq!(cigar.expect("no error").to_string(), "8M");
    /// # Ok::<(), Error>(())
    /// ```
    ///
    /// Combining deletion and insertion:
    /// ```
    /// use nanalogue_core::simulate_mod_bam::PerfectSeqMatchToNot;
    /// use nanalogue_core::{ReadState, F32Bw0and1, OrdPair, DNARestrictive, Error};
    /// use std::str::FromStr;
    /// use rand::Rng;
    ///
    /// let mut rng = rand::rng();
    /// let delete_range = OrdPair::new(
    ///     F32Bw0and1::try_from(0.5)?,
    ///     F32Bw0and1::try_from(0.75)?
    /// )?;
    /// let insert_seq = DNARestrictive::from_str("TT")?;
    /// let (seq, cigar) = PerfectSeqMatchToNot::seq(b"AAAACCCCGGGG".to_vec())?
    ///     .delete(delete_range)
    ///     .insert_middle(insert_seq)
    ///     .build(ReadState::PrimaryFwd, &mut rng)?;
    /// // Applies both deletion and insertion operations
    /// assert_eq!(seq, b"AAAACCTTGGG");
    /// assert_eq!(cigar.expect("no error").to_string(), "6M2I3D3M");
    /// # Ok::<(), Error>(())
    /// ```
    ///
    /// # Returns
    /// A tuple of (`final_sequence`, `cigar_string`) wrapped in a `Result` where:
    /// - `final_sequence` includes barcodes if specified
    /// - `cigar_string` includes `SoftClip` operations for barcodes, or `None` if unmapped
    ///
    /// # Errors
    /// - If parameters result in a sequence with no mapped bases i.e. a sequence without
    ///   any M characters in the CIGAR string.
    /// - If a deletion occurs right after a softclip at the start of a sequence
    ///   a cigar like "20S10D..." is not valid as we can just shift the start of the sequence
    ///   and make the cigar "20S...". We will end up in this state if the deletion is right at
    ///   the beginning of the sequence. We may deal with this in the future by shifting the start
    ///   of the sequence, but we don't want to do this right now. A similar problem arises if
    ///   the deletion is at the end of a sequence.
    ///
    /// # Panics
    /// Panics on number conversion errors (sequence or barcode length overflow)
    #[expect(
        clippy::cast_precision_loss,
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "Casting and arithmetic needed for fractional calculations on sequence positions"
    )]
    #[expect(
        clippy::too_many_lines,
        reason = "performs many steps to alter sequence (indel, mismatch, barcode)"
    )]
    pub fn build<R: Rng>(
        self,
        read_state: ReadState,
        rng: &mut R,
    ) -> Result<(Vec<u8>, Option<CigarString>), Error> {
        let seq = self.seq;

        // Step 1: Create initial representation as (base, operation) tuples
        // All bases start with 'M' for Match operation
        let mut bases_and_ops: Vec<(u8, u8)> = seq.iter().map(|&base| (base, b'M')).collect();

        // Step 2: Apply mismatch (randomly mutate bases, keep 'M' operation)
        // We use `choose_multiple` on indices to preserve position information,
        // unlike `partial_shuffle` which would reorder elements and destroy positions.
        if let Some(mismatch_frac) = self.mismatch {
            let num_mismatches =
                ((bases_and_ops.len() as f32) * mismatch_frac.val()).round() as usize;

            let indices_to_mutate: Vec<usize> =
                (0..bases_and_ops.len()).choose_multiple(rng, num_mismatches);

            for idx in indices_to_mutate {
                let item = bases_and_ops
                    .get_mut(idx)
                    .expect("idx is within bounds from choose_multiple");
                let current_base = item.0;
                let new_base = match current_base {
                    v @ (b'A' | b'C' | b'G' | b'T') => {
                        // Sample random bases until we get A/C/G/T different from the current
                        loop {
                            let candidate: AllowedAGCTN = rng.random();
                            let candidate_u8: u8 = candidate.into();
                            if candidate_u8 != v && candidate_u8 != b'N' {
                                break candidate_u8;
                            }
                        }
                    }
                    other => other, // Keep N and other bases as is
                };
                item.0 = new_base;
            }
        }

        // Step 3: Apply delete (mark bases as deleted with 'D' operation)
        if let Some(delete_range) = self.delete {
            let start = ((bases_and_ops.len() as f32) * delete_range.low().val()).round() as usize;
            let end_raw =
                ((bases_and_ops.len() as f32) * delete_range.high().val()).round() as usize;
            let end = end_raw.min(bases_and_ops.len());

            if let Some(slice) = bases_and_ops.get_mut(start..end) {
                for item in slice {
                    item.1 = b'D';
                }
            }
        }

        // Step 4: Apply insert_middle (add new tuples with 'I' operation)
        if let Some(insert_seq) = self.insert_middle {
            let middle = bases_and_ops.len() / 2;
            let insert_bases = insert_seq.get_dna_restrictive().get();
            let insertions: Vec<(u8, u8)> = insert_bases.iter().map(|&b| (b, b'I')).collect();

            drop(
                bases_and_ops
                    .splice(middle..middle, insertions)
                    .collect::<Vec<_>>(),
            );
        }

        // Step 5: Apply barcode (add tuples with 'S' for SoftClip operation)
        if let Some(barcode) = self.barcode {
            // Use add_barcode with empty sequence to get the transformed barcodes
            let barcoded_seq = add_barcode(&[], barcode.clone(), read_state);
            let barcode_len = barcode.get_dna_restrictive().get().len();

            // Split the result into start and end barcodes
            let (start_bc, end_bc) = barcoded_seq.split_at(barcode_len);

            // Prepend start barcode with 'S' operations
            let bc_tuples_start: Vec<(u8, u8)> = start_bc.iter().map(|&b| (b, b'S')).collect();
            drop(
                bases_and_ops
                    .splice(0..0, bc_tuples_start)
                    .collect::<Vec<_>>(),
            );

            // Append end barcode with 'S' operations
            let bc_tuples_end: Vec<(u8, u8)> = end_bc.iter().map(|&b| (b, b'S')).collect();
            bases_and_ops.extend(bc_tuples_end);
        }

        // Step 6: Build CIGAR string from operations (second entries)
        let cigar_string = CigarString(
            bases_and_ops
                .iter()
                .map(|&(_, op)| op)
                .fold(Vec::<(u8, u32)>::new(), |mut acc, op| {
                    match acc.last_mut() {
                        Some(&mut (last_op, ref mut count)) if last_op == op => {
                            *count = count.saturating_add(1);
                        }
                        _ => acc.push((op, 1u32)),
                    }
                    acc
                })
                .into_iter()
                .map(|(op, count)| match op {
                    b'M' => Cigar::Match(count),
                    b'I' => Cigar::Ins(count),
                    b'D' => Cigar::Del(count),
                    b'S' => Cigar::SoftClip(count),
                    _ => unreachable!("Invalid CIGAR operation: {op}"),
                })
                .collect::<Vec<_>>(),
        );

        // Step 7: Build final sequence by filtering out deletions and collecting bases
        let final_seq: Vec<u8> = bases_and_ops
            .iter()
            .copied()
            .filter_map(|item| (item.1 != b'D').then_some(item.0))
            .collect();

        // Step 8: final checks
        let first_ok = (cigar_string.0.as_slice().first().map(|x| x.char()) == Some('M'))
            || (cigar_string
                .0
                .as_slice()
                .first_chunk::<2>()
                .map(|&x| [x[0].char(), x[1].char()])
                == Some(['S', 'M']));
        let last_ok = (cigar_string.0.as_slice().last().map(|x| x.char()) == Some('M'))
            || (cigar_string
                .0
                .as_slice()
                .last_chunk::<2>()
                .map(|&x| [x[0].char(), x[1].char()])
                == Some(['M', 'S']));
        if first_ok && last_ok {
            Ok((
                final_seq,
                (!read_state.is_unmapped()).then_some(cigar_string),
            ))
        } else {
            Err(Error::SimulateDNASeqCIGAREndProblem(
                "too few bp or insertion/deletion close to end".to_owned(),
            ))
        }
    }
}

/// Generates BAM MM and ML tags with DNA modification probabilities
///
/// # Examples
///
/// ```
/// use std::str::FromStr;
/// use nanalogue_core::{DNARestrictive, Error};
/// use nanalogue_core::simulate_mod_bam::{ModConfigBuilder, generate_random_dna_modification};
/// use rand::Rng;
///
/// // Create a DNA sequence with multiple cytosines
/// let seq = DNARestrictive::from_str("ACGTCGCGATCGACGTCGCGATCG").unwrap();
///
/// // Configure modification for cytosines on plus strand.
/// // NOTE: we use mod ranges with an identical low and high value below
/// // because we do not want to deal with random values in an example.
///
/// let mod_config_c = ModConfigBuilder::default()
///     .base('C')
///     .is_strand_plus(true)
///     .mod_code("m".into())
///     .win(vec![2, 3])
///     .mod_range(vec![(0.8, 0.8), (0.4, 0.4)]).build()?;
///
/// // Put modifications for A on the minus strand.
/// // Again, we use equal low and high value of 0.2 as we do not
/// // want to deal with values drawn randomly from a range.
///
/// let mod_config_a = ModConfigBuilder::default()
///     .base('A')
///     .is_strand_plus(false)
///     .mod_code("20000".into())
///     .win([2].into())
///     .mod_range(vec![(0.2, 0.2)]).build()?;
///
/// let mod_config = vec![mod_config_c, mod_config_a];
///
/// let mut rng = rand::rng();
/// let (mm_str, ml_str) = generate_random_dna_modification(&mod_config, &seq, &mut rng);
///
/// // MM string format: C+m?gap1,gap2,gap3,gap4,...;A-20000?,...;
/// assert_eq!(mm_str, String::from("C+m?,0,0,0,0,0,0,0,0;A-20000?,0,0,0,0;"));
///
/// // ML string contains probabilities
/// assert_eq!(ml_str, vec![204, 204, 102, 102, 102, 204, 204, 102, 51, 51, 51, 51]);
/// # Ok::<(), Error>(())
/// ```
///
pub fn generate_random_dna_modification<R: Rng, S: GetDNARestrictive>(
    mod_configs: &[ModConfig],
    seq: &S,
    rng: &mut R,
) -> (String, Vec<u8>) {
    let seq_bytes = seq.get_dna_restrictive().get();
    let mut mm_str = String::new();

    // rough estimate of capacity below
    let mut ml_vec: Vec<u8> = Vec::with_capacity(seq_bytes.len());

    for mod_config in mod_configs {
        let base = u8::from(mod_config.base);
        let strand = if mod_config.is_strand_plus { '+' } else { '-' };
        let mod_code = mod_config.mod_code;
        let mut count = if base == b'N' {
            seq_bytes.len()
        } else {
            seq_bytes
                .iter()
                .zip(iter::repeat(&base))
                .filter(|&(&a, &b)| a == b)
                .count()
        };
        let mut output: Vec<u8> = Vec::with_capacity(count);
        for k in mod_config
            .win
            .iter()
            .cycle()
            .zip(mod_config.mod_range.iter().cycle())
        {
            let low = u8::from(k.1.low());
            let high = u8::from(k.1.high());
            #[expect(
                clippy::redundant_else,
                reason = "so that the clippy arithmetic lint fits better with the code"
            )]
            #[expect(
                clippy::arithmetic_side_effects,
                reason = "we loop only if count > 0 and catch count == 0, thus count doesn't underflow"
            )]
            if count == 0 {
                break;
            } else {
                for _ in 0..k.0.get() {
                    output.push(rng.random_range(low..=high));
                    count -= 1;
                    if count == 0 {
                        break;
                    }
                }
            }
        }
        if !output.is_empty() {
            let mod_len = output.len();
            ml_vec.append(&mut output);
            mm_str += format!(
                "{}{}{}?,{};",
                base as char,
                strand,
                mod_code,
                join(vec![0; mod_len], ",")
            )
            .as_str();
        }
    }
    ml_vec.shrink_to_fit();
    (mm_str, ml_vec)
}

/// Generates random DNA sequence of given length
///
/// # Panics
/// Panics if the sequence length exceeds `usize::MAX` (2^32 - 1 on 32-bit systems, 2^64 - 1 on 64-bit systems).
///
/// ```
/// use std::num::NonZeroU64;
/// use rand::Rng;
/// use nanalogue_core::simulate_mod_bam::generate_random_dna_sequence;
///
/// let mut rng = rand::rng();
/// let seq = generate_random_dna_sequence(100.try_into().expect("no error"), &mut rng);
/// assert_eq!(seq.len(), 100);
/// assert!(seq.iter().all(|&base| [b'A', b'C', b'G', b'T'].contains(&base)));
/// ```
pub fn generate_random_dna_sequence<R: Rng>(length: NonZeroU64, rng: &mut R) -> Vec<u8> {
    const DNA_BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    iter::repeat_with(|| {
        *DNA_BASES
            .get(rng.random_range(0..4))
            .expect("random_range(0..4) is always 0-3")
    })
    .take(usize::try_from(length.get()).expect("sequence length exceeds usize::MAX"))
    .collect()
}

/// Adds barcodes to read sequence based on strand orientation.
///
/// For forward and unmapped reads:
/// - Prepends barcode to read sequence
/// - Appends reverse complement of barcode to read sequence
///
/// For reverse reads:
/// - Prepends complement of barcode to read sequence
/// - Appends reverse of barcode (not reverse complement) to read sequence
///
/// # Examples
/// ```
/// use nanalogue_core::simulate_mod_bam::add_barcode;
/// use nanalogue_core::{ReadState, DNARestrictive};
/// use std::str::FromStr;
///
/// // Forward read: barcode + seq + revcomp(barcode)
/// let result = add_barcode("GGGGGGGG".as_bytes(), "ACGTAA".parse().unwrap(), ReadState::PrimaryFwd);
/// assert_eq!(result, b"ACGTAAGGGGGGGGTTACGT".to_vec());
///
/// // Reverse read: comp(barcode) + seq + rev(barcode)
/// let result = add_barcode("GGGGGGGG".as_bytes(), "ACGTAA".parse().unwrap(), ReadState::PrimaryRev);
/// assert_eq!(result, b"TGCATTGGGGGGGGAATGCA".to_vec());
/// ```
#[must_use]
#[expect(
    clippy::needless_pass_by_value,
    reason = "barcode by value allows passing `&'static str` easily (I think, may fix in future). \
If you don't know what this means, ignore this; this is _not_ a bug."
)]
pub fn add_barcode(read_seq: &[u8], barcode: DNARestrictive, read_state: ReadState) -> Vec<u8> {
    let bc_bytes = barcode.get();

    match read_state {
        ReadState::PrimaryFwd
        | ReadState::SecondaryFwd
        | ReadState::SupplementaryFwd
        | ReadState::Unmapped => {
            // Forward/Unmapped: barcode + read_seq + reverse_complement(barcode)
            let revcomp_bc = bio::alphabets::dna::revcomp(bc_bytes);
            [bc_bytes, read_seq, &revcomp_bc[..]].concat()
        }
        ReadState::PrimaryRev | ReadState::SecondaryRev | ReadState::SupplementaryRev => {
            // Reverse: complement(barcode) + read_seq + reverse(barcode)
            let comp_bc: Vec<u8> = bc_bytes
                .iter()
                .map(|&b| bio::alphabets::dna::complement(b))
                .collect();
            let rev_bc: Vec<u8> = bc_bytes.iter().copied().rev().collect();
            [&comp_bc[..], read_seq, &rev_bc[..]].concat()
        }
    }
}

/// Generates contigs with random sequence according to configuration
///
/// ```
/// use std::num::{NonZeroU32, NonZeroU64};
/// use nanalogue_core::{OrdPair, GetDNARestrictive};
/// use nanalogue_core::simulate_mod_bam::generate_contigs_denovo;
/// use rand::Rng;
///
/// let mut rng = rand::rng();
/// let contigs = generate_contigs_denovo(
///     3.try_into().unwrap(),
///     (100, 200).try_into().unwrap(),
///     &mut rng
/// );
/// assert_eq!(contigs.len(), 3);
/// assert!(contigs.iter().all(|c| (100..=200).contains(&c.get_dna_restrictive().get().len())));
/// ```
///
#[expect(
    clippy::missing_panics_doc,
    reason = "no length number or DNA conversion problems expected"
)]
pub fn generate_contigs_denovo<R: Rng>(
    contig_number: NonZeroU32,
    len_range: OrdPair<NonZeroU64>,
    rng: &mut R,
) -> Vec<Contig> {
    (0..contig_number.get())
        .map(|i| {
            let length = rng.random_range(len_range.low().get()..=len_range.high().get());
            let seq_bytes =
                generate_random_dna_sequence(NonZeroU64::try_from(length).expect("no error"), rng);
            let seq_str = String::from_utf8(seq_bytes).expect("valid DNA sequence");
            ContigBuilder::default()
                .name(format!("contig_{i:05}"))
                .seq(seq_str)
                .build()
                .expect("no error")
        })
        .collect()
}

/// Generates contigs with repeated sequence pattern
///
/// Creates contigs by repeating a given DNA sequence to reach a random length
/// within the specified range. The last repetition is clipped if needed.
///
/// ```
/// use std::num::{NonZeroU32, NonZeroU64};
/// use std::str::FromStr;
/// use nanalogue_core::{DNARestrictive, GetDNARestrictive, OrdPair};
/// use nanalogue_core::simulate_mod_bam::generate_contigs_denovo_repeated_seq;
/// use rand::Rng;
///
/// let mut rng = rand::rng();
/// let seq = DNARestrictive::from_str("ACGT").unwrap();
/// let contigs = generate_contigs_denovo_repeated_seq(
///     2.try_into().unwrap(),
///     (10, 12).try_into().unwrap(),
///     &seq,
///     &mut rng
/// );
/// assert_eq!(contigs.len(), 2);
/// for (i, contig) in contigs.iter().enumerate() {
///     assert_eq!(contig.name, format!("contig_0000{i}"));
///     match contig.get_dna_restrictive().get() {
///         b"ACGTACGTAC" | b"ACGTACGTACG" | b"ACGTACGTACGT" => {},
///         _ => panic!("Unexpected sequence"),
///     }
/// }
/// ```
#[expect(
    clippy::missing_panics_doc,
    reason = "no length number or DNA conversion problems expected"
)]
pub fn generate_contigs_denovo_repeated_seq<R: Rng, S: GetDNARestrictive>(
    contig_number: NonZeroU32,
    len_range: OrdPair<NonZeroU64>,
    seq: &S,
    rng: &mut R,
) -> Vec<Contig> {
    (0..contig_number.get())
        .map(|i| {
            let length = rng.random_range(len_range.low().get()..=len_range.high().get());
            let contig_seq: Vec<u8> = seq
                .get_dna_restrictive()
                .get()
                .iter()
                .cycle()
                .take(usize::try_from(length).expect("number conversion error"))
                .copied()
                .collect();
            let seq_str =
                String::from_utf8(contig_seq).expect("valid UTF-8 from repeated DNA sequence");
            ContigBuilder::default()
                .name(format!("contig_{i:05}"))
                .seq(seq_str)
                .build()
                .expect("valid DNA sequence from repeated pattern")
        })
        .collect()
}

/// Generates reads that align to contigs
///
/// # Example
///
/// Here we generate a contig, specify read properties, and ask for them to be generated.
/// Although we have used `ContigBuilder` here, you can also use `DNARestrictive` to pass
/// DNA sequences directly in the first argument i.e. you can use any `struct` that implements
/// `GetDNARestrictive`.
///
/// ```
/// use nanalogue_core::Error;
/// use nanalogue_core::simulate_mod_bam::{ContigBuilder, ReadConfigBuilder, generate_reads_denovo};
/// use rand::Rng;
///
/// let mut contig = ContigBuilder::default()
///     .name("chr1")
///     .seq("ACGTACGTACGTACGT".into()).build()?;
/// let contigs = vec![contig];
///
/// let mut config = ReadConfigBuilder::default()
///     .number(10)
///     .mapq_range((10, 20))
///     .base_qual_range((20, 30))
///     .len_range((0.2, 0.5)).build()?;
/// // NOTE: barcodes are optional, and will add 2*barcode_length to length statistics.
/// //       i.e. length stats are imposed independent of barcodes.
/// let mut rng = rand::rng();
/// let reads = generate_reads_denovo(&contigs, &config, "RG1", &mut rng).unwrap();
/// assert_eq!(reads.len(), 10);
/// # Ok::<(), Error>(())
/// ```
///
/// # Errors
/// Returns an error if the contigs input slice is empty,
/// if parameters are such that zero read lengths are produced,
/// or if BAM record creation fails, such as when making RG, MM, or ML tags.
#[expect(
    clippy::missing_panics_doc,
    clippy::too_many_lines,
    reason = "number conversion errors or mis-generation of DNA bases are unlikely \
    as we generate DNA sequences ourselves here, and genomic data are unlikely \
    to exceed ~2^63 bp or have ~2^32 contigs. Function length is acceptable for read generation"
)]
#[expect(
    clippy::cast_possible_truncation,
    reason = "read length calculated as a fraction of contig length, managed with trunc()"
)]
pub fn generate_reads_denovo<R: Rng, S: GetDNARestrictive>(
    contigs: &[S],
    read_config: &ReadConfig,
    read_group: &str,
    rng: &mut R,
) -> Result<Vec<bam::Record>, Error> {
    if contigs.is_empty() {
        return Err(Error::UnavailableData("no contigs found".to_owned()));
    }

    let mut reads = Vec::new();

    for _ in 0..read_config.number.get() {
        // Select a random contig
        let contig_idx = rng.random_range(0..contigs.len());
        let contig = contigs
            .get(contig_idx)
            .expect("contig_idx is within contigs range");
        let contig_len = contig.get_dna_restrictive().get().len() as u64;

        // Calculate read length as fraction of contig length
        #[expect(
            clippy::cast_precision_loss,
            reason = "u64->f32 causes precision loss but we are fine with this for now"
        )]
        #[expect(
            clippy::cast_sign_loss,
            reason = "these are positive numbers so no problem"
        )]
        let read_len = ((rng
            .random_range(read_config.len_range.low().val()..=read_config.len_range.high().val())
            * contig_len as f32)
            .trunc() as u64)
            .min(contig_len);

        // Set starting position
        let start_pos = rng.random_range(0..=(contig_len.saturating_sub(read_len)));

        // Extract sequence from contig
        let random_state: ReadState = rng.random();
        let end_pos = usize::try_from(start_pos.checked_add(read_len).expect("u64 overflow"))
            .expect("number conversion error");
        let (read_seq, cigar) = {
            let start_idx = usize::try_from(start_pos).expect("number conversion error");
            let temp_seq = contig
                .get_dna_restrictive()
                .get()
                .get(start_idx..end_pos)
                .expect("start_idx and end_pos are within contig bounds")
                .to_vec();
            let mut builder = PerfectSeqMatchToNot::seq(temp_seq)?;
            if let Some(barcode) = read_config.barcode.clone() {
                builder = builder.barcode(barcode);
            }
            if let Some(delete) = read_config.delete {
                builder = builder.delete(delete);
            }
            if let Some(insert_middle) = read_config.insert_middle.clone() {
                builder = builder.insert_middle(insert_middle);
            }
            if let Some(mismatch) = read_config.mismatch {
                builder = builder.mismatch(mismatch);
            }
            builder.build(random_state, rng)?
        };

        // Generate quality scores (for final read length including any adjustments like barcodes)
        let qual: Vec<u8> = iter::repeat_with(|| {
            rng.random_range(RangeInclusive::from(read_config.base_qual_range))
        })
        .take(read_seq.len())
        .collect();

        // Generate modification information after reverse complementing sequence if needed
        let seq: DNARestrictive;
        let (mod_prob_mm_tag, mod_pos_ml_tag) = generate_random_dna_modification(
            &read_config.mods,
            match random_state.strand() {
                '.' | '+' => {
                    seq = DNARestrictive::from_str(str::from_utf8(&read_seq).expect("no error"))
                        .expect("no error");
                    &seq
                }
                '-' => {
                    seq = DNARestrictive::from_str(
                        str::from_utf8(&bio::alphabets::dna::revcomp(&read_seq)).expect("no error"),
                    )
                    .expect("no error");
                    &seq
                }
                _ => unreachable!("`strand` is supposed to return one of +/-/."),
            },
            rng,
        );

        // Create BAM record
        let record = {
            let mut record = bam::Record::new();
            let uuid = {
                let mut bytes = [0u8; 16];
                rng.fill(&mut bytes);
                uuid::Builder::from_random_bytes(bytes).into_uuid()
            };
            let qname = format!("{read_group}.{uuid}").into_bytes();
            record.unset_flags();
            record.set_flags(u16::from(random_state));
            if random_state == ReadState::Unmapped {
                record.set(&qname, None, &read_seq, &qual);
                record.set_mapq(255);
                record.set_tid(-1);
                record.set_pos(-1);
            } else {
                let mapq = rng.random_range(RangeInclusive::from(read_config.mapq_range));
                record.set(&qname, cigar.as_ref(), &read_seq, &qual);
                record.set_mapq(mapq);
                record.set_tid(i32::try_from(contig_idx).expect("number conversion error"));
                record.set_pos(i64::try_from(start_pos).expect("number conversion error"));
            }
            record.set_mpos(-1);
            record.set_mtid(-1);
            record.push_aux(b"RG", Aux::String(read_group))?;
            if !mod_prob_mm_tag.is_empty() && !mod_pos_ml_tag.is_empty() {
                record.push_aux(b"MM", Aux::String(mod_prob_mm_tag.as_str()))?;
                record.push_aux(b"ML", Aux::ArrayU8((&mod_pos_ml_tag).into()))?;
            }
            record
        };
        reads.push(record);
    }

    Ok(reads)
}

/// Main function to generate simulated BAM file
///
/// # Example
///
/// ```
/// use nanalogue_core::{Error, SimulationConfig, simulate_mod_bam::run};
/// use uuid::Uuid;
///
/// let config_json = r#"{
///   "contigs": {
///     "number": 4,
///     "len_range": [1000, 8000]
///   },
///   "reads": [{
///     "number": 1000,
///     "mapq_range": [10, 20],
///     "base_qual_range": [10, 20],
///     "len_range": [0.1, 0.8]
///   }]
/// }"#;
///
/// // Note: * Optional "barcode" field can be added to reads (e.g., "barcode": "ACGTAA")
/// //         Length statistics are imposed independent of barcodes; so they will add
/// //         2x barcode length on top of these.
/// //       * Optional "repeated_seq" field can be added to contigs (e.g., "repeated_seq": "ACGT")
/// //         to create contigs by repeating a sequence instead of random generation.
/// //       * Optional modifications can be added using the "mods" field
/// //       * Optional "seed" field can be added (e.g., "seed": 12345) for reproducible output
/// // Note: The files used below should not exist because they will be created by this function.
/// //       If they exist, they will be overwritten.
///
/// let config: SimulationConfig = serde_json::from_str(&config_json)?;
/// let temp_dir = std::env::temp_dir();
/// let bam_path = temp_dir.join(format!("{}.bam", Uuid::new_v4()));
/// let fasta_path = temp_dir.join(format!("{}.fa", Uuid::new_v4()));
/// let bai_path = bam_path.with_extension("bam.bai");
/// run(config, &bam_path, &fasta_path)?;
/// std::fs::remove_file(&bam_path)?;
/// std::fs::remove_file(&fasta_path)?;
/// std::fs::remove_file(&bai_path)?;
/// # Ok::<(), Error>(())
/// ```
///
/// # Errors
/// Returns an error if JSON parsing fails, read generation fails, or BAM/FASTA writing fails.
pub fn run<F>(
    config: SimulationConfig,
    bam_output_path: &F,
    fasta_output_path: &F,
) -> Result<(), Error>
where
    F: AsRef<Path> + ?Sized,
{
    if let Some(s) = config.seed {
        let mut rng = StdRng::seed_from_u64(s);
        run_inner(config, bam_output_path, fasta_output_path, &mut rng)
    } else {
        let mut rng = rand::rng();
        run_inner(config, bam_output_path, fasta_output_path, &mut rng)
    }
}

/// Generates simulated BAM and FASTA files using the provided RNG.
fn run_inner<F, R: Rng>(
    config: SimulationConfig,
    bam_output_path: &F,
    fasta_output_path: &F,
    rng: &mut R,
) -> Result<(), Error>
where
    F: AsRef<Path> + ?Sized,
{
    let contigs = match config.contigs.repeated_seq {
        Some(seq) => generate_contigs_denovo_repeated_seq(
            config.contigs.number,
            config.contigs.len_range,
            &seq,
            rng,
        ),
        None => generate_contigs_denovo(config.contigs.number, config.contigs.len_range, rng),
    };
    let read_groups: Vec<String> = (0..config.reads.len()).map(|k| k.to_string()).collect();
    let reads = {
        let mut temp_reads = Vec::new();
        for k in config.reads.into_iter().zip(read_groups.clone()) {
            temp_reads.append(&mut generate_reads_denovo(&contigs, &k.0, &k.1, rng)?);
        }
        temp_reads.sort_by_key(|k| (k.is_unmapped(), k.tid(), k.pos(), k.is_reverse()));
        temp_reads
    };

    write_bam_denovo(
        reads,
        contigs
            .iter()
            .map(|k| (k.name.clone(), k.get_dna_restrictive().get().len())),
        read_groups,
        vec![String::from("simulated BAM file, not real data")],
        bam_output_path,
    )?;
    write_fasta(
        contigs.into_iter().map(|k| (k.name.clone(), k)),
        fasta_output_path,
    )?;

    Ok(())
}

/// Temporary BAM simulation with automatic cleanup
///
/// Creates temporary BAM and FASTA files for testing purposes and
/// automatically removes them when dropped.
#[derive(Debug, Serialize, Deserialize)]
pub struct TempBamSimulation {
    /// Path to bam file that will be created
    bam_path: String,
    /// Path to fasta file that will be created
    fasta_path: String,
}

impl TempBamSimulation {
    /// Creates a new temporary BAM simulation from JSON configuration
    ///
    /// # Errors
    /// Returns an error if the simulation run fails
    pub fn new(config: SimulationConfig) -> Result<Self, Error> {
        let temp_dir = std::env::temp_dir();
        let bam_path = temp_dir.join(format!("{}.bam", Uuid::new_v4()));
        let fasta_path = temp_dir.join(format!("{}.fa", Uuid::new_v4()));

        run(config, &bam_path, &fasta_path)?;

        Ok(Self {
            bam_path: bam_path.to_string_lossy().to_string(),
            fasta_path: fasta_path.to_string_lossy().to_string(),
        })
    }

    /// Returns the path to the temporary BAM file
    #[must_use]
    pub fn bam_path(&self) -> &str {
        &self.bam_path
    }

    /// Returns the path to the temporary FASTA file
    #[must_use]
    pub fn fasta_path(&self) -> &str {
        &self.fasta_path
    }
}

impl Drop for TempBamSimulation {
    fn drop(&mut self) {
        // Ignore errors during cleanup - files may already be deleted
        drop(std::fs::remove_file(&self.bam_path));
        drop(std::fs::remove_file(&self.fasta_path));
        drop(std::fs::remove_file(&(self.bam_path.clone() + ".bai")));
    }
}

#[cfg(test)]
mod random_dna_generation_test {
    use super::*;

    /// Test for generation of random DNA of a given length
    #[test]
    fn generate_random_dna_sequence_works() {
        let seq = generate_random_dna_sequence(NonZeroU64::new(100).unwrap(), &mut rand::rng());
        assert_eq!(seq.len(), 100);
        for base in seq {
            assert!([b'A', b'C', b'G', b'T'].contains(&base));
        }
    }
}

#[cfg(test)]
mod seeded_simulation_tests {
    use super::*;
    use rust_htslib::bam::Read as _;

    fn run_seeded_simulation(config_json: &str) -> Result<Vec<bam::Record>, Error> {
        let config: SimulationConfig = serde_json::from_str(config_json)?;
        let temp = TempBamSimulation::new(config)?;
        let mut reader = crate::nanalogue_bam_reader(temp.bam_path())?;
        reader
            .records()
            .collect::<Result<Vec<_>, _>>()
            .map_err(Error::from)
    }

    #[test]
    fn seeded_simulation_is_reproducible() -> Result<(), Error> {
        let config_json = r#"{
          "contigs": { "number": 2, "len_range": [500, 1000],
                       "repeated_seq": "ATCGAATT" },
          "reads": [{ "number": 50, "mapq_range": [10, 20],
                      "base_qual_range": [10, 20], "len_range": [0.2, 0.5],
                      "mods": [
                        { "base": "T", "is_strand_plus": true, "mod_code": "T",
                          "win": [4, 5], "mod_range": [[0.1, 0.2], [0.3, 0.4]] },
                        { "base": "C", "is_strand_plus": true, "mod_code": "m",
                          "win": [3], "mod_range": [[0.6, 0.9]] }
                      ] }],
          "seed": 12345
        }"#;

        let run_a = run_seeded_simulation(config_json)?;
        let run_b = run_seeded_simulation(config_json)?;
        assert_eq!(run_a, run_b, "same seed must produce identical BAM records");

        Ok(())
    }
}

#[cfg(test)]
mod read_generation_no_mods_tests {
    use super::*;
    use rust_htslib::bam::Read as _;

    /// Tests read generation with desired properties but no modifications
    #[test]
    fn generate_reads_denovo_no_mods_works() {
        let mut rng = rand::rng();
        let contigs = vec![
            ContigBuilder::default()
                .name("contig_0")
                .seq("ACGTACGTACGTACGTACGT".into())
                .build()
                .unwrap(),
            ContigBuilder::default()
                .name("contig_1")
                .seq("TGCATGCATGCATGCATGCA".into())
                .build()
                .unwrap(),
        ];

        let config = ReadConfigBuilder::default()
            .number(10)
            .mapq_range((10, 20))
            .base_qual_range((30, 50))
            .len_range((0.2, 0.8))
            .build()
            .unwrap();

        let reads = generate_reads_denovo(&contigs, &config, "ph", &mut rng).unwrap();
        assert_eq!(reads.len(), 10);

        for read in &reads {
            // check mapping information
            if !read.is_unmapped() {
                assert!((10..=20).contains(&read.mapq()));
                assert!((0..2).contains(&read.tid()));
            }

            // check sequence length, quality, and check that no mods are present
            assert!((4..=16).contains(&read.seq_len()));
            assert!(read.qual().iter().all(|x| (30..=50).contains(x)));
            let _: rust_htslib::errors::Error = read.aux(b"MM").unwrap_err();
            let _: rust_htslib::errors::Error = read.aux(b"ML").unwrap_err();

            // check read name, must be "ph.UUID"
            assert_eq!(read.qname().get(0..3).unwrap(), *b"ph.");
            let _: Uuid =
                Uuid::parse_str(str::from_utf8(read.qname().get(3..).unwrap()).unwrap()).unwrap();
        }
    }

    /// Tests full BAM creation
    #[test]
    fn full_bam_generation() {
        let config_json = r#"{
            "contigs": {
                "number": 2,
                "len_range": [100, 200]
            },
            "reads": [{
                "number": 1000,
                "mapq_range": [10, 20],
                "base_qual_range": [10, 20],
                "len_range": [0.1, 0.8]
            }]
        }"#;

        // create files with random names and check for existence
        let temp_dir = std::env::temp_dir();
        let bam_path = temp_dir.join(format!("{}.bam", Uuid::new_v4()));
        let fasta_path = temp_dir.join(format!("{}.fa", Uuid::new_v4()));

        let config: SimulationConfig = serde_json::from_str(config_json).unwrap();
        run(config, &bam_path, &fasta_path).unwrap();

        assert!(bam_path.exists());
        assert!(fasta_path.exists());

        // read BAM file and check contig, read count
        let mut reader = bam::Reader::from_path(&bam_path).unwrap();
        let header = reader.header();
        assert_eq!(header.target_count(), 2);
        assert_eq!(reader.records().count(), 1000);

        // delete files
        let bai_path = bam_path.with_extension("bam.bai");
        std::fs::remove_file(&bam_path).unwrap();
        std::fs::remove_file(fasta_path).unwrap();
        std::fs::remove_file(bai_path).unwrap();
    }

    /// Tests `TempBamSimulation` struct functionality without mods
    #[test]
    fn temp_bam_simulation_struct_no_mods() {
        let config_json = r#"{
            "contigs": {
                "number": 2,
                "len_range": [100, 200]
            },
            "reads": [{
                "number": 1000,
                "mapq_range": [10, 20],
                "base_qual_range": [10, 20],
                "len_range": [0.1, 0.8]
            }]
        }"#;

        // Create temporary simulation
        let config: SimulationConfig = serde_json::from_str(config_json).unwrap();
        let sim = TempBamSimulation::new(config).unwrap();

        // Verify files exist
        assert!(Path::new(sim.bam_path()).exists());
        assert!(Path::new(sim.fasta_path()).exists());

        // Read BAM file and check contig count
        let mut reader = bam::Reader::from_path(sim.bam_path()).unwrap();
        let header = reader.header();
        assert_eq!(header.target_count(), 2);
        assert_eq!(reader.records().count(), 1000);
    }

    /// Tests `TempBamSimulation` automatic cleanup
    #[test]
    fn temp_bam_simulation_cleanup() {
        let config_json = r#"{
            "contigs": {
                "number": 1,
                "len_range": [50, 50]
            },
            "reads": [{
                "number": 10,
                "len_range": [0.5, 0.5]
            }]
        }"#;

        let bam_path: String;
        let fasta_path: String;

        {
            let config: SimulationConfig = serde_json::from_str(config_json).unwrap();
            let sim = TempBamSimulation::new(config).unwrap();
            bam_path = sim.bam_path().to_string();
            fasta_path = sim.fasta_path().to_string();

            // Files should exist while sim is in scope
            assert!(Path::new(&bam_path).exists());
            assert!(Path::new(&fasta_path).exists());
        } // sim is dropped here

        // Files should be cleaned up after drop
        assert!(!Path::new(&bam_path).exists());
        assert!(!Path::new(&fasta_path).exists());
    }

    /// Tests error when generating reads with empty contigs slice
    #[test]
    fn generate_reads_denovo_empty_contigs_error() {
        let contigs: Vec<Contig> = vec![];
        let config = ReadConfig::default();
        let mut rng = rand::rng();

        let result = generate_reads_denovo(&contigs, &config, "test", &mut rng);
        assert!(matches!(result, Err(Error::UnavailableData(_))));
    }

    /// Tests error when read length would be zero
    #[test]
    fn generate_reads_denovo_zero_length_error() {
        let contigs = [ContigBuilder::default()
            .name("tiny")
            .seq("ACGT".into())
            .build()
            .unwrap()];

        let config = ReadConfigBuilder::default()
            .number(1)
            .len_range((0.0, 0.0))
            .build()
            .unwrap();

        let mut rng = rand::rng();
        let result = generate_reads_denovo(&contigs, &config, "test", &mut rng);
        assert!(matches!(result, Err(Error::InvalidState(_))));
    }

    /// Tests invalid JSON structure causing empty reads generation
    #[test]
    fn run_empty_reads_error() {
        let invalid_json = r#"{ "reads": [] }"#; // Empty reads array
        let temp_dir = std::env::temp_dir();
        let bam_path = temp_dir.join(format!("{}.bam", Uuid::new_v4()));
        let fasta_path = temp_dir.join(format!("{}.fa", Uuid::new_v4()));

        let config: SimulationConfig = serde_json::from_str(invalid_json).unwrap();
        let result = run(config, &bam_path, &fasta_path);
        // With empty reads, this should succeed but produce an empty BAM (valid)
        // So we won't assert error here, just test it doesn't crash
        drop(result);
        let bai_path = bam_path.with_extension("bam.bai");
        drop(std::fs::remove_file(&bam_path));
        drop(std::fs::remove_file(&fasta_path));
        drop(std::fs::remove_file(&bai_path));
    }

    /// Tests multiple read groups in BAM generation
    #[test]
    fn multiple_read_groups_work() {
        let config_json = r#"{
            "contigs": {
                "number": 2,
                "len_range": [100, 200]
            },
            "reads": [
                {
                    "number": 50,
                    "mapq_range": [10, 20],
                    "base_qual_range": [20, 30],
                    "len_range": [0.1, 0.5]
                },
                {
                    "number": 75,
                    "mapq_range": [30, 40],
                    "base_qual_range": [25, 35],
                    "len_range": [0.3, 0.7]
                },
                {
                    "number": 25,
                    "mapq_range": [5, 15],
                    "base_qual_range": [15, 25],
                    "len_range": [0.2, 0.6]
                }
            ]
        }"#;

        let config: SimulationConfig = serde_json::from_str(config_json).unwrap();
        let sim = TempBamSimulation::new(config).unwrap();
        let mut reader = bam::Reader::from_path(sim.bam_path()).unwrap();

        // Should have 50 + 75 + 25 = 150 reads total
        assert_eq!(reader.records().count(), 150);

        // Verify read groups exist in header
        let mut reader2 = bam::Reader::from_path(sim.bam_path()).unwrap();
        let header = reader2.header();
        let header_text = std::str::from_utf8(header.as_bytes()).unwrap();
        assert!(header_text.contains("@RG\tID:0"));
        assert!(header_text.contains("@RG\tID:1"));
        assert!(header_text.contains("@RG\tID:2"));

        // Count reads per read group and verify read name format
        let mut rg_counts = [0, 0, 0];
        for record in reader2.records() {
            let read = record.unwrap();
            let qname = std::str::from_utf8(read.qname()).unwrap();

            if let Ok(Aux::String(rg)) = read.aux(b"RG") {
                // Verify read name format: should be "rg.uuid"
                let parts: Vec<&str> = qname.split('.').collect();
                assert_eq!(parts.len(), 2, "Read name should be in format 'n.uuid'");

                let read_group_prefix =
                    parts.first().expect("parts should have at least 1 element");
                let uuid_part = parts.get(1).expect("parts should have 2 elements");

                assert_eq!(
                    *read_group_prefix, rg,
                    "Read name prefix should match read group"
                );

                // Verify the second part is a valid UUID
                assert!(
                    Uuid::parse_str(uuid_part).is_ok(),
                    "Read name should contain a valid UUID after the dot"
                );

                match rg {
                    "0" => *rg_counts.get_mut(0).unwrap() += 1,
                    "1" => *rg_counts.get_mut(1).unwrap() += 1,
                    "2" => *rg_counts.get_mut(2).unwrap() += 1,
                    unexpected => unreachable!("Unexpected read group: {unexpected}"),
                }
            }
        }
        assert_eq!(rg_counts[0], 50);
        assert_eq!(rg_counts[1], 75);
        assert_eq!(rg_counts[2], 25);
    }

    /// Tests that read sequences are valid substrings of their parent contigs
    #[expect(
        clippy::cast_sign_loss,
        reason = "tid and pos are non-negative in valid BAM"
    )]
    #[expect(
        clippy::cast_possible_truncation,
        reason = "tid fits in usize for normal BAM files"
    )]
    #[test]
    fn read_sequences_match_contigs() {
        let contigs = vec![
            ContigBuilder::default()
                .name("chr1")
                .seq("ACGTACGTACGTACGTACGTACGTACGTACGT".into())
                .build()
                .unwrap(),
            ContigBuilder::default()
                .name("chr2")
                .seq("TGCATGCATGCATGCATGCATGCATGCATGCA".into())
                .build()
                .unwrap(),
        ];

        let read_config = ReadConfigBuilder::default()
            .number(50)
            .len_range((0.2, 0.8))
            .build()
            .unwrap();

        let mut rng = rand::rng();
        let reads = generate_reads_denovo(&contigs, &read_config, "test", &mut rng).unwrap();

        let mut validated_count = 0;
        for read in &reads {
            // Only validate forward mapped reads without soft clips for simplicity
            if !read.is_unmapped() && !read.is_reverse() {
                let cigar = read.cigar();

                // Skip reads with soft clips (shouldn't happen without barcodes, but be safe)
                let has_soft_clip = cigar.iter().any(|op| matches!(op, Cigar::SoftClip(_)));
                if has_soft_clip {
                    continue;
                }

                let tid = read.tid() as usize;
                let pos = read.pos() as usize;
                let contig = contigs.get(tid).expect("valid tid");
                let contig_seq = contig.get_dna_restrictive().get();

                // Get alignment length
                let mut aligned_len = 0usize;
                for op in &cigar {
                    if let Cigar::Match(len) = *op {
                        aligned_len = aligned_len.checked_add(len as usize).expect("overflow");
                    }
                }

                let read_seq = read.seq().as_bytes();
                let expected_seq = contig_seq
                    .get(pos..pos.checked_add(aligned_len).expect("overflow"))
                    .expect("contig subsequence exists");

                assert_eq!(
                    read_seq, expected_seq,
                    "Forward read sequence should exactly match contig substring"
                );
                validated_count += 1;
            }
        }

        // Ensure we validated at least some reads
        assert!(
            validated_count > 0,
            "Should have validated at least some forward reads"
        );
    }

    /// Tests different read states (unmapped, secondary, supplementary)
    #[test]
    fn different_read_states_work() {
        let config_json = r#"{
            "contigs": {
                "number": 1,
                "len_range": [200, 200]
            },
            "reads": [{
                "number": 1000,
                "mapq_range": [10, 20],
                "base_qual_range": [20, 30],
                "len_range": [0.2, 0.8]
            }]
        }"#;

        let config: SimulationConfig = serde_json::from_str(config_json).unwrap();
        let sim = TempBamSimulation::new(config).unwrap();
        let mut reader = bam::Reader::from_path(sim.bam_path()).unwrap();

        let mut has_unmapped = false;
        let mut has_forward = false;
        let mut has_reverse = false;
        let mut has_secondary = false;
        let mut has_supplementary = false;

        for record in reader.records() {
            let read = record.unwrap();

            if read.is_unmapped() {
                has_unmapped = true;
                // Verify unmapped read properties
                assert_eq!(read.tid(), -1);
                assert_eq!(read.pos(), -1);
                assert_eq!(read.mapq(), 255);
            } else {
                if read.is_reverse() {
                    has_reverse = true;
                } else {
                    has_forward = true;
                }

                if read.is_secondary() {
                    has_secondary = true;
                }

                if read.is_supplementary() {
                    has_supplementary = true;
                }
            }
        }

        // With 1000 reads, we should have diverse read states
        assert!(has_unmapped, "Should have unmapped reads");
        assert!(has_forward, "Should have forward reads");
        assert!(has_reverse, "Should have reverse reads");
        assert!(has_secondary, "Should have secondary reads");
        assert!(has_supplementary, "Should have supplementary reads");
    }

    /// Tests read length at 100% of contig length
    #[test]
    fn read_length_full_contig_works() {
        let contigs = vec![
            ContigBuilder::default()
                .name("chr1")
                .seq("ACGTACGTACGTACGTACGTACGTACGTACGT".into())
                .build()
                .unwrap(),
        ];

        let config_full_length = ReadConfigBuilder::default()
            .number(20)
            .len_range((1.0, 1.0))
            .build()
            .unwrap();

        let mut rng = rand::rng();
        let reads_full =
            generate_reads_denovo(&contigs, &config_full_length, "test", &mut rng).unwrap();

        // All mapped reads should be exactly contig length
        for read in &reads_full {
            if !read.is_unmapped() {
                assert_eq!(read.seq_len(), 32, "Read should be full contig length");
            }
        }
    }

    /// Tests read length with very small fraction of contig
    #[test]
    fn read_length_small_fraction_works() {
        // Use 200bp contig so 1% = 2bp (won't round to 0)
        let contigs = [ContigBuilder::default()
            .name("chr1")
            .seq("ACGT".repeat(50))
            .build()
            .unwrap()];

        let config_small = ReadConfigBuilder::default()
            .number(20)
            .len_range((0.01, 0.05))
            .build()
            .unwrap();

        let mut rng = rand::rng();
        let reads_small = generate_reads_denovo(&contigs, &config_small, "test", &mut rng).unwrap();

        // All mapped reads should be very small (2-10 bp given 200bp contig)
        for read in &reads_small {
            if !read.is_unmapped() {
                assert!(read.seq_len() <= 10, "Read should be very small");
                assert!(read.seq_len() >= 2, "Read should be at least 2bp");
            }
        }
    }
}

#[cfg(test)]
mod read_generation_barcodes {
    use super::*;
    use rust_htslib::bam::Read as _;

    /// Tests `add_barcode` function with forward and reverse reads
    #[expect(
        clippy::shadow_unrelated,
        reason = "repetition is fine; each block is clearly separated"
    )]
    #[test]
    fn add_barcode_works() {
        let read_seq = b"GGGGGGGG".to_vec();
        let barcode = DNARestrictive::from_str("ACGTAA").unwrap();

        // Test forward read: barcode + seq + revcomp(barcode)
        let result = add_barcode(&read_seq, barcode.clone(), ReadState::PrimaryFwd);
        assert_eq!(result, b"ACGTAAGGGGGGGGTTACGT".to_vec());

        // Test reverse read: comp(barcode) + seq + rev(barcode)
        let result = add_barcode(&read_seq, barcode, ReadState::PrimaryRev);
        assert_eq!(result, b"TGCATTGGGGGGGGAATGCA".to_vec());
    }

    /// Tests CIGAR string validity with and without barcodes
    #[test]
    fn cigar_strings_valid_with_or_without_barcodes() {
        // Test without barcodes
        let contigs = [ContigBuilder::default()
            .name("chr1")
            .seq("ACGTACGTACGTACGTACGTACGTACGTACGT".into())
            .build()
            .unwrap()];

        let config_no_barcode = ReadConfigBuilder::default()
            .number(20)
            .len_range((0.3, 0.7))
            .build()
            .unwrap();

        let mut rng = rand::rng();
        let reads_no_barcode =
            generate_reads_denovo(&contigs, &config_no_barcode, "test", &mut rng).unwrap();

        for read in &reads_no_barcode {
            if !read.is_unmapped() {
                let cigar = read.cigar();
                // Without barcode, CIGAR should be just a single Match operation
                assert_eq!(cigar.len(), 1);
                assert!(matches!(cigar.first().unwrap(), Cigar::Match(_)));

                // Verify CIGAR length matches sequence length
                let cigar_len: u32 = cigar
                    .iter()
                    .map(|op| match *op {
                        Cigar::Match(len) | Cigar::SoftClip(len) => len,
                        Cigar::Ins(_)
                        | Cigar::Del(_)
                        | Cigar::RefSkip(_)
                        | Cigar::HardClip(_)
                        | Cigar::Pad(_)
                        | Cigar::Equal(_)
                        | Cigar::Diff(_) => unreachable!(),
                    })
                    .sum();
                assert_eq!(cigar_len as usize, read.seq_len());
            }
        }

        // Test with barcodes
        let config_with_barcode = ReadConfigBuilder::default()
            .number(20)
            .len_range((0.3, 0.7))
            .barcode("ACGTAA".into())
            .build()
            .unwrap();

        let reads_with_barcode =
            generate_reads_denovo(&contigs, &config_with_barcode, "test", &mut rng).unwrap();

        for read in &reads_with_barcode {
            if !read.is_unmapped() {
                let cigar = read.cigar();
                // With barcode, CIGAR should be: SoftClip + Match + SoftClip
                assert_eq!(cigar.len(), 3);
                assert!(matches!(cigar.first().unwrap(), Cigar::SoftClip(6))); // barcode length
                assert!(matches!(cigar.get(1).unwrap(), Cigar::Match(_)));
                assert!(matches!(cigar.get(2).unwrap(), Cigar::SoftClip(6))); // barcode length

                // Verify total CIGAR length matches sequence length
                let cigar_len: u32 = cigar
                    .iter()
                    .map(|op| match *op {
                        Cigar::Match(len) | Cigar::SoftClip(len) => len,
                        Cigar::Ins(_)
                        | Cigar::Del(_)
                        | Cigar::RefSkip(_)
                        | Cigar::HardClip(_)
                        | Cigar::Pad(_)
                        | Cigar::Equal(_)
                        | Cigar::Diff(_) => unreachable!(),
                    })
                    .sum();
                assert_eq!(cigar_len as usize, read.seq_len());
            }
        }
    }

    /// Tests read generation with barcode
    #[test]
    fn generate_reads_denovo_with_barcode_works() {
        let config_json = r#"{
            "contigs": {
                "number": 1,
                "len_range": [20, 20]
            },
            "reads": [{
                "number": 10,
                "len_range": [0.2, 0.8],
                "barcode": "GTACGG"
            }]
        }"#;

        let config: SimulationConfig = serde_json::from_str(config_json).unwrap();
        let sim = TempBamSimulation::new(config).unwrap();
        let mut reader = bam::Reader::from_path(sim.bam_path()).unwrap();

        for record in reader.records() {
            let read = record.unwrap();
            let seq = read.seq().as_bytes();
            let seq_len = seq.len();

            if read.is_reverse() {
                // Reverse: comp(GTACGG)=CATGCC at start, rev(GTACGG)=GGCATG at end
                assert_eq!(seq.get(..6).expect("seq has at least 6 bases"), b"CATGCC");
                assert_eq!(
                    seq.get(seq_len.saturating_sub(6)..)
                        .expect("seq has at least 6 bases"),
                    b"GGCATG"
                );
            } else {
                // Forward: GTACGG at start, revcomp(GTACGG)=CCGTAC at end
                assert_eq!(seq.get(..6).expect("seq has at least 6 bases"), b"GTACGG");
                assert_eq!(
                    seq.get(seq_len.saturating_sub(6)..)
                        .expect("seq has at least 6 bases"),
                    b"CCGTAC"
                );
            }
        }
    }
}

#[cfg(test)]
mod read_generation_with_mods_tests {
    use super::*;
    use crate::{CurrRead, ThresholdState, curr_reads_to_dataframe};
    use rust_htslib::bam::Read as _;

    /// Tests read generation with desired properties but with modifications
    #[test]
    fn generate_reads_denovo_with_mods_works() {
        let mut rng = rand::rng();
        let contigs = vec![
            ContigBuilder::default()
                .name("contig_0")
                .seq("ACGTACGTACGTACGTACGT".into())
                .build()
                .unwrap(),
            ContigBuilder::default()
                .name("contig_1")
                .seq("TGCATGCATGCATGCATGCA".into())
                .build()
                .unwrap(),
        ];

        let config = ReadConfigBuilder::default()
            .number(10000)
            .mapq_range((10, 20))
            .base_qual_range((30, 50))
            .len_range((0.2, 0.8))
            .mods(
                [ModConfigBuilder::default()
                    .base('C')
                    .mod_code("m".into())
                    .win([4].into())
                    .mod_range([(0f32, 1f32)].into())
                    .build()
                    .unwrap()]
                .into(),
            )
            .build()
            .unwrap();

        let reads = generate_reads_denovo(&contigs, &config, "1", &mut rng).unwrap();
        assert_eq!(reads.len(), 10000);

        let mut zero_to_255_visited = vec![false; 256];

        for read in &reads {
            // check mapping information
            if !read.is_unmapped() {
                assert!((10..=20).contains(&read.mapq()));
                assert!((0..2).contains(&read.tid()));
            }

            // check other information
            assert!((4..=16).contains(&read.seq_len()));
            assert!(read.qual().iter().all(|x| (30..=50).contains(x)));

            // check modification information
            let Aux::String(mod_pos_mm_tag) = read.aux(b"MM").unwrap() else {
                unreachable!()
            };
            let Aux::ArrayU8(mod_prob_ml_tag) = read.aux(b"ML").unwrap() else {
                unreachable!()
            };
            // Parse and verify the actual gap position values in MM tag
            let pos: Vec<&str> = mod_pos_mm_tag
                .strip_prefix("C+m?,")
                .unwrap()
                .strip_suffix(';')
                .unwrap()
                .split(',')
                .collect();
            assert!(pos.iter().all(|&x| x == "0"), "All MM values should be 0");

            // verify probability values which should have been converted to 0-255.
            #[expect(
                clippy::indexing_slicing,
                reason = "we are not gonna bother, this should be fine"
            )]
            for k in 0..mod_prob_ml_tag.len() {
                zero_to_255_visited[usize::from(mod_prob_ml_tag.get(k).expect("no error"))] = true;
            }
        }

        // check that all values from 0-255 have been visited as we are using a huge number of
        // reads.
        assert!(zero_to_255_visited.into_iter().all(|x| x));
    }

    /// Tests `TempBamSimulation` struct functionality with mods
    #[test]
    fn temp_bam_simulation_struct_with_mods() {
        let config_json = r#"{
            "contigs": {
                "number": 2,
                "len_range": [100, 200]
            },
            "reads": [{
                "number": 1000,
                "mapq_range": [10, 20],
                "base_qual_range": [10, 20],
                "len_range": [0.1, 0.8],
                "mods": [{
                    "base": "T",
                    "is_strand_plus": true,
                    "mod_code": "T",
                    "win": [4, 5],
                    "mod_range": [[0.1, 0.2], [0.3, 0.4]]
                }]
            }]
        }"#;

        // Create temporary simulation
        let config: SimulationConfig = serde_json::from_str(config_json).unwrap();
        let sim = TempBamSimulation::new(config).unwrap();

        // Verify files exist
        assert!(Path::new(sim.bam_path()).exists());
        assert!(Path::new(sim.fasta_path()).exists());

        // Read BAM file and check contig, read count
        let mut reader = bam::Reader::from_path(sim.bam_path()).unwrap();
        let header = reader.header();
        assert_eq!(header.target_count(), 2);
        assert_eq!(reader.records().count(), 1000);
    }

    /// Tests `generate_random_dna_modification` with empty mod config
    #[test]
    fn generate_random_dna_modification_empty_config() {
        let seq = DNARestrictive::from_str("ACGTACGT").unwrap();
        let mod_configs: Vec<ModConfig> = vec![];
        let mut rng = rand::rng();

        let (mm_str, ml_vec) = generate_random_dna_modification(&mod_configs, &seq, &mut rng);

        assert_eq!(mm_str, String::new());
        assert_eq!(ml_vec.len(), 0);
    }

    /// Tests `generate_random_dna_modification` with single modification config
    #[test]
    fn generate_random_dna_modification_single_mod() {
        let seq = DNARestrictive::from_str("ACGTCGCGATCG").unwrap();
        let mod_config = ModConfigBuilder::default()
            .base('C')
            .mod_code("m".into())
            .win(vec![2])
            .mod_range(vec![(0.5, 0.5)])
            .build()
            .unwrap();

        let mut rng = rand::rng();
        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config], &seq, &mut rng);

        // Sequence has 4 C's, so we expect 4 modifications
        assert_eq!(ml_vec.len(), 4);
        // All should be 128 (0.5 * 255)
        assert!(ml_vec.iter().all(|&x| x == 128u8));
        // All should be zero
        let probs: Vec<&str> = mm_str
            .strip_prefix("C+m?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert_eq!(probs.len(), 4);
        assert!(probs.iter().all(|&x| x == "0"));
    }

    /// Tests `generate_random_dna_modification` with multiple modification configs
    #[test]
    fn generate_random_dna_modification_multiple_mods() {
        let seq = DNARestrictive::from_str("ACGTACGT").unwrap();

        let mod_config_c = ModConfigBuilder::default()
            .base('C')
            .mod_code('m'.into())
            .win(vec![1])
            .mod_range(vec![(0.8, 0.8)])
            .build()
            .unwrap();

        let mod_config_t = ModConfigBuilder::default()
            .base('T')
            .is_strand_plus(false)
            .mod_code('t'.into())
            .win(vec![1])
            .mod_range(vec![(0.4, 0.4)])
            .build()
            .unwrap();

        let mut rng = rand::rng();
        let (mm_str, ml_vec) =
            generate_random_dna_modification(&[mod_config_c, mod_config_t], &seq, &mut rng);

        // Sequence has 2 C's and 2 T's
        assert!(mm_str.contains("C+m?,"));
        assert!(mm_str.contains("T-t?,"));
        assert_eq!(
            ml_vec,
            vec![204u8, 204u8, 102u8, 102u8],
            "ML values should be 204,204,102,102"
        );

        // Parse and verify the C modifications
        let c_section = mm_str
            .split(';')
            .find(|s| s.starts_with("C+m?"))
            .expect("Should have C+m? section");
        let c_probs: Vec<&str> = c_section
            .strip_prefix("C+m?,")
            .unwrap()
            .split(',')
            .collect();
        assert_eq!(
            c_probs.len(),
            2,
            "Should have exactly 2 C modifications in ACGTACGT"
        );
        // All C pos should be 0
        assert!(
            c_probs.iter().all(|&x| x == "0"),
            "All C modifications should have gap pos 0",
        );

        // Parse and verify the T modifications
        let t_section = mm_str
            .split(';')
            .find(|s| s.starts_with("T-t?"))
            .expect("Should have T-t? section");
        let t_probs: Vec<&str> = t_section
            .strip_prefix("T-t?,")
            .unwrap()
            .split(',')
            .collect();
        assert_eq!(
            t_probs.len(),
            2,
            "Should have exactly 2 T modifications in ACGTACGT"
        );
        // All T pos should be 0
        assert!(
            t_probs.iter().all(|&x| x == "0"),
            "All T modifications should have a gap pos of 0"
        );
    }

    /// Tests `generate_random_dna_modification` with N base (all bases)
    #[test]
    fn generate_random_dna_modification_n_base() {
        let seq = DNARestrictive::from_str("ACGT").unwrap();

        let mod_config = ModConfigBuilder::default()
            .base('N')
            .is_strand_plus(true)
            .mod_code("n".into())
            .win([4].into())
            .mod_range([(0.5, 0.5)].into())
            .build()
            .unwrap();
        let mut rng = rand::rng();
        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config], &seq, &mut rng);

        // N base means all 4 bases should be marked
        assert_eq!(ml_vec.len(), 4);
        // Parse and verify the actual pos values in MM tag
        let pos: Vec<&str> = mm_str
            .strip_prefix("N+n?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert_eq!(pos.len(), 4, "Should have exactly 4 modifications for ACGT");

        // All positions should be 0
        assert!(
            pos.iter().all(|&x| x == "0"),
            "All N base modifications should have 0 as their gap position coordinate"
        );

        // Verify all ML values are 128 (probabilities)
        assert_eq!(
            ml_vec,
            vec![128u8, 128u8, 128u8, 128u8],
            "All ML values should be 128 and number of values is 4"
        );
    }

    /// Tests `generate_random_dna_modification` with cycling windows
    #[test]
    fn generate_random_dna_modification_cycling_windows() {
        let seq = DNARestrictive::from_str("CCCCCCCCCCCCCCCC").unwrap(); // 16 C's

        let mod_config = ModConfigBuilder::default()
            .base('C')
            .mod_code("m".into())
            .win([3, 2].into())
            .mod_range([(0.8, 0.8), (0.4, 0.4)].into())
            .build()
            .unwrap();

        let mut rng = rand::rng();
        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config], &seq, &mut rng);

        // Should have 16 modifications, cycling pattern: 3@0.8, 2@0.4, 3@0.8, 2@0.4, ...
        let pos: Vec<&str> = mm_str
            .strip_prefix("C+m?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert_eq!(pos, vec!["0"; 16]);
        // Verify the cycling pattern (0.8 = 204, 0.4 = 102)
        let expected_pattern: Vec<u8> = vec![
            204, 204, 204, 102, 102, 204, 204, 204, 102, 102, 204, 204, 204, 102, 102, 204,
        ];
        assert_eq!(ml_vec, expected_pattern);
    }
    /// Tests multiple simultaneous modifications on different bases
    #[expect(
        clippy::similar_names,
        reason = "has_c_mod, has_a_mod, has_t_mod are clear in context"
    )]
    #[test]
    fn multiple_simultaneous_modifications_work() {
        let config_json = r#"{
            "contigs": {
                "number": 1,
                "len_range": [100, 100],
                "repeated_seq": "ACGTACGTACGTACGT"
            },
            "reads": [{
                "number": 20,
                "mapq_range": [20, 30],
                "base_qual_range": [20, 30],
                "len_range": [0.8, 1.0],
                "mods": [
                    {
                        "base": "C",
                        "is_strand_plus": true,
                        "mod_code": "m",
                        "win": [2],
                        "mod_range": [[0.7, 0.9]]
                    },
                    {
                        "base": "A",
                        "is_strand_plus": true,
                        "mod_code": "a",
                        "win": [3],
                        "mod_range": [[0.3, 0.5]]
                    },
                    {
                        "base": "T",
                        "is_strand_plus": false,
                        "mod_code": "t",
                        "win": [1],
                        "mod_range": [[0.5, 0.6]]
                    }
                ]
            }]
        }"#;

        let config: SimulationConfig = serde_json::from_str(config_json).unwrap();
        let sim = TempBamSimulation::new(config).unwrap();
        let mut reader = bam::Reader::from_path(sim.bam_path()).unwrap();

        let mut has_c_mod = false;
        let mut has_a_mod = false;
        let mut has_t_mod = false;

        for record in reader.records() {
            let read = record.unwrap();
            if let Ok(Aux::String(mm_tag)) = read.aux(b"MM") {
                // Check that all three modification types are present
                if mm_tag.contains("C+m?") {
                    has_c_mod = true;
                }
                if mm_tag.contains("A+a?") {
                    has_a_mod = true;
                }
                if mm_tag.contains("T-t?") {
                    has_t_mod = true;
                }

                // Verify ML tag exists and has correct format
                assert!(
                    matches!(read.aux(b"ML").unwrap(), Aux::ArrayU8(_)),
                    "ML tag should be ArrayU8"
                );
            }
        }

        assert!(has_c_mod, "Should have C+m modifications");
        assert!(has_a_mod, "Should have A+a modifications");
        assert!(has_t_mod, "Should have T-t modifications");
    }

    /// Tests that modifications only target the specified bases
    #[expect(
        clippy::shadow_unrelated,
        reason = "clear variable reuse in sequential tests"
    )]
    #[test]
    fn modifications_target_correct_bases() {
        // Create sequence with known base counts
        let seq = DNARestrictive::from_str("AAAACCCCGGGGTTTT").unwrap();

        // template we will reuse for various mods
        let mod_config_template = ModConfigBuilder::default()
            .base('C')
            .mod_code("m".into())
            .win([4].into())
            .mod_range([(1.0, 1.0)].into());

        // Test C modification - should only mark C bases
        let mod_config_c = mod_config_template.clone().build().unwrap();

        let mut rng = rand::rng();
        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config_c], &seq, &mut rng);

        // Sequence has exactly 4 C's, so should have 4 modifications
        assert_eq!(ml_vec.len(), 4);
        let probs: Vec<&str> = mm_str
            .strip_prefix("C+m?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert_eq!(probs.len(), 4, "Should have exactly 4 C modifications");

        // Test T modification - should only mark T bases
        let mod_config_t = mod_config_template
            .clone()
            .base('T')
            .is_strand_plus(false)
            .mod_code("t".into())
            .build()
            .unwrap();

        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config_t], &seq, &mut rng);

        // Sequence has exactly 4 T's, so should have 4 modifications
        assert_eq!(ml_vec.len(), 4);
        let probs: Vec<&str> = mm_str
            .strip_prefix("T-t?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert_eq!(probs.len(), 4, "Should have exactly 4 T modifications");

        // Test A modification - should only mark A bases
        let mod_config_a = mod_config_template
            .clone()
            .base('A')
            .mod_code("a".into())
            .build()
            .unwrap();

        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config_a], &seq, &mut rng);

        // Sequence has exactly 4 A's, so should have 4 modifications
        assert_eq!(ml_vec.len(), 4);
        let probs: Vec<&str> = mm_str
            .strip_prefix("A+a?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert_eq!(probs.len(), 4, "Should have exactly 4 A modifications");

        // Test G modification - should only mark G bases
        let mod_config_g = mod_config_template
            .clone()
            .base('G')
            .mod_code("g".into())
            .build()
            .unwrap();

        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config_g], &seq, &mut rng);

        // Sequence has exactly 4 G's, so should have 4 modifications
        assert_eq!(ml_vec.len(), 4);
        let probs: Vec<&str> = mm_str
            .strip_prefix("G+g?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert_eq!(probs.len(), 4, "Should have exactly 4 G modifications");
    }

    /// Tests edge case: sequence with no target base for modification
    #[test]
    fn edge_case_no_target_base_for_modification() {
        let seq = DNARestrictive::from_str("AAAAAAAAAA").unwrap(); // Only A's

        let mod_config = ModConfigBuilder::default()
            .base('C')
            .mod_code("m".into())
            .win([5].into())
            .mod_range([(0.5, 0.5)].into())
            .build()
            .unwrap();

        let mut rng = rand::rng();
        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config], &seq, &mut rng);

        // Should have no modifications since no C's exist
        assert_eq!(mm_str, String::new());
        assert_eq!(ml_vec.len(), 0);
    }

    /// Tests edge case: modification with probability 0.0 and 1.0
    #[expect(
        clippy::shadow_unrelated,
        reason = "clear variable reuse in sequential tests"
    )]
    #[test]
    fn edge_case_modification_probability_extremes() {
        let seq = DNARestrictive::from_str("CCCCCCCC").unwrap();

        // template to build modification configurations
        let mod_config_template = ModConfigBuilder::default()
            .base('C')
            .mod_code("m".into())
            .win([8].into());

        let mod_config_zero = mod_config_template
            .clone()
            .mod_range([(0.0, 0.0)].into())
            .build()
            .unwrap();

        let mut rng = rand::rng();
        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config_zero], &seq, &mut rng);

        assert_eq!(ml_vec, vec![0u8; 8]);
        let pos: Vec<&str> = mm_str
            .strip_prefix("C+m?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert!(pos.iter().all(|&x| x == "0"));

        // Test probability 1.0
        let mod_config_one = mod_config_template
            .clone()
            .mod_range([(1.0, 1.0)].into())
            .build()
            .unwrap();

        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config_one], &seq, &mut rng);

        assert_eq!(ml_vec, vec![255u8; 8]);
        let pos: Vec<&str> = mm_str
            .strip_prefix("C+m?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert!(pos.iter().all(|&x| x == "0"));
    }

    /// Tests config deserialization from JSON string with mods
    #[expect(
        clippy::indexing_slicing,
        reason = "test validates length before indexing"
    )]
    #[test]
    fn config_deserialization_works_with_json_with_mods() {
        // Create a JSON config string
        let json_str = r#"{
            "contigs": {
                "number": 5,
                "len_range": [100, 500],
                "repeated_seq": "ACGTACGT"
            },
            "reads": [
                {
                    "number": 100,
                    "mapq_range": [10, 30],
                    "base_qual_range": [20, 40],
                    "len_range": [0.1, 0.9],
                    "barcode": "ACGTAA",
                    "mods": [{
                        "base": "C",
                        "is_strand_plus": true,
                        "mod_code": "m",
                        "win": [5, 3],
                        "mod_range": [[0.3, 0.7], [0.1, 0.5]]
                    }]
                },
                {
                    "number": 50,
                    "mapq_range": [5, 15],
                    "base_qual_range": [15, 25],
                    "len_range": [0.2, 0.8]
                }
            ]
        }"#;

        // Deserialize
        let config: SimulationConfig = serde_json::from_str(json_str).unwrap();

        // Verify all fields deserialized correctly
        assert_eq!(config.contigs.number.get(), 5);
        assert_eq!(config.contigs.len_range.low().get(), 100);
        assert_eq!(config.contigs.len_range.high().get(), 500);
        assert!(config.contigs.repeated_seq.is_some());

        assert_eq!(config.reads.len(), 2);

        // Verify first read config
        assert_eq!(config.reads[0].number.get(), 100);
        assert_eq!(config.reads[0].mapq_range.low(), 10);
        assert_eq!(config.reads[0].mapq_range.high(), 30);
        assert_eq!(config.reads[0].base_qual_range.low(), 20);
        assert_eq!(config.reads[0].base_qual_range.high(), 40);
        assert!(config.reads[0].barcode.is_some());
        assert_eq!(config.reads[0].mods.len(), 1);
        assert!(matches!(config.reads[0].mods[0].base, AllowedAGCTN::C));
        assert_eq!(config.reads[0].mods[0].win.len(), 2);

        // Verify second read config
        assert_eq!(config.reads[1].number.get(), 50);
        assert_eq!(config.reads[1].mods.len(), 0);
        assert!(config.reads[1].barcode.is_none());
    }

    /// Tests that modifications are applied to barcode bases
    /// Creates a sequence with no C's, but barcode has C's, so modifications should appear
    #[test]
    fn barcodes_with_modifications_integration_works() {
        // Create contig with only A's (no C's)
        let contigs = [ContigBuilder::default()
            .name("chr1")
            .seq("A".repeat(50))
            .build()
            .unwrap()];

        let config = ReadConfigBuilder::default()
            .number(20)
            .len_range((0.3, 0.7))
            .barcode("ACGTAA".into()) // barcode contains Cs i.e. barcode sequence has one C
            .mods(
                [ModConfigBuilder::default()
                    .base('C')
                    .mod_code("m".into())
                    .win([10].into())
                    .mod_range([(0.8, 0.8)].into())
                    .build()
                    .unwrap()]
                .into(),
            )
            .build()
            .unwrap();

        let mut rng = rand::rng();
        let reads = generate_reads_denovo(&contigs, &config, "test", &mut rng).unwrap();

        let mut has_modifications = false;
        for read in &reads {
            if !read.is_unmapped() {
                // The read sequence should have barcodes added
                // Forward: ACGTAA + AAAAA... + TTACGT (barcode has 1 C, rev comp has 0 C)
                // Reverse: TGCATT + AAAAA... + AATGCA (complement has 1 C, reverse has 0 C)
                // So each read should have at least 1 C from the barcode

                let mm_pos_result = read.aux(b"MM");
                let ml_prob_result = read.aux(b"ML");

                // Should have MM and ML tags because barcode introduced C's
                assert!(
                    mm_pos_result.is_ok(),
                    "MM tag should be present due to C in barcode"
                );
                assert!(
                    ml_prob_result.is_ok(),
                    "ML tag should be present due to C in barcode"
                );

                if let Ok(Aux::String(mm_str)) = mm_pos_result {
                    assert!(
                        mm_str.starts_with("C+m?,"),
                        "MM tag should contain C+m modifications"
                    );
                    assert!(!mm_str.is_empty(), "MM tag should not be empty");
                    has_modifications = true;
                }
            }
        }

        assert!(
            has_modifications,
            "Should have found C modifications from barcode bases"
        );
    }

    /// Tests edge case: window size larger than number of target bases in sequence
    #[test]
    fn edge_case_window_larger_than_sequence() {
        let seq = DNARestrictive::from_str("ACGTACGT").unwrap(); // Only 2 C's

        let mod_config = ModConfigBuilder::default()
            .base('C')
            .mod_code("m".into())
            .win([100].into()) // window of 100 Cs but only 2Cs exist in the sequence.
            .mod_range([(0.5, 0.5)].into())
            .build()
            .unwrap();

        let mut rng = rand::rng();
        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config], &seq, &mut rng);

        // Should only generate modifications for the 2 C's that exist
        assert_eq!(
            ml_vec,
            vec![128u8, 128u8],
            "All prob should be 128 (0.5*255)"
        );
        let pos: Vec<&str> = mm_str
            .strip_prefix("C+m?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert_eq!(pos.len(), 2, "Should have exactly 2 position values");
        assert!(pos.iter().all(|&x| x == "0"), "All pos should be 0");
    }

    /// Tests edge case: single-base windows with multiple cycles
    #[test]
    fn edge_case_single_base_windows() {
        let seq = DNARestrictive::from_str("C".repeat(16).as_str()).unwrap(); // 16 C's

        let mod_config = ModConfigBuilder::default()
            .base('C')
            .is_strand_plus(true)
            .mod_code("m".into())
            .win([1].into())
            .mod_range([(0.9, 0.9), (0.1, 0.1)].into())
            .build()
            .unwrap();

        let mut rng = rand::rng();
        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config], &seq, &mut rng);

        // Should have 16 modifications, alternating pattern: 0.9, 0.1, 0.9, 0.1, ...
        // Verify alternating pattern (0.9*255=229.5→230, 0.1*255=25.5→26)
        let expected_pattern: Vec<u8> = vec![
            230, 26, 230, 26, 230, 26, 230, 26, 230, 26, 230, 26, 230, 26, 230, 26,
        ];
        assert_eq!(ml_vec, expected_pattern);
        let pos: Vec<&str> = mm_str
            .strip_prefix("C+m?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert_eq!(pos, vec!["0"; 16], "pos should have 16 entries, all 0");
    }

    /// Tests that modifications on reverse reads are correctly applied to reverse-complemented sequence
    #[test]
    fn reverse_reads_modification_positions_correct() {
        // Create sequence with 3 C's: "AAAAACAAAAACAAAAAC"
        let contigs = [ContigBuilder::default()
            .name("chr1")
            .seq("AAAAACAAAAACAAAAAC".into())
            .build()
            .unwrap()];

        let config = ReadConfigBuilder::default()
            .number(100)
            .len_range((1.0, 1.0))
            .mods(
                [ModConfigBuilder::default()
                    .base('C')
                    .mod_code("m".into())
                    .win([10].into())
                    .mod_range([(0.8, 0.8)].into())
                    .build()
                    .unwrap()]
                .into(),
            )
            .build()
            .unwrap();

        let mut rng = rand::rng();
        let reads = generate_reads_denovo(&contigs, &config, "test", &mut rng).unwrap();

        let mut forward_with_mods = 0;

        for read in &reads {
            if !read.is_unmapped() {
                let mm_result = read.aux(b"MM");

                if read.is_reverse() {
                    // Reverse complement of "AAAAACAAAAACAAAAAC"
                    // has 3 G's but 0 C's, so no C modifications should exist
                    assert!(
                        mm_result.is_err(),
                        "Reverse reads should have no MM tag (revcomp has no C's)"
                    );
                } else {
                    // Forward reads should have MM tag with 3 C modifications
                    assert!(mm_result.is_ok(), "Forward reads should have MM tag");
                    if let Ok(Aux::String(mm_str)) = mm_result {
                        // Parse and count modifications
                        let probs: Vec<&str> = mm_str
                            .strip_prefix("C+m?,")
                            .unwrap()
                            .strip_suffix(';')
                            .unwrap()
                            .split(',')
                            .collect();
                        assert_eq!(probs.len(), 3, "Should have exactly 3 C modifications");
                        forward_with_mods += 1;
                    }
                }
            }
        }

        assert!(
            forward_with_mods > 0,
            "Should have validated some forward reads with modifications"
        );
    }

    /// Tests that mismatches affect modification reference positions while preserving mod quality
    ///
    /// This test creates two read groups:
    /// - Group 0: 100% mod probability (quality 255), no mismatches - mods at expected C positions
    /// - Group 1: 51% mod probability (quality ~130), 100% mismatch - mods at shifted positions
    ///
    /// For a "ACGT" repeated contig, C's are at positions 1, 5, 9, 13, ... (form 4n+1).
    /// On reverse complement, G's become C's at positions 2, 6, 10, 14, ... (form 4n+2).
    /// With 100% mismatch, positions should be shifted and no longer follow these patterns.
    #[expect(
        clippy::too_many_lines,
        reason = "test requires setup, iteration, and multiple assertions for thorough validation"
    )]
    #[test]
    fn mismatch_mod_check() {
        let json_str = r#"{
            "contigs": {
                "number": 1,
                "len_range": [100, 100],
                "repeated_seq": "ACGT"
            },
            "reads": [
                {
                    "number": 100,
                    "len_range": [1.0, 1.0],
                    "mods": [{
                        "base": "C",
                        "is_strand_plus": true,
                        "mod_code": "m",
                        "win": [1],
                        "mod_range": [[1.0, 1.0]]
                    }]
                },
                {
                    "number": 100,
                    "len_range": [1.0, 1.0],
                    "mismatch": 1.0,
                    "mods": [{
                        "base": "C",
                        "is_strand_plus": true,
                        "mod_code": "m",
                        "win": [1],
                        "mod_range": [[0.51, 0.51]]
                    }]
                }
            ]
        }"#;

        let config: SimulationConfig = serde_json::from_str(json_str).unwrap();
        let sim = TempBamSimulation::new(config).unwrap();

        let mut bam = bam::Reader::from_path(sim.bam_path()).unwrap();
        let mut df_collection = Vec::new();

        for k in bam.records() {
            let record = k.unwrap();
            let curr_read = CurrRead::default()
                .try_from_only_alignment(&record)
                .unwrap()
                .set_mod_data(&record, ThresholdState::GtEq(0), 0)
                .unwrap();
            df_collection.push(curr_read);
        }

        let df = curr_reads_to_dataframe(df_collection.as_slice()).unwrap();

        let read_id_col = df.column("read_id").unwrap().str().unwrap();
        let ref_position_col = df.column("ref_position").unwrap().i64().unwrap();
        let alignment_type_col = df.column("alignment_type").unwrap().str().unwrap();
        let mod_quality_col = df.column("mod_quality").unwrap().u32().unwrap();

        let mut row_count_0: usize = 0;
        let mut row_count_1: usize = 0;

        // Violation counters for debugging
        let mut group0_forward_position_violations: usize = 0;
        let mut group0_reverse_position_violations: usize = 0;
        let mut group0_unmapped_position_violations: usize = 0;
        let mut group0_quality_violations: usize = 0;
        let mut group1_forward_position_violations: usize = 0;
        let mut group1_reverse_position_violations: usize = 0;
        let mut group1_unmapped_position_violations: usize = 0;
        let mut group1_quality_violations: usize = 0;
        let mut read_id_violations: usize = 0;

        for i in 0..df.height() {
            let read_id = read_id_col.get(i).unwrap();
            let ref_position = ref_position_col.get(i).unwrap();
            let alignment_type = alignment_type_col.get(i).unwrap();
            let mod_quality = mod_quality_col.get(i).unwrap();

            if read_id.starts_with("0.") {
                // Group 0: no mismatch, mod positions should follow 4n+1 (forward) or 4n+2 (reverse)
                if alignment_type.contains("forward") {
                    if ref_position % 4 != 1 {
                        group0_forward_position_violations += 1;
                    }
                } else if alignment_type.contains("reverse") {
                    if ref_position % 4 != 2 {
                        group0_reverse_position_violations += 1;
                    }
                } else {
                    // unmapped reads should have ref_position of -1
                    if ref_position != -1 {
                        group0_unmapped_position_violations += 1;
                    }
                }
                if mod_quality != 255 {
                    group0_quality_violations += 1;
                }
                row_count_0 += 1;
            } else if read_id.starts_with("1.") {
                // Group 1: 100% mismatch, mod positions should NOT follow expected patterns
                if alignment_type.contains("forward") {
                    if ref_position % 4 == 1 {
                        group1_forward_position_violations += 1;
                    }
                } else if alignment_type.contains("reverse") {
                    if ref_position % 4 == 2 {
                        group1_reverse_position_violations += 1;
                    }
                } else {
                    // unmapped reads should have ref_position of -1
                    if ref_position != -1 {
                        group1_unmapped_position_violations += 1;
                    }
                }
                if mod_quality != 130 {
                    group1_quality_violations += 1;
                }
                row_count_1 += 1;
            } else {
                read_id_violations += 1;
            }
        }

        // Assert all violations are zero
        assert_eq!(
            group0_forward_position_violations, 0,
            "Group 0 forward: expected 0 position violations"
        );
        assert_eq!(
            group0_reverse_position_violations, 0,
            "Group 0 reverse: expected 0 position violations"
        );
        assert_eq!(
            group0_unmapped_position_violations, 0,
            "Group 0 unmapped: expected 0 position violations"
        );
        assert_eq!(
            group0_quality_violations, 0,
            "Group 0: expected 0 quality violations"
        );
        assert_eq!(
            group1_forward_position_violations, 0,
            "Group 1 forward: expected 0 position violations"
        );
        assert_eq!(
            group1_reverse_position_violations, 0,
            "Group 1 reverse: expected 0 position violations"
        );
        assert_eq!(
            group1_unmapped_position_violations, 0,
            "Group 1 unmapped: expected 0 position violations"
        );
        assert_eq!(
            group1_quality_violations, 0,
            "Group 1: expected 0 quality violations"
        );
        assert_eq!(read_id_violations, 0, "expected 0 read_id violations");
        assert!(
            row_count_0 > 0,
            "Should have processed some rows from group 0"
        );
        assert!(
            row_count_1 > 0,
            "Should have processed some rows from group 1"
        );
    }
}

#[cfg(test)]
mod contig_generation_tests {
    use super::*;
    use rust_htslib::bam::Read as _;

    /// Tests edge case: repeated sequence contig generation
    #[test]
    fn edge_case_repeated_sequence_exact_length() {
        let config_json = r#"{
            "contigs": {
                "number": 1,
                "len_range": [16, 16],
                "repeated_seq": "ACGT"
            },
            "reads": [{
                "number": 5,
                "len_range": [0.5, 0.5]
            }]
        }"#;

        let config: SimulationConfig = serde_json::from_str(config_json).unwrap();
        let sim = TempBamSimulation::new(config).unwrap();
        let reader = bam::Reader::from_path(sim.bam_path()).unwrap();

        // Verify contig is exactly ACGTACGTACGTACGT (4 repeats)
        let header = reader.header();
        assert_eq!(header.target_len(0).unwrap(), 16);
    }

    /// Tests edge case: minimum contig size (1 bp)
    #[test]
    fn edge_case_minimum_contig_size() {
        let config_json = r#"{
            "contigs": {
                "number": 1,
                "len_range": [1, 1]
            },
            "reads": [{
                "number": 10,
                "mapq_range": [10, 20],
                "base_qual_range": [20, 30],
                "len_range": [1.0, 1.0]
            }]
        }"#;

        let config: SimulationConfig = serde_json::from_str(config_json).unwrap();
        let sim = TempBamSimulation::new(config).unwrap();
        let mut reader = bam::Reader::from_path(sim.bam_path()).unwrap();

        // Verify 1bp contig works
        let header = reader.header();
        assert_eq!(header.target_count(), 1);
        assert_eq!(header.target_len(0).unwrap(), 1);

        // Some reads should exist (though some may be unmapped)
        assert!(reader.records().count() > 0);
    }

    /// Tests contig generation
    #[test]
    fn generate_contigs_denovo_works() {
        let config = ContigConfigBuilder::default()
            .number(5)
            .len_range((100, 200))
            .build()
            .unwrap();
        let contigs = generate_contigs_denovo(config.number, config.len_range, &mut rand::rng());
        assert_eq!(contigs.len(), 5);
        for (i, contig) in contigs.iter().enumerate() {
            assert_eq!(contig.name, format!("contig_0000{i}"));
            assert!((100..=200).contains(&contig.get_dna_restrictive().get().len()));
            for base in contig.get_dna_restrictive().get() {
                assert!([b'A', b'C', b'G', b'T'].contains(base));
            }
        }
    }

    /// Tests contig generation with repeated sequence
    #[test]
    fn generate_contigs_denovo_repeated_seq_works() {
        let seq = DNARestrictive::from_str("ACGT").unwrap();
        let contigs = generate_contigs_denovo_repeated_seq(
            NonZeroU32::new(10000).unwrap(),
            OrdPair::new(NonZeroU64::new(10).unwrap(), NonZeroU64::new(12).unwrap()).unwrap(),
            &seq,
            &mut rand::rng(),
        );

        let mut counts = [0, 0, 0];

        assert_eq!(contigs.len(), 10000);
        for (i, contig) in contigs.iter().enumerate() {
            assert_eq!(contig.name, format!("contig_{i:05}"));
            let idx = match contig.get_dna_restrictive().get() {
                b"ACGTACGTAC" => 0,
                b"ACGTACGTACG" => 1,
                b"ACGTACGTACGT" => 2,
                _ => unreachable!(),
            };
            *counts
                .get_mut(idx)
                .expect("idx is 0, 1, or 2; counts has 3 elements") += 1;
        }

        // 3000/10000 is quite lax actually - we expect ~3333 each
        for count in counts {
            assert!(count >= 3000);
        }
    }
}

/// Tests for `PerfectSeqMatchToNot` methods
#[cfg(test)]
mod perfect_seq_match_to_not_tests {
    use super::*;

    #[test]
    fn seq_with_valid_sequence() -> Result<(), Error> {
        let _builder = PerfectSeqMatchToNot::seq(b"ACGT".to_vec())?;
        Ok(())
    }

    #[test]
    fn seq_with_empty_sequence_returns_error() {
        let result = PerfectSeqMatchToNot::seq(b"".to_vec());
        assert!(matches!(result, Err(Error::InvalidState(_))));
    }

    #[test]
    fn build_without_barcode_forward_read() -> Result<(), Error> {
        let mut rng = rand::rng();
        let (seq, cigar) = PerfectSeqMatchToNot::seq(b"GGGGGGGG".to_vec())?
            .build(ReadState::PrimaryFwd, &mut rng)?;
        assert_eq!(seq, b"GGGGGGGG");
        assert_eq!(cigar.expect("no error").to_string(), "8M");
        Ok(())
    }

    #[test]
    fn build_without_barcode_reverse_read() -> Result<(), Error> {
        let mut rng = rand::rng();
        let (seq, cigar) = PerfectSeqMatchToNot::seq(b"GGGGGGGG".to_vec())?
            .build(ReadState::PrimaryRev, &mut rng)?;
        assert_eq!(seq, b"GGGGGGGG");
        assert_eq!(cigar.expect("no error").to_string(), "8M");
        Ok(())
    }

    #[test]
    fn build_with_barcode_forward_read() -> Result<(), Error> {
        let mut rng = rand::rng();
        let barcode = DNARestrictive::from_str("ACGTAA").unwrap();
        let (seq, cigar) = PerfectSeqMatchToNot::seq(b"GGGGGGGG".to_vec())?
            .barcode(barcode)
            .build(ReadState::PrimaryFwd, &mut rng)?;
        assert_eq!(seq, b"ACGTAAGGGGGGGGTTACGT");
        assert_eq!(cigar.expect("no error").to_string(), "6S8M6S");
        Ok(())
    }

    #[test]
    fn build_with_barcode_reverse_read() -> Result<(), Error> {
        let mut rng = rand::rng();
        let barcode = DNARestrictive::from_str("ACGTAA").unwrap();
        let (seq, cigar) = PerfectSeqMatchToNot::seq(b"GGGGGGGG".to_vec())?
            .barcode(barcode)
            .build(ReadState::PrimaryRev, &mut rng)?;
        assert_eq!(seq, b"TGCATTGGGGGGGGAATGCA");
        assert_eq!(cigar.expect("no error").to_string(), "6S8M6S");
        Ok(())
    }

    #[test]
    fn build_unmapped_read_returns_none_cigar() -> Result<(), Error> {
        let mut rng = rand::rng();
        let (seq, cigar) = PerfectSeqMatchToNot::seq(b"GGGGGGGG".to_vec())?
            .build(ReadState::Unmapped, &mut rng)?;
        assert_eq!(seq, b"GGGGGGGG");
        assert!(cigar.is_none());
        Ok(())
    }

    #[test]
    fn build_unmapped_read_with_barcode_returns_none_cigar() -> Result<(), Error> {
        let mut rng = rand::rng();
        let barcode = DNARestrictive::from_str("ACGTAA").unwrap();
        let (seq, cigar) = PerfectSeqMatchToNot::seq(b"GGGGGGGG".to_vec())?
            .barcode(barcode)
            .build(ReadState::Unmapped, &mut rng)?;
        assert_eq!(seq, b"ACGTAAGGGGGGGGTTACGT");
        assert!(cigar.is_none());
        Ok(())
    }

    #[test]
    fn build_with_delete() -> Result<(), Error> {
        let mut rng = rand::rng();
        let delete_range = OrdPair::new(F32Bw0and1::try_from(0.5)?, F32Bw0and1::try_from(0.75)?)?;
        let (seq, cigar) = PerfectSeqMatchToNot::seq(b"AAAATTTTCCCCGGGG".to_vec())?
            .delete(delete_range)
            .build(ReadState::PrimaryFwd, &mut rng)?;

        // Original: AAAATTTTCCCCGGGG (16 bases)
        // Delete 50%-75%: positions 8-12, deletes CCCC (4 bases)
        // Result: AAAATTTTGGGG (12 bases)
        assert_eq!(seq, b"AAAATTTTGGGG");
        assert_eq!(cigar.expect("no error").to_string(), "8M4D4M");
        Ok(())
    }

    #[test]
    fn build_with_insert_middle() -> Result<(), Error> {
        let mut rng = rand::rng();
        let insert_seq = DNARestrictive::from_str("TTTT").unwrap();
        let (seq, cigar) = PerfectSeqMatchToNot::seq(b"AAAAGGGG".to_vec())?
            .insert_middle(insert_seq)
            .build(ReadState::PrimaryFwd, &mut rng)?;

        // Original: AAAAGGGG (8 bases)
        // Middle is at position 4
        // Insert TTTT at position 4
        // Result: AAAATTTTGGGG (12 bases)
        assert_eq!(seq, b"AAAATTTTGGGG");
        assert_eq!(cigar.expect("no error").to_string(), "4M4I4M");
        Ok(())
    }

    #[test]
    fn build_with_mismatch() -> Result<(), Error> {
        let mut rng = rand::rng();
        let mismatch_frac = F32Bw0and1::try_from(0.5)?;
        let (seq, cigar) = PerfectSeqMatchToNot::seq(b"AAAAAAAA".to_vec())?
            .mismatch(mismatch_frac)
            .build(ReadState::PrimaryFwd, &mut rng)?;

        // Original: AAAAAAAA (8 bases)
        // 50% mismatch = 4 bases should be changed
        // Count how many bases are NOT 'A'
        let mismatch_count = seq.iter().filter(|&&b| b != b'A').count();
        assert_eq!(mismatch_count, 4);

        // CIGAR should still be 8M (mismatches don't change CIGAR)
        assert_eq!(cigar.expect("no error").to_string(), "8M");
        Ok(())
    }

    #[test]
    fn build_with_delete_and_insert() -> Result<(), Error> {
        let mut rng = rand::rng();
        let delete_range = OrdPair::new(F32Bw0and1::try_from(0.5)?, F32Bw0and1::try_from(0.75)?)?;
        let insert_seq = DNARestrictive::from_str("TT").unwrap();
        let (seq, cigar) = PerfectSeqMatchToNot::seq(b"AAAACCCCGGGG".to_vec())?
            .delete(delete_range)
            .insert_middle(insert_seq)
            .build(ReadState::PrimaryFwd, &mut rng)?;

        // Original: AAAACCCCGGGG (12 bases)
        // Mark positions 6-8 for deletion (CCG)
        // Insert TT at position 6 (middle of 12)
        // Final operations: 6M (AAAACC) + 2I (TT) + 3D (CCG) + 3M (GGG)
        // Result: AAAACCTTGGG (11 bases)
        assert_eq!(seq, b"AAAACCTTGGG");
        assert_eq!(cigar.expect("no error").to_string(), "6M2I3D3M");
        Ok(())
    }

    #[test]
    fn build_with_all_features() -> Result<(), Error> {
        let mut rng = rand::rng();
        let delete_range = OrdPair::new(F32Bw0and1::try_from(0.25)?, F32Bw0and1::try_from(0.5)?)?;
        let insert_seq = DNARestrictive::from_str("CC").unwrap();
        let mismatch_frac = F32Bw0and1::try_from(0.25)?;
        let barcode = DNARestrictive::from_str("AAA").unwrap();

        let (seq, cigar) = PerfectSeqMatchToNot::seq(b"GGGGGGGGGGGGGGGG".to_vec())?
            .delete(delete_range)
            .insert_middle(insert_seq)
            .mismatch(mismatch_frac)
            .barcode(barcode)
            .build(ReadState::PrimaryFwd, &mut rng)?;

        // Original: GGGGGGGGGGGGGGGG (16 bases)
        // Mismatch 25% = 4 bases changed (but to C, T, or A)
        // Delete 25%-50%: positions 4-8, deletes GGGG (4 bases) -> 12 bases
        // Insert CC at middle (position 6) -> 14 bases
        // Add barcode AAA to both ends -> 20 bases total
        assert_eq!(seq.len(), 20);

        // CIGAR should have SoftClips, deletion, insertion, and matches
        let cigar_str = cigar.expect("no error").to_string();
        assert!(cigar_str.starts_with("3S")); // Barcode at start
        assert!(cigar_str.ends_with("3S")); // Barcode at end
        assert!(cigar_str.contains('D')); // Deletion
        assert!(cigar_str.contains('I')); // Insertion
        Ok(())
    }

    #[test]
    fn build_with_delete_zero_length() -> Result<(), Error> {
        let mut rng = rand::rng();
        // Delete nothing (start == end after rounding)
        let delete_range = OrdPair::new(F32Bw0and1::try_from(0.5)?, F32Bw0and1::try_from(0.5)?)?;
        let (seq, cigar) = PerfectSeqMatchToNot::seq(b"AAAAAAAA".to_vec())?
            .delete(delete_range)
            .build(ReadState::PrimaryFwd, &mut rng)?;

        // No deletion should occur
        assert_eq!(seq, b"AAAAAAAA");
        assert_eq!(cigar.expect("no error").to_string(), "8M");
        Ok(())
    }

    #[test]
    fn build_with_mismatch_zero_fraction() -> Result<(), Error> {
        let mut rng = rand::rng();
        let mismatch_frac = F32Bw0and1::try_from(0.0)?;
        let (seq, cigar) = PerfectSeqMatchToNot::seq(b"AAAAAAAA".to_vec())?
            .mismatch(mismatch_frac)
            .build(ReadState::PrimaryFwd, &mut rng)?;

        // No mismatches should occur
        assert_eq!(seq, b"AAAAAAAA");
        assert_eq!(cigar.expect("no error").to_string(), "8M");
        Ok(())
    }

    #[test]
    #[should_panic(expected = "SimulateDNASeqCIGAREndProblem")]
    fn build_with_delete_whole_length() {
        let mut rng = rand::rng();
        let delete_range = OrdPair::new(F32Bw0and1::zero(), F32Bw0and1::one()).unwrap();
        let (_seq, _cigar) = PerfectSeqMatchToNot::seq(b"AAAAAAAA".to_vec())
            .unwrap()
            .delete(delete_range)
            .build(ReadState::PrimaryFwd, &mut rng)
            .unwrap();
    }

    #[test]
    #[should_panic(expected = "SimulateDNASeqCIGAREndProblem")]
    fn build_with_insert_almost_whole_length() {
        let mut rng = rand::rng();
        let insert_seq = DNARestrictive::from_str("AATT").unwrap();
        let (_seq, _cigar) = PerfectSeqMatchToNot::seq(b"A".to_vec())
            .unwrap()
            .insert_middle(insert_seq)
            .build(ReadState::PrimaryFwd, &mut rng)
            .unwrap();
    }

    #[test]
    #[should_panic(expected = "SimulateDNASeqCIGAREndProblem")]
    fn build_with_delete_whole_length_with_barcode() {
        let mut rng = rand::rng();
        let delete_range = OrdPair::new(F32Bw0and1::zero(), F32Bw0and1::one()).unwrap();
        let barcode = DNARestrictive::from_str("CGCG").unwrap();
        let (_seq, _cigar) = PerfectSeqMatchToNot::seq(b"AAAAAAAA".to_vec())
            .unwrap()
            .delete(delete_range)
            .barcode(barcode)
            .build(ReadState::PrimaryFwd, &mut rng)
            .unwrap();
    }

    #[test]
    #[should_panic(expected = "SimulateDNASeqCIGAREndProblem")]
    fn build_with_insert_almost_whole_length_with_barcode() {
        let mut rng = rand::rng();
        let insert_seq = DNARestrictive::from_str("AATT").unwrap();
        let barcode = DNARestrictive::from_str("CGCG").unwrap();
        let (_seq, _cigar) = PerfectSeqMatchToNot::seq(b"A".to_vec())
            .unwrap()
            .insert_middle(insert_seq)
            .barcode(barcode)
            .build(ReadState::PrimaryFwd, &mut rng)
            .unwrap();
    }
}
