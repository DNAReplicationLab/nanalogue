//! # Write Simulated Mod BAM
//! Generates simulated BAM files with base modifications for testing purposes.
//! Accepts JSON configuration to produce both BAM and FASTA reference files.
//! Please note that both BAM files and FASTA files are created from scratch,
//! so please do not specify pre-existing BAM or FASTA files in the output
//! path; if so, they will be overwritten.
//!
//! ## Example Usage
//!
//! ```no_run
//! use nanalogue_core::simulate_mod_bam::run;
//!
//! let config_json = r#"{
//!   "contigs": {
//!     "number": 4,
//!     "len_range": [1000, 8000]
//!   },
//!   "reads": [{
//!     "number": 1000,
//!     "mapq_range": [10, 20],
//!     "base_qual_range": [10, 20],
//!     "len_range": [0.1, 0.8],
//!     "barcode": "ACGTAA",
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
//! // Note: * "barcode" field is optional and can be omitted; barcodes added to both ends.
//! //         As sequences and lengths are generated independently of barcodes,
//! //         a barcode will _add_ 2 times so many bp to sequence length statistics.
//! //       * "repeated_seq" field is optional in contigs; if set, contigs are made by
//! //         repeating this sequence instead of generating random sequences.
//! //       * "mods" are optional as well, omitting them will create a BAM file
//! //          with no modification information.
//! //       * mod information is specified using windows i.e. in window of given size,
//! //         we enforce that mod probabilities per base are in the given range.
//! //         Multiple windows can be specified, and they repeat along the
//! //         read until it ends.
//!
//! run(
//!     config_json,
//!     "output.bam",
//!     "reference.fasta"
//! ).unwrap();
//! // Paths used here must not exist already as these are files created anew.
//! ```

use crate::{
    AllowedAGCTN, DNARestrictive, Error, F32Bw0and1, GetDNARestrictive, ModChar, OrdPair, ReadState,
};
use crate::{write_bam_denovo, write_fasta};
use itertools::join;
use rand::{Rng, random};
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use serde::{Deserialize, Serialize};
use std::iter;
use std::num::{NonZeroU32, NonZeroU64};
use std::ops::RangeInclusive;
use std::path::Path;
use std::str::FromStr as _;
use uuid::Uuid;

/// Creates an `OrdPair` of `F32Bw0and1` values
macro_rules! f32_ord_pair_bw_0_and_1 {
    ($low:expr, $high:expr) => {
        OrdPair::new(
            F32Bw0and1::new($low).unwrap(),
            F32Bw0and1::new($high).unwrap(),
        )
        .unwrap()
    };
}

/// Main configuration struct for simulation
#[derive(Debug, Default, Clone, Serialize, Deserialize)]
#[serde(default)]
#[non_exhaustive]
pub struct SimulationConfig {
    /// Configuration for contig generation
    pub contigs: ContigConfig,
    /// Configuration for read generation
    pub reads: Vec<ReadConfig>,
}

/// Configuration for contig generation
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
#[non_exhaustive]
pub struct ContigConfig {
    /// Number of contigs to generate
    pub number: NonZeroU32,
    /// Contig length range in bp [min, max]
    pub len_range: OrdPair<NonZeroU64>,
    /// Optional repeated sequence to use for contigs instead of random generation
    pub repeated_seq: Option<DNARestrictive>,
}

/// Configuration for read generation
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
#[non_exhaustive]
pub struct ReadConfig {
    /// Total number of reads to generate
    pub number: NonZeroU32,
    /// Mapping quality range [min, max]
    pub mapq_range: OrdPair<u8>,
    /// Base quality score range [min, max]
    pub base_qual_range: OrdPair<u8>,
    /// Read length range as fraction of contig [min, max] (e.g.: [0.1, 0.8] = 10% to 80%)
    pub len_range: OrdPair<F32Bw0and1>,
    /// Optional barcode DNA sequence to add to read ends
    pub barcode: Option<DNARestrictive>,
    /// Modification configurations
    pub mods: Vec<ModConfig>,
}

/// Configuration for modification generation
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
#[non_exhaustive]
pub struct ModConfig {
    /// Base that is modified (A, C, G, T, or N)
    pub base: AllowedAGCTN,
    /// Whether this is on the plus strand
    pub is_strand_plus: bool,
    /// Modification code (character or numeric)
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
    pub win: Vec<NonZeroU32>,
    /// Vector of modification density range e.g. [[0.4, 0.6], [0.1, 0.2]].
    /// Also see description of win above.
    pub mod_range: Vec<OrdPair<F32Bw0and1>>,
}

/// Represents a contig with name and sequence
#[derive(Debug, Clone, Serialize, Deserialize)]
#[non_exhaustive]
pub struct Contig {
    /// Contig name
    pub name: String,
    /// Contig sequence (A, C, G, T)
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
            len_range: f32_ord_pair_bw_0_and_1!(0.0, 0.0),
            barcode: None,
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
            mod_range: vec![f32_ord_pair_bw_0_and_1!(0.0, 1.0)],
        }
    }
}

/// Generates iterator with DNA modification probabilities
///
/// # Examples
///
/// ```
/// use std::num::NonZeroU32;
/// use std::str::FromStr;
/// use nanalogue_core::{AllowedAGCTN, DNARestrictive, F32Bw0and1, ModChar, OrdPair};
/// use nanalogue_core::simulate_mod_bam::{ModConfig, generate_random_dna_modification};
/// use rand::Rng;
///
/// // Create a DNA sequence with multiple cytosines
/// let seq = DNARestrictive::from_str("ACGTCGCGATCGACGTCGCGATCG").unwrap();
///
/// // Configure modification for cytosines on plus strand.
/// // NOTE: we use mod ranges with an identical low and high value below
/// // because we do not want to deal with random values in an example.
///
/// let mut mod_config_c = ModConfig::default();
/// mod_config_c.base = AllowedAGCTN::C;
/// mod_config_c.is_strand_plus = true;
/// mod_config_c.mod_code = ModChar::new('m');
/// mod_config_c.win = vec![NonZeroU32::new(2).unwrap(), NonZeroU32::new(3).unwrap()];
/// mod_config_c.mod_range = vec![
///     OrdPair::new(F32Bw0and1::new(0.8).unwrap(), F32Bw0and1::new(0.8).unwrap()).unwrap(),
///     OrdPair::new(F32Bw0and1::new(0.4).unwrap(), F32Bw0and1::new(0.4).unwrap()).unwrap()
/// ];
///
/// let mut mod_config_a = ModConfig::default();
/// mod_config_a.base = AllowedAGCTN::A;
/// mod_config_a.is_strand_plus = false;
/// mod_config_a.mod_code = ModChar::new('a');
/// mod_config_a.win = vec![NonZeroU32::new(2).unwrap()];
/// mod_config_a.mod_range = vec![
///     OrdPair::new(F32Bw0and1::new(0.2).unwrap(), F32Bw0and1::new(0.2).unwrap()).unwrap()
/// ];
///
/// let mod_config = vec![mod_config_c, mod_config_a];
///
/// let mut rng = rand::rng();
/// let (mm_str, ml_str) = generate_random_dna_modification(&mod_config, &seq, &mut rng);
///
/// // MM string format: C+m?prob1,prob2,prob3,prob4,...;A-a?,...;
/// assert_eq!(mm_str, String::from("C+m?,204,204,102,102,102,204,204,102;A-a?,51,51,51,51;"));
///
/// // ML string contains gap coordinates
/// assert_eq!(ml_str, vec![0; 12]);
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
            let low = u8::from(k.1.get_low());
            let high = u8::from(k.1.get_high());
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
            ml_vec.append(&mut vec![0; output.len()]);
            mm_str += format!(
                "{}{}{}?,{};",
                base as char,
                strand,
                mod_code,
                join(output, ",")
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
/// let seq = generate_random_dna_sequence(NonZeroU64::new(100).unwrap(), &mut rng);
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
/// let read_seq = b"GGGGGGGG".to_vec();
/// let barcode = DNARestrictive::from_str("ACGTAA").unwrap();
///
/// // Forward read: barcode + seq + revcomp(barcode)
/// let result = add_barcode(&read_seq, &barcode, ReadState::PrimaryFwd);
/// assert_eq!(result, b"ACGTAAGGGGGGGGTTACGT".to_vec());
///
/// // Reverse read: comp(barcode) + seq + rev(barcode)
/// let result = add_barcode(&read_seq, &barcode, ReadState::PrimaryRev);
/// assert_eq!(result, b"TGCATTGGGGGGGGAATGCA".to_vec());
/// ```
pub fn add_barcode<S: GetDNARestrictive>(
    read_seq: &[u8],
    barcode: &S,
    read_state: ReadState,
) -> Vec<u8> {
    let bc_bytes = barcode.get_dna_restrictive().get();

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
///     NonZeroU32::new(3).unwrap(),
///     OrdPair::new(NonZeroU64::new(100).unwrap(), NonZeroU64::new(200).unwrap()).unwrap(),
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
            let length = rng.random_range(len_range.get_low().get()..=len_range.get_high().get());
            let seq_bytes =
                generate_random_dna_sequence(NonZeroU64::try_from(length).expect("no error"), rng);
            let seq_str = std::str::from_utf8(&seq_bytes).expect("valid DNA sequence");
            Contig {
                name: format!("contig_{i:05}"),
                seq: DNARestrictive::from_str(seq_str).expect("valid DNA sequence"),
            }
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
///     NonZeroU32::new(2).unwrap(),
///     OrdPair::new(NonZeroU64::new(10).unwrap(), NonZeroU64::new(12).unwrap()).unwrap(),
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
    let seq_bytes = seq.get_dna_restrictive().get();

    (0..contig_number.get())
        .map(|i| {
            let length = rng.random_range(len_range.get_low().get()..=len_range.get_high().get());
            let contig_seq: Vec<u8> = seq_bytes
                .iter()
                .cycle()
                .take(usize::try_from(length).expect("number conversion error"))
                .copied()
                .collect();
            let seq_str =
                std::str::from_utf8(&contig_seq).expect("valid UTF-8 from repeated DNA sequence");
            Contig {
                name: format!("contig_{i:05}"),
                seq: DNARestrictive::from_str(seq_str)
                    .expect("valid DNA sequence from repeated pattern"),
            }
        })
        .collect()
}

/// Generates reads that align to contigs
///
/// ```
/// use std::num::NonZeroU32;
/// use std::str::FromStr;
/// use nanalogue_core::{DNARestrictive, OrdPair, F32Bw0and1};
/// use nanalogue_core::simulate_mod_bam::{Contig, ReadConfig, generate_reads_denovo};
/// use rand::Rng;
///
/// let mut contig = Contig::default();
/// contig.name = "chr1".to_string();
/// contig.seq = DNARestrictive::from_str("ACGTACGTACGTACGT").unwrap();
/// let contigs = vec![contig];
///
/// let mut config = ReadConfig::default();
/// config.number = NonZeroU32::new(10).unwrap();
/// config.mapq_range = OrdPair::new(10, 20).unwrap();
/// config.base_qual_range = OrdPair::new(20, 30).unwrap();
/// config.len_range = OrdPair::new(F32Bw0and1::new(0.2).unwrap(),
///     F32Bw0and1::new(0.5).unwrap()).unwrap();
/// config.barcode = None;
/// config.mods = vec![];
/// // NOTE: barcodes are optional, and will add 2*barcode_length to length statistics.
/// //       i.e. length stats are imposed independent of barcodes.
/// let mut rng = rand::rng();
/// let reads = generate_reads_denovo(&contigs, &config, "RG1", &mut rng).unwrap();
/// assert_eq!(reads.len(), 10);
/// ```
///
/// # Errors
/// Returns an error if the contigs input slice is empty,
/// if parameters are such that zero read lengths are produced,
/// or if BAM record creation fails, such as when making RG, MM, or ML tags.
#[expect(
    clippy::missing_panics_doc,
    reason = "number conversion errors or mis-generation of DNA bases are unlikely \
    as we generate DNA sequences ourselves here, and genomic data are unlikely \
    to exceed ~2^63 bp or have ~2^32 contigs"
)]
#[expect(
    clippy::cast_possible_truncation,
    reason = "read length calculated as a fraction of contig length, managed with trunc()"
)]
#[expect(
    clippy::too_many_lines,
    reason = "Complex read generation with multiple configuration options"
)]
pub fn generate_reads_denovo<R: Rng, S: GetDNARestrictive>(
    contigs: &[S],
    read_config: &ReadConfig,
    read_group: &str,
    rng: &mut R,
) -> Result<Vec<bam::Record>, Error> {
    if contigs.is_empty() {
        return Err(Error::UnavailableData);
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
        let read_len = ((rng.random_range(
            read_config.len_range.get_low().val()..=read_config.len_range.get_high().val(),
        ) * contig_len as f32)
            .trunc() as u64)
            .min(contig_len);

        // Ensure read length is at least 1; this also checks if contig_len is non-zero
        if read_len == 0 {
            return Err(Error::InvalidState(
                "Read length calculated as 0; increase len_range or contig size".into(),
            ));
        }

        // Set starting position
        let start_pos = rng.random_range(0..=(contig_len.saturating_sub(read_len)));

        // Extract sequence from contig
        let random_state: ReadState = random();
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
            let cigar_ops = vec![Cigar::Match(
                u32::try_from(read_len).expect("number conversion error"),
            )];

            // Add barcode if specified
            match read_config.barcode.as_ref() {
                Some(barcode) => {
                    let barcode_len = u32::try_from(barcode.get_dna_restrictive().get().len())
                        .expect("number conversion error");
                    (
                        add_barcode(&temp_seq, barcode, random_state),
                        CigarString(
                            [
                                &[Cigar::SoftClip(barcode_len)],
                                &cigar_ops[..],
                                &[Cigar::SoftClip(barcode_len)],
                            ]
                            .concat(),
                        ),
                    )
                }
                None => (temp_seq, CigarString(cigar_ops)),
            }
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
            match random_state {
                ReadState::Unmapped
                | ReadState::PrimaryFwd
                | ReadState::SecondaryFwd
                | ReadState::SupplementaryFwd => {
                    seq = DNARestrictive::from_str(str::from_utf8(&read_seq).expect("no error"))
                        .expect("no error");
                    &seq
                }
                ReadState::PrimaryRev | ReadState::SecondaryRev | ReadState::SupplementaryRev => {
                    seq = DNARestrictive::from_str(
                        str::from_utf8(&bio::alphabets::dna::revcomp(&read_seq)).expect("no error"),
                    )
                    .expect("no error");
                    &seq
                }
            },
            rng,
        );

        // Create BAM record
        let record = {
            let mut record = bam::Record::new();
            let qname = format!("{}.{}", read_group, Uuid::new_v4()).into_bytes();
            record.unset_flags();
            record.set_flags(u16::from(random_state));
            if random_state == ReadState::Unmapped {
                record.set(&qname, None, &read_seq, &qual);
                record.set_mapq(255);
                record.set_tid(-1);
                record.set_pos(-1);
            } else {
                let mapq = rng.random_range(RangeInclusive::from(read_config.mapq_range));
                record.set(&qname, Some(&cigar), &read_seq, &qual);
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
/// ```no_run
/// use nanalogue_core::simulate_mod_bam::run;
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
///     "len_range": [0.1, 0.8],
///   }]
/// }"#;
///
/// // Note: Optional "barcode" field can be added to reads (e.g., "barcode": "ACGTAA")
/// //       Length statistics are imposed independent of barcodes; so they will add
/// //       2x barcode length on top of these.
/// //       Optional "repeated_seq" field can be added to contigs (e.g., "repeated_seq": "ACGT")
/// //       to create contigs by repeating a sequence instead of random generation.
/// run(config_json, "output.bam", "reference.fasta").unwrap();
/// ```
///
/// # Errors
/// Returns an error if JSON parsing fails, read generation fails, or BAM/FASTA writing fails.
pub fn run<F>(config_json: &str, bam_output_path: &F, fasta_output_path: &F) -> Result<(), Error>
where
    F: AsRef<Path> + ?Sized,
{
    let config: SimulationConfig = serde_json::from_str(config_json)?;

    let mut rng = rand::rng();

    let contigs = match config.contigs.repeated_seq {
        Some(seq) => generate_contigs_denovo_repeated_seq(
            config.contigs.number,
            config.contigs.len_range,
            &seq,
            &mut rng,
        ),
        None => generate_contigs_denovo(config.contigs.number, config.contigs.len_range, &mut rng),
    };
    let read_groups: Vec<String> = (0..config.reads.len()).map(|k| k.to_string()).collect();
    let reads = {
        let mut temp_reads = Vec::new();
        for k in config.reads.into_iter().zip(read_groups.clone()) {
            temp_reads.append(&mut generate_reads_denovo(&contigs, &k.0, &k.1, &mut rng)?);
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
        contigs.into_iter().map(|k| {
            let seq_vec = k.get_dna_restrictive().get().to_vec();
            (k.name, seq_vec)
        }),
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
    pub fn new(config_json: &str) -> Result<Self, Error> {
        let temp_dir = std::env::temp_dir();
        let bam_path = temp_dir.join(format!("{}.bam", Uuid::new_v4()));
        let fasta_path = temp_dir.join(format!("{}.fa", Uuid::new_v4()));

        run(config_json, &bam_path, &fasta_path)?;

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
mod read_generation_no_mods_tests {
    use super::*;
    use rust_htslib::bam::Read as _;

    /// Tests read generation with desired properties but no modifications
    #[test]
    fn generate_reads_denovo_no_mods_works() {
        let mut rng = rand::rng();
        let contigs = vec![
            Contig {
                name: "contig_0".to_string(),
                seq: DNARestrictive::from_str("ACGTACGTACGTACGTACGT").unwrap(),
            },
            Contig {
                name: "contig_1".to_string(),
                seq: DNARestrictive::from_str("TGCATGCATGCATGCATGCA").unwrap(),
            },
        ];

        let config = ReadConfig {
            number: NonZeroU32::new(10).unwrap(),
            mapq_range: OrdPair::new(10, 20).unwrap(),
            base_qual_range: OrdPair::new(30, 50).unwrap(),
            len_range: f32_ord_pair_bw_0_and_1!(0.2, 0.8),
            barcode: None,
            mods: vec![],
        };

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
        run(config_json, &bam_path, &fasta_path).unwrap();
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
        let sim = TempBamSimulation::new(config_json).unwrap();

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
            let sim = TempBamSimulation::new(config_json).unwrap();
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
        assert!(matches!(result, Err(Error::UnavailableData)));
    }

    /// Tests error when read length would be zero
    #[test]
    fn generate_reads_denovo_zero_length_error() {
        let contigs = vec![Contig {
            name: "tiny".to_string(),
            seq: DNARestrictive::from_str("ACGT").unwrap(),
        }];

        let config = ReadConfig {
            number: NonZeroU32::new(1).unwrap(),
            len_range: f32_ord_pair_bw_0_and_1!(0.0, 0.0),
            ..Default::default()
        };

        let mut rng = rand::rng();
        let result = generate_reads_denovo(&contigs, &config, "test", &mut rng);
        assert!(matches!(result, Err(Error::InvalidState(_))));
    }

    /// Tests error when JSON deserialization fails
    #[test]
    fn run_bad_json_error() {
        let bad_json = r"{ invalid json }";
        let temp_dir = std::env::temp_dir();
        let bam_path = temp_dir.join(format!("{}.bam", Uuid::new_v4()));
        let fasta_path = temp_dir.join(format!("{}.fa", Uuid::new_v4()));

        let result = run(bad_json, &bam_path, &fasta_path);
        assert!(result.is_err());
    }

    /// Tests invalid JSON structure causing empty reads generation
    #[test]
    fn run_empty_reads_error() {
        let invalid_json = r#"{ "reads": [] }"#; // Empty reads array
        let temp_dir = std::env::temp_dir();
        let bam_path = temp_dir.join(format!("{}.bam", Uuid::new_v4()));
        let fasta_path = temp_dir.join(format!("{}.fa", Uuid::new_v4()));

        let result = run(invalid_json, &bam_path, &fasta_path);
        // With empty reads, this should succeed but produce an empty BAM (valid)
        // So we won't assert error here, just test it doesn't crash
        drop(result);
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

        let sim = TempBamSimulation::new(config_json).unwrap();
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

        // Count reads per read group
        let mut rg_counts = [0, 0, 0];
        for record in reader2.records() {
            let read = record.unwrap();
            if let Ok(Aux::String(rg)) = read.aux(b"RG") {
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
            Contig {
                name: "chr1".to_string(),
                seq: DNARestrictive::from_str("ACGTACGTACGTACGTACGTACGTACGTACGT").unwrap(),
            },
            Contig {
                name: "chr2".to_string(),
                seq: DNARestrictive::from_str("TGCATGCATGCATGCATGCATGCATGCATGCA").unwrap(),
            },
        ];

        let read_config = ReadConfig {
            number: NonZeroU32::new(50).unwrap(),
            mapq_range: OrdPair::new(10, 20).unwrap(),
            base_qual_range: OrdPair::new(20, 30).unwrap(),
            len_range: f32_ord_pair_bw_0_and_1!(0.2, 0.8),
            barcode: None, // No barcodes to simplify validation
            mods: vec![],
            ..Default::default()
        };

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

        let sim = TempBamSimulation::new(config_json).unwrap();
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
        let result = add_barcode(&read_seq, &barcode, ReadState::PrimaryFwd);
        assert_eq!(result, b"ACGTAAGGGGGGGGTTACGT".to_vec());

        // Test reverse read: comp(barcode) + seq + rev(barcode)
        let result = add_barcode(&read_seq, &barcode, ReadState::PrimaryRev);
        assert_eq!(result, b"TGCATTGGGGGGGGAATGCA".to_vec());
    }

    /// Tests CIGAR string validity with and without barcodes
    #[test]
    fn cigar_strings_valid_with_or_without_barcodes() {
        // Test without barcodes
        let contigs = vec![Contig {
            name: "chr1".to_string(),
            seq: DNARestrictive::from_str("ACGTACGTACGTACGTACGTACGTACGTACGT").unwrap(),
        }];

        let config_no_barcode = ReadConfig {
            number: NonZeroU32::new(20).unwrap(),
            mapq_range: OrdPair::new(10, 20).unwrap(),
            base_qual_range: OrdPair::new(20, 30).unwrap(),
            len_range: f32_ord_pair_bw_0_and_1!(0.3, 0.7),
            barcode: None,
            mods: vec![],
        };

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
        let config_with_barcode = ReadConfig {
            number: NonZeroU32::new(20).unwrap(),
            mapq_range: OrdPair::new(10, 20).unwrap(),
            base_qual_range: OrdPair::new(20, 30).unwrap(),
            len_range: f32_ord_pair_bw_0_and_1!(0.3, 0.7),
            barcode: Some(DNARestrictive::from_str("ACGTAA").unwrap()),
            mods: vec![],
        };

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

        let sim = TempBamSimulation::new(config_json).unwrap();
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
    use rust_htslib::bam::Read as _;

    /// Tests read generation with desired properties but with modifications
    #[test]
    fn generate_reads_denovo_with_mods_works() {
        let mut rng = rand::rng();
        let contigs = vec![
            Contig {
                name: "contig_0".to_string(),
                seq: DNARestrictive::from_str("ACGTACGTACGTACGTACGT").unwrap(),
            },
            Contig {
                name: "contig_1".to_string(),
                seq: DNARestrictive::from_str("TGCATGCATGCATGCATGCA").unwrap(),
            },
        ];

        let config = ReadConfig {
            number: NonZeroU32::new(10).unwrap(),
            mapq_range: OrdPair::new(10, 20).unwrap(),
            base_qual_range: OrdPair::new(30, 50).unwrap(),
            len_range: f32_ord_pair_bw_0_and_1!(0.2, 0.8),
            barcode: None,
            mods: vec![ModConfig {
                base: AllowedAGCTN::C,
                is_strand_plus: true,
                mod_code: ModChar::new('m'),
                win: vec![NonZeroU32::new(4).unwrap()],
                mod_range: vec![
                    OrdPair::<F32Bw0and1>::new(F32Bw0and1::zero(), F32Bw0and1::one()).unwrap(),
                ],
            }],
        };

        let reads = generate_reads_denovo(&contigs, &config, "1", &mut rng).unwrap();
        assert_eq!(reads.len(), 10);

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
            assert!(mod_pos_mm_tag.starts_with("C+m?,"));
            assert!(mod_pos_mm_tag.ends_with(';'));
            for k in 0..mod_prob_ml_tag.len() {
                assert_eq!(mod_prob_ml_tag.get(k).unwrap(), 0u8);
            }
        }
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
        let sim = TempBamSimulation::new(config_json).unwrap();

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
        let mod_config = ModConfig {
            base: AllowedAGCTN::C,
            is_strand_plus: true,
            mod_code: ModChar::new('m'),
            win: vec![NonZeroU32::new(2).unwrap()],
            mod_range: vec![f32_ord_pair_bw_0_and_1!(0.5, 0.5)],
            ..Default::default()
        };

        let mut rng = rand::rng();
        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config], &seq, &mut rng);

        // Sequence has 4 C's, so we expect 4 modifications
        assert!(mm_str.starts_with("C+m?,"));
        assert!(mm_str.ends_with(';'));
        assert_eq!(ml_vec.len(), 4);
        // All should be gap coordinate 0
        assert!(ml_vec.iter().all(|&x| x == 0));
        // All should be probability 128 (0.5 * 255)
        let probs: Vec<&str> = mm_str
            .strip_prefix("C+m?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert_eq!(probs.len(), 4);
        assert!(probs.iter().all(|&x| x == "128"));
    }

    /// Tests `generate_random_dna_modification` with multiple modification configs
    #[test]
    fn generate_random_dna_modification_multiple_mods() {
        let seq = DNARestrictive::from_str("ACGTACGT").unwrap();

        let mod_config_c = ModConfig {
            base: AllowedAGCTN::C,
            is_strand_plus: true,
            mod_code: ModChar::new('m'),
            win: vec![NonZeroU32::new(1).unwrap()],
            mod_range: vec![f32_ord_pair_bw_0_and_1!(0.8, 0.8)],
            ..Default::default()
        };

        let mod_config_t = ModConfig {
            base: AllowedAGCTN::T,
            is_strand_plus: false,
            mod_code: ModChar::new('t'),
            win: vec![NonZeroU32::new(1).unwrap()],
            mod_range: vec![f32_ord_pair_bw_0_and_1!(0.4, 0.4)],
            ..Default::default()
        };

        let mut rng = rand::rng();
        let (mm_str, ml_vec) =
            generate_random_dna_modification(&[mod_config_c, mod_config_t], &seq, &mut rng);

        // Sequence has 2 C's and 2 T's
        assert!(mm_str.contains("C+m?,"));
        assert!(mm_str.contains("T-t?,"));
        assert_eq!(ml_vec.len(), 4);
        assert!(ml_vec.iter().all(|&x| x == 0), "All ML values should be 0");

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
        // All C probabilities should be 204 (0.8 * 255)
        assert!(
            c_probs.iter().all(|&x| x == "204"),
            "All C modifications should have probability 204 (0.8)"
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
        // All T probabilities should be 102 (0.4 * 255)
        assert!(
            t_probs.iter().all(|&x| x == "102"),
            "All T modifications should have probability 102 (0.4)"
        );
    }

    /// Tests `generate_random_dna_modification` with N base (all bases)
    #[test]
    fn generate_random_dna_modification_n_base() {
        let seq = DNARestrictive::from_str("ACGT").unwrap();

        let mod_config = ModConfig {
            base: AllowedAGCTN::N,
            is_strand_plus: true,
            mod_code: ModChar::new('n'),
            win: vec![NonZeroU32::new(4).unwrap()],
            mod_range: vec![f32_ord_pair_bw_0_and_1!(0.5, 0.5)],
            ..Default::default()
        };

        let mut rng = rand::rng();
        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config], &seq, &mut rng);

        // N base means all 4 bases should be marked
        assert!(mm_str.starts_with("N+n?,"));
        assert!(mm_str.ends_with(';'));
        assert_eq!(ml_vec.len(), 4);

        // Parse and verify the actual probability values in MM tag
        let probs: Vec<&str> = mm_str
            .strip_prefix("N+n?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert_eq!(
            probs.len(),
            4,
            "Should have exactly 4 modifications for ACGT"
        );

        // All probabilities should be 128 (0.5 * 255)
        assert!(
            probs.iter().all(|&x| x == "128"),
            "All N base modifications should have probability 128 (0.5)"
        );

        // Verify all ML values are 0 (gap coordinates)
        assert!(ml_vec.iter().all(|&x| x == 0), "All ML values should be 0");
    }

    /// Tests `generate_random_dna_modification` with cycling windows
    #[test]
    fn generate_random_dna_modification_cycling_windows() {
        let seq = DNARestrictive::from_str("CCCCCCCCCCCCCCCC").unwrap(); // 16 C's

        let mod_config = ModConfig {
            base: AllowedAGCTN::C,
            is_strand_plus: true,
            mod_code: ModChar::new('m'),
            win: vec![NonZeroU32::new(3).unwrap(), NonZeroU32::new(2).unwrap()],
            mod_range: vec![
                f32_ord_pair_bw_0_and_1!(0.8, 0.8),
                f32_ord_pair_bw_0_and_1!(0.4, 0.4),
            ],
            ..Default::default()
        };

        let mut rng = rand::rng();
        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config], &seq, &mut rng);

        // Should have 16 modifications, cycling pattern: 3@0.8, 2@0.4, 3@0.8, 2@0.4, ...
        assert_eq!(ml_vec.len(), 16);
        let probs: Vec<&str> = mm_str
            .strip_prefix("C+m?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert_eq!(probs.len(), 16);

        // Verify the cycling pattern (0.8 = 204, 0.4 = 102)
        let expected_pattern = vec![
            "204", "204", "204", "102", "102", "204", "204", "204", "102", "102", "204", "204",
            "204", "102", "102", "204",
        ];
        assert_eq!(probs, expected_pattern);
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

        let sim = TempBamSimulation::new(config_json).unwrap();
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

        // Test C modification - should only mark C bases
        let mod_config_c = ModConfig {
            base: AllowedAGCTN::C,
            is_strand_plus: true,
            mod_code: ModChar::new('m'),
            win: vec![NonZeroU32::new(4).unwrap()],
            mod_range: vec![f32_ord_pair_bw_0_and_1!(1.0, 1.0)],
            ..Default::default()
        };

        let mut rng = rand::rng();
        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config_c], &seq, &mut rng);

        // Sequence has exactly 4 C's, so should have 4 modifications
        assert_eq!(ml_vec.len(), 4);
        assert!(mm_str.starts_with("C+m?,"));
        let probs: Vec<&str> = mm_str
            .strip_prefix("C+m?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert_eq!(probs.len(), 4, "Should have exactly 4 C modifications");

        // Test T modification - should only mark T bases
        let mod_config_t = ModConfig {
            base: AllowedAGCTN::T,
            is_strand_plus: false,
            mod_code: ModChar::new('t'),
            win: vec![NonZeroU32::new(4).unwrap()],
            mod_range: vec![f32_ord_pair_bw_0_and_1!(1.0, 1.0)],
            ..Default::default()
        };

        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config_t], &seq, &mut rng);

        // Sequence has exactly 4 T's, so should have 4 modifications
        assert_eq!(ml_vec.len(), 4);
        assert!(mm_str.starts_with("T-t?,"));
        let probs: Vec<&str> = mm_str
            .strip_prefix("T-t?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert_eq!(probs.len(), 4, "Should have exactly 4 T modifications");

        // Test A modification - should only mark A bases
        let mod_config_a = ModConfig {
            base: AllowedAGCTN::A,
            is_strand_plus: true,
            mod_code: ModChar::new('a'),
            win: vec![NonZeroU32::new(4).unwrap()],
            mod_range: vec![f32_ord_pair_bw_0_and_1!(1.0, 1.0)],
            ..Default::default()
        };

        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config_a], &seq, &mut rng);

        // Sequence has exactly 4 A's, so should have 4 modifications
        assert_eq!(ml_vec.len(), 4);
        assert!(mm_str.starts_with("A+a?,"));
        let probs: Vec<&str> = mm_str
            .strip_prefix("A+a?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert_eq!(probs.len(), 4, "Should have exactly 4 A modifications");

        // Test G modification - should only mark G bases
        let mod_config_g = ModConfig {
            base: AllowedAGCTN::G,
            is_strand_plus: true,
            mod_code: ModChar::new('g'),
            win: vec![NonZeroU32::new(4).unwrap()],
            mod_range: vec![f32_ord_pair_bw_0_and_1!(1.0, 1.0)],
            ..Default::default()
        };

        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config_g], &seq, &mut rng);

        // Sequence has exactly 4 G's, so should have 4 modifications
        assert_eq!(ml_vec.len(), 4);
        assert!(mm_str.starts_with("G+g?,"));
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

        let mod_config = ModConfig {
            base: AllowedAGCTN::C, // Looking for C's
            is_strand_plus: true,
            mod_code: ModChar::new('m'),
            win: vec![NonZeroU32::new(5).unwrap()],
            mod_range: vec![f32_ord_pair_bw_0_and_1!(0.5, 0.5)],
            ..Default::default()
        };

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

        // Test probability 0.0
        let mod_config_zero = ModConfig {
            base: AllowedAGCTN::C,
            is_strand_plus: true,
            mod_code: ModChar::new('m'),
            win: vec![NonZeroU32::new(8).unwrap()],
            mod_range: vec![f32_ord_pair_bw_0_and_1!(0.0, 0.0)],
            ..Default::default()
        };

        let mut rng = rand::rng();
        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config_zero], &seq, &mut rng);

        assert_eq!(ml_vec.len(), 8);
        let probs: Vec<&str> = mm_str
            .strip_prefix("C+m?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert!(probs.iter().all(|&x| x == "0"));

        // Test probability 1.0
        let mod_config_one = ModConfig {
            base: AllowedAGCTN::C,
            is_strand_plus: true,
            mod_code: ModChar::new('m'),
            win: vec![NonZeroU32::new(8).unwrap()],
            mod_range: vec![f32_ord_pair_bw_0_and_1!(1.0, 1.0)],
            ..Default::default()
        };

        let (mm_str, ml_vec) = generate_random_dna_modification(&[mod_config_one], &seq, &mut rng);

        assert_eq!(ml_vec.len(), 8);
        let probs: Vec<&str> = mm_str
            .strip_prefix("C+m?,")
            .unwrap()
            .strip_suffix(';')
            .unwrap()
            .split(',')
            .collect();
        assert!(probs.iter().all(|&x| x == "255"));
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
        assert_eq!(config.contigs.len_range.get_low().get(), 100);
        assert_eq!(config.contigs.len_range.get_high().get(), 500);
        assert!(config.contigs.repeated_seq.is_some());

        assert_eq!(config.reads.len(), 2);

        // Verify first read config
        assert_eq!(config.reads[0].number.get(), 100);
        assert_eq!(config.reads[0].mapq_range.get_low(), 10);
        assert_eq!(config.reads[0].mapq_range.get_high(), 30);
        assert_eq!(config.reads[0].base_qual_range.get_low(), 20);
        assert_eq!(config.reads[0].base_qual_range.get_high(), 40);
        assert!(config.reads[0].barcode.is_some());
        assert_eq!(config.reads[0].mods.len(), 1);
        assert!(matches!(config.reads[0].mods[0].base, AllowedAGCTN::C));
        assert_eq!(config.reads[0].mods[0].win.len(), 2);

        // Verify second read config
        assert_eq!(config.reads[1].number.get(), 50);
        assert_eq!(config.reads[1].mods.len(), 0);
        assert!(config.reads[1].barcode.is_none());
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

        let sim = TempBamSimulation::new(config_json).unwrap();
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

        let sim = TempBamSimulation::new(config_json).unwrap();
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
        let config = ContigConfig {
            number: NonZeroU32::new(5).unwrap(),
            len_range: OrdPair::new(NonZeroU64::new(100).unwrap(), NonZeroU64::new(200).unwrap())
                .unwrap(),
            repeated_seq: None,
        };
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
