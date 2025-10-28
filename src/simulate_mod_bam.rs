//! # Write Simulated Mod BAM
//! Generates simulated BAM files with base modifications for testing purposes.
//! Accepts JSON configuration to produce both BAM and FASTA reference files.
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
//!     "mods": []
//!   }]
//! }"#;
//!
//! // Note: * "barcode" field is optional and can be omitted; barcodes added to both ends.
//! //         As sequences and lengths are generated independently of barcodes,
//! //         a barcode will _add_ 2 times so many bp to sequence length statistics.
//! //       * "repeated_seq" field is optional in contigs; if set, contigs are made by
//! //         repeating this sequence instead of generating random sequences.
//!
//! run(
//!     config_json,
//!     "output.bam",
//!     "reference.fasta"
//! ).unwrap();
//! ```

use crate::{DNARestrictive, Error, F32Bw0and1, GetDNARestrictive, ModChar, OrdPair, ReadState};
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
use std::str::FromStr;
use uuid::Uuid;

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
    /// Base that is modified (A, C, G, T, etc.)
    pub base: char,
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
            len_range: OrdPair::new(F32Bw0and1::new(0.0).unwrap(), F32Bw0and1::new(0.0).unwrap())
                .unwrap(),
            barcode: None,
            mods: Vec::new(),
        }
    }
}

impl Default for ModConfig {
    fn default() -> Self {
        Self {
            base: 'C',
            is_strand_plus: true,
            mod_code: ModChar::new('m'),
            win: vec![NonZeroU32::new(1).unwrap()],
            mod_range: vec![
                OrdPair::new(F32Bw0and1::new(0.0).unwrap(), F32Bw0and1::new(1.0).unwrap()).unwrap(),
            ],
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
/// use nanalogue_core::{DNARestrictive, F32Bw0and1, ModChar, OrdPair};
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
/// mod_config_c.base = 'C';
/// mod_config_c.is_strand_plus = true;
/// mod_config_c.mod_code = ModChar::new('m');
/// mod_config_c.win = vec![NonZeroU32::new(2).unwrap(), NonZeroU32::new(3).unwrap()];
/// mod_config_c.mod_range = vec![
///     OrdPair::new(F32Bw0and1::new(0.8).unwrap(), F32Bw0and1::new(0.8).unwrap()).unwrap(),
///     OrdPair::new(F32Bw0and1::new(0.4).unwrap(), F32Bw0and1::new(0.4).unwrap()).unwrap()
/// ];
///
/// let mut mod_config_a = ModConfig::default();
/// mod_config_a.base = 'A';
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
/// # Panics
/// Panics if the `mod_config` contains a base that is not A, G, C, T, or N.
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
        let base = match mod_config.base {
            'A' => b'A',
            'G' => b'G',
            'C' => b'C',
            'T' => b'T',
            'N' => b'N',
            _ => panic!("wrong base"),
        };
        let strand = if mod_config.is_strand_plus { '+' } else { '-' };
        let mod_code = mod_config.mod_code;
        let mut count = if base == b'N' {
            seq_bytes.len()
        } else {
            seq_bytes
                .iter()
                .zip(iter::repeat(&base))
                .filter(|(a, b)| a == b)
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
    iter::repeat_with(|| DNA_BASES[rng.random_range(0..4)])
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
        let contig = &contigs[contig_idx];
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
            let temp_seq = contig.get_dna_restrictive().get()
                [usize::try_from(start_pos).expect("number conversion error")..end_pos]
                .to_vec();
            let cigar_ops = vec![Cigar::Match(
                u32::try_from(read_len).expect("number conversion error"),
            )];

            // Add barcode if specified
            match &read_config.barcode {
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
    bam_path: String,
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
mod tests {
    use super::*;
    use rust_htslib::bam::Read;

    /// Test for generation of random DNA of a given length
    #[test]
    fn generate_random_dna_sequence_works() {
        let seq = generate_random_dna_sequence(NonZeroU64::new(100).unwrap(), &mut rand::rng());
        assert_eq!(seq.len(), 100);
        for base in seq {
            assert!([b'A', b'C', b'G', b'T'].contains(&base));
        }
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

    /// Tests read generation with desired properties.
    #[test]
    fn generate_reads_denovo_works() {
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
            len_range: OrdPair::new(F32Bw0and1::new(0.2).unwrap(), F32Bw0and1::new(0.8).unwrap())
                .unwrap(),
            barcode: None,
            mods: vec![],
        };

        let reads = generate_reads_denovo(&contigs, &config, "1", &mut rng).unwrap();
        assert_eq!(reads.len(), 10);

        for read in &reads {
            if !read.is_unmapped() {
                assert!((4..=16).contains(&read.seq_len()));
                assert!((10..=20).contains(&read.mapq()));
                assert!((0..2).contains(&read.tid()));
                assert!(read.qual().iter().all(|x| (30..=50).contains(x)));
            }
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
                "len_range": [0.1, 0.8],
                "mods": []
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

    /// Tests `TempBamSimulation` struct functionality
    #[test]
    fn temp_bam_simulation_struct() {
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
                "mods": []
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

        // Files will be automatically cleaned up when sim is dropped
    }

    /// Tests `add_barcode` function with forward and reverse reads
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
                assert_eq!(&seq[..6], b"CATGCC");
                assert_eq!(&seq[seq_len - 6..], b"GGCATG");
            } else {
                // Forward: GTACGG at start, revcomp(GTACGG)=CCGTAC at end
                assert_eq!(&seq[..6], b"GTACGG");
                assert_eq!(&seq[seq_len - 6..], b"CCGTAC");
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
                _ => panic!("Unexpected sequence"),
            };
            counts[idx] += 1;
        }

        // 3000/10000 is quite lax actually - we expect ~3333 each
        for count in counts {
            assert!(count >= 3000);
        }
    }
}
