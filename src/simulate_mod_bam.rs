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
//! // Note: "barcode" field is optional and can be omitted; barcodes added to both ends.
//! //       As sequences and lengths are generated independently of barcodes,
//! //       a barcode will _add_ 2 times so many bp to sequence length statistics.
//!
//! run(
//!     config_json,
//!     "output.bam",
//!     "reference.fasta"
//! ).unwrap();
//! ```

use crate::{Error, F32Bw0and1, ModChar, OrdPair, ReadState};
use crate::{write_bam_denovo, write_fasta};
use rand::{Rng, random};
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use serde::{Deserialize, Serialize};
use std::num::{NonZeroU32, NonZeroU64};
use std::ops::RangeInclusive;
use std::path::Path;
use uuid::Uuid;

/// Main configuration struct for simulation
#[derive(Debug, Default, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct SimulationConfig {
    /// Configuration for contig generation
    pub contigs: ContigConfig,
    /// Configuration for read generation
    pub reads: Vec<ReadConfig>,
}

/// Configuration for contig generation
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(default)]
pub struct ContigConfig {
    /// Number of contigs to generate
    pub number: NonZeroU32,
    /// Contig length range in bp [min, max]
    pub len_range: OrdPair<NonZeroU64>,
}

/// Configuration for read generation
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
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
    pub barcode: Option<String>,
    /// Modification configurations
    pub mods: Vec<ModConfig>,
}

/// Configuration for modification generation
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(default)]
pub struct ModConfig {
    /// Base that is modified (A, C, G, T, etc.)
    pub base: char,
    /// Whether this is on the plus strand
    pub is_strand_plus: bool,
    /// Modification code (character or numeric)
    pub mod_code: ModChar,
    /// Window size for modification density variation
    pub win: NonZeroU32,
    /// Modification density range e.g. [0.2, 0.8]
    pub mod_range: OrdPair<F32Bw0and1>,
}

/// Represents a contig with name and sequence
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct Contig {
    /// Contig name
    pub name: String,
    /// Contig sequence (A, C, G, T)
    pub seq: Vec<u8>,
}

impl Default for ContigConfig {
    fn default() -> Self {
        Self {
            number: NonZeroU32::new(1).unwrap(),
            len_range: OrdPair::new(NonZeroU64::new(1).unwrap(), NonZeroU64::new(1).unwrap())
                .unwrap(),
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
            win: NonZeroU32::new(1).unwrap(),
            mod_range: OrdPair::new(F32Bw0and1::new(0.0).unwrap(), F32Bw0and1::new(0.0).unwrap())
                .unwrap(),
        }
    }
}

/// Generates random DNA sequence of given length
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
    (0..length.get())
        .map(|_| DNA_BASES[rng.random_range(0..4)])
        .collect()
}

/// Validates that a DNA sequence contains only valid bases (A, C, G, T).
/// Does not accept ambiguous bases like 'N'.
///
/// # Examples
/// ```
/// use nanalogue_core::simulate_mod_bam::is_valid_dna_restrictive;
///
/// assert!(is_valid_dna_restrictive("ACGT"));
/// assert!(is_valid_dna_restrictive("acgt"));
/// assert!(!is_valid_dna_restrictive("ACGTN"));
/// assert!(!is_valid_dna_restrictive(""));
/// ```
pub fn is_valid_dna_restrictive(seq: &str) -> bool {
    if seq.is_empty() {
        return false;
    }
    seq.bytes()
        .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
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
/// use nanalogue_core::ReadState;
///
/// let read_seq = b"GGGGGGGG".to_vec();
/// let barcode = "ACGTAA".to_string();
///
/// // Forward read: barcode + seq + revcomp(barcode)
/// let result = add_barcode(read_seq.clone(), barcode.clone(), ReadState::PrimaryFwd).unwrap();
/// assert_eq!(result, b"ACGTAAGGGGGGGGTTACGT".to_vec());
///
/// // Reverse read: comp(barcode) + seq + rev(barcode)
/// let result = add_barcode(read_seq.clone(), barcode.clone(), ReadState::PrimaryRev).unwrap();
/// assert_eq!(result, b"TGCATTGGGGGGGGAATGCA".to_vec());
/// ```
pub fn add_barcode(
    read_seq: Vec<u8>,
    barcode: String,
    read_state: ReadState,
) -> Result<Vec<u8>, Error> {
    // Validate barcode
    if !is_valid_dna_restrictive(&barcode) {
        return Err(Error::InvalidSeq);
    }

    let bc_bytes: Vec<u8> = barcode.bytes().map(|b| b.to_ascii_uppercase()).collect();

    let result = match read_state {
        ReadState::PrimaryFwd
        | ReadState::SecondaryFwd
        | ReadState::SupplementaryFwd
        | ReadState::Unmapped => {
            // Forward/Unmapped: barcode + read_seq + reverse_complement(barcode)
            let revcomp_bc = bio::alphabets::dna::revcomp(&bc_bytes);
            [bc_bytes, read_seq, revcomp_bc].concat()
        }
        ReadState::PrimaryRev | ReadState::SecondaryRev | ReadState::SupplementaryRev => {
            // Reverse: complement(barcode) + read_seq + reverse(barcode)
            let comp_bc: Vec<u8> = bc_bytes
                .iter()
                .map(|&b| bio::alphabets::dna::complement(b))
                .collect();
            let rev_bc: Vec<u8> = bc_bytes.iter().copied().rev().collect();
            [comp_bc, read_seq, rev_bc].concat()
        }
    };

    Ok(result)
}

/// Generates contigs according to configuration
///
/// ```
/// use std::num::{NonZeroU32, NonZeroU64};
/// use nanalogue_core::OrdPair;
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
/// assert!(contigs.iter().all(|c| (100..=200).contains(&c.seq.len())));
/// ```
pub fn generate_contigs_denovo<R: Rng>(
    contig_number: NonZeroU32,
    len_range: OrdPair<NonZeroU64>,
    rng: &mut R,
) -> Vec<Contig> {
    (0..contig_number.get())
        .map(|i| {
            let length = rng.random_range(len_range.get_low().get()..=len_range.get_high().get());
            Contig {
                name: format!("contig_{}", i),
                seq: generate_random_dna_sequence(
                    NonZeroU64::try_from(length).expect("no error"),
                    rng,
                ),
            }
        })
        .collect()
}

/// Generates reads that align to contigs
///
/// ```
/// use std::num::NonZeroU32;
/// use nanalogue_core::{OrdPair, F32Bw0and1};
/// use nanalogue_core::simulate_mod_bam::{Contig, ReadConfig, generate_reads_denovo};
/// use rand::Rng;
///
/// let contigs = vec![Contig {
///     name: "chr1".to_string(),
///     seq: b"ACGTACGTACGTACGT".to_vec(),
/// }];
/// let config = ReadConfig {
///     number: NonZeroU32::new(10).unwrap(),
///     mapq_range: OrdPair::new(10, 20).unwrap(),
///     base_qual_range: OrdPair::new(20, 30).unwrap(),
///     len_range: OrdPair::new(F32Bw0and1::new(0.2).unwrap(),
///         F32Bw0and1::new(0.5).unwrap()).unwrap(),
///     barcode: None,
///     mods: vec![],
/// };
/// // NOTE: barcodes are optional, and will add 2*barcode_length to length statistics.
/// //       i.e. length stats are imposed independent of barcodes.
/// let mut rng = rand::rng();
/// let reads = generate_reads_denovo(&contigs, config, "RG1", &mut rng).unwrap();
/// assert_eq!(reads.len(), 10);
/// ```
pub fn generate_reads_denovo<R: Rng>(
    contigs: &[Contig],
    config: ReadConfig,
    read_group: &str,
    rng: &mut R,
) -> Result<Vec<bam::Record>, Error> {
    let mut reads = Vec::new();

    for _ in 0..config.number.get() {
        // Select a random contig
        let contig_idx = rng.random_range(0..contigs.len());
        let contig = &contigs[contig_idx];
        let contig_len = contig.seq.len() as u64;

        // Calculate read length as fraction of contig length
        let min_read_len = (config.len_range.get_low().val() * contig_len as f32) as u64;
        let max_read_len = (config.len_range.get_high().val() * contig_len as f32) as u64;
        let read_len = rng.random_range(min_read_len..=max_read_len);

        // Set starting position
        let start_pos = {
            if read_len > contig_len {
                unreachable!()
            } else {
                rng.random_range(0..=(contig_len - read_len))
            }
        };

        // Extract sequence from contig
        let random_state: ReadState = random();
        let end_pos = (start_pos + read_len) as usize;
        let read_seq = {
            let temp_seq = contig.seq[start_pos as usize..end_pos].to_vec();
            // Add barcode if specified
            match &config.barcode {
                Some(barcode) => add_barcode(temp_seq, barcode.clone(), random_state)?,
                None => temp_seq,
            }
        };

        // Generate quality scores (for final read length including barcodes)
        let qual: Vec<u8> = (0..read_seq.len())
            .map(|_| rng.random_range(RangeInclusive::from(config.base_qual_range)))
            .collect();

        // Generate MAPQ
        let mapq = rng.random_range(RangeInclusive::from(config.mapq_range));

        // Create BAM record
        let record = {
            let mut record = bam::Record::new();
            let qname = format!("{}", Uuid::new_v4()).into_bytes();
            record.unset_flags();
            record.set_flags(u16::from(random_state));
            match random_state {
                ReadState::Unmapped => {
                    record.set(&qname, None, &read_seq, &qual);
                    record.set_mapq(255);
                    record.set_tid(-1);
                    record.set_pos(-1);
                }
                _ => {
                    let cigar = match &config.barcode {
                        Some(barcode) => {
                            let barcode_len = barcode.len() as u32;
                            CigarString(vec![
                                Cigar::SoftClip(barcode_len),
                                Cigar::Match(read_len as u32),
                                Cigar::SoftClip(barcode_len),
                            ])
                        }
                        None => CigarString(vec![Cigar::Match(read_len as u32)]),
                    };
                    record.set(&qname, Some(&cigar), &read_seq, &qual);
                    record.set_mapq(mapq);
                    record.set_tid(contig_idx as i32);
                    record.set_pos(start_pos as i64);
                }
            }
            record.set_mpos(-1);
            record.set_mtid(-1);
            record.push_aux(b"RG", Aux::String(read_group))?;
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
/// run(config_json, "output.bam", "reference.fasta").unwrap();
/// ```
pub fn run<F>(config_json: &str, bam_output_path: &F, fasta_output_path: &F) -> Result<bool, Error>
where
    F: AsRef<Path> + ?Sized,
{
    let config: SimulationConfig = serde_json::from_str(config_json)?;

    let mut rng = rand::rng();

    let contigs =
        generate_contigs_denovo(config.contigs.number, config.contigs.len_range, &mut rng);
    let read_groups: Vec<String> = (0..config.reads.len()).map(|k| k.to_string()).collect();
    let reads = {
        let mut temp_reads = Vec::new();
        for k in config.reads.into_iter().zip(read_groups.clone()) {
            temp_reads.append(&mut generate_reads_denovo(&contigs, k.0, &k.1, &mut rng)?);
        }
        temp_reads.sort_by_key(|k| (k.is_unmapped(), k.tid(), k.pos(), k.is_reverse()));
        temp_reads
    };

    write_bam_denovo(
        reads,
        contigs.iter().map(|k| (k.name.clone(), k.seq.len())),
        read_groups,
        vec![String::from("simulated BAM file, not real data")],
        bam_output_path,
    )?;
    write_fasta(
        contigs.into_iter().map(|k| (k.name, k.seq)),
        fasta_output_path,
    )?;

    Ok(true)
}

/// Temporary BAM simulation with automatic cleanup
///
/// Creates temporary BAM and FASTA files for testing purposes and
/// automatically removes them when dropped.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct TempBamSimulation {
    bam_path: String,
    fasta_path: String,
}

impl TempBamSimulation {
    /// Creates a new temporary BAM simulation from JSON configuration
    pub fn new(config_json: &str) -> Result<Self, Error> {
        let bam_path = format!("/tmp/{}.bam", Uuid::new_v4());
        let fasta_path = format!("/tmp/{}.fa", Uuid::new_v4());

        if !(run(config_json, &bam_path, &fasta_path)?) {
            return Err(Error::UnknownError);
        }

        Ok(Self {
            bam_path,
            fasta_path,
        })
    }

    /// Returns the path to the temporary BAM file
    pub fn bam_path(&self) -> &str {
        &self.bam_path
    }

    /// Returns the path to the temporary FASTA file
    pub fn fasta_path(&self) -> &str {
        &self.fasta_path
    }
}

impl Drop for TempBamSimulation {
    fn drop(&mut self) {
        // Ignore errors during cleanup - files may already be deleted
        let _ = std::fs::remove_file(&self.bam_path);
        let _ = std::fs::remove_file(&self.fasta_path);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::Read;
    use std::path::Path;

    /// Test for generation of random DNA of a given length
    #[test]
    fn test_generate_random_dna_sequence() {
        let seq = generate_random_dna_sequence(NonZeroU64::new(100).unwrap(), &mut rand::rng());
        assert_eq!(seq.len(), 100);
        for base in seq {
            assert!([b'A', b'C', b'G', b'T'].contains(&base));
        }
    }

    /// Tests contig generation
    #[test]
    fn test_generate_contigs_denovo() {
        let config = ContigConfig {
            number: NonZeroU32::new(5).unwrap(),
            len_range: OrdPair::new(NonZeroU64::new(100).unwrap(), NonZeroU64::new(200).unwrap())
                .unwrap(),
        };
        let contigs = generate_contigs_denovo(config.number, config.len_range, &mut rand::rng());
        assert_eq!(contigs.len(), 5);
        for (i, contig) in contigs.iter().enumerate() {
            assert_eq!(contig.name, format!("contig_{}", i));
            assert!((100..=200).contains(&contig.seq.len()));
            for base in contig.seq.clone() {
                assert!([b'A', b'C', b'G', b'T'].contains(&base));
            }
        }
    }

    /// Tests read generation with desired properties.
    #[test]
    fn test_generate_reads_denovo() {
        let mut rng = rand::rng();
        let contigs = vec![
            Contig {
                name: "contig_0".to_string(),
                seq: b"ACGTACGTACGTACGTACGT".to_vec(),
            },
            Contig {
                name: "contig_1".to_string(),
                seq: b"TGCATGCATGCATGCATGCA".to_vec(),
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

        let reads = generate_reads_denovo(&contigs, config, "1", &mut rng).unwrap();
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
    fn test_full_bam_generation() {
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
        let bam_path = format!("/tmp/{}.bam", Uuid::new_v4());
        let fasta_path = format!("/tmp/{}.fa", Uuid::new_v4());
        assert!(run(config_json, &bam_path, &fasta_path).unwrap());
        assert!(Path::new(&bam_path).exists());
        assert!(Path::new(&fasta_path).exists());

        // read BAM file and check contig, read count
        let mut reader = bam::Reader::from_path(&bam_path).unwrap();
        let header = reader.header();
        assert_eq!(header.target_count(), 2);
        assert_eq!(reader.records().count(), 1000);

        // delete files
        std::fs::remove_file(bam_path).unwrap();
        std::fs::remove_file(fasta_path).unwrap();
    }

    /// Tests TempBamSimulation struct functionality
    #[test]
    fn test_temp_bam_simulation_struct() {
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
}
