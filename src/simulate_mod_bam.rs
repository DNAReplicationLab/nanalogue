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
//!     "mods": []
//!   }]
//! }"#;
//!
//! run(
//!     config_json,
//!     "output.bam",
//!     "reference.fasta"
//! ).unwrap();
//! ```

use crate::{Error, F32Bw0and1, ModChar, OrdPair, ReadState};
use crate::{write_denovo_bam, write_fasta};
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

/// Generates contigs according to configuration
///
/// ```
/// use std::num::{NonZeroU32, NonZeroU64};
/// use nanalogue_core::OrdPair;
/// use nanalogue_core::simulate_mod_bam::generate_contigs;
/// use rand::Rng;
///
/// let mut rng = rand::rng();
/// let contigs = generate_contigs(
///     NonZeroU32::new(3).unwrap(),
///     OrdPair::new(NonZeroU64::new(100).unwrap(), NonZeroU64::new(200).unwrap()).unwrap(),
///     &mut rng
/// );
/// assert_eq!(contigs.len(), 3);
/// assert!(contigs.iter().all(|c| (100..=200).contains(&c.seq.len())));
/// ```
pub fn generate_contigs<R: Rng>(
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
/// use nanalogue_core::simulate_mod_bam::{Contig, ReadConfig, generate_reads};
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
///     mods: vec![],
/// };
/// let mut rng = rand::rng();
/// let reads = generate_reads(&contigs, config, "RG1", &mut rng).unwrap();
/// assert_eq!(reads.len(), 10);
/// ```
pub fn generate_reads<R: Rng>(
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
        let end_pos = (start_pos + read_len) as usize;
        let read_seq = contig.seq[start_pos as usize..end_pos].to_vec();

        // Generate quality scores
        let qual: Vec<u8> = (0..read_len)
            .map(|_| rng.random_range(RangeInclusive::from(config.base_qual_range)))
            .collect();

        // Generate MAPQ
        let mapq = rng.random_range(RangeInclusive::from(config.mapq_range));

        // Create BAM record
        let record = {
            let mut record = bam::Record::new();
            let qname = format!("{}", Uuid::new_v4()).into_bytes();
            record.unset_flags();
            let random_state: ReadState = random();
            record.set_flags(u16::from(random_state));
            match random_state {
                ReadState::Unmapped => {
                    record.set(&qname, None, &read_seq, &qual);
                    record.set_mapq(255);
                    record.set_tid(-1);
                    record.set_pos(-1);
                }
                _ => {
                    let cigar = CigarString(vec![Cigar::Match(read_len as u32)]);
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
///     "mods": []
///   }]
/// }"#;
///
/// run(config_json, "output.bam", "reference.fasta").unwrap();
/// ```
pub fn run<F>(config_json: &str, bam_output_path: &F, fasta_output_path: &F) -> Result<bool, Error>
where
    F: AsRef<Path> + ?Sized,
{
    let config: SimulationConfig = serde_json::from_str(config_json)?;

    let mut rng = rand::rng();

    let contigs = generate_contigs(config.contigs.number, config.contigs.len_range, &mut rng);
    let read_groups: Vec<String> = (0..config.reads.len()).map(|k| k.to_string()).collect();
    let reads = {
        let mut temp_reads = Vec::new();
        for k in config.reads.into_iter().zip(read_groups.clone()) {
            temp_reads.append(&mut generate_reads(&contigs, k.0, &k.1, &mut rng)?);
        }
        temp_reads.sort_by_key(|k| (k.is_unmapped(), k.tid(), k.pos(), k.is_reverse()));
        temp_reads
    };

    write_denovo_bam(
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
    fn test_generate_contigs() {
        let config = ContigConfig {
            number: NonZeroU32::new(5).unwrap(),
            len_range: OrdPair::new(NonZeroU64::new(100).unwrap(), NonZeroU64::new(200).unwrap())
                .unwrap(),
        };
        let contigs = generate_contigs(config.number, config.len_range, &mut rand::rng());
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
    fn test_generate_reads() {
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
            mods: vec![],
        };

        let reads = generate_reads(&contigs, config, "1", &mut rng).unwrap();
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
}
