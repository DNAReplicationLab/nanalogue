//! # Write Simulated Mod BAM
//! Generates simulated BAM files with base modifications for testing purposes.
//! Accepts JSON configuration to produce both BAM and FASTA reference files.
//!
//! ## Example Usage
//!
//! ```no_run
//! use nanalogue_core::write_simulated_mod_bam::run;
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
use rand::{Rng, random};
use rust_htslib::bam;
use rust_htslib::bam::Header;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::Write;
use std::num::{NonZeroU32, NonZeroU64};
use std::ops::RangeInclusive;
use uuid::Uuid;

/// Main configuration struct for simulation
#[derive(Debug, Clone, Serialize, Deserialize)]
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
    /// Use implicit notation, i.e. "." or "" to say missing bases are unmodified.
    /// Defaults to using "?" i.e. missing bases are missing.
    pub implicit: bool,
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
            implicit: false,
        }
    }
}

impl Default for SimulationConfig {
    fn default() -> Self {
        Self {
            contigs: ContigConfig::default(),
            reads: Vec::new(),
        }
    }
}

/// Generates random DNA sequence of given length
fn generate_random_dna_sequence<R: Rng>(length: NonZeroU64, rng: &mut R) -> Vec<u8> {
    const DNA_BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
    (0..length.get())
        .map(|_| DNA_BASES[rng.random_range(0..4)])
        .collect()
}

/// Generates contigs according to configuration
fn generate_contigs<R: Rng>(
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

/// Writes contigs to a FASTA file
fn write_fasta(contigs: &[Contig], output_path: &str) -> Result<(), Error> {
    let mut file = File::create(output_path)?;
    for contig in contigs {
        writeln!(file, ">{}", contig.name)?;
        writeln!(file, "{}", String::from_utf8_lossy(&contig.seq))?;
    }
    Ok(())
}

/// Generates reads that align to contigs
fn generate_reads<R: Rng>(
    contigs: &[Contig],
    config: &ReadConfig,
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
            record.set_flags({
                let random_state: ReadState = random();
                u16::try_from(random_state)?
            });
            if !record.is_unmapped() {
                let cigar = CigarString(vec![Cigar::Match(read_len as u32)]);
                record.set(&qname, Some(&cigar), &read_seq, &qual);
                record.set_mapq(mapq);
                record.set_tid(contig_idx as i32);
                record.set_pos(start_pos as i64);
            } else {
                record.set(&qname, None, &read_seq, &qual);
            }
            record.push_aux(b"RG", Aux::U8(1))?;
            record
        };
        reads.push(record);
    }

    Ok(reads)
}

/// Writes BAM file with reads
fn write_bam(reads: &[bam::Record], contigs: &[Contig], output_path: &str) -> Result<(), Error> {
    let header = {
        let mut header = Header::new();
        for contig in contigs {
            let _ = header.push_record(
                &bam::header::HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", &contig.name)
                    .push_tag(b"LN", contig.seq.len()),
            );
        }
        let _ = header.push_record(
            &bam::header::HeaderRecord::new(b"RG")
                .push_tag(b"ID", 1)
                .push_tag(b"PL", "ULTIMA")
                .push_tag(b"LB", "blank")
                .push_tag(b"SM", "blank")
                .push_tag(b"PU", "blank"),
        );
        header
    };

    let mut writer = bam::Writer::from_path(output_path, &header, bam::Format::Bam)?;

    for read in reads {
        writer.write(read)?;
    }

    Ok(())
}

/// Main function to generate simulated BAM file
///
/// # Example
///
/// ```no_run
/// use nanalogue_core::write_simulated_mod_bam::run;
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
pub fn run(
    config_json: &str,
    bam_output_path: &str,
    fasta_output_path: &str,
) -> Result<bool, Error> {
    let config: SimulationConfig = serde_json::from_str(config_json)?;

    let mut rng = rand::rng();

    let contigs = generate_contigs(config.contigs.number, config.contigs.len_range, &mut rng);
    let reads = {
        let mut temp_reads = Vec::new();
        for k in config.reads {
            temp_reads.append(&mut generate_reads(&contigs, &k, &mut rng)?);
        }
        temp_reads
    };

    write_fasta(&contigs, fasta_output_path)?;
    write_bam(&reads, &contigs, bam_output_path)?;

    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::Read;

    #[test]
    fn test_generate_random_dna_sequence() {
        let seq = generate_random_dna_sequence(NonZeroU64::new(100).unwrap(), &mut rand::rng());
        assert_eq!(seq.len(), 100);
        for base in seq {
            assert!([b'A', b'C', b'G', b'T'].contains(&base));
        }
    }

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
            assert!(contig.seq.len() >= 100 && contig.seq.len() <= 200);
            for base in contig.seq.clone() {
                assert!([b'A', b'C', b'G', b'T'].contains(&base));
            }
        }
    }

    #[test]
    fn test_write_fasta() {
        let contigs = vec![
            Contig {
                name: "test_contig_0".to_string(),
                seq: b"ACGT".to_vec(),
            },
            Contig {
                name: "test_contig_1".to_string(),
                seq: b"TGCA".to_vec(),
            },
        ];

        let temp_path = format!("/tmp/{}.fa", Uuid::new_v4());
        write_fasta(&contigs, &temp_path).expect("no error");

        let content = std::fs::read_to_string(temp_path.clone()).expect("no error");
        assert_eq!(content, ">test_contig_0\nACGT\n>test_contig_1\nTGCA\n");

        std::fs::remove_file(&temp_path).expect("no error");
    }

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
            base_qual_range: OrdPair::new(10, 20).unwrap(),
            len_range: OrdPair::new(F32Bw0and1::new(0.2).unwrap(), F32Bw0and1::new(0.8).unwrap())
                .unwrap(),
            mods: vec![],
        };

        let reads = generate_reads(&contigs, &config, &mut rng).unwrap();
        assert!(reads.len() <= 10);

        for read in &reads {
            if !read.is_unmapped() {
                assert!(read.mapq() >= 10 && read.mapq() <= 20);
                assert!(read.tid() >= 0 && read.tid() < 2);
            }
        }
    }

    #[test]
    fn test_full_bam_generation() {
        let config_json = r#"{
            "contigs": {
                "number": 2,
                "len_range": [100, 200]
            },
            "reads": [{
                "number": 5,
                "mapq_range": [10, 20],
                "base_qual_range": [10, 20],
                "len_range": [0.1, 0.8],
                "mods": []
            }]
        }"#;

        let bam_path = format!("/tmp/{}.bam", Uuid::new_v4());
        let fasta_path = format!("/tmp/{}.fa", Uuid::new_v4());

        assert!(run(config_json, &bam_path, &fasta_path).unwrap());

        assert!(std::path::Path::new(&bam_path).exists());
        assert!(std::path::Path::new(&fasta_path).exists());

        let mut reader = bam::Reader::from_path(&bam_path).unwrap();
        let header = reader.header();
        assert_eq!(header.target_count(), 2);

        let read_count = {
            let mut cnt = 0;
            for _ in reader.records() {
                cnt += 1;
            }
            cnt
        };
        assert!(read_count <= 5);

        std::fs::remove_file(bam_path).unwrap();
        std::fs::remove_file(fasta_path).unwrap();
    }
}
