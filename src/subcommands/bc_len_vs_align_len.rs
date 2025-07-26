use crate::{nanalogue_bam_reader, nanalogue_mm_ml_parser};
use csv::ReaderBuilder;
use fibertools_rs::utils::basemods::BaseMods;
use rust_htslib::{bam::Read, bam::ext::BamRecordExtensions};
use serde::Deserialize;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{self, Write};
use std::str;

// Enum has three states:
// - only BAM file alignment information is available
// - only sequencing summary file information is available
// - both are available, in which case we retain only the sequence
//   length from the sequencing summary file and the alignment length
//   from the BAM file, discarding the sequence length from the BAM file.
enum ReadLenState {
    OnlyAlign { align_len: u64, seq_len: u64 },
    OnlyBc(u64),
    BothAlignBc { align_len: u64, bc_len: u64 },
}

// Implement a structure representing the above state
// and other information in the read
struct ReadLen {
    state: ReadLenState,
    mod_count: Option<u64>,
}

impl ReadLen {
    fn new_bc_len(l: u64) -> Self {
        Self {
            state: ReadLenState::OnlyBc(l),
            mod_count: None,
        }
    }
    fn new_align_len(align_len: u64, seq_len: u64, mod_count: Option<u64>) -> Self {
        Self {
            state: ReadLenState::OnlyAlign { align_len, seq_len },
            mod_count,
        }
    }
    fn add_align_len(&mut self, align_len: u64, mod_count: Option<u64>) {
        match &self.state {
            ReadLenState::OnlyAlign {
                align_len: _,
                seq_len: _,
            }
            | ReadLenState::BothAlignBc {
                align_len: _,
                bc_len: _,
            } => {
                eprintln!("Alignment statistics duplicate detected!");
                std::process::exit(1);
            }
            ReadLenState::OnlyBc(bl) => {
                self.state = ReadLenState::BothAlignBc {
                    align_len,
                    bc_len: *bl,
                };
            }
        }
        self.mod_count = mod_count;
    }

    fn add_bc_len(&mut self, bc_len: u64) {
        match &self.state {
            ReadLenState::OnlyBc(_)
            | ReadLenState::BothAlignBc {
                align_len: _,
                bc_len: _,
            } => {
                eprintln!("Basecalled length duplicate detected!");
                std::process::exit(1);
            }
            ReadLenState::OnlyAlign {
                align_len: al,
                seq_len: _,
            } => {
                self.state = ReadLenState::BothAlignBc {
                    align_len: *al,
                    bc_len,
                };
            }
        }
    }
}

/// Represents a record with the columns of interest from the TSV file.
#[derive(Debug, Deserialize)]
struct TSVRecord {
    read_id: String,
    sequence_length_template: u64,
}

/// Opens a TSV file, extracts 'read_id' and 'sequence_length_template' columns,
/// and builds a HashMap.
///
/// # Arguments
///
/// * `file_path` - A string slice that holds the path to the TSV file.
///
/// # Returns
///
/// A `Result` which is either a `HashMap` with `read_id` as the key and
/// a `Read` as the value, or an error.
fn process_tsv(file_path: &str) -> Result<HashMap<String, ReadLen>, Box<dyn Error>> {
    let mut data_map = HashMap::<String, ReadLen>::new();

    match file_path {
        "" => {}
        fp => {
            let file = File::open(fp)?;

            let mut rdr = ReaderBuilder::new()
                .comment(Some(b'#'))
                .delimiter(b'\t')
                .from_reader(file);

            for result in rdr.deserialize() {
                let record: TSVRecord = result?;
                data_map
                    .entry(record.read_id)
                    .and_modify(|entry| entry.add_bc_len(record.sequence_length_template))
                    .or_insert(ReadLen::new_bc_len(record.sequence_length_template));
            }
        }
    };

    Ok(data_map)
}

pub fn run(bam_path: &str, seq_summ_path: &str, is_mod_count: bool) -> Result<(), Box<dyn Error>> {
    // read TSV file and convert into hashmap
    let mut data_map = match process_tsv(seq_summ_path) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("An error occurred: {e}");
            std::process::exit(1);
        }
    };

    // set up a flag to check if sequencing summary file has data
    let is_seq_summ_data: bool = !data_map.is_empty();

    // open BAM file
    let mut bam = nanalogue_bam_reader(bam_path);

    // Go record by record in the BAM file,
    // for primary alignments, get the read id and the alignment length,
    // and put it in the hashmap
    for r in bam.records() {
        // read primary records
        let record = match r {
            Ok(v) => {
                if v.is_unmapped() || v.is_secondary() || v.is_supplementary() {
                    continue;
                } else {
                    v
                }
            }
            Err(e) => {
                eprintln!("Some error while reading records {e}");
                std::process::exit(1)
            }
        };

        // get read id
        let qname: String = match str::from_utf8(record.qname()) {
            Ok(v) => v.to_string(),
            Err(e) => {
                eprintln!("Invalid UTF-8 sequence: {e}");
                std::process::exit(1)
            }
        };

        // get length of alignment
        let align_len: u64 = (record.reference_end() - record.pos()).try_into().unwrap_or_else(|_| {
                eprintln!("Problem getting alignment length");
                std::process::exit(1);
        });

        // get length of sequence
        let seq_len: u64 = record.seq_len().try_into().unwrap_or_else(|_| {
                eprintln!("Error while getting sequencing length!");
                std::process::exit(1);
        });

        // get modification information
        let mod_count: Option<u64> = match is_mod_count {
            false => None,
            true => {
                let BaseMods { base_mods: v } = nanalogue_mm_ml_parser(&record, 128);
                Some({
                    let mut cnt: u64 = 0;
                    for k in v {
                        cnt += k.ranges.qual.len() as u64;
                    }
                    cnt
                })
            }
        };

        // add data depending on whether an entry is already present
        // in the hashmap from the sequencing summary file
        data_map
            .entry(qname)
            .and_modify(|entry| entry.add_align_len(align_len, mod_count))
            .or_insert(ReadLen::new_align_len(align_len, seq_len, mod_count));
    }

    // set up an output header string
    let mut output_header = "# bam file: ".to_owned() + bam_path + "\n";
    if is_seq_summ_data {
        output_header = output_header + "# seq summ file: " + seq_summ_path + "\n"
    }
    if is_mod_count {
        output_header += "read_id\talign_length\tsequence_length_template\tmod_count";
    } else {
        output_header += "read_id\talign_length\tsequence_length_template";
    }

    // This apparently helps writing to the terminal faster,
    // according to https://rust-cli.github.io/book/tutorial/output.html
    let stdout = io::stdout();
    let mut handle = io::BufWriter::new(stdout);

    // print the output header
    match writeln!(handle, "{output_header}") {
        Ok(_) => {}
        Err(_) => {
            eprintln!("Error while writing output!");
            std::process::exit(1);
        }
    };

    // print output tsv data
    // If both seq summ and BAM file are available, then the length in the seq
    // summ file takes precedence. If only the BAM file is available, then
    // we use the sequence length in the BAM file as the basecalled sequence length.
    for (key, val) in data_map.iter() {
        match (&val, is_seq_summ_data, is_mod_count) {
            (
                ReadLen {
                    state:
                        ReadLenState::BothAlignBc {
                            align_len: al,
                            bc_len: bl,
                        },
                    mod_count: Some(ml),
                },
                true,
                true,
            ) => match writeln!(handle, "{key}\t{al}\t{bl}\t{ml}") {
                Ok(_) => {}
                Err(_) => {
                    eprintln!("Error while writing output!");
                    std::process::exit(1);
                }
            },
            (
                ReadLen {
                    state:
                        ReadLenState::OnlyAlign {
                            align_len: al,
                            seq_len: sl,
                        },
                    mod_count: Some(ml),
                },
                false,
                true,
            ) => match writeln!(handle, "{key}\t{al}\t{sl}\t{ml}") {
                Ok(_) => {}
                Err(_) => {
                    eprintln!("Error while writing output!");
                    std::process::exit(1);
                }
            },
            (
                ReadLen {
                    state:
                        ReadLenState::BothAlignBc {
                            align_len: al,
                            bc_len: bl,
                        },
                    mod_count: _,
                },
                true,
                false,
            ) => match writeln!(handle, "{key}\t{al}\t{bl}") {
                Ok(_) => {}
                Err(_) => {
                    eprintln!("Error while writing output!");
                    std::process::exit(1);
                }
            },
            (
                ReadLen {
                    state:
                        ReadLenState::OnlyAlign {
                            align_len: al,
                            seq_len: sl,
                        },
                    mod_count: _,
                },
                false,
                false,
            ) => match writeln!(handle, "{key}\t{al}\t{sl}") {
                Ok(_) => {}
                Err(_) => {
                    eprintln!("Error while writing output!");
                    std::process::exit(1);
                }
            },
            _ => {}
        }
    }

    Ok(())
}
