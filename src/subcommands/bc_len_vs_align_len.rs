use crate::{nanalogue_bam_reader, CurrRead, ReadState, Error};
use csv::ReaderBuilder;
use rust_htslib::{bam::Read};
use serde::Deserialize;
use std::collections::HashMap;
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
    OnlyAlign { align_len: u64, seq_len: u64, mod_count: Option<String>, read_state: ReadState },
    OnlyBc(u64),
    BothAlignBc { align_len: u64, bc_len: u64, mod_count: Option<String>, read_state: ReadState },
}

// Implement a structure representing the above state
struct ReadLen {
    state: ReadLenState,
}

impl ReadLen {
    fn new_bc_len(l: u64) -> Self {
        Self {
            state: ReadLenState::OnlyBc(l),
        }
    }

    fn new_align_len(align_len: u64, seq_len: u64, mod_count: Option<String>, read_state: ReadState) -> Self {
        Self {
            state: ReadLenState::OnlyAlign { align_len, seq_len, mod_count, read_state },
        }
    }

    fn add_align_len(&mut self, align_len: u64, mod_count: Option<String>,
        read_state: ReadState) {
        match &self.state {
            ReadLenState::OnlyAlign {
                align_len: _,
                seq_len: _,
                mod_count: _,
                read_state: _,
            }
            | ReadLenState::BothAlignBc {
                align_len: _,
                bc_len: _,
                mod_count: _,
                read_state: _,
            } => {
                eprintln!("error, duplicates detected!");
                std::process::exit(1);
            },
            ReadLenState::OnlyBc(bl) => {
                self.state = ReadLenState::BothAlignBc {
                    align_len,
                    bc_len: *bl,
                    mod_count,
                    read_state
                };
            },
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
fn process_tsv(file_path: &str) -> Result<HashMap<String, ReadLen>, Error> {
    let mut data_map = HashMap::<String, ReadLen>::new();

    match file_path {
        "" => {},
        fp => {
            let file = File::open(fp)?;
            let mut rdr = ReaderBuilder::new()
                .comment(Some(b'#'))
                .delimiter(b'\t')
                .from_reader(file);

            for result in rdr.deserialize() {
                let record: TSVRecord = result?;
                if data_map.insert(record.read_id, ReadLen::new_bc_len(record.sequence_length_template)).is_some() {
                    Err(Error::InvalidDuplicates(file_path.to_string()))?;
                };
            }
        }
    };

    Ok(data_map)
}

pub fn run(bam_path: &str, seq_summ_path: &str, is_mod_count: bool) -> Result<bool, Error> {
    // read TSV file and convert into hashmap
    let mut data_map = process_tsv(seq_summ_path)?;

    // set up a flag to check if sequencing summary file has data
    let is_seq_summ_data: bool = !data_map.is_empty();

    // open BAM file
    let mut bam = nanalogue_bam_reader(bam_path);

    // Go record by record in the BAM file,
    // get the read id and the alignment length, and put it in the hashmap
    for r in bam.records() {

        // read records
        let record = r?;

        // get information of current read
        let mut curr_read_state = CurrRead::new();
        let qname :String = match curr_read_state.set_read_id(&record){
            Ok(false) | Err(_) => continue,
            Ok(true) => {
                let Some(ref v) = curr_read_state.read_id else {
                    continue;
                };
                v.to_string()
            },
        };
        curr_read_state.set_read_state(&record)?;
        let read_state = match curr_read_state.state {
            ReadState::PrimaryFwd | ReadState::PrimaryRev => curr_read_state.state,
            _ => continue,
        };
        let Some(align_len) = curr_read_state.set_align_len(&record)? else {
            continue;
        };
        let Some(seq_len) = curr_read_state.set_seq_len(&record)? else {
            continue;
        };

        // get modification information
        let mod_count: Option<String> = match is_mod_count {
            false => None,
            true => {
                curr_read_state.set_mod_data(&record, 128);
                let mut output_string = String::from("");
                match curr_read_state.mod_count_per_mod() {
                    None => {output_string += "0;"},
                    Some(v) => {
                        for (key,value) in v.into_iter(){
                            output_string = output_string + &format!("{}:{};",
                                match key {
                                    'A'..='Z' | 'a'..='z' => key.to_string(),
                                    _ => format!("{}", key as u32),
                                }, value).to_string();
                        }
                    }
                }
                output_string.pop();
                Some(output_string)
            }
        };

        // add data depending on whether an entry is already present
        // in the hashmap from the sequencing summary file
        data_map
            .entry(qname)
            .and_modify(|entry| entry.add_align_len(align_len, mod_count.clone(), read_state))
            .or_insert(ReadLen::new_align_len(align_len, seq_len, mod_count, read_state));
    }

    // set up an output header string
    let mut output_header = "# bam file: ".to_owned() + bam_path + "\n";
    if is_seq_summ_data {
        output_header = output_header + "# seq summ file: " + seq_summ_path + "\n"
    }
    output_header += "read_id\talign_length\tsequence_length_template\talignment_type";
    if is_mod_count {
        output_header += "\tmod_count";
    }

    // This apparently helps writing to the terminal faster,
    // according to https://rust-cli.github.io/book/tutorial/output.html
    let stdout = io::stdout();
    let mut handle = io::BufWriter::new(stdout);

    // print the output header
    writeln!(handle, "{output_header}")?;

    // print output tsv data
    // If both seq summ and BAM file are available, then the length in the seq
    // summ file takes precedence. If only the BAM file is available, then
    // we use the sequence length in the BAM file as the basecalled sequence length.
    for (key, val) in data_map.iter() {
        match (&val, is_seq_summ_data) {
            (
                ReadLen {
                    state:
                        ReadLenState::BothAlignBc {
                            align_len: al,
                            bc_len: bl,
                            mod_count: Some(ml),
                            read_state: rs,
                        },
                },
                true,
            ) => writeln!(handle, "{key}\t{al}\t{bl}\t{rs}\t{ml}")?, 
            (
                ReadLen {
                    state:
                        ReadLenState::OnlyAlign {
                            align_len: al,
                            seq_len: sl,
                            mod_count: Some(ml),
                            read_state: rs,
                        },
                },
                false,
            ) => writeln!(handle, "{key}\t{al}\t{sl}\t{rs}\t{ml}")?,
            (
                ReadLen {
                    state:
                        ReadLenState::BothAlignBc {
                            align_len: al,
                            bc_len: bl,
                            mod_count: None,
                            read_state: rs,
                        },
                },
                true,
            ) => writeln!(handle, "{key}\t{al}\t{bl}\t{rs}")?,
            (
                ReadLen {
                    state:
                        ReadLenState::OnlyAlign {
                            align_len: al,
                            seq_len: sl,
                            mod_count: None,
                            read_state: rs,
                        },
                },
                false,
            ) => writeln!(handle, "{key}\t{al}\t{sl}\t{rs}")?,
            (
                ReadLen {
                    state: ReadLenState::OnlyBc(_)
                },
                _,
            ) |
            (
                ReadLen {
                    state: ReadLenState::OnlyAlign {
                        align_len: _,
                        seq_len: _,
                        mod_count: _,
                        read_state: _,
                    }
                },
                true,
            ) => {},
            (
                ReadLen {
                    state: ReadLenState::BothAlignBc {
                        align_len: _,
                        bc_len: _,
                        mod_count: _,
                        read_state: _,
                    }
                },
                false,
            ) => Err(Error::InvalidState("invalid state while writing output".to_string()))?,
        }
    }

    Ok(true)
}
