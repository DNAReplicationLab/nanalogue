use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use csv::ReaderBuilder;
use serde::Deserialize;
use std::str;
use rust_htslib::{bam, bam::Read, bam::ext::BamRecordExtensions};
use std::io::{self, Write};


// Enum has three states:
// - only BAM file alignment information is available
// - only sequencing summary file information is available
// - both are available, in which case we retain only the sequence
//   length from the sequencing summary file and the alignment length
//   from the BAM file, discarding the sequence length from the BAM file.
enum ReadLenState {
    OnlyAlignLen {
        align_len: u64,
        seq_len: u64,
    },
    OnlyBcLen(u64),
    BothAlignBcLen {
        align_len: u64,
        bc_len: u64,
    },
}

// Implement a structure representing the above state
struct ReadLen {
    state: ReadLenState,
}

impl ReadLen {
    fn add_align_len(&mut self, align_len: u64) {
        match &self.state {
            ReadLenState::OnlyAlignLen { align_len: _, seq_len: _ } | 
                ReadLenState::BothAlignBcLen { align_len: _, bc_len: _ } => {
                eprintln!("Alignment statistics duplicate detected!");
                std::process::exit(1);
            },
            ReadLenState::OnlyBcLen(bl) => {
                self.state = ReadLenState::BothAlignBcLen{
                    align_len: align_len,
                    bc_len: *bl,
                };
            },
        }
    }
    fn add_bc_len(&mut self, bc_len: u64) {
        match &self.state {
            ReadLenState::OnlyBcLen(_) | ReadLenState::BothAlignBcLen { align_len: _, bc_len: _ } => {
                eprintln!("Basecalled length duplicate detected!");
                std::process::exit(1);
            },
            ReadLenState::OnlyAlignLen {align_len: al, seq_len: _} => {
                self.state = ReadLenState::BothAlignBcLen {
                    align_len: *al,
                    bc_len: bc_len,
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
fn process_tsv(file_path: &str) -> Result<HashMap<String, ReadLen>, Box<dyn Error>> {

    let mut data_map = HashMap::<String, ReadLen>::new();

    match file_path{
        "" => {},
        fp => { 
            let file = File::open(fp)?;

            let mut rdr = ReaderBuilder::new()
                .comment(Some(b'#'))
                .delimiter(b'\t')
                .from_reader(file);


            for result in rdr.deserialize() {
                let record: TSVRecord = result?;
                data_map.entry(record.read_id)
                    .and_modify( |entry| entry.add_bc_len(record.sequence_length_template))
                    .or_insert( ReadLen { state: ReadLenState::OnlyBcLen(record.sequence_length_template) });
            }
        },
    };

    Ok(data_map)
}

pub fn run(bam_path: &str, seq_summ_path: &str) -> Result<(), Box<dyn Error>> {


    // read TSV file and convert into hashmap
    let mut data_map = match process_tsv(seq_summ_path) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("An error occurred: {}", e);
            std::process::exit(1);
        },
    };

    // set up a flag to check if sequencing summary file has data
    let is_seq_summ_data :bool = data_map.len() > 0;

    // open BAM file
    let mut bam = match bam::Reader::from_path(bam_path) {
        Ok(v) => {
            v
        },
        Err(e) => {
            eprintln!("Problem opening file, error: {}", e);
            std::process::exit(1)
        },
    };

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
            },
            Err(e) => {
                eprintln!("Some error while reading records {}", e);
                std::process::exit(1)
            },
        };

        // get read id
        let qname :String = match str::from_utf8(record.qname()) {
            Ok(v) => v.to_string(),
            Err(e) => {
                eprintln!("Invalid UTF-8 sequence: {}", e);
                std::process::exit(1)
            },
        };

        // get length of alignment
        let align_len :u64 = match (record.reference_end() - record.pos()).try_into() {
            Ok(v) => v,
            Err(_) => {
                eprintln!("Problem getting alignment length");
                std::process::exit(1);
            },
        };

        // get length of sequence
        let seq_len :u64 = match record.seq_len().try_into() {
            Ok(v) => v,
            Err(_) => {
                eprintln!("Error while getting sequencing length!");
                std::process::exit(1);
            },
        };

        // add data depending on whether an entry is already present
        // in the hashmap from the sequencing summary file
        data_map.entry(qname)
            .and_modify( |entry| entry.add_align_len(align_len))
            .or_insert( ReadLen { 
                state: ReadLenState::OnlyAlignLen {
                    align_len: align_len,
                    seq_len: seq_len
                }});

    };

    // set up an output header string
    let mut output_header = "# bam file: ".to_owned() + bam_path + "\n";
    match is_seq_summ_data {
        true => output_header = output_header + "# seq summ file: " + seq_summ_path + "\n",
        false => {},
    }
    output_header = output_header + "read_id\talign_length\tsequence_length_template";

    // This apparently helps writing to the terminal faster,
    // according to https://rust-cli.github.io/book/tutorial/output.html
    let stdout = io::stdout(); 
    let mut handle = io::BufWriter::new(stdout); 

    // print the output header
    match writeln!(handle, "{}", output_header){
        Ok(_) => {},
        Err(_) => {
            eprintln!("Error while writing output!");
            std::process::exit(1);
        },
    };

    // print output tsv data
    // If both seq summ and BAM file are available, then the length in the seq
    // summ file takes precedence. If only the BAM file is available, then
    // we use the sequence length in the BAM file as the basecalled sequence length.
    for (key, val) in data_map.iter() {
        match (&val.state, is_seq_summ_data) {
            (ReadLenState::BothAlignBcLen {align_len: al, bc_len: bl}, true) => 
                match writeln!(handle, "{key}\t{0}\t{1}", al, bl){
                    Ok(_) => {},
                    Err(_) => {
                        eprintln!("Error while writing output!");
                        std::process::exit(1);
                    },
                },
            (ReadLenState::OnlyAlignLen {align_len: al, seq_len: sl }, false) => 
                match writeln!(handle, "{key}\t{0}\t{1}", al, sl){
                    Ok(_) => {},
                    Err(_) => {
                        eprintln!("Error while writing output!");
                        std::process::exit(1);
                    },
                },
            _ => {},
        }
    };

    Ok(())

}
