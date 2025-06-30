use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use csv::ReaderBuilder;
use serde::Deserialize;
use std::str;
use rust_htslib::{bam, bam::Read, bam::ext::BamRecordExtensions};
use clap::Parser;
use std::io::{self, Write};


#[derive(Parser)]
struct Cli {
    command: String,
    bam_filename: String,
    seq_summ_file: String
}

/// Represents a record with the columns of interest from the TSV file.
#[derive(Debug, Deserialize)]
struct Record {
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
/// a tuple containing `sequence_length_template` as the value, or an error.
fn process_tsv(file_path: &str) -> Result<HashMap<String, (u64, u64)>, Box<dyn Error>> {
    let file = File::open(file_path)?;
    let mut rdr = ReaderBuilder::new()
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .from_reader(file);

    let mut data_map = HashMap::<String, (u64, u64)>::new();

    for result in rdr.deserialize() {
        let record: Record = result?;
        data_map.insert(record.read_id, (record.sequence_length_template, 0));
    }

    Ok(data_map)
}

fn main() {

    // parse arguments
    let args = Cli::parse();

    // only one command implemented for now
    match args.command.as_ref() {
        "bcLen_vs_alignLen" => "bcLen_vs_alignLen",
        _ => {
            eprintln!("Only bcLen_vs_alignLen implemented for now!");
            std::process::exit(1);
        },
    };

    // set up an output header string
    let mut output_header = String::from("read_id\talign_length\tsequence_length_template");

    // read TSV file and convert into hashmap
    let mut data_map = match process_tsv(&args.seq_summ_file) {
        Ok(v) => {
            output_header = "# seq summ file: ".to_owned() + &args.seq_summ_file + "\n" + &output_header;
            v
        }
        Err(e) => {
            eprintln!("An error occurred: {}", e);
            std::process::exit(1);
        }
    };

    // add BAM file to output header
    output_header = "# bam file: ".to_owned() + &args.bam_filename + "\n" + &output_header;

    // open BAM file
    let mut bam = match bam::Reader::from_path(args.bam_filename) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("Problem opening file, error: {}", e);
            std::process::exit(1)
        },
    };

    // Go record by record in the BAM file,
    // for primary alignments, get the read id and the alignment length,
    // and put it in the hashmap
    for r in bam.records() {

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

        let qname :String = match str::from_utf8(record.qname()) {
            Ok(v) => v.to_string(),
            Err(e) => {
                eprintln!("Invalid UTF-8 sequence: {}", e);
                std::process::exit(1)
            },
        };

        let align_len :u64 = match (record.reference_end() - record.pos()).try_into() {
            Ok(v) => v,
            Err(_) => {
                eprintln!("Problem getting alignment length");
                std::process::exit(1);
            },
        };

        data_map.entry(qname).and_modify(|entry| *entry = (entry.0, align_len));

    };

    // This apparently helps writing to the terminal faster,
    // according to https://rust-cli.github.io/book/tutorial/output.html
    let stdout = io::stdout(); 
    let mut handle = io::BufWriter::new(stdout); 

    // print the output
    match writeln!(handle, "{}", output_header){
        Ok(v) => v,
        Err(_) => {
            eprintln!("Error while writing output!");
            std::process::exit(1);
        },
    };
    for (key, val) in data_map.iter() {
        if val.1 > 0 {
            match writeln!(handle, "{key}\t{0}\t{1}", val.1, val.0){
                Ok(v) => v,
                Err(_) => {
                    eprintln!("Error while writing output!");
                    std::process::exit(1);
                },
            };
        }
    }

}
