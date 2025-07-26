use crate::{CurrRead, ReadState, ReadTransition};
use crate::{nanalogue_bam_reader, nanalogue_mm_ml_parser};
use fibertools_rs::utils::basemods::BaseMods;
use rust_htslib::{bam::Read, bam::ext::BamRecordExtensions};
use std::error::Error;

pub fn run(bam_path: &str, read_id: &str) -> Result<(), Box<dyn Error>> {
    // open BAM file
    let mut bam = nanalogue_bam_reader(bam_path);

    let mut output_string = String::from("");

    // Go record by record in the BAM file,
    for r in bam.records() {
        // read records
        let record = match r {
            Ok(v) => {
                // get read id
                let qname: String = str::from_utf8(v.qname()).to_string().unwrap_or_else(|e| {
                    eprintln!("Invalid UTF-8 sequence: {e}");
                    std::process::exit(1);
                });
                if qname == read_id {
                    v
                } else {
                    continue;
                }
            },
            Err(e) => {
                eprintln!("Some error while reading records {e}");
                std::process::exit(1)
            },
        };

        // set the read state using the type of alignment
        let mut curr_read_state = CurrRead::new();
        curr_read_state.read_id = Some(read_id.to_string());
        if record.is_unmapped() {
            curr_read_state.transition(ReadTransition::PrimaryToUnmapped);
        }
        if record.is_secondary() {
            curr_read_state.transition(ReadTransition::PrimaryToSecondary);
        }
        if record.is_supplementary() {
            curr_read_state.transition(ReadTransition::PrimaryToSupplementary);
        }

        // get length of sequence
        let seq_len: u64 = record.seq_len().try_into().unwrap_or_else( |_| {
            eprintln!("Error while getting sequencing length!");
            std::process::exit(1);
        });
        curr_read_state.seq_len = Some(seq_len);

        // get length of alignment
        match curr_read_state.state {
            ReadState::Unmapped => {}
            _ => {
                let align_len: u64 = (record.reference_end() - record.pos()).try_into()
                    .unwrap_or_else(|_| {
                        eprintln!("Problem getting alignment length");
                        std::process::exit(1);
                    });
                curr_read_state.align_len = Some(align_len);
            }
        };

        // get modification information
        let BaseMods { base_mods: v } = nanalogue_mm_ml_parser(&record, 128);
        curr_read_state.mods = Some(v);

        output_string = output_string + &curr_read_state.to_string() + "\n";
    }

    if !output_string.is_empty() {
        println!("{}", CurrRead::header_string());
        print!("{output_string}");
    }

    Ok(())
}
