use crate::{CurrRead};
use crate::{nanalogue_bam_reader, nanalogue_mm_ml_parser};
use fibertools_rs::utils::basemods::BaseMods;
use rust_htslib::{bam::Read};
use std::error::Error;

pub fn run(bam_path: &str, read_id: &str) -> Result<(), Box<dyn Error>> {
    // open BAM file
    let mut bam = nanalogue_bam_reader(bam_path);

    let mut output_string = String::from("");

    // Go record by record in the BAM file,
    for r in bam.records() {
        // read records
        let mut curr_read_state = CurrRead::new();
        match r {
            Ok(v) => {
                match curr_read_state.set_read_id(&v){
                    Ok(true) => {
                        if let Some(ref w) = curr_read_state.read_id {
                            if w != read_id {
                                continue;
                            } else {
                                // set the read state
                                curr_read_state.set_read_state(&v).expect("failed to set read state!");
                                curr_read_state.set_seq_len(&v).expect("failed to set sequence length!");
                                curr_read_state.set_align_len(&v).expect("failed to set alignment length!");
                                let BaseMods { base_mods: x } = nanalogue_mm_ml_parser(&v, 128);
                                curr_read_state.mods = Some(x);
                            }
                        }
                    },
                    Ok(false) => continue,
                    Err(e) => {
                        eprintln!("Some error while reading records {e}");
                        std::process::exit(1)
                    }
                }
            },
            Err(e) => {
                eprintln!("Some error while reading records {e}");
                std::process::exit(1)
            },
        };
        output_string = output_string + &curr_read_state.to_string() + "\n";
    }

    if !output_string.is_empty() {
        println!("{}", CurrRead::header_string());
        print!("{output_string}");
    }

    Ok(())
}
