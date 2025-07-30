use crate::{nanalogue_bam_reader, CurrRead, Error};
use rust_htslib::bam::Read;

pub fn run(bam_path: &str, read_id: &str) -> Result<bool, Error> {

    // open BAM file
    let mut bam = nanalogue_bam_reader(bam_path)?;

    // initialize output string
    let mut output_string = String::from("");

    // Go record by record in the BAM file,
    for r in bam.records() {
        // read records
        let mut curr_read_state = CurrRead::new();
        let record = r?;
        match curr_read_state.set_read_id(&record){
            Ok(true) => {
                if curr_read_state.get_read_id()? != read_id {
                    continue;
                } else {
                    // set the read state
                    curr_read_state.set_read_state(&record)?;
                    curr_read_state.set_seq_len(&record)?;
                    curr_read_state.set_align_len(&record)?;
                    curr_read_state.set_mod_data(&record, 128);
                    curr_read_state.set_contig_and_start(&record)?;
                }
            },
            Ok(false) => continue,
            Err(e) => Err(e)?,
        }
        output_string = output_string + &curr_read_state.to_string() + "\n";
    }

    if !output_string.is_empty() {
        print!("{output_string}");
    }

    Ok(true)
}
