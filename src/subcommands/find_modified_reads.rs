use crate::{nanalogue_bam_reader, CurrRead, Error, process_mod_type, 
    F32Bw0and1, OrdPair, InputBam};
use rust_htslib::bam::Read;
use std::num::NonZeroU32;

pub fn run(bam_options: &mut InputBam, tag: &str, win: NonZeroU32, slide: NonZeroU32, 
    dens_limits: OrdPair<F32Bw0and1>, invert: bool) -> Result<bool, Error> {

    // open BAM file
    let mut bam = nanalogue_bam_reader(bam_options)?;

    // prepare the tag
    let tag_char = process_mod_type(tag)?;

    // get density limits
    let dens_min = dens_limits.get_low();
    let dens_max = dens_limits.get_high();

    // prepare output string
    let mut output_string = String::from("");

    // Go record by record in the BAM file,
    for r in bam.records() {
        // read records
        let mut curr_read_state = CurrRead::new();
        let record = r?;
        curr_read_state.set_read_id(&record)?;

        // set the modified read state
        curr_read_state.set_mod_data_one_tag(&record, 0, tag_char);

        // catch if one window meets our criterion,
        // and react accordingly using invert's state
        if match curr_read_state.windowed_mod_data(win.get().try_into()?, slide.get().try_into()?, tag_char)?{
            Some(v) => ! (v.iter().any(|k| *k > dens_max || *k < dens_min) ^ invert),
            None => false,
        } {
            output_string = output_string + curr_read_state.get_read_id()? + "\n";
        }
    }

    if !output_string.is_empty() {
        println!("# input bam {}", bam_options.bam_path);
        print!("{output_string}");
    }

    Ok(true)
}
