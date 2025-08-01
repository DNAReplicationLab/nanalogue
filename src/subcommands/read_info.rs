use crate::{nanalogue_bam_reader, CurrRead, Error, InputBam};
use rust_htslib::bam::{Read,Record};

pub fn run(bam_options: &mut InputBam, read_id: &str) -> Result<bool, Error> {

    // open BAM file
    let mut bam = nanalogue_bam_reader(bam_options)?;

    // convert read id into bytes
    let read_id_bytes = read_id.as_bytes();

    // initialize output string
    let mut output_string = String::from("");

    // Go record by record in the BAM file,
    // and collect entries that match our read id
    for k in bam.records().filter(|r|{
        match r {
            Ok(v) => v.qname() == read_id_bytes,
            Err(_) => true,
        }
    }).collect::<Result<Vec<Record>, _>>()?{
        output_string = output_string + &CurrRead::try_from(k)?.to_string() + "\n";
    }
    if !output_string.is_empty() {
        print!("{output_string}");
    }
    Ok(true)
}
