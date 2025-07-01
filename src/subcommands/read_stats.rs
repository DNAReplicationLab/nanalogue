use rust_htslib::{bam};
use std::error::Error;

pub fn run(bam_path: &str) -> Result<(), Box<dyn Error>> {
    println!("Executing read_stats on: {}", bam_path);

    // We just open the file for now, do nothing else
    let _ = bam::Reader::from_path(bam_path)?;

    Ok(())
}

