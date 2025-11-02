//! # Nanalogue Simulate BAM
//!
//! Companion tool to nanalogue which creates artificial BAM or mod BAM files
//! for developers wishing to test BAM parsing or BAM modification data parsing.
use clap::Parser;
use nanalogue_core::{Error, simulate_mod_bam};

/// Main command line parsing struct that gets paths to files to be created.
#[derive(Parser, Debug)]
#[command(author, version, 
    about = "Simulate BAM with or without modifications. Aimed at developers who wish to test their BAM parsers",
    long_about = None)]
struct Cli {
    /// Input JSON file path
    json: String,
    /// Output mod BAM file path; if pre-existing, the file will be overwritten.
    bam: String,
    /// Output fasta file path; if pre-existing, the file will be overwritten.
    fasta: String,
}

/// Main function, run the program. All business logic handled by `run`
///
/// This separation of function between `main` and `run` is so that we
/// can test the functionality of `run` easily through our code without
/// actually running the program on the command line like an external user.
fn main() {
    // Parse command line options
    let cli = Cli::parse();

    // call the run function and get the result
    match run(&cli) {
        Ok(()) => {}
        Err(e) => {
            eprintln!("Error during execution: {e}");
            std::process::exit(1);
        }
    }
}

/// Simple wrapper around `simulate_mod_bam`.
///
/// # Errors
/// Returns errors from simulating BAM files
fn run(cli: &Cli) -> Result<(), Error> {
    let json_str = std::fs::read_to_string(&cli.json)?;
    simulate_mod_bam::run(&json_str, &cli.bam, &cli.fasta)
}
