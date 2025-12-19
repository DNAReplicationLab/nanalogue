//! # Nanalogue (Nucleic Acid Analogue)
//!
//! We process and calculate data associated with DNA/RNA molecules, their alignments to
//! reference genomes, modification information on them, and other miscellaneous
//! information from BAM files.
use clap::Parser as _;
use nanalogue_core::commands;
use rust_htslib::htslib;
use std::io;

/// Main function, run the program. All business logic handled by [`commands::run`].
///
/// This separation of function between `main` and [`commands::run`] is so that we
/// can test the functionality of `run` easily through our code without
/// actually running the program on the command line like an external user.
fn main() {
    // Initialize SSL certificates for HTTPS support
    nanalogue_core::init_ssl_certificates();

    // we do not want to print `htslib` errors but want to deal with
    // them using our own error handling.
    unsafe {
        htslib::hts_set_log_level(0);
    }

    // Parse command line options
    let cli = commands::Cli::parse();

    // set up writers to stdout
    let stdout = io::stdout();
    let handle = io::BufWriter::new(stdout);

    // call the run function and get the result
    match commands::run(cli, handle) {
        Ok(()) => {}
        Err(e) => {
            eprintln!("Error during execution: {e}");
            std::process::exit(1);
        }
    }
}
