use clap::{Parser, Subcommand};
use nanalogue_core::{self, subcommands}; 
use std::error::Error;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Correlates basecalled length with alignment length.
    BcLenVsAlignLen {
        /// Input BAM file
        bam_file: String,
        #[clap(default_value_t = String::from(""))]
        /// Input sequence summary file from Guppy/Dorado (optional)
        seq_summ_file: String,
    },
    /// Calculates various statistics on reads.
    ReadStats {
        /// Input BAM file
        bam_file: String,
    },
    /// Print information about a read
    ReadInfo {
        /// Input BAM file
        bam_file: String,
        /// Input read id
        read_id: String,
    },
}

fn main() -> Result<(), Box<dyn Error>> {
    let cli = Cli::parse();

    // Match on the subcommand and call the corresponding logic from the library
    match cli.command {
        Commands::BcLenVsAlignLen { bam_file, seq_summ_file } => {
            subcommands::bc_len_vs_align_len::run(&bam_file, &seq_summ_file)?;
        }
        Commands::ReadStats { bam_file } => {
            subcommands::read_stats::run(&bam_file)?;
        }
        Commands::ReadInfo { bam_file, read_id } => {
            subcommands::read_info::run(&bam_file, &read_id)?;
        }
    }

    Ok(())
}
