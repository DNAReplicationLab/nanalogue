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
    /// Print basecalled len, align len, mod count per molecule
    BcLenAlignLenModCount {
        /// Turn off mod counts
        #[arg(short, long, default_value_t = false)]
        no_mod_count: bool,
        /// Input BAM file
        bam_file: String,
        #[clap(default_value_t = String::from(""))]
        /// Input sequence summary file from Guppy/Dorado (optional)
        seq_summ_file: String,
    },
    /// Calculates various summary statistics on all reads.
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
        Commands::BcLenAlignLenModCount { bam_file, seq_summ_file, no_mod_count } => {
            subcommands::bc_len_vs_align_len::run(&bam_file, &seq_summ_file, !no_mod_count)?;
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
