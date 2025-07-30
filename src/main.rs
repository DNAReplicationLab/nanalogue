use clap::{Parser, Subcommand};
use nanalogue_core::{self, subcommands, Error};
use std::num::NonZeroU32;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Print basecalled len, align len, mod count per molecule
    ReadsTableWithMods {
        /// Input BAM file
        bam_file: String,
        /// Input sequence summary file from Guppy/Dorado (optional)
        #[clap(default_value_t = String::from(""))]
        seq_summ_file: String,
    },
    ReadsTableNoMods {
        /// Input BAM file
        bam_file: String,
        /// Input sequence summary file from Guppy/Dorado (optional)
        #[clap(default_value_t = String::from(""))]
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
    /// Find reads with user-supplied level of modifications
    FindModifiedReads {
        /// Input BAM file
        bam_file: String,
        /// modified tag
        #[clap(long)]
        tag: String,
        /// window size
        #[clap(long)]
        win: NonZeroU32,
        /// slide
        #[clap(long)]
        slide: NonZeroU32,
        /// maximum density
        #[clap(long)]
        dens_max: f32,
        /// invert above filter
        #[clap(long, default_value = "false")]
        invert: bool,
    }
}

fn main() -> Result<(), Error> {
    let cli = Cli::parse();

    // Match on the subcommand and call the corresponding logic from the library
    let result = match cli.command {
        Commands::ReadsTableWithMods {
            bam_file,
            seq_summ_file,
        } => {
            subcommands::bc_len_vs_align_len::run(&bam_file, &seq_summ_file, true)
        },
        Commands::ReadsTableNoMods {
            bam_file,
            seq_summ_file,
        } => {
            subcommands::bc_len_vs_align_len::run(&bam_file, &seq_summ_file, false)
        },
        Commands::ReadStats { bam_file } => {
            subcommands::read_stats::run(&bam_file)
        },
        Commands::ReadInfo { bam_file, read_id } => {
            subcommands::read_info::run(&bam_file, &read_id)
        },
        Commands::FindModifiedReads { bam_file, tag, win, slide, dens_max, invert } => {
            subcommands::find_modified_reads::run(&bam_file, &tag, win, slide, dens_max, invert)
        },
    };

    match result {
        Ok(true) | Ok(false) => Ok(()),
        Err(e) => {
            eprintln!("Error during execution: {e}");
            std::process::exit(1);
        },
    }
}
