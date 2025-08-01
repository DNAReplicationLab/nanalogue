use clap::{Parser, Subcommand};
use nanalogue_core::{self, subcommands, Error, F32Bw0and1, OrdPair, InputBam, ModChar};
use std::num::NonZeroU32;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Prints basecalled len, align len, mod count per molecule
    ReadsTableWithMods {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
        /// Input sequence summary file from Guppy/Dorado (optional)
        #[clap(default_value_t = String::from(""))]
        seq_summ_file: String,
    },
    /// Prints basecalled len, align len per molecule
    ReadsTableNoMods {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
        /// Input sequence summary file from Guppy/Dorado (optional)
        #[clap(default_value_t = String::from(""))]
        seq_summ_file: String,
    },
    /// Calculates various summary statistics on all reads.
    ReadStats {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
    },
    /// Print information about a read
    ReadInfo {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
        /// Input read id
        read_id: String,
    },
    /// Find reads with user-supplied limits of modifications
    FindModifiedReads {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
        /// modified tag
        #[clap(long)]
        tag: ModChar,
        /// window size
        #[clap(long)]
        win: NonZeroU32,
        /// slide
        #[clap(long)]
        slide: NonZeroU32,
        /// density limits, specify as min,max where these are 2 numbers between 0 and 1
        #[clap(long)]
        dens_limits: OrdPair<F32Bw0and1>,
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
            mut bam,
            seq_summ_file,
        } => {
            subcommands::bc_len_vs_align_len::run(&mut bam, &seq_summ_file, true)
        },
        Commands::ReadsTableNoMods {
            mut bam,
            seq_summ_file,
        } => {
            subcommands::bc_len_vs_align_len::run(&mut bam, &seq_summ_file, false)
        },
        Commands::ReadStats { mut bam } => {
            subcommands::read_stats::run(&mut bam)
        },
        Commands::ReadInfo { mut bam, read_id } => {
            subcommands::read_info::run(&mut bam, &read_id)
        },
        Commands::FindModifiedReads { mut bam, tag, win, slide, dens_limits, invert } => {
            subcommands::find_modified_reads::run(&mut bam, tag, win, slide, dens_limits, invert)
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
