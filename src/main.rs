use clap::{Parser, Subcommand};
use nanalogue_core::{
    self, BamPreFilt, BamRcRecords, Error, F32Bw0and1, InputBam, ModChar, OrdPair,
    nanalogue_bam_reader, subcommands,
};
use std::num::NonZeroU32;
use std::ops::RangeInclusive;

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
        /// only find reads such that all windowed density values are in these limits.
        /// specify as low,high e.g. 0.2,0.7 so that the condition is that all windowed_values
        /// satisfy low <= windowed_value <= high.
        #[clap(long)]
        all_dens_between: OrdPair<F32Bw0and1>,
        /// invert above filter, so the filter is now: atleast one windowed_value is < low or
        /// at least one windowed_value is > high.
        #[clap(long, default_value = "false")]
        invert: bool,
    },
}

fn main() -> Result<(), Error> {
    let cli = Cli::parse();

    /// pre filtering the BAM file according to input options
    macro_rules! pre_filt {
        ( $b : expr, $c : expr) => {
            $b.rc_records.filter(|r| match r {
                Ok(v) => v.pre_filt($c),
                Err(_) => true,
            })
        };
    }

    // Match on the subcommand and call the corresponding logic from the library
    let result = match cli.command {
        Commands::ReadsTableWithMods { bam, seq_summ_file } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader)?;
            subcommands::bc_len_vs_align_len::run(
                pre_filt!(bam_rc_records, &bam),
                &seq_summ_file,
                true,
            )
        }
        Commands::ReadsTableNoMods { bam, seq_summ_file } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader)?;
            subcommands::bc_len_vs_align_len::run(
                pre_filt!(bam_rc_records, &bam),
                &seq_summ_file,
                false,
            )
        }
        Commands::ReadStats { bam } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader)?;
            subcommands::read_stats::run(pre_filt!(bam_rc_records, &bam))
        }
        Commands::ReadInfo { bam, read_id } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader)?;
            let contig_names = bam_rc_records.contig_names.clone();
            subcommands::read_info::run(
                pre_filt!(bam_rc_records, &bam),
                &read_id,
                Some(contig_names),
            )
        }
        Commands::FindModifiedReads {
            bam,
            tag,
            win,
            slide,
            all_dens_between,
            invert,
        } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader)?;
            subcommands::find_modified_reads::run(
                pre_filt!(bam_rc_records, &bam),
                tag,
                win,
                slide,
                |x| RangeInclusive::from(all_dens_between).contains(x),
                invert,
            )
        }
    };

    match result {
        Ok(true) | Ok(false) => Ok(()),
        Err(e) => {
            eprintln!("Error during execution: {e}");
            std::process::exit(1);
        }
    }
}
