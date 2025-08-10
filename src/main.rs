use clap::{Parser, Subcommand};
use nanalogue_core::{
    self, BamPreFilt, BamRcRecords, Contains, Error, F32Bw0and1, InputBam, ModChar, OrdPair,
    RestrictModCalledStrand, ThresholdState, nanalogue_bam_reader, subcommands,
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
    /// Find names of modified reads through criteria specified by sub commands
    FindModifiedReads {
        #[command(subcommand)]
        command: FindModReadsCommands,
    },
}

#[derive(Subcommand, Debug)]
enum FindModReadsCommands {
    /// Find reads with all windowed modification densities within specified limits
    AllDensBetween {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
        /// modified tag
        #[clap(long)]
        tag: ModChar,
        /// modified strand
        #[clap(long)]
        strand: Option<RestrictModCalledStrand>,
        /// window size
        #[clap(long)]
        win: NonZeroU32,
        /// step window by this size
        #[clap(long)]
        step: NonZeroU32,
        /// only find reads such that all windowed density values are in these limits.
        /// specify as low,high e.g. 0.2,0.7 so that the condition is that all windowed_values
        /// satisfy low <= windowed_value <= high.
        #[clap(long)]
        dens_limits: OrdPair<F32Bw0and1>,
    },
    /// Find reads with windowed modification density such that at least one window is
    /// at or above the high value.
    AnyDensAbove {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
        /// modified tag
        #[clap(long)]
        tag: ModChar,
        /// modified strand
        #[clap(long)]
        strand: Option<RestrictModCalledStrand>,
        /// window size
        #[clap(long)]
        win: NonZeroU32,
        /// step window by this size
        #[clap(long)]
        step: NonZeroU32,
        /// high value of criterion i.e. at least one window >= high
        #[clap(long)]
        high: F32Bw0and1,
    },
    /// Find reads with windowed modification density such that at least one window is
    /// at or below the low value
    AnyDensBelow {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
        /// modified tag
        #[clap(long)]
        tag: ModChar,
        /// modified strand
        #[clap(long)]
        strand: Option<RestrictModCalledStrand>,
        /// window size
        #[clap(long)]
        win: NonZeroU32,
        /// step window by this size
        #[clap(long)]
        step: NonZeroU32,
        /// low value of criterion i.e. at least one window <= low
        #[clap(long)]
        low: F32Bw0and1,
    },
    /// Find reads with windowed modification density such that at least one window is
    /// at or below the low value and at least one window is at or above the high value.
    /// This operation may enrich for reads with spatial gradients in modification density.
    AnyDensBelowAndAnyDensAbove {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
        /// modified tag
        #[clap(long)]
        tag: ModChar,
        /// modified strand
        #[clap(long)]
        strand: Option<RestrictModCalledStrand>,
        /// window size
        #[clap(long)]
        win: NonZeroU32,
        /// step window by this size
        #[clap(long)]
        step: NonZeroU32,
        /// low criterion i.e. at least one window <= low,
        /// anded with the high criterion
        #[clap(long)]
        low: F32Bw0and1,
        /// high criterion i.e. at least one window >= high,
        /// anded with the low criterion
        #[clap(long)]
        high: F32Bw0and1,
    },
    /// Find reads with windowed modification density such that max of all densities
    /// minus min of all densities is at least the value specified.
    /// This operation may enrich for reads with spatial gradients in modification density.
    DensRangeAbove {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
        /// modified tag
        #[clap(long)]
        tag: ModChar,
        /// modified strand
        #[clap(long)]
        strand: Option<RestrictModCalledStrand>,
        /// window size
        #[clap(long)]
        win: NonZeroU32,
        /// step window by this size
        #[clap(long)]
        step: NonZeroU32,
        /// max(windowed densities) - min(windowed densities)
        /// is at least this value
        #[clap(long)]
        min_range: F32Bw0and1,
    },
    /// Find reads such that absolute value of gradient in modification density
    /// measured in windows is at least the value specified.
    AnyAbsGradAbove {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
        /// modified tag
        #[clap(long)]
        tag: ModChar,
        /// modified strand
        #[clap(long)]
        strand: Option<RestrictModCalledStrand>,
        /// window size
        #[clap(long)]
        win: NonZeroU32,
        /// step window by this size
        #[clap(long)]
        step: NonZeroU32,
        /// gradient is at least this value. e.g. a gradient of 0.005 with a window
        /// size of 100 means you expect a variation of 0.005 * 100 = 0.5 over at least
        /// one window i.e. greater than 0.5 or smaller than -0.5. For your guidance,
        /// the theoretical maximum gradient in a window of size W for our data is approximately
        /// 3/(2W) (when W is much larger than 1), so you can use this to think what
        /// gradient level you want to set.
        #[clap(long)]
        min_grad: F32Bw0and1,
    },
}

fn main() -> Result<(), Error> {
    let cli = Cli::parse();

    // threshold and calculate mean modification density per window
    let threshold_and_mean =
        |mod_list: &[u8], _: &[Option<i64>], _: &[Option<i64>]| -> Result<F32Bw0and1, Error> {
            let win_size: usize = mod_list.len();
            let count_mod: usize = mod_list
                .iter()
                .filter(|x| ThresholdState::GtEq(128).contains(x))
                .count();
            F32Bw0and1::new(count_mod as f32 / win_size as f32)
        };

    // threshold and calculate gradient of mod density per window.
    // NOTE: to calculate gradient, we need (x, y) data where x is the
    // coordinate along the read and y is 0 or 1 depending on whether the
    // base is modified or not. But, in the calculation below, we use only
    // the y values i.e. the modified values and assume even spacing
    // along the x direction i.e. along the read coordinate. This assumption
    // will break down if modifications occur in bursts, or if positions
    // are missing along the read, or in other scenarios. As our goal
    // here is to use a simple method to calculate gradients and select
    // interesting reads, we are o.k. with an approximate calculation.
    let threshold_and_abs_gradient =
        |mod_list: &[u8], _: &[Option<i64>], _: &[Option<i64>]| -> Result<F32Bw0and1, Error> {
            let win_size = mod_list.len();
            let x_mean: f32 = (win_size as f32 + 1.0) / 2.0;
            let numerator: f32 = mod_list
                .iter()
                .enumerate()
                .map(|(i, x)| {
                    if ThresholdState::GtEq(128).contains(x) {
                        i as f32 + 1.0 - x_mean
                    } else {
                        0.0
                    }
                })
                .sum();
            let denominator: f32 = (1..=win_size)
                .map(|x| (x as f32 - x_mean) * (x as f32 - x_mean))
                .sum();
            F32Bw0and1::new(f32::abs(numerator) / denominator)
        };

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
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &bam)?;
            subcommands::bc_len_vs_align_len::run(
                pre_filt!(bam_rc_records, &bam),
                &seq_summ_file,
                true,
            )
        }
        Commands::ReadsTableNoMods { bam, seq_summ_file } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &bam)?;
            subcommands::bc_len_vs_align_len::run(
                pre_filt!(bam_rc_records, &bam),
                &seq_summ_file,
                false,
            )
        }
        Commands::ReadStats { bam } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &bam)?;
            subcommands::read_stats::run(pre_filt!(bam_rc_records, &bam))
        }
        Commands::ReadInfo { bam, read_id } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &bam)?;
            let contig_names = bam_rc_records.contig_names.clone();
            subcommands::read_info::run(
                pre_filt!(bam_rc_records, &bam),
                &read_id,
                Some(contig_names),
            )
        }
        Commands::FindModifiedReads {
            command:
                FindModReadsCommands::AllDensBetween {
                    bam,
                    tag,
                    strand,
                    win,
                    step,
                    dens_limits,
                },
        } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &bam)?;
            subcommands::find_modified_reads::run(
                pre_filt!(bam_rc_records, &bam),
                tag,
                strand,
                win,
                step,
                threshold_and_mean,
                |x| {
                    x.iter()
                        .all(|r| RangeInclusive::from(dens_limits).contains(r))
                },
            )
        }
        Commands::FindModifiedReads {
            command:
                FindModReadsCommands::AnyDensAbove {
                    bam,
                    tag,
                    strand,
                    win,
                    step,
                    high,
                },
        } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &bam)?;
            let interval_high_to_1 = OrdPair::new(high, F32Bw0and1::one())?;
            subcommands::find_modified_reads::run(
                pre_filt!(bam_rc_records, &bam),
                tag,
                strand,
                win,
                step,
                threshold_and_mean,
                |x| {
                    x.iter()
                        .any(|r| RangeInclusive::from(interval_high_to_1).contains(r))
                },
            )
        }
        Commands::FindModifiedReads {
            command:
                FindModReadsCommands::AnyDensBelow {
                    bam,
                    tag,
                    strand,
                    win,
                    step,
                    low,
                },
        } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &bam)?;
            let interval_0_to_low = OrdPair::new(F32Bw0and1::zero(), low)?;
            subcommands::find_modified_reads::run(
                pre_filt!(bam_rc_records, &bam),
                tag,
                strand,
                win,
                step,
                threshold_and_mean,
                |x| {
                    x.iter()
                        .any(|r| RangeInclusive::from(interval_0_to_low).contains(r))
                },
            )
        }
        Commands::FindModifiedReads {
            command:
                FindModReadsCommands::AnyDensBelowAndAnyDensAbove {
                    bam,
                    tag,
                    strand,
                    win,
                    step,
                    low,
                    high,
                },
        } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &bam)?;
            let interval_0_to_low = OrdPair::new(F32Bw0and1::zero(), low)?;
            let interval_high_to_1 = OrdPair::new(high, F32Bw0and1::one())?;
            subcommands::find_modified_reads::run(
                pre_filt!(bam_rc_records, &bam),
                tag,
                strand,
                win,
                step,
                threshold_and_mean,
                |x| {
                    x.iter()
                        .any(|r| RangeInclusive::from(interval_0_to_low).contains(r))
                        && x.iter()
                            .any(|r| RangeInclusive::from(interval_high_to_1).contains(r))
                },
            )
        }
        Commands::FindModifiedReads {
            command:
                FindModReadsCommands::DensRangeAbove {
                    bam,
                    tag,
                    strand,
                    win,
                    step,
                    min_range,
                },
        } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &bam)?;
            subcommands::find_modified_reads::run(
                pre_filt!(bam_rc_records, &bam),
                tag,
                strand,
                win,
                step,
                threshold_and_mean,
                |x| {
                    x.iter()
                        .map(|r| r.get_val())
                        .reduce(f32::max)
                        .unwrap_or(0.0)
                        - x.iter()
                            .map(|r| r.get_val())
                            .reduce(f32::min)
                            .unwrap_or(0.0)
                        >= min_range.get_val()
                },
            )
        }
        Commands::FindModifiedReads {
            command:
                FindModReadsCommands::AnyAbsGradAbove {
                    bam,
                    tag,
                    strand,
                    win,
                    step,
                    min_grad,
                },
        } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &bam)?;
            let interval_min_grad_to_1 = OrdPair::new(min_grad, F32Bw0and1::one())?;
            subcommands::find_modified_reads::run(
                pre_filt!(bam_rc_records, &bam),
                tag,
                strand,
                win,
                step,
                threshold_and_abs_gradient,
                |x| {
                    x.iter()
                        .any(|r| RangeInclusive::from(interval_min_grad_to_1).contains(r))
                },
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
