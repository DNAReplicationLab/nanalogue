#![deny(
    missing_copy_implementations,
    missing_debug_implementations,
    missing_docs,
    trivial_casts,
    trivial_numeric_casts,
    unused_extern_crates,
    unused_import_braces,
    unused_qualifications,
    unused_results
)]
//! # Nanalogue (Nucleic Acid Analogue)
//!
//! We process and calculate data associated with DNA molecules, their alignments to
//! reference genomes, modification information on them, and other miscellaneous
//! information.
//!
use bedrs::{Bed3, Coordinates};
use clap::{Parser, Subcommand};
use nanalogue_core::{
    self, BamPreFilt, BamRcRecords, Error, F32Bw0and1, GenomicRegion, InputBam, InputMods,
    InputWindowing, OptionalTag, OrdPair, RequiredTag, analysis, nanalogue_bam_reader,
    find_modified_reads, read_info, read_stats, reads_table, window_reads, simulate_mod_bam,
};
use std::io;
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
    ReadTableShowMods {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
        /// Input modification options
        #[clap(flatten)]
        mods: InputMods<OptionalTag>,
        /// Genomic region from which basecalled sequences are displayed (optional)
        #[clap(long)]
        seq_region: Option<GenomicRegion>,
        /// Displays entire basecalled sequence (optional)
        #[clap(long, conflicts_with = "seq_region")]
        seq_full: bool,
        /// Input sequence summary file from Guppy/Dorado (optional)
        #[clap(default_value_t = String::from(""))]
        seq_summ_file: String,
    },
    /// Prints basecalled len, align len per molecule
    ReadTableHideMods {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
        /// Genomic region from which basecalled sequences are displayed (optional)
        #[clap(long)]
        seq_region: Option<GenomicRegion>,
        /// Displays entire basecalled sequence (optional)
        #[clap(long, conflicts_with = "seq_region")]
        seq_full: bool,
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
    /// Prints information about reads
    ReadInfo {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
    },
    /// Find names of modified reads through criteria specified by sub commands
    FindModifiedReads {
        #[command(subcommand)]
        command: FindModReadsCommands,
    },
    /// Output windowed densities of all reads
    WindowDens {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
        /// Input windowing options
        #[clap(flatten)]
        win: InputWindowing,
        /// Input modification options
        #[clap(flatten)]
        mods: InputMods<OptionalTag>,
    },
    /// Output windowed gradients of all reads
    WindowGrad {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
        /// Input windowing options
        #[clap(flatten)]
        win: InputWindowing,
        /// Input modification options
        #[clap(flatten)]
        mods: InputMods<OptionalTag>,
    },
    /// Write a simulated mod BAM file according to json options
    WriteSimulatedModBAM {
        /// Input JSON file path
        json: String,
        /// Output mod BAM file path
        bam: String,
        /// Output fasta file path
        fasta: String,
    },
}

#[derive(Subcommand, Debug)]
enum FindModReadsCommands {
    /// Find reads with all windowed modification densities within specified limits
    AllDensBetween {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
        /// Input windowing options
        #[clap(flatten)]
        win: InputWindowing,
        /// Input modification options
        #[clap(flatten)]
        mods: InputMods<RequiredTag>,
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
        /// Input windowing options
        #[clap(flatten)]
        win: InputWindowing,
        /// Input modification options
        #[clap(flatten)]
        mods: InputMods<RequiredTag>,
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
        /// Input windowing options
        #[clap(flatten)]
        win: InputWindowing,
        /// Input modification options
        #[clap(flatten)]
        mods: InputMods<RequiredTag>,
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
        /// Input windowing options
        #[clap(flatten)]
        win: InputWindowing,
        /// Input modification options
        #[clap(flatten)]
        mods: InputMods<RequiredTag>,
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
        /// Input windowing options
        #[clap(flatten)]
        win: InputWindowing,
        /// Input modification options
        #[clap(flatten)]
        mods: InputMods<RequiredTag>,
        /// max(windowed densities) - min(windowed densities)
        /// is at least this value
        #[clap(long)]
        min_range: F32Bw0and1,
    },
    /// Find reads such that absolute value of gradient in modification density
    /// measured in windows is at least the value specified.
    /// This operation enriches for reads with spatial gradients in modification density.
    AnyAbsGradAbove {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
        /// Input windowing options
        #[clap(flatten)]
        win: InputWindowing,
        /// Input modification options
        #[clap(flatten)]
        mods: InputMods<RequiredTag>,
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

    // set up writers to stdout
    let stdout = io::stdout();
    let mut handle = io::BufWriter::new(stdout);

    /// pre filtering the BAM file according to input options
    macro_rules! pre_filt {
        ( $b : expr, $c : expr) => {
            $b.rc_records.filter(|r| match r {
                Ok(v) => v.pre_filt($c),
                Err(_) => true,
            })
        };
    }

    /// decides whether to display sequence from a region,
    /// the full basecalled sequence, or neither
    macro_rules! seq_display {
        ( $b : expr, $c : expr, $d : expr) => {
            match ($b, $c) {
                (None, false) => None,
                (Some(v), false) => Some(v.try_to_bed3($d)?),
                (None, true) => Some(Bed3::<i32, u64>::empty()),
                (Some(_), true) => {
                    return Err(Error::NotImplementedError(
                        "cannot call seq-region and seq-full together".to_string(),
                    ));
                }
            }
        };
    }

    // Match on the subcommand and call the corresponding logic from the library
    let result = match cli.command {
        Commands::ReadTableShowMods {
            mut bam,
            mut mods,
            seq_region,
            seq_full,
            seq_summ_file,
        } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            reads_table::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                Some(mods),
                seq_display!(seq_region, seq_full, bam_rc_records.header),
                &seq_summ_file,
            )
        }
        Commands::ReadTableHideMods {
            mut bam,
            seq_region,
            seq_full,
            seq_summ_file,
        } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records =
                BamRcRecords::new(&mut bam_reader, &mut bam, &mut InputMods::default())?;
            reads_table::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                None,
                seq_display!(seq_region, seq_full, bam_rc_records.header),
                &seq_summ_file,
            )
        }
        Commands::ReadStats { mut bam } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records =
                BamRcRecords::new(&mut bam_reader, &mut bam, &mut InputMods::default())?;
            read_stats::run(&mut handle, pre_filt!(bam_rc_records, &bam))
        }
        Commands::ReadInfo { mut bam } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records =
                BamRcRecords::new(&mut bam_reader, &mut bam, &mut InputMods::default())?;
            read_info::run(&mut handle, pre_filt!(bam_rc_records, &bam))
        }
        Commands::FindModifiedReads {
            command:
                FindModReadsCommands::AllDensBetween {
                    mut bam,
                    win,
                    mut mods,
                    dens_limits,
                },
        } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            find_modified_reads::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                win,
                mods,
                analysis::threshold_and_mean,
                |x| {
                    x.iter()
                        .all(|r| RangeInclusive::from(dens_limits).contains(r))
                },
            )
        }
        Commands::FindModifiedReads {
            command:
                FindModReadsCommands::AnyDensAbove {
                    mut bam,
                    win,
                    mut mods,
                    high,
                },
        } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            let interval_high_to_1 = OrdPair::new(high, F32Bw0and1::one())?;
            find_modified_reads::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                win,
                mods,
                analysis::threshold_and_mean,
                |x| {
                    x.iter()
                        .any(|r| RangeInclusive::from(interval_high_to_1).contains(r))
                },
            )
        }
        Commands::FindModifiedReads {
            command:
                FindModReadsCommands::AnyDensBelow {
                    mut bam,
                    win,
                    mut mods,
                    low,
                },
        } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            let interval_0_to_low = OrdPair::new(F32Bw0and1::zero(), low)?;
            find_modified_reads::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                win,
                mods,
                analysis::threshold_and_mean,
                |x| {
                    x.iter()
                        .any(|r| RangeInclusive::from(interval_0_to_low).contains(r))
                },
            )
        }
        Commands::FindModifiedReads {
            command:
                FindModReadsCommands::AnyDensBelowAndAnyDensAbove {
                    mut bam,
                    win,
                    mut mods,
                    low,
                    high,
                },
        } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            let interval_0_to_low = OrdPair::new(F32Bw0and1::zero(), low)?;
            let interval_high_to_1 = OrdPair::new(high, F32Bw0and1::one())?;
            find_modified_reads::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                win,
                mods,
                analysis::threshold_and_mean,
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
                    mut bam,
                    win,
                    mut mods,
                    min_range,
                },
        } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            find_modified_reads::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                win,
                mods,
                analysis::threshold_and_mean,
                |x| {
                    x.iter().map(|r| r.val()).reduce(f32::max).unwrap_or(0.0)
                        - x.iter().map(|r| r.val()).reduce(f32::min).unwrap_or(0.0)
                        >= min_range.val()
                },
            )
        }
        Commands::FindModifiedReads {
            command:
                FindModReadsCommands::AnyAbsGradAbove {
                    mut bam,
                    win,
                    mut mods,
                    min_grad,
                },
        } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            let interval_min_grad_to_1 = OrdPair::new(min_grad, F32Bw0and1::one())?;
            find_modified_reads::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                win,
                mods,
                |x| {
                    Ok(F32Bw0and1::abs_f32_abs_val_below_1(
                        analysis::threshold_and_gradient(x)?,
                    ))
                },
                |x| {
                    x.iter()
                        .any(|r| RangeInclusive::from(interval_min_grad_to_1).contains(r))
                },
            )
        }
        Commands::WindowDens {
            mut bam,
            win,
            mut mods,
        } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            window_reads::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                win,
                mods,
                |x| analysis::threshold_and_mean(x).map(|y| y.into()),
            )
        }
        Commands::WindowGrad {
            mut bam,
            win,
            mut mods,
        } => {
            let mut bam_reader = nanalogue_bam_reader(&bam.bam_path)?;
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            window_reads::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                win,
                mods,
                analysis::threshold_and_gradient,
            )
        }
        Commands::WriteSimulatedModBAM { json, bam, fasta } => {
            let json_str = std::fs::read_to_string(&json)?;
            simulate_mod_bam::run(&json_str, &bam, &fasta)
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

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn verify_cli() {
        Cli::command().debug_assert();
    }
}
