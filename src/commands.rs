//! # Commands run in `main.rs`
//!
//! We set up the commands and their code in this file.
use crate::{
    BamPreFilt as _, BamRcRecords, Error, F32Bw0and1, GenomicRegion, InputBam, InputMods,
    InputWindowing, OptionalTag, OrdPair, PathOrURLOrStdin, RequiredTag, SeqDisplayOptions,
    analysis, find_modified_reads, nanalogue_bam_reader, nanalogue_bam_reader_from_stdin,
    nanalogue_bam_reader_from_url, nanalogue_indexed_bam_reader,
    nanalogue_indexed_bam_reader_from_url, peek, read_info, read_stats, reads_table, window_reads,
};
use clap::{Parser, Subcommand};
use derive_builder::Builder;
use rust_htslib::bam;
use rust_htslib::errors::Error as RHError;
use std::io;
use std::ops::RangeInclusive;

/// Main command line parsing struct
#[derive(Builder, Parser, Debug, Clone)]
#[command(author, version, about, long_about = None)]
#[non_exhaustive]
pub struct Cli {
    /// Our subcommands
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug, Clone)]
#[non_exhaustive]
/// Commands of nanalogue CLI used in `main.rs`
pub enum Commands {
    /// Prints basecalled len, align len, mod count per molecule
    ReadTableShowMods {
        /// Input BAM file
        #[clap(flatten)]
        bam: InputBam,
        /// Input modification options
        #[clap(flatten)]
        mods: InputMods<OptionalTag>,
        /// Genomic region from which basecalled sequences are displayed (optional)
        #[clap(long, conflicts_with = "seq_full", group = "seq")]
        seq_region: Option<GenomicRegion>,
        /// Displays entire basecalled sequence (optional)
        #[clap(long, conflicts_with = "seq_region", group = "seq")]
        seq_full: bool,
        /// Displays basecalling qualities (optional)
        #[clap(long, requires = "seq")]
        show_base_qual: bool,
        /// Show insertions in lower case
        #[clap(long, requires = "seq_region")]
        show_ins_lowercase: bool,
        /// Shows modified bases as Z (or z depending on other options)
        #[clap(long, requires = "seq_region")]
        show_mod_z: bool,
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
        #[clap(long, conflicts_with = "seq_full", group = "seq")]
        seq_region: Option<GenomicRegion>,
        /// Displays entire basecalled sequence (optional)
        #[clap(long, conflicts_with = "seq_region", group = "seq")]
        seq_full: bool,
        /// Displays basecalling qualities (optional)
        #[clap(long, requires = "seq")]
        show_base_qual: bool,
        /// Show insertions in lower case
        #[clap(long, requires = "seq_region")]
        show_ins_lowercase: bool,
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
        /// Input modification options
        #[clap(flatten)]
        mods: InputMods<OptionalTag>,
        /// Print detailed modification data (JSON)
        #[clap(long, conflicts_with = "detailed_pretty")]
        detailed: bool,
        /// Pretty-print detailed modification data (JSON)
        #[clap(long, conflicts_with = "detailed")]
        detailed_pretty: bool,
    },
    /// Find names of modified reads through criteria specified by sub commands
    FindModifiedReads {
        /// Subcommands corresponding to finding modified reads
        #[command(subcommand)]
        command: FindModReadsSubcommands,
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
    /// Display BAM file contigs, contig lengths, and mod types from a "peek" at the header and first 100 records
    Peek {
        /// Input BAM file (path, URL, or '-' for stdin)
        bam: PathOrURLOrStdin,
    },
}

#[derive(Subcommand, Debug, Clone)]
#[non_exhaustive]
/// Subcommands to find single-molecules that fit some criterion.
pub enum FindModReadsSubcommands {
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
        /// specify as low,high e.g. 0.2,0.7 so that the condition is that all `windowed_values`
        /// satisfy low <= `windowed_value` <= high.
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

impl FindModReadsSubcommands {
    /// Get the BAM command line options alone
    #[must_use]
    pub fn bam(self) -> InputBam {
        match self {
            FindModReadsSubcommands::AllDensBetween { bam, .. }
            | FindModReadsSubcommands::AnyDensAbove { bam, .. }
            | FindModReadsSubcommands::AnyDensBelow { bam, .. }
            | FindModReadsSubcommands::AnyDensBelowAndAnyDensAbove { bam, .. }
            | FindModReadsSubcommands::DensRangeAbove { bam, .. }
            | FindModReadsSubcommands::AnyAbsGradAbove { bam, .. } => bam,
        }
    }
}

impl Commands {
    /// Get the BAM command line options alone
    #[must_use]
    pub fn bam(self) -> InputBam {
        match self {
            Commands::FindModifiedReads { command } => command.bam(),
            Commands::ReadTableShowMods { bam, .. }
            | Commands::ReadTableHideMods { bam, .. }
            | Commands::ReadStats { bam, .. }
            | Commands::ReadInfo { bam, .. }
            | Commands::WindowDens { bam, .. }
            | Commands::WindowGrad { bam, .. } => bam,
            Commands::Peek { bam } => bam.into(),
        }
    }
}

/// Subcommands are handed off through this function, which accepts
/// command line options and a writeable handle as input.
///
/// # Errors
/// Returns errors associated with the subcommands or if command line
/// options are problematic, or if BAM file opening fails.
pub fn run<W>(cli: Cli, handle: W) -> Result<(), Error>
where
    W: io::Write,
{
    // Open BAM file/path/data and call on subcommands.
    // NOTE: we filter by region twice, once using the index, and
    // once when we read the records obtained. So, even if we do not have
    // an index, we can still filter by region, but this operation could be
    // slower or much slower. So, we are fine ignoring the region at this step
    // if the input is from stdin or if the index cannot be found.
    let bam_cli = cli.command.clone().bam();
    match (bam_cli.bam_path, bam_cli.region) {
        (PathOrURLOrStdin::Stdin, _) => run_on_bam(cli, handle, nanalogue_bam_reader_from_stdin()?),
        (PathOrURLOrStdin::Path(v), None) => run_on_bam(cli, handle, nanalogue_bam_reader(&v)?),
        (PathOrURLOrStdin::URL(v), None) => {
            run_on_bam(cli, handle, nanalogue_bam_reader_from_url(&v)?)
        }
        (PathOrURLOrStdin::Path(v), Some(w)) => {
            match nanalogue_indexed_bam_reader(&v, (&w).try_into()?) {
                Err(Error::RustHtslibError(RHError::BamInvalidIndex { .. })) => {
                    println!("# cannot find index file. region retrieval could be slower.");
                    run_on_bam(cli, handle, nanalogue_bam_reader(&v)?)
                }
                Err(e) => Err(e),
                Ok(x) => run_on_bam(cli, handle, x),
            }
        }
        (PathOrURLOrStdin::URL(v), Some(w)) => {
            match nanalogue_indexed_bam_reader_from_url(&v, (&w).try_into()?) {
                Err(Error::RustHtslibError(RHError::BamInvalidIndex { .. })) => {
                    println!("# cannot find index file. region retrieval could be slower.");
                    run_on_bam(cli, handle, nanalogue_bam_reader_from_url(&v)?)
                }
                Err(e) => Err(e),
                Ok(x) => run_on_bam(cli, handle, x),
            }
        }
    }
}

/// Subcommands are executed through this function, which accepts
/// command line options, a writeable handle, and a readable BAM as input.
///
/// # Errors
/// Returns errors associated with the subcommands or if command line
/// options are problematic
#[expect(clippy::too_many_lines, reason = "Comprehensive CLI command routing")]
pub fn run_on_bam<R, W>(cli: Cli, mut handle: W, mut bam_reader: R) -> Result<(), Error>
where
    W: io::Write,
    R: bam::Read,
{
    /// pre filtering the BAM file according to input options
    macro_rules! pre_filt {
        ( $b : expr, $c : expr) => {
            $b.rc_records
                .filter(|r| r.as_ref().map_or(true, |v| v.pre_filt($c)))
        };
    }

    /// decides whether to display sequence from a region,
    /// the full basecalled sequence, or neither, and optionally the basecalled qualities.
    /// Note: clap's `conflicts_with` ensures `seq_region` and `seq_full` are mutually exclusive
    macro_rules! seq_display {
        ( $b : expr, $c : expr, $d : expr, $e : expr, $f : expr, $g : expr) => {
            match ($b, $c, $d, $e, $f) {
                (None, false, _, _, _) => SeqDisplayOptions::No,
                (Some(v), false, flag_qual, flag_ins, flag_mods) => SeqDisplayOptions::Region {
                    show_base_qual: flag_qual,
                    show_ins_lowercase: flag_ins,
                    region: v.try_to_bed3(&$g)?,
                    show_mod_z: flag_mods,
                },
                (None, true, flag_qual, false, false) => SeqDisplayOptions::Full {
                    show_base_qual: flag_qual,
                },
                (None, true, _, true, _) | (None, true, _, _, true) | (Some(_), true, _, _, _) => {
                    unreachable!(
                        "clap prevents seq_region and seq_full, or setting insert \
pos retrieval/mod colouring without seq_region"
                    )
                }
            }
        };
    }

    // Match on the subcommand and call the corresponding logic from the library
    match cli.command {
        Commands::ReadTableShowMods {
            mut bam,
            mut mods,
            seq_region,
            seq_full,
            show_base_qual,
            show_ins_lowercase,
            show_mod_z,
            seq_summ_file,
        } => {
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            reads_table::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                Some(mods),
                seq_display!(
                    seq_region,
                    seq_full,
                    show_base_qual,
                    show_ins_lowercase,
                    show_mod_z,
                    bam_rc_records.header
                ),
                &seq_summ_file,
            )
        }
        Commands::ReadTableHideMods {
            mut bam,
            seq_region,
            seq_full,
            show_base_qual,
            show_ins_lowercase,
            seq_summ_file,
        } => {
            let bam_rc_records = BamRcRecords::new(
                &mut bam_reader,
                &mut bam,
                &mut InputMods::<OptionalTag>::default(),
            )?;
            reads_table::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                None,
                seq_display!(
                    seq_region,
                    seq_full,
                    show_base_qual,
                    show_ins_lowercase,
                    false,
                    bam_rc_records.header
                ),
                &seq_summ_file,
            )
        }
        Commands::ReadStats { mut bam } => {
            let bam_rc_records = BamRcRecords::new(
                &mut bam_reader,
                &mut bam,
                &mut InputMods::<OptionalTag>::default(),
            )?;
            read_stats::run(&mut handle, pre_filt!(bam_rc_records, &bam))
        }
        Commands::ReadInfo {
            mut bam,
            mut mods,
            detailed,
            detailed_pretty,
        } => {
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            read_info::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                mods,
                (detailed || detailed_pretty).then_some(detailed_pretty),
            )
        }
        Commands::FindModifiedReads {
            command:
                FindModReadsSubcommands::AllDensBetween {
                    mut bam,
                    win,
                    mut mods,
                    dens_limits,
                },
        } => {
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            find_modified_reads::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                win,
                &mods,
                analysis::threshold_and_mean,
                |x| {
                    x.iter()
                        .all(|r| RangeInclusive::from(dens_limits).contains(r))
                },
            )
        }
        Commands::FindModifiedReads {
            command:
                FindModReadsSubcommands::AnyDensAbove {
                    mut bam,
                    win,
                    mut mods,
                    high,
                },
        } => {
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            let interval_high_to_1 = OrdPair::new(high, F32Bw0and1::one())?;
            find_modified_reads::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                win,
                &mods,
                analysis::threshold_and_mean,
                |x| {
                    x.iter()
                        .any(|r| RangeInclusive::from(interval_high_to_1).contains(r))
                },
            )
        }
        Commands::FindModifiedReads {
            command:
                FindModReadsSubcommands::AnyDensBelow {
                    mut bam,
                    win,
                    mut mods,
                    low,
                },
        } => {
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            let interval_0_to_low = OrdPair::new(F32Bw0and1::zero(), low)?;
            find_modified_reads::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                win,
                &mods,
                analysis::threshold_and_mean,
                |x| {
                    x.iter()
                        .any(|r| RangeInclusive::from(interval_0_to_low).contains(r))
                },
            )
        }
        Commands::FindModifiedReads {
            command:
                FindModReadsSubcommands::AnyDensBelowAndAnyDensAbove {
                    mut bam,
                    win,
                    mut mods,
                    low,
                    high,
                },
        } => {
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            let interval_0_to_low = OrdPair::new(F32Bw0and1::zero(), low)?;
            let interval_high_to_1 = OrdPair::new(high, F32Bw0and1::one())?;
            find_modified_reads::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                win,
                &mods,
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
                FindModReadsSubcommands::DensRangeAbove {
                    mut bam,
                    win,
                    mut mods,
                    min_range,
                },
        } => {
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            find_modified_reads::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                win,
                &mods,
                analysis::threshold_and_mean,
                |x| {
                    x.iter()
                        .map(F32Bw0and1::val)
                        .reduce(f32::max)
                        .unwrap_or(0.0)
                        - x.iter()
                            .map(F32Bw0and1::val)
                            .reduce(f32::min)
                            .unwrap_or(0.0)
                        >= min_range.val()
                },
            )
        }
        Commands::FindModifiedReads {
            command:
                FindModReadsSubcommands::AnyAbsGradAbove {
                    mut bam,
                    win,
                    mut mods,
                    min_grad,
                },
        } => {
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            let interval_min_grad_to_1 = OrdPair::new(min_grad, F32Bw0and1::one())?;
            find_modified_reads::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                win,
                &mods,
                |x| {
                    Ok(F32Bw0and1::abs_f32_abs_val_at_most_1(
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
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            window_reads::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                win,
                &mods,
                |x| analysis::threshold_and_mean(x).map(Into::into),
            )
        }
        Commands::WindowGrad {
            mut bam,
            win,
            mut mods,
        } => {
            let bam_rc_records = BamRcRecords::new(&mut bam_reader, &mut bam, &mut mods)?;
            window_reads::run(
                &mut handle,
                pre_filt!(bam_rc_records, &bam),
                win,
                &mods,
                analysis::threshold_and_gradient,
            )
        }
        Commands::Peek { bam } => {
            let mut input_bam: InputBam = bam.into();
            let bam_rc_records = BamRcRecords::new(
                &mut bam_reader,
                &mut input_bam,
                &mut InputMods::<OptionalTag>::default(),
            )?;
            peek::run(
                &mut handle,
                &bam_rc_records.header,
                bam_rc_records.rc_records.take(100),
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory as _;

    #[test]
    fn verify_cli() {
        Cli::command().debug_assert();
    }
}
