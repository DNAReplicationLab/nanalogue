#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Replays the BAM viewer's record fetches against one indexed BAM.
//!
//! The viewer calls [`RegionSequenceReader::sequences`] (table mode) or
//! [`RegionSequenceReader::profiles`] (individual mode) when it opens, and again on
//! every horizontal move, goto, and width-changing resize. This benchmark performs the
//! same sequence of calls a user session would: open the BAM, fetch the initial window,
//! press `l` repeatedly, then jump to a distant position with `g`.
//!
//! Every invocation is a fresh process, so run it after dropping the page cache to
//! measure cold-start behaviour, for example:
//!
//! ```sh
//! sync; echo 3 > /proc/sys/vm/drop_caches
//! cargo bench -q --features bam-viewer --bench bam_viewer_fetch -- \
//!     in.bam contig_00000 100000 table 101 20
//! ```
//!
//! Arguments: `BAM CONTIG START MODE WINDOW STEPS [MOD_CODE] [WIN]`, where `MODE` is
//! `table` or `individual`. `MOD_CODE` defaults to `m` and `WIN` (individual mode's
//! modified-base window) defaults to 300. Timings in milliseconds and an output digest,
//! used to confirm that optimisations leave the viewer's data unchanged, are printed as
//! one JSON line.

#[cfg(feature = "bam-viewer")]
use nanalogue_core::region_sequences::RegionSequenceReader;
#[cfg(feature = "bam-viewer")]
use nanalogue_core::{Error, ModChar};
#[cfg(feature = "bam-viewer")]
use std::hash::{DefaultHasher, Hash as _, Hasher as _};
#[cfg(feature = "bam-viewer")]
use std::num::NonZeroU32;
#[cfg(feature = "bam-viewer")]
use std::str::FromStr as _;
#[cfg(feature = "bam-viewer")]
use std::time::Instant;

/// Explains that viewer fetches are only available with the viewer feature.
#[cfg(not(feature = "bam-viewer"))]
#[expect(
    clippy::print_stderr,
    reason = "the benchmark explains why it did nothing"
)]
fn main() {
    eprintln!("bam_viewer_fetch requires the `bam-viewer` feature");
}

#[cfg(feature = "bam-viewer")]
/// Parsed command-line options.
struct Options {
    /// Indexed BAM path.
    bam: String,
    /// Reference name.
    contig: String,
    /// Zero-based first visible position.
    start: u32,
    /// Whether to fetch individual-mode profiles instead of table sequences.
    individual: bool,
    /// Window length, as the viewer derives it from the terminal width.
    window: u32,
    /// Number of `l` presses to replay.
    steps: u32,
    /// Modification type highlighted or plotted.
    mod_type: ModChar,
    /// Modified bases per individual-mode window.
    win: NonZeroU32,
}

#[cfg(feature = "bam-viewer")]
/// Parses positional arguments.
fn parse_options() -> Result<Options, Error> {
    let args: Vec<String> = std::env::args()
        .skip(1)
        .filter(|a| a != "--bench")
        .collect();
    let arg = |index: usize| -> Result<&str, Error> {
        args.get(index).map(String::as_str).ok_or_else(|| {
            Error::InvalidState(String::from(
                "usage: BAM CONTIG START MODE WINDOW STEPS [MOD_CODE] [WIN]",
            ))
        })
    };
    let individual = match arg(3)? {
        "table" => false,
        "individual" => true,
        other => return Err(Error::InvalidState(format!("unknown mode {other}"))),
    };
    Ok(Options {
        bam: arg(0)?.to_owned(),
        contig: arg(1)?.to_owned(),
        start: arg(2)?.parse()?,
        individual,
        window: arg(4)?.parse()?,
        steps: arg(5)?.parse()?,
        mod_type: ModChar::from_str(args.get(6).map_or("m", String::as_str))?,
        win: NonZeroU32::new(args.get(7).map_or(Ok(300), |w| w.parse())?)
            .ok_or_else(|| Error::InvalidState(String::from("WIN must be positive")))?,
    })
}

#[cfg(feature = "bam-viewer")]
/// Fetches one viewer window and folds its content into the digest.
///
/// Returns the number of reads fetched.
fn fetch(
    reader: &mut RegionSequenceReader,
    options: &Options,
    tid: u32,
    start: u32,
    digest: &mut DefaultHasher,
) -> Result<usize, Error> {
    let contig_len = reader
        .target_len(tid)
        .ok_or_else(|| Error::InvalidState(String::from("bad tid")))?;
    let end = start.saturating_add(options.window).min(contig_len);
    if options.individual {
        let profiles = reader.profiles(tid, start, end, options.mod_type, options.win)?;
        for profile in &profiles {
            profile.read_id().hash(digest);
            profile.is_reverse().hash(digest);
            profile.align_start().hash(digest);
            profile.align_end().hash(digest);
            profile.calls().hash(digest);
            for &(win_start, win_end, value) in profile.windows() {
                (win_start, win_end, value.val().to_bits()).hash(digest);
            }
            profile.window_series_starts().hash(digest);
        }
        Ok(profiles.len())
    } else {
        let rows = reader.sequences(tid, start, end, Some(options.mod_type))?;
        for row in &rows {
            row.read_id().hash(digest);
            row.region_offset().hash(digest);
            row.sequence().hash(digest);
            row.sequence_with_insertions().hash(digest);
            row.modifications().hash(digest);
            row.modifications_with_insertions().hash(digest);
            row.is_reverse().hash(digest);
        }
        Ok(rows.len())
    }
}

#[cfg(feature = "bam-viewer")]
/// Peak resident set size in KiB, read from Linux's `/proc/self/status` (0 elsewhere).
#[cfg(feature = "bam-viewer")]
fn peak_rss_kib() -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|status| {
            status
                .lines()
                .find_map(|line| line.strip_prefix("VmHWM:"))
                .and_then(|value| value.trim().trim_end_matches("kB").trim().parse().ok())
        })
        .unwrap_or(0)
}

/// Milliseconds elapsed since an instant.
#[cfg(feature = "bam-viewer")]
fn millis(since: Instant) -> f64 {
    since.elapsed().as_secs_f64() * 1000.0
}

#[cfg(feature = "bam-viewer")]
#[expect(clippy::print_stdout, reason = "benchmark results are printed")]
fn main() -> Result<(), Error> {
    let options = parse_options()?;
    let mut digest = DefaultHasher::new();

    let total = Instant::now();
    let reader_start = Instant::now();
    let mut reader = RegionSequenceReader::from_path(&options.bam)?;
    let tid = reader
        .target_id(&options.contig)
        .ok_or_else(|| Error::InvalidState(String::from("contig not in header")))?;
    let contig_len = reader
        .target_len(tid)
        .ok_or_else(|| Error::InvalidState(String::from("bad tid")))?;
    let open_ms = millis(reader_start);

    let first = Instant::now();
    let first_reads = fetch(&mut reader, &options, tid, options.start, &mut digest)?;
    let first_fetch_ms = millis(first);

    let steps = Instant::now();
    let mut start = options.start;
    let mut step_fetches = 0u32;
    for _ in 0..options.steps {
        let next_start = start
            .saturating_add(options.window)
            .min(contig_len.saturating_sub(options.window))
            .max(start);
        if next_start == start {
            break;
        }
        start = next_start;
        let _reads = fetch(&mut reader, &options, tid, start, &mut digest)?;
        step_fetches = step_fetches.saturating_add(1);
    }
    let steps_ms = millis(steps);

    let goto = Instant::now();
    let goto_start = options
        .start
        .saturating_add(contig_len.checked_div(4).unwrap_or(0))
        .min(contig_len.saturating_sub(options.window));
    let goto_reads = fetch(&mut reader, &options, tid, goto_start, &mut digest)?;
    let goto_ms = millis(goto);

    println!(
        "{{\"mode\":\"{}\",\"open_ms\":{open_ms:.3},\"first_fetch_ms\":{first_fetch_ms:.3},\
         \"first_reads\":{first_reads},\"steps\":{},\"step_fetches\":{step_fetches},\
         \"steps_total_ms\":{steps_ms:.3},\
         \"goto_ms\":{goto_ms:.3},\"goto_reads\":{goto_reads},\"total_ms\":{:.3},\
         \"peak_rss_kib\":{},\"digest\":\"{:016x}\"}}",
        if options.individual {
            "individual"
        } else {
            "table"
        },
        options.steps,
        millis(total),
        peak_rss_kib(),
        digest.finish(),
    );
    Ok(())
}
