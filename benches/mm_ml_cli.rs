#![cfg_attr(coverage_nightly, feature(coverage_attribute))]
#![expect(
    clippy::print_stdout,
    reason = "the manual benchmark reports progress to the terminal"
)]

//! Manual end-to-end benchmark for MM/ML-related CLI commands.
//!
//! The benchmark generates deterministic `ModBAM` fixtures on demand, optionally compares a
//! baseline executable with the current checkout, alternates execution order between runs, and
//! writes descriptive JSON without enforcing a fixed performance threshold.
//!
//! Run the representative warm-cache workload with:
//! `cargo bench --bench mm_ml_cli`
//!
//! Compare another binary with verified rootless cold caches with:
//! `cargo bench --bench mm_ml_cli -- --baseline-bin /path/to/old/nanalogue --cache cold`
//!
//! Run the large `window-dens` workload with:
//! `cargo bench --bench mm_ml_cli -- --fixture window-dens-e2e --command window-dens --cache cold`

use clap::{Parser, ValueEnum};
use serde_json::{Map, Value, json};
use std::collections::BTreeMap;
use std::fs;
use std::io::{self, Read as _};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

/// Default number of timed runs per binary.
const DEFAULT_RUNS: usize = 5;
/// Default worker count passed to `nanalogue`.
const DEFAULT_THREADS: usize = 2;

/// Command-line options for the manual benchmark.
#[derive(Debug, Parser)]
#[command(
    about = "Generate deterministic ModBAM fixtures and benchmark nanalogue CLI commands",
    after_long_help = "Cold mode uses POSIX_FADV_DONTNEED and mincore on 64-bit x86 or ARM Linux to verify that the BAM, index, and executables have zero resident pages before every timed run. It does not require root, but the files must be owned by the current user.\n\nOrdinary window-dens workloads use --win 1 --step 1000 so every fixture produces output without flooding dense workloads. The window-dens-e2e fixture uses --win 300 --step 300."
)]
struct Cli {
    /// Marker passed automatically by Cargo when running a harness-free benchmark.
    #[arg(long = "bench", hide = true)]
    _cargo_bench: bool,
    /// Fixture to run; repeat as needed. Defaults to a representative set.
    #[arg(long, value_enum)]
    fixture: Vec<FixtureName>,
    /// CLI command to run; repeat as needed. Defaults to read-info.
    #[arg(long, value_enum)]
    command: Vec<BenchCommand>,
    /// Timed runs per binary.
    #[arg(long, default_value_t = NonZeroUsize::new(DEFAULT_RUNS).expect("default is nonzero"))]
    runs: NonZeroUsize,
    /// Worker threads passed to nanalogue.
    #[arg(long, default_value_t = NonZeroUsize::new(DEFAULT_THREADS).expect("default is nonzero"))]
    threads: NonZeroUsize,
    /// File-cache state before timed runs.
    #[arg(long, value_enum, default_value_t = CacheMode::Warm)]
    cache: CacheMode,
    /// Binary to benchmark. Defaults to a release build of the current checkout.
    #[arg(long)]
    candidate_bin: Option<PathBuf>,
    /// Optional comparison binary. Speedup is reported as baseline divided by candidate.
    #[arg(long)]
    baseline_bin: Option<PathBuf>,
    /// Simulator executable. Defaults to a release build of the current checkout.
    #[arg(long)]
    simulator_bin: Option<PathBuf>,
    /// Fixture and result directory. Defaults beneath target/.
    #[arg(long)]
    work_dir: Option<PathBuf>,
    /// JSON result path. Defaults to `WORK_DIR/results.json`.
    #[arg(long)]
    output: Option<PathBuf>,
    /// Regenerate fixtures even when matching complete files already exist.
    #[arg(long)]
    regenerate: bool,
}

/// Deterministic fixture choices.
#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
#[non_exhaustive]
enum FixtureName {
    /// One dense C+m group.
    DenseSingle,
    /// One sparse explicit C+m group.
    SparseSingle,
    /// One sparse implicit C+m group.
    ImplicitSparse,
    /// Dense C+m and C+h groups with valid aggregate probabilities.
    SameBaseMulti,
    /// Dense groups on A, C, G, and T.
    FourBaseDense,
    /// Many reads between 100 and 200 bases.
    ShortMany,
    /// Ultra-long reads with sparse N+n calls.
    UltraSparseN,
    /// One billion bases for the dense `window-dens` benchmark.
    WindowDensE2e,
}

impl FixtureName {
    /// Stable directory/result name.
    const fn name(self) -> &'static str {
        match self {
            Self::DenseSingle => "dense_single",
            Self::SparseSingle => "sparse_single",
            Self::ImplicitSparse => "implicit_sparse",
            Self::SameBaseMulti => "same_base_multi",
            Self::FourBaseDense => "four_base_dense",
            Self::ShortMany => "short_many",
            Self::UltraSparseN => "ultra_sparse_n",
            Self::WindowDensE2e => "window_dens_e2e",
        }
    }

    /// Deterministic simulator configuration.
    fn config(self) -> Value {
        let (contigs, reads, seed) = match self {
            Self::DenseSingle => (
                json!({"number": 1, "len_range": [200_000, 200_000]}),
                vec![read_group(
                    2_500,
                    [0.05, 0.1],
                    &[modification("C", "m", 1, None, "?", [0.0, 1.0])],
                )],
                301,
            ),
            Self::SparseSingle => (
                json!({"number": 1, "len_range": [200_000, 200_000]}),
                vec![read_group(
                    2_500,
                    [0.05, 0.1],
                    &[modification("C", "m", 1_000, Some(999), "?", [0.0, 1.0])],
                )],
                302,
            ),
            Self::ImplicitSparse => (
                json!({"number": 1, "len_range": [200_000, 200_000]}),
                vec![read_group(
                    2_500,
                    [0.05, 0.1],
                    &[modification("C", "m", 1_000, Some(999), ".", [0.0, 1.0])],
                )],
                303,
            ),
            Self::SameBaseMulti => (
                json!({"number": 1, "len_range": [200_000, 200_000]}),
                vec![read_group(
                    2_500,
                    [0.05, 0.1],
                    &[
                        modification("C", "m", 1, None, "?", [0.35, 0.45]),
                        modification("C", "h", 1, None, "?", [0.35, 0.45]),
                    ],
                )],
                304,
            ),
            Self::FourBaseDense => (
                json!({"number": 1, "len_range": [200_000, 200_000]}),
                vec![read_group(
                    1_250,
                    [0.05, 0.1],
                    &[
                        modification("A", "a", 1, None, "?", [0.0, 1.0]),
                        modification("C", "m", 1, None, "?", [0.0, 1.0]),
                        modification("G", "o", 1, None, "?", [0.0, 1.0]),
                        modification("T", "T", 1, None, "?", [0.0, 1.0]),
                    ],
                )],
                305,
            ),
            Self::ShortMany => (
                json!({"number": 1, "len_range": [10_000, 10_000]}),
                vec![read_group(
                    50_000,
                    [0.01, 0.02],
                    &[modification("C", "m", 1, None, "?", [0.0, 1.0])],
                )],
                306,
            ),
            Self::UltraSparseN => (
                json!({
                    "number": 1,
                    "len_range": [1_000_000, 1_000_000],
                    "repeated_seq": "ACGT"
                }),
                vec![read_group(
                    50,
                    [0.5, 0.8],
                    &[modification("N", "n", 1_000, Some(999), "?", [0.0, 1.0])],
                )],
                307,
            ),
            Self::WindowDensE2e => (
                json!({
                    "number": 12,
                    "len_range": [1_000_000, 1_000_000],
                    "repeated_seq": "ACGT"
                }),
                vec![json!({
                    "number": 100_000,
                    "mapq_range": [60, 60],
                    "base_qual_range": [30, 30],
                    "len_range": [0.01, 0.01],
                    "mods": [modification("T", "T", 2_500, None, "?", [0.8, 1.0])]
                })],
                42,
            ),
        };
        json!({"contigs": contigs, "reads": reads, "seed": seed})
    }

    /// Window and step sizes used by `window-dens` for this fixture.
    const fn window_dens_args(self) -> (&'static str, &'static str) {
        if matches!(self, Self::WindowDensE2e) {
            ("300", "300")
        } else {
            ("1", "1000")
        }
    }

    /// Expected output lines for workloads that require an exact validation pass.
    const fn expected_output_lines(self, command: BenchCommand) -> Option<u64> {
        if matches!(self, Self::WindowDensE2e) && matches!(command, BenchCommand::WindowDens) {
            Some(800_001)
        } else {
            None
        }
    }

    /// Input bases used to report throughput for the large end-to-end workload.
    const fn benchmark_input_bases(self, command: BenchCommand) -> Option<u64> {
        if matches!(self, Self::WindowDensE2e) && matches!(command, BenchCommand::WindowDens) {
            Some(1_000_000_000)
        } else {
            None
        }
    }
}

/// Commands that exercise MM/ML parsing.
#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
#[non_exhaustive]
enum BenchCommand {
    /// Print read information.
    ReadInfo,
    /// Print the read table with modification counts.
    ReadTableShowMods,
    /// Calculate sparse density windows.
    WindowDens,
}

impl BenchCommand {
    /// CLI subcommand spelling.
    const fn name(self) -> &'static str {
        match self {
            Self::ReadInfo => "read-info",
            Self::ReadTableShowMods => "read-table-show-mods",
            Self::WindowDens => "window-dens",
        }
    }
}

/// Cache state requested for measured runs.
#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
#[non_exhaustive]
enum CacheMode {
    /// Warm files with an unmeasured command before timing.
    Warm,
    /// Evict and verify every watched file before each timed command.
    Cold,
}

impl CacheMode {
    /// Stable output spelling.
    const fn name(self) -> &'static str {
        match self {
            Self::Warm => "warm",
            Self::Cold => "cold",
        }
    }
}

/// Descriptive statistics for one binary.
#[derive(Clone, Copy, Debug)]
struct Summary {
    /// Arithmetic mean in seconds.
    mean: f64,
    /// Sample standard deviation in seconds.
    standard_deviation: f64,
    /// Median in seconds.
    median: f64,
    /// Minimum in seconds.
    minimum: f64,
    /// Maximum in seconds.
    maximum: f64,
}

impl Summary {
    /// Convert statistics to JSON.
    fn json(self) -> Value {
        json!({
            "mean_s": self.mean,
            "sd_s": self.standard_deviation,
            "median_s": self.median,
            "min_s": self.minimum,
            "max_s": self.maximum,
        })
    }
}

/// Build one modification configuration.
fn modification(
    base: &str,
    code: &str,
    window: u32,
    configured_drop: Option<u32>,
    suffix: &str,
    probability: [f64; 2],
) -> Value {
    let mut value = json!({
        "base": base,
        "is_strand_plus": true,
        "mod_code": code,
        "mm_suffix": suffix,
        "win": [window],
        "mod_range": [probability],
    });
    if let Some(drop_count) = configured_drop {
        let _previous = value
            .as_object_mut()
            .expect("literal is an object")
            .insert("drop".to_owned(), json!([drop_count]));
    }
    value
}

/// Build common read simulation settings.
fn read_group(number: u32, lengths: [f64; 2], modifications: &[Value]) -> Value {
    json!({
        "number": number,
        "mapq_range": [10, 60],
        "base_qual_range": [10, 40],
        "len_range": lengths,
        "mismatch": 0.02,
        "mods": modifications,
    })
}

/// Return items in first-seen order without duplicates.
fn unique<T: Copy + Eq>(items: &[T]) -> Vec<T> {
    let mut result = Vec::new();
    for &item in items {
        if !result.contains(&item) {
            result.push(item);
        }
    }
    result
}

/// Repository root known when Cargo compiles this benchmark.
fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Run a command and return an error on failure.
fn run_checked(command: &mut Command, description: &str) -> Result<(), String> {
    let status = command
        .status()
        .map_err(|error| format!("failed to run {description}: {error}"))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("{description} exited with {status}"))
    }
}

/// Find the host Rust target triple.
fn host_target(root: &Path) -> Result<String, String> {
    let output = Command::new("rustc")
        .arg("-vV")
        .current_dir(root)
        .output()
        .map_err(|error| format!("failed to run rustc -vV: {error}"))?;
    if !output.status.success() {
        return Err(format!("rustc -vV exited with {}", output.status));
    }
    let stdout = String::from_utf8(output.stdout)
        .map_err(|error| format!("rustc -vV returned non-UTF-8 output: {error}"))?;
    stdout
        .lines()
        .find_map(|line| line.strip_prefix("host: ").map(str::to_owned))
        .ok_or_else(|| "could not determine the host Rust target".to_owned())
}

/// Paths used for current-checkout release binaries.
fn release_binary_paths(root: &Path) -> Result<(String, PathBuf, PathBuf), String> {
    let host = host_target(root)?;
    let release_dir = root
        .join("target/nanalogue-mm-ml-cli-build")
        .join(&host)
        .join("release");
    Ok((
        host,
        release_dir.join("nanalogue"),
        release_dir.join("nanalogue_sim_bam"),
    ))
}

/// Build only missing benchmark executables in a dedicated target directory.
fn build_release_binaries(
    root: &Path,
    build_candidate: bool,
    build_simulator: bool,
) -> Result<(Option<PathBuf>, Option<PathBuf>), String> {
    let (host, candidate, simulator) = release_binary_paths(root)?;
    let target_dir = root.join("target/nanalogue-mm-ml-cli-build");
    let mut command = Command::new("cargo");
    let _: &mut Command = command.args([
        "build",
        "--locked",
        "--quiet",
        "--release",
        "--target",
        &host,
        "--target-dir",
    ]);
    let _: &mut Command = command.arg(&target_dir).current_dir(root);
    if build_candidate {
        let _: &mut Command = command.args(["--bin", "nanalogue"]);
    }
    if build_simulator {
        let _: &mut Command = command.args(["--bin", "nanalogue_sim_bam"]);
    }
    run_checked(&mut command, "release binary build")?;
    Ok((
        build_candidate.then_some(candidate),
        build_simulator.then_some(simulator),
    ))
}

/// Resolve and validate a file used as an executable.
fn checked_executable(path: &Path, description: &str) -> Result<PathBuf, String> {
    let resolved = fs::canonicalize(path)
        .map_err(|error| format!("cannot resolve {description} {}: {error}", path.display()))?;
    if resolved.is_file() {
        Ok(resolved)
    } else {
        Err(format!(
            "{description} is not a file: {}",
            resolved.display()
        ))
    }
}

/// Reject a build output that aliases an executable supplied for another role.
fn reject_build_collision(output: &Path, supplied: &[(&str, Option<&Path>)]) -> Result<(), String> {
    let comparable_output = fs::canonicalize(output).unwrap_or_else(|_| output.to_path_buf());
    for &(description, supplied_path) in supplied {
        if supplied_path.is_some_and(|path| {
            fs::canonicalize(path).unwrap_or_else(|_| path.to_path_buf()) == comparable_output
        }) {
            return Err(format!(
                "refusing to overwrite the supplied {description} at {}; copy it outside target/nanalogue-mm-ml-cli-build",
                output.display()
            ));
        }
    }
    Ok(())
}

/// Return whether all fixture outputs exist and are non-empty.
fn fixture_complete(paths: &[&Path]) -> bool {
    paths.iter().all(|path| {
        path.metadata()
            .is_ok_and(|metadata| metadata.is_file() && metadata.len() > 0)
    })
}

/// Generate or reuse a fixture whose successful manifest matches its configuration.
fn prepare_fixture(
    fixture: FixtureName,
    simulator: &Path,
    work_dir: &Path,
    regenerate: bool,
) -> Result<PathBuf, String> {
    let fixture_dir = work_dir.join(fixture.name());
    fs::create_dir_all(&fixture_dir)
        .map_err(|error| format!("cannot create {}: {error}", fixture_dir.display()))?;
    let manifest = fixture_dir.join("config.json");
    let simulator_input = fixture_dir.join("simulator-input.json");
    let bam = fixture_dir.join("simulated.bam");
    let bai = PathBuf::from(format!("{}.bai", bam.display()));
    let fasta = fixture_dir.join("reference.fasta");
    let rendered = format!(
        "{}\n",
        serde_json::to_string_pretty(&fixture.config())
            .map_err(|error| format!("cannot serialize fixture: {error}"))?
    );
    let config_matches = fs::read_to_string(&manifest).is_ok_and(|existing| existing == rendered);
    let complete = fixture_complete(&[&bam, &bai, &fasta]);

    if regenerate || !config_matches || !complete {
        if let Err(error) = fs::remove_file(&manifest)
            && error.kind() != io::ErrorKind::NotFound
        {
            return Err(format!("cannot invalidate {}: {error}", manifest.display()));
        }
        fs::write(&simulator_input, &rendered)
            .map_err(|error| format!("cannot write {}: {error}", simulator_input.display()))?;
        println!("Generating fixture {}...", fixture.name());
        run_checked(
            Command::new(simulator).args([
                simulator_input.as_os_str(),
                bam.as_os_str(),
                fasta.as_os_str(),
            ]),
            "fixture simulator",
        )?;
        if !fixture_complete(&[&bam, &bai, &fasta]) {
            return Err(format!(
                "simulator did not create a complete fixture in {}",
                fixture_dir.display()
            ));
        }
        fs::write(&manifest, rendered)
            .map_err(|error| format!("cannot publish {}: {error}", manifest.display()))?;
    } else {
        println!("Reusing fixture {}: {}", fixture.name(), bam.display());
    }
    Ok(bam)
}

/// Capture UTF-8 stdout from a successful command.
fn captured_stdout(command: &mut Command, description: &str) -> Result<String, String> {
    let output = command
        .output()
        .map_err(|error| format!("failed to run {description}: {error}"))?;
    if !output.status.success() {
        return Err(format!("{description} exited with {}", output.status));
    }
    String::from_utf8(output.stdout)
        .map_err(|error| format!("{description} returned non-UTF-8 output: {error}"))
}

/// Validate read counts and lengths reported for the `window-dens-e2e` fixture.
fn validate_window_dens_e2e_read_stats(stdout: &str) -> Result<(), String> {
    let values: BTreeMap<&str, &str> = stdout
        .lines()
        .skip(1)
        .filter_map(|line| line.split_once('\t'))
        .collect();
    let read_count = [
        "n_primary_alignments",
        "n_secondary_alignments",
        "n_supplementary_alignments",
        "n_unmapped_reads",
    ]
    .iter()
    .try_fold(0u64, |total, key| {
        let value = values
            .get(key)
            .ok_or_else(|| format!("read-stats output is missing {key}"))?
            .parse::<u64>()
            .map_err(|error| format!("invalid read-stats value for {key}: {error}"))?;
        total
            .checked_add(value)
            .ok_or_else(|| "read count overflowed u64".to_owned())
    })?;
    if read_count != 100_000 {
        return Err(format!(
            "window-dens-e2e fixture contains {read_count} reads; expected 100000"
        ));
    }
    for key in ["seq_len_mean", "seq_len_min", "seq_len_max"] {
        let observed = values
            .get(key)
            .ok_or_else(|| format!("read-stats output is missing {key}"))?;
        if *observed != "10000" {
            return Err(format!(
                "window-dens-e2e fixture has {key}={observed}; expected 10000"
            ));
        }
    }
    Ok(())
}

/// Validate contigs and modification types reported for the `window-dens-e2e` fixture.
fn validate_window_dens_e2e_peek(stdout: &str) -> Result<(), String> {
    let remainder = stdout
        .strip_prefix("contigs_and_lengths:\n")
        .ok_or_else(|| "peek output is missing the contig heading".to_owned())?;
    let (contigs, modifications) = remainder
        .split_once("\n\nmodifications:\n")
        .ok_or_else(|| "peek output is missing the modification heading".to_owned())?;
    let valid_contigs = contigs.lines().filter(|line| {
        line.split_once('\t').is_some_and(|(name, length)| {
            name.strip_prefix("contig_").is_some_and(|suffix| {
                !suffix.is_empty() && suffix.bytes().all(|byte| byte.is_ascii_digit())
            }) && length == "1000000"
        })
    });
    if valid_contigs.count() != 12 || contigs.lines().count() != 12 {
        return Err("window-dens-e2e fixture must contain twelve 1 Mb contigs".to_owned());
    }
    if !modifications.lines().any(|line| line == "T+T") {
        return Err("window-dens-e2e fixture does not contain T+T modifications".to_owned());
    }
    Ok(())
}

/// Validate the large `window-dens-e2e` fixture before timing it.
fn validate_window_dens_e2e_fixture(
    fixture: FixtureName,
    binary: &Path,
    threads: usize,
    bam: &Path,
) -> Result<(), String> {
    if fixture != FixtureName::WindowDensE2e {
        return Ok(());
    }
    let thread_count = threads.to_string();
    let read_stats = captured_stdout(
        Command::new(binary)
            .arg("read-stats")
            .args(["--threads", &thread_count])
            .arg(bam),
        "window-dens-e2e fixture read-stats validation",
    )?;
    validate_window_dens_e2e_read_stats(&read_stats)?;
    let peek = captured_stdout(
        Command::new(binary).arg("peek").arg(bam),
        "window-dens-e2e fixture peek validation",
    )?;
    validate_window_dens_e2e_peek(&peek)
}

/// Build one benchmarked nanalogue invocation.
fn invocation(
    binary: &Path,
    fixture: FixtureName,
    command: BenchCommand,
    threads: usize,
    bam: &Path,
) -> Command {
    let mut invocation = Command::new(binary);
    let thread_count = threads.to_string();
    let _: &mut Command = invocation
        .arg(command.name())
        .args(["--threads", &thread_count]);
    if command == BenchCommand::WindowDens {
        let (window, step) = fixture.window_dens_args();
        let _: &mut Command = invocation.args(["--win", window, "--step", step]);
    }
    let _: &mut Command = invocation.arg(bam).stdout(Stdio::null());
    invocation
}

/// Run one measured CLI command and return elapsed seconds.
fn run_once(
    binary: &Path,
    fixture: FixtureName,
    command: BenchCommand,
    threads: usize,
    bam: &Path,
) -> Result<f64, String> {
    let start = Instant::now();
    run_checked(
        &mut invocation(binary, fixture, command, threads, bam),
        "benchmarked nanalogue command",
    )?;
    Ok(start.elapsed().as_secs_f64())
}

/// Count output lines from one untimed validation invocation.
#[expect(
    clippy::naive_bytecount,
    reason = "untimed output validation does not justify another dependency"
)]
fn output_line_count(
    binary: &Path,
    fixture: FixtureName,
    command: BenchCommand,
    threads: usize,
    bam: &Path,
) -> Result<u64, String> {
    let mut child = invocation(binary, fixture, command, threads, bam)
        .stdout(Stdio::piped())
        .spawn()
        .map_err(|error| format!("failed to start output validation: {error}"))?;
    let mut stdout = child
        .stdout
        .take()
        .ok_or_else(|| "output validation did not capture stdout".to_owned())?;
    let mut buffer = [0u8; 0x4000];
    let mut lines = 0u64;
    loop {
        let bytes_read = stdout
            .read(&mut buffer)
            .map_err(|error| format!("failed to read validation output: {error}"))?;
        if bytes_read == 0 {
            break;
        }
        let newlines = buffer
            .get(..bytes_read)
            .expect("read byte count cannot exceed the supplied buffer length")
            .iter()
            .filter(|byte| **byte == b'\n')
            .count();
        lines = lines
            .checked_add(u64::try_from(newlines).map_err(|error| error.to_string())?)
            .ok_or_else(|| "validation output line count overflowed u64".to_owned())?;
    }
    let status = child
        .wait()
        .map_err(|error| format!("failed to wait for output validation: {error}"))?;
    if status.success() {
        Ok(lines)
    } else {
        Err(format!("output validation exited with {status}"))
    }
}

/// Calculate descriptive statistics.
fn summarize(values: &[f64]) -> Result<Summary, String> {
    let count_u32 = u32::try_from(values.len()).map_err(|error| error.to_string())?;
    let count = f64::from(count_u32);
    let mean = values.iter().sum::<f64>() / count;
    let standard_deviation = if values.len() > 1 {
        let denominator = f64::from(
            count_u32
                .checked_sub(1)
                .expect("more than one timing has a nonzero denominator"),
        );
        (values
            .iter()
            .map(|value| (value - mean).powi(2))
            .sum::<f64>()
            / denominator)
            .sqrt()
    } else {
        0.0
    };
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let midpoint = sorted
        .len()
        .checked_div(2)
        .expect("division by the nonzero constant two succeeds");
    let median = if sorted.len().is_multiple_of(2) {
        let lower_index = midpoint
            .checked_sub(1)
            .expect("an even nonempty timing list has a lower median");
        let lower = *sorted
            .get(lower_index)
            .expect("lower median index is in bounds");
        let upper = *sorted
            .get(midpoint)
            .expect("upper median index is in bounds");
        f64::midpoint(lower, upper)
    } else {
        *sorted.get(midpoint).expect("median index is in bounds")
    };
    Ok(Summary {
        mean,
        standard_deviation,
        median,
        minimum: *sorted.first().expect("at least one run is required"),
        maximum: *sorted.last().expect("at least one run is required"),
    })
}

/// One labelled executable participating in a benchmark.
#[derive(Clone, Copy, Debug)]
struct LabelledBinary<'a> {
    /// Label stored in output.
    label: &'static str,
    /// Executable path.
    path: &'a Path,
}

/// Benchmark one fixture/command combination.
fn benchmark_case(
    binaries: &[LabelledBinary<'_>],
    fixture: FixtureName,
    command: BenchCommand,
    threads: usize,
    bam: &Path,
    runs: usize,
    cache: CacheMode,
) -> Result<Value, String> {
    let bai = PathBuf::from(format!("{}.bai", bam.display()));
    if let Some(expected_lines) = fixture.expected_output_lines(command) {
        for binary in binaries {
            let observed_lines = output_line_count(binary.path, fixture, command, threads, bam)?;
            if observed_lines != expected_lines {
                return Err(format!(
                    "{} produced {observed_lines} lines; expected {expected_lines}",
                    binary.label
                ));
            }
        }
    }
    for binary in binaries {
        let _warmup_seconds = run_once(binary.path, fixture, command, threads, bam)?;
    }
    let mut timings: BTreeMap<&str, Vec<f64>> = binaries
        .iter()
        .map(|binary| (binary.label, Vec::with_capacity(runs)))
        .collect();

    for run_index in 0..runs {
        let mut order = binaries.to_vec();
        if run_index & 1 == 1 {
            order.reverse();
        }
        for binary in order {
            if cache == CacheMode::Cold {
                let mut watched = vec![bam, bai.as_path()];
                watched.extend(binaries.iter().map(|item| item.path));
                cold_cache::evict_files(&watched)?;
            }
            let elapsed = run_once(binary.path, fixture, command, threads, bam)?;
            timings
                .get_mut(binary.label)
                .expect("every binary has a timing vector")
                .push(elapsed);
            println!(
                "  run {}/{} {:<9} {:8.3} s",
                run_index.saturating_add(1),
                runs,
                binary.label,
                elapsed
            );
        }
    }

    let mut timing_json = Map::new();
    let mut summary_json = Map::new();
    let mut summaries = BTreeMap::new();
    for (label, values) in &timings {
        let summary = summarize(values)?;
        let _: Option<Value> = timing_json.insert((*label).to_owned(), json!(values));
        let _: Option<Value> = summary_json.insert((*label).to_owned(), summary.json());
        let _: Option<Summary> = summaries.insert(*label, summary);
    }
    let mut result = Map::new();
    let _: Option<Value> = result.insert(
        "bam_bytes".to_owned(),
        json!(bam.metadata().map_err(|error| error.to_string())?.len()),
    );
    let _: Option<Value> = result.insert("timings_s".to_owned(), Value::Object(timing_json));
    let _: Option<Value> = result.insert("summaries".to_owned(), Value::Object(summary_json));
    if let (Some(baseline), Some(candidate)) =
        (summaries.get("baseline"), summaries.get("candidate"))
    {
        let _: Option<Value> = result.insert(
            "speedup_baseline_over_candidate".to_owned(),
            json!(baseline.mean / candidate.mean),
        );
    }
    if let Some(input_bases) = fixture.benchmark_input_bases(command) {
        let input_megabases = f64::from(
            u32::try_from(input_bases).expect("benchmark input-base count is bounded by u32::MAX"),
        ) / 1_000_000.0;
        let throughput: Map<String, Value> = summaries
            .iter()
            .map(|(label, summary)| ((*label).to_owned(), json!(input_megabases / summary.median)))
            .collect();
        let _: Option<Value> = result.insert(
            "median_throughput_mb_per_s".to_owned(),
            Value::Object(throughput),
        );
    }
    Ok(Value::Object(result))
}

/// Linux page-cache eviction and residency verification.
#[cfg(all(
    target_os = "linux",
    any(target_arch = "x86_64", target_arch = "aarch64")
))]
mod cold_cache {
    use std::ffi::{c_int, c_long, c_void};
    use std::fs::File;
    use std::os::fd::AsRawFd as _;
    use std::path::Path;
    use std::{io, ptr};

    /// `mmap` read protection.
    const PROT_READ: c_int = 1;
    /// Shared `mmap` mapping.
    const MAP_SHARED: c_int = 1;
    /// Linux `POSIX_FADV_DONTNEED` value.
    const POSIX_FADV_DONTNEED: c_int = 4;
    /// Linux `_SC_PAGESIZE` value.
    const SC_PAGESIZE: c_int = 30;

    unsafe extern "C" {
        /// Advise the kernel about file access.
        fn posix_fadvise(fd: c_int, offset: i64, length: i64, advice: c_int) -> c_int;
        /// Map file pages without reading them.
        fn mmap(
            address: *mut c_void,
            length: usize,
            protection: c_int,
            flags: c_int,
            fd: c_int,
            offset: i64,
        ) -> *mut c_void;
        /// Query mapped-page residency.
        fn mincore(address: *mut c_void, length: usize, vector: *mut u8) -> c_int;
        /// Release a mapping.
        fn munmap(address: *mut c_void, length: usize) -> c_int;
        /// Read a system configuration value.
        fn sysconf(name: c_int) -> c_long;
        /// Flush pending filesystem writes before page-cache eviction.
        fn sync();
    }

    /// Return a file's resident bytes without faulting its pages into memory.
    fn resident_bytes(path: &Path) -> Result<u64, String> {
        let file =
            File::open(path).map_err(|error| format!("cannot open {}: {error}", path.display()))?;
        let size_u64 = file
            .metadata()
            .map_err(|error| format!("cannot stat {}: {error}", path.display()))?
            .len();
        if size_u64 == 0 {
            return Ok(0);
        }
        let size = usize::try_from(size_u64)
            .map_err(|error| format!("{} is too large to map: {error}", path.display()))?;
        // SAFETY: `sysconf` has no pointer arguments or memory-safety preconditions.
        let page_size_raw = unsafe { sysconf(SC_PAGESIZE) };
        let page_size = usize::try_from(page_size_raw)
            .map_err(|error| format!("invalid system page size {page_size_raw}: {error}"))?;
        // SAFETY: The file descriptor remains open for the mapping call, the length is the file
        // length, and no memory is accessed through the returned pointer.
        let address = unsafe {
            mmap(
                ptr::null_mut(),
                size,
                PROT_READ,
                MAP_SHARED,
                file.as_raw_fd(),
                0,
            )
        };
        if address as usize == usize::MAX {
            return Err(format!(
                "mmap failed for {}: {}",
                path.display(),
                io::Error::last_os_error()
            ));
        }
        let mut vector = vec![0u8; size.div_ceil(page_size)];
        // SAFETY: `address` is a successful mapping of `size` bytes and `vector` has one byte per
        // mapped page, as required by `mincore`.
        let mincore_result = unsafe { mincore(address, size, vector.as_mut_ptr()) };
        let mincore_error = (mincore_result != 0).then(io::Error::last_os_error);
        // SAFETY: `address` and `size` identify the successful mapping created above.
        let unmap_result = unsafe { munmap(address, size) };
        if let Some(error) = mincore_error {
            return Err(format!("mincore failed for {}: {error}", path.display()));
        }
        if unmap_result != 0 {
            return Err(format!(
                "munmap failed for {}: {}",
                path.display(),
                io::Error::last_os_error()
            ));
        }
        let resident_pages = vector.iter().filter(|value| **value & 1 == 1).count();
        u64::try_from(resident_pages)
            .ok()
            .and_then(|pages| pages.checked_mul(u64::try_from(page_size).ok()?))
            .ok_or_else(|| "resident-byte count overflowed u64".to_owned())
    }

    /// Evict files and require zero resident pages afterwards.
    pub(super) fn evict_files(paths: &[&Path]) -> Result<(), String> {
        // SAFETY: `sync` has no arguments or memory-safety preconditions.
        unsafe {
            sync();
        }
        for path in paths {
            let file = File::open(path)
                .map_err(|error| format!("cannot open {}: {error}", path.display()))?;
            // SAFETY: The descriptor is valid for this call, offsets cover the whole file, and
            // `POSIX_FADV_DONTNEED` is a valid Linux advice constant.
            let result = unsafe { posix_fadvise(file.as_raw_fd(), 0, 0, POSIX_FADV_DONTNEED) };
            if result != 0 {
                return Err(format!(
                    "POSIX_FADV_DONTNEED failed for {}: {}",
                    path.display(),
                    io::Error::from_raw_os_error(result)
                ));
            }
        }
        for path in paths {
            let resident = resident_bytes(path)?;
            if resident != 0 {
                return Err(format!(
                    "cache eviction failed: {} has {resident} resident bytes; cold mode requires files owned by the current user",
                    path.display(),
                ));
            }
        }
        Ok(())
    }
}

/// Cold-cache mode is available only on Linux.
#[cfg(not(all(
    target_os = "linux",
    any(target_arch = "x86_64", target_arch = "aarch64")
)))]
mod cold_cache {
    use std::path::Path;

    /// Return a platform support error.
    pub(super) fn evict_files(_paths: &[&Path]) -> Result<(), String> {
        Err("cold cache mode supports only 64-bit x86 or ARM Linux".to_owned())
    }
}

/// Program implementation separated from process-level error reporting.
#[expect(
    clippy::too_many_lines,
    reason = "this top-level benchmark workflow is clearest in execution order"
)]
fn run(cli: Cli) -> Result<PathBuf, String> {
    let root = repo_root();
    let work_dir = cli
        .work_dir
        .unwrap_or_else(|| root.join("target/nanalogue-mm-ml-cli-benchmark"));
    fs::create_dir_all(&work_dir)
        .map_err(|error| format!("cannot create {}: {error}", work_dir.display()))?;
    let supplied_candidate = cli
        .candidate_bin
        .as_deref()
        .map(|path| checked_executable(path, "candidate binary"))
        .transpose()?;
    let supplied_simulator = cli
        .simulator_bin
        .as_deref()
        .map(|path| checked_executable(path, "simulator binary"))
        .transpose()?;
    let baseline = cli
        .baseline_bin
        .as_deref()
        .map(|path| checked_executable(path, "baseline binary"))
        .transpose()?;
    let needs_candidate = supplied_candidate.is_none();
    let needs_simulator = supplied_simulator.is_none();
    let (generated_candidate, generated_simulator) = if needs_candidate || needs_simulator {
        let (_, candidate_output, simulator_output) = release_binary_paths(&root)?;
        let supplied_paths = [
            ("candidate binary", supplied_candidate.as_deref()),
            ("simulator binary", supplied_simulator.as_deref()),
            ("baseline binary", baseline.as_deref()),
        ];
        if needs_candidate {
            reject_build_collision(&candidate_output, &supplied_paths)?;
        }
        if needs_simulator {
            reject_build_collision(&simulator_output, &supplied_paths)?;
        }
        println!("Building missing current-checkout release binaries...");
        build_release_binaries(&root, needs_candidate, needs_simulator)?
    } else {
        (None, None)
    };
    let candidate = supplied_candidate
        .or(generated_candidate)
        .ok_or_else(|| "candidate binary was not supplied or built".to_owned())?;
    let simulator = supplied_simulator
        .or(generated_simulator)
        .ok_or_else(|| "simulator binary was not supplied or built".to_owned())?;

    let fixtures = if cli.fixture.is_empty() {
        vec![
            FixtureName::DenseSingle,
            FixtureName::SparseSingle,
            FixtureName::ShortMany,
            FixtureName::UltraSparseN,
        ]
    } else {
        unique(&cli.fixture)
    };
    let commands = if cli.command.is_empty() {
        vec![BenchCommand::ReadInfo]
    } else {
        unique(&cli.command)
    };
    let fixture_paths: Vec<(FixtureName, PathBuf)> = fixtures
        .iter()
        .map(|fixture| {
            prepare_fixture(*fixture, &simulator, &work_dir, cli.regenerate)
                .map(|path| (*fixture, path))
        })
        .collect::<Result<_, _>>()?;
    for fixture_path in &fixture_paths {
        validate_window_dens_e2e_fixture(
            fixture_path.0,
            &candidate,
            cli.threads.get(),
            &fixture_path.1,
        )?;
    }
    let mut binaries = Vec::with_capacity(2);
    if let Some(path) = baseline.as_deref() {
        binaries.push(LabelledBinary {
            label: "baseline",
            path,
        });
    }
    binaries.push(LabelledBinary {
        label: "candidate",
        path: &candidate,
    });

    let mut benchmark_results = Map::new();
    for (fixture, bam) in fixture_paths {
        for &command in &commands {
            let key = format!("{}/{}", fixture.name(), command.name());
            println!("Benchmarking {key} ({} cache)...", cli.cache.name());
            let result = benchmark_case(
                &binaries,
                fixture,
                command,
                cli.threads.get(),
                &bam,
                cli.runs.get(),
                cli.cache,
            )?;
            let candidate_summary = result
                .get("summaries")
                .and_then(|value| value.get("candidate"))
                .ok_or_else(|| "candidate summary is missing".to_owned())?;
            println!(
                "  candidate mean={:.3} s sd={:.3} median={:.3}",
                candidate_summary["mean_s"].as_f64().unwrap_or_default(),
                candidate_summary["sd_s"].as_f64().unwrap_or_default(),
                candidate_summary["median_s"].as_f64().unwrap_or_default(),
            );
            if let Some(throughput) = result
                .get("median_throughput_mb_per_s")
                .and_then(|value| value.get("candidate"))
                .and_then(Value::as_f64)
            {
                println!("  candidate median throughput={throughput:.2} Mb/s");
            }
            let _: Option<Value> = benchmark_results.insert(key, result);
        }
    }

    let created = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|error| format!("system clock precedes the Unix epoch: {error}"))?
        .as_secs();
    let output_value = json!({
        "created_unix_seconds": created,
        "cache": cli.cache.name(),
        "runs": cli.runs.get(),
        "threads": cli.threads.get(),
        "binaries": {
            "baseline": baseline.as_ref().map(|path| path.display().to_string()),
            "candidate": candidate.display().to_string(),
            "simulator": simulator.display().to_string(),
        },
        "fixtures": fixtures.iter().map(|fixture| (fixture.name(), fixture.config())).collect::<BTreeMap<_, _>>(),
        "benchmarks": benchmark_results,
    });
    let output = cli.output.unwrap_or_else(|| work_dir.join("results.json"));
    if let Some(parent) = output.parent() {
        fs::create_dir_all(parent)
            .map_err(|error| format!("cannot create {}: {error}", parent.display()))?;
    }
    fs::write(
        &output,
        format!(
            "{}\n",
            serde_json::to_string_pretty(&output_value)
                .map_err(|error| format!("cannot serialize results: {error}"))?
        ),
    )
    .map_err(|error| format!("cannot write {}: {error}", output.display()))?;
    Ok(output)
}

#[expect(
    clippy::print_stdout,
    clippy::print_stderr,
    reason = "the manual benchmark reports progress and errors to the terminal"
)]
fn main() {
    match run(Cli::parse()) {
        Ok(output) => println!("Results: {}", output.display()),
        Err(error) => {
            eprintln!("error: {error}");
            std::process::exit(1);
        }
    }
}

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
/// Regression tests for deterministic benchmark workloads.
mod tests {
    #[test]
    fn accepts_cargo_bench_marker() {
        let cli = <super::Cli as clap::Parser>::try_parse_from(["mm_ml_cli", "--bench"])
            .expect("Cargo's benchmark marker should be accepted");

        assert!(cli._cargo_bench);

        for help_flag in ["-h", "--help"] {
            let help = <super::Cli as clap::Parser>::try_parse_from(["mm_ml_cli", help_flag])
                .expect_err("help should exit before running the benchmark");
            assert_eq!(help.kind(), clap::error::ErrorKind::DisplayHelp);
            assert!(!help.to_string().contains("--bench"));
        }
    }

    #[test]
    fn window_dens_e2e_fixture_preserves_workload() {
        let fixture = super::FixtureName::WindowDensE2e;
        assert_eq!(
            fixture.window_dens_args(),
            ("300", "300"),
            "window-dens-e2e window and step sizes changed"
        );
        assert_eq!(
            fixture.expected_output_lines(super::BenchCommand::WindowDens),
            Some(800_001),
            "window-dens-e2e output line count changed"
        );
        assert_eq!(
            fixture.benchmark_input_bases(super::BenchCommand::WindowDens),
            Some(1_000_000_000),
            "window-dens-e2e input-base count changed"
        );

        let expected_config = serde_json::json!({
            "contigs": {
                "number": 12,
                "len_range": [1_000_000, 1_000_000],
                "repeated_seq": "ACGT"
            },
            "reads": [{
                "number": 100_000,
                "mapq_range": [60, 60],
                "base_qual_range": [30, 30],
                "len_range": [0.01, 0.01],
                "mods": [super::modification("T", "T", 2_500, None, "?", [0.8, 1.0])]
            }],
            "seed": 42
        });
        assert_eq!(
            fixture.config(),
            expected_config,
            "window-dens-e2e simulator configuration changed"
        );
    }

    #[test]
    fn window_dens_e2e_fixture_rejects_wrong_read_lengths() {
        let valid = "key\tvalue\n\
n_primary_alignments\t100000\n\
n_secondary_alignments\t0\n\
n_supplementary_alignments\t0\n\
n_unmapped_reads\t0\n\
seq_len_mean\t10000\n\
seq_len_min\t10000\n\
seq_len_max\t10000\n";
        super::validate_window_dens_e2e_read_stats(valid)
            .expect("window-dens-e2e read statistics should validate");

        let wrong_lengths = valid.replace("\t10000\n", "\t9600\n");
        assert_eq!(
            super::validate_window_dens_e2e_read_stats(&wrong_lengths),
            Err("window-dens-e2e fixture has seq_len_mean=9600; expected 10000".to_owned()),
            "incorrect read lengths must fail the length validation"
        );
    }
}
