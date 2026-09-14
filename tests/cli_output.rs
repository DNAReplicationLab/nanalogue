//! Executable-level tests for CLI output stream contracts.

/// Tests that invoke the compiled CLI executable.
#[cfg(test)]
mod tests {
    use std::{path::Path, process::Command};

    const MAIN_HELP_CASES: &[(&[&str], &str)] = &[
        (&[], "nanalogue.txt"),
        (
            &["read-table-show-mods"],
            "nanalogue_read-table-show-mods.txt",
        ),
        (
            &["read-table-hide-mods"],
            "nanalogue_read-table-hide-mods.txt",
        ),
        (&["read-stats"], "nanalogue_read-stats.txt"),
        (&["read-info"], "nanalogue_read-info.txt"),
        (
            &["find-modified-reads"],
            "nanalogue_find-modified-reads.txt",
        ),
        (
            &["find-modified-reads", "all-dens-between"],
            "nanalogue_find-modified-reads_all-dens-between.txt",
        ),
        (
            &["find-modified-reads", "any-dens-above"],
            "nanalogue_find-modified-reads_any-dens-above.txt",
        ),
        (
            &["find-modified-reads", "any-dens-below"],
            "nanalogue_find-modified-reads_any-dens-below.txt",
        ),
        (
            &["find-modified-reads", "any-dens-below-and-any-dens-above"],
            "nanalogue_find-modified-reads_any-dens-below-and-any-dens-above.txt",
        ),
        (
            &["find-modified-reads", "dens-range-above"],
            "nanalogue_find-modified-reads_dens-range-above.txt",
        ),
        (
            &["find-modified-reads", "any-abs-grad-above"],
            "nanalogue_find-modified-reads_any-abs-grad-above.txt",
        ),
        (&["window-dens"], "nanalogue_window-dens.txt"),
        (&["window-grad"], "nanalogue_window-grad.txt"),
        (&["peek"], "nanalogue_peek.txt"),
    ];

    fn help_output(executable: &str, args: &[&str]) -> std::process::Output {
        Command::new(executable)
            .args(args)
            .arg("--help")
            .env("TERM", "dumb")
            .env("NO_COLOR", "1")
            .env("CLICOLOR", "0")
            .env("CLICOLOR_FORCE", "0")
            .env("COLUMNS", "100")
            .output()
            .expect("CLI executable should run")
    }

    fn assert_help_matches_golden(executable: &str, args: &[&str], golden_name: &str) {
        let output = help_output(executable, args);
        assert!(
            output.status.success(),
            "{} --help failed: {}",
            args.join(" "),
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            output.stderr.is_empty(),
            "{} --help should not write to stderr",
            args.join(" ")
        );

        let golden_path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/golden/cli_help")
            .join(golden_name);
        let expected =
            std::fs::read_to_string(&golden_path).expect("golden help output should be readable");
        let actual = String::from_utf8(output.stdout).expect("help output should be UTF-8");
        assert_eq!(actual, expected, "help differs from {golden_name}");
    }

    fn run_main(args: &[&str]) -> std::process::Output {
        Command::new(env!("CARGO_BIN_EXE_nanalogue"))
            .args(args)
            .env("TERM", "dumb")
            .env("NO_COLOR", "1")
            .env("CLICOLOR", "0")
            .env("CLICOLOR_FORCE", "0")
            .env("COLUMNS", "100")
            .output()
            .expect("nanalogue executable should run")
    }

    fn assert_main_error(args: &[&str], expected_fragments: &[&str]) {
        let output = run_main(args);
        assert_eq!(
            output.status.code(),
            Some(2),
            "invalid arguments should use Clap's exit code: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            output.stdout.is_empty(),
            "argument errors should not write to stdout"
        );
        let stderr = String::from_utf8(output.stderr).expect("diagnostic should be UTF-8");
        assert!(stderr.contains("error:"));
        assert!(stderr.contains("For more information, try '--help'."));
        for fragment in expected_fragments {
            assert!(
                stderr.contains(fragment),
                "diagnostic for {} should contain {fragment}: {stderr}",
                args.join(" ")
            );
        }
    }

    /// Every reachable main-CLI command has deterministic long help captured in a golden file.
    #[test]
    fn main_help_matches_recursive_goldens() {
        for &(args, golden) in MAIN_HELP_CASES {
            assert_help_matches_golden(env!("CARGO_BIN_EXE_nanalogue"), args, golden);
        }
    }

    /// The simulator is also Clap-based and its long help is captured in full.
    #[test]
    fn simulator_help_matches_golden() {
        assert_help_matches_golden(
            env!("CARGO_BIN_EXE_nanalogue_sim_bam"),
            &[],
            "nanalogue_sim_bam.txt",
        );
    }

    /// Defaults and representative flattened arguments parse through the real executable.
    #[test]
    fn representative_flattened_arguments_are_accepted() {
        let bam = concat!(env!("CARGO_MANIFEST_DIR"), "/examples/example_1.bam");
        for args in [
            &["read-info", bam][..],
            &[
                "read-info",
                "--threads",
                "1",
                "--min-seq-len",
                "0",
                "--sample-fraction",
                "1",
                "--tag",
                "m",
                "--mod-prob-filter",
                "0.4,0.6",
                "--trim-read-ends-mod",
                "0",
                bam,
            ][..],
        ] {
            let output = run_main(args);
            assert!(
                output.status.success(),
                "{} should succeed: {}",
                args.join(" "),
                String::from_utf8_lossy(&output.stderr)
            );
        }
    }

    /// A nested subcommand accepts its required flattened arguments through the real executable.
    #[test]
    fn nested_subcommand_arguments_are_accepted() {
        let bam = concat!(env!("CARGO_MANIFEST_DIR"), "/examples/example_1.bam");
        let args = [
            "find-modified-reads",
            "any-dens-above",
            "--win",
            "2",
            "--step",
            "1",
            "--tag",
            "m",
            "--high",
            "0.7",
            bam,
        ];
        let output = run_main(&args);
        assert!(
            output.status.success(),
            "nested command should succeed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    /// Conflicting options are rejected before command execution.
    #[test]
    fn conflicting_arguments_have_diagnostics() {
        for (args, expected) in [
            (
                &["read-info", "--detailed", "--detailed-pretty", "reads.bam"][..],
                &["--detailed", "--detailed-pretty"][..],
            ),
            (
                &[
                    "read-stats",
                    "--read-id",
                    "one",
                    "--read-id-list",
                    "ids.txt",
                    "reads.bam",
                ][..],
                &["--read-id", "--read-id-list"][..],
            ),
            (
                &[
                    "read-table-hide-mods",
                    "--seq-region",
                    "chr1:1-10",
                    "--seq-full",
                    "reads.bam",
                ][..],
                &["--seq-region", "--seq-full"][..],
            ),
        ] {
            assert_main_error(args, expected);
        }
    }

    /// Direct requirements, required flattened options, and groups are enforced.
    #[test]
    fn missing_required_arguments_have_diagnostics() {
        for (args, expected) in [
            (
                &["read-table-hide-mods", "--show-base-qual", "reads.bam"][..],
                &["--seq-region", "--seq-full"][..],
            ),
            (
                &["read-stats", "--full-region", "reads.bam"][..],
                &["--region"][..],
            ),
            (
                &[
                    "find-modified-reads",
                    "any-dens-above",
                    "--win",
                    "20",
                    "--step",
                    "5",
                    "--high",
                    "0.7",
                    "reads.bam",
                ][..],
                &["--tag"][..],
            ),
            (
                &["window-dens", "--step", "5", "reads.bam"][..],
                &["--win"][..],
            ),
            (
                &["window-dens", "--win", "20", "reads.bam"][..],
                &["--step"][..],
            ),
        ] {
            assert_main_error(args, expected);
        }

        let bam = concat!(env!("CARGO_MANIFEST_DIR"), "/examples/example_1.bam");
        let output = run_main(&[
            "read-table-hide-mods",
            "--show-base-qual",
            "--seq-full",
            bam,
        ]);
        assert!(
            output.status.success(),
            "seq-full should satisfy the seq group: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    /// Invalid values and unknown subcommands produce actionable diagnostics.
    #[test]
    fn invalid_inputs_have_diagnostics() {
        for (args, expected) in [
            (
                &["read-stats", "--sample-fraction", "1.5", "reads.bam"][..],
                &["--sample-fraction", "1.5"][..],
            ),
            (
                &["read-info", "--mod-prob-filter", "0.6,0.4", "reads.bam"][..],
                &["--mod-prob-filter", "0.6,0.4"][..],
            ),
            (
                &["read-info", "--threads", "0", "reads.bam"][..],
                &["invalid value", "--threads", "0"][..],
            ),
            (&["not-a-command"][..], &["not-a-command", "Usage:"][..]),
        ] {
            assert_main_error(args, expected);
        }
    }

    /// Version requests succeed for both Clap binaries and use stdout only.
    #[test]
    fn version_uses_the_cli_output_contract() {
        for (executable, name) in [
            (env!("CARGO_BIN_EXE_nanalogue"), "nanalogue"),
            (env!("CARGO_BIN_EXE_nanalogue_sim_bam"), "nanalogue"),
        ] {
            let output = Command::new(executable)
                .arg("--version")
                .output()
                .expect("CLI executable should run");

            assert!(output.status.success(), "version request should succeed");
            let stdout = String::from_utf8(output.stdout).expect("version output should be UTF-8");
            assert_eq!(stdout, format!("{name} {}\n", env!("CARGO_PKG_VERSION")));
            assert!(
                output.stderr.is_empty(),
                "version should not write to stderr"
            );
        }
    }

    /// Detailed read information remains valid JSON when region lookup falls back
    /// to reading an unindexed local BAM sequentially.
    #[test]
    fn detailed_read_info_keeps_index_warning_off_stdout() {
        let bam_path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/examples/example_1_copy_no_index.bam"
        );
        let output = Command::new(env!("CARGO_BIN_EXE_nanalogue"))
            .args([
                "read-info",
                "--detailed",
                "--region",
                "dummyI:1-22",
                bam_path,
            ])
            .output()
            .expect("nanalogue executable should run");

        assert!(
            output.status.success(),
            "command failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let _json = serde_json::from_slice::<serde_json::Value>(&output.stdout)
            .expect("detailed stdout should contain only valid JSON");

        let warning = "# cannot find index file. region retrieval could be slower.";
        assert!(
            !String::from_utf8_lossy(&output.stdout).contains(warning),
            "warning should not be written to stdout"
        );
        assert!(
            String::from_utf8_lossy(&output.stderr).contains(warning),
            "warning should be written to stderr"
        );
    }
}

/// Combinatorial black-box tests for the simulator executable.
#[cfg(test)]
mod simulator_cli_combinatorial_tests {
    use nanalogue_core::uuid;
    use rand::{RngExt as _, SeedableRng as _, rngs::StdRng};
    use rust_htslib::bam::{self, Read as _};
    use std::{
        collections::HashSet,
        ffi::{OsStr, OsString},
        fs,
        path::{Path, PathBuf},
        process::{Command, Output},
    };

    /// Fixed seed for reproducible positional and unknown-flag generation.
    const CLI_MATRIX_SEED: u64 = 0x05EE_DC11;
    /// Sanity limit for committed help and manifest inputs.
    const MAX_SOURCE_BYTES: usize = 50_000;
    /// Sanity limit for parsing the package manifest.
    const MAX_MANIFEST_LINES: usize = 50_000;
    /// Sanity limit for paths, arguments, diagnostics, and command output.
    const MAX_TEXT_BYTES: usize = 10_000;
    /// Sanity limit for one CLI invocation's argument count.
    const MAX_ARGUMENTS: usize = 40;
    /// Resource limit for exhaustive argument permutation.
    const MAX_PERMUTATION_ARGUMENTS: usize = 10;
    /// Maximum iterations permitted while generating argument permutations.
    const MAX_PERMUTATION_LOOPS: usize = 3_628_800;
    /// Sanity limit for the compiled simulator executable.
    const MAX_EXECUTABLE_BYTES: u64 = 1_000_000_000;
    /// Expected FASTA content derived directly from the committed fixture.
    const EXPECTED_FASTA: &[u8] = b">contig_00000\n\
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";

    /// Exact display output accepted from simulator CLI actions.
    struct ExpectedOutputs {
        /// Concise `-h` output.
        short_help: Vec<u8>,
        /// Detailed `--help` output.
        long_help: Vec<u8>,
        /// Package name and version output.
        version: Vec<u8>,
    }

    impl ExpectedOutputs {
        /// Load exact short-help, long-help, and version output expectations.
        fn load() -> Self {
            let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
            let golden_dir = manifest_dir.join("tests/golden/cli_help");
            let short_help = fs::read(golden_dir.join("nanalogue_sim_bam_short.txt"))
                .expect("simulator short-help golden should be readable");
            let long_help = fs::read(golden_dir.join("nanalogue_sim_bam.txt"))
                .expect("simulator long-help golden should be readable");
            let manifest = fs::read_to_string(manifest_dir.join("Cargo.toml"))
                .expect("Cargo.toml should be readable");

            // Overflow guards for the small committed inputs parsed by this test.
            for (label, byte_len) in [
                ("short-help golden", short_help.len()),
                ("long-help golden", long_help.len()),
                ("Cargo.toml", manifest.len()),
            ] {
                assert!(byte_len > 0, "{label} should be nonempty");
                assert!(
                    byte_len < MAX_SOURCE_BYTES,
                    "{label} should be smaller than {MAX_SOURCE_BYTES} bytes"
                );
            }

            let mut in_package = false;
            let mut name = None;
            let mut version = None;

            for (line_index, raw_line) in manifest.lines().enumerate() {
                assert!(
                    line_index < MAX_MANIFEST_LINES,
                    "Cargo.toml should contain fewer than {MAX_MANIFEST_LINES} lines"
                );
                let line = raw_line.trim();
                if line.starts_with('[') {
                    if in_package {
                        break;
                    }
                    in_package = line == "[package]";
                    continue;
                }
                if !in_package {
                    continue;
                }
                let Some((key, raw_value)) = line.split_once('=') else {
                    continue;
                };
                let value = raw_value.trim();
                let Some(unquoted) = value
                    .strip_prefix('"')
                    .and_then(|without_prefix| without_prefix.strip_suffix('"'))
                else {
                    continue;
                };
                match key.trim() {
                    "name" => assert!(
                        name.replace(unquoted.to_owned()).is_none(),
                        "duplicate package name"
                    ),
                    "version" => assert!(
                        version.replace(unquoted.to_owned()).is_none(),
                        "duplicate package version"
                    ),
                    _ => {}
                }
            }

            let package_name = name.expect("Cargo.toml [package] must contain a quoted name");
            let package_version =
                version.expect("Cargo.toml [package] must contain a quoted version");
            Self {
                short_help,
                long_help,
                version: format!("{package_name} {package_version}\n").into_bytes(),
            }
        }
    }

    /// Temporary directory removed when the owning scope exits, including after panic.
    #[derive(Debug)]
    struct TempDir {
        /// Unique directory path.
        path: PathBuf,
    }

    impl TempDir {
        /// Create a uniquely named temporary directory.
        fn new(label: &str) -> Self {
            let path =
                std::env::temp_dir().join(format!("nanalogue_cli_{label}_{}", uuid::v4_random()));
            let path_len = path.as_os_str().as_encoded_bytes().len();
            assert!(path_len > 0, "temporary path should be nonempty");
            assert!(
                path_len < MAX_TEXT_BYTES,
                "temporary path should be shorter than {MAX_TEXT_BYTES} bytes"
            );
            fs::create_dir_all(&path).expect("temporary directory should be created");
            Self { path }
        }

        /// Return the temporary directory path.
        fn path(&self) -> &Path {
            &self.path
        }
    }

    impl Drop for TempDir {
        fn drop(&mut self) {
            let path_len = self.path.as_os_str().as_encoded_bytes().len();
            assert!(path_len > 0, "temporary path should be nonempty");
            assert!(
                path_len < MAX_TEXT_BYTES,
                "temporary path should be shorter than {MAX_TEXT_BYTES} bytes"
            );
            drop(fs::remove_dir_all(&self.path));
        }
    }

    /// Paths used by a real three-positional simulator invocation.
    struct SimulationPaths {
        /// Input JSON copied from the committed fixture.
        input: PathBuf,
        /// BAM output path.
        bam: PathBuf,
        /// BAM index output path.
        index: PathBuf,
        /// FASTA output path.
        fasta: PathBuf,
    }

    impl SimulationPaths {
        /// Copy the seeded fixture and construct all output paths.
        fn new(temp: &TempDir) -> Self {
            let temp_path_len = temp.path().as_os_str().as_encoded_bytes().len();
            assert!(temp_path_len > 0, "temporary path should be nonempty");
            assert!(
                temp_path_len < MAX_TEXT_BYTES,
                "temporary path should be shorter than {MAX_TEXT_BYTES} bytes"
            );
            let paths = Self {
                input: temp.path().join("input.json"),
                bam: temp.path().join("output.bam"),
                index: temp.path().join("output.bam.bai"),
                fasta: temp.path().join("output.fasta"),
            };
            let copied_bytes = fs::copy(
                Path::new(env!("CARGO_MANIFEST_DIR"))
                    .join("tests/fixtures/sim_bam/seeded_repeated_sequence.json"),
                &paths.input,
            )
            .expect("seeded fixture should copy into temporary directory");
            assert!(copied_bytes > 0, "seeded fixture should be nonempty");
            assert!(
                copied_bytes < MAX_TEXT_BYTES as u64,
                "seeded fixture should be smaller than {MAX_TEXT_BYTES} bytes"
            );
            paths
        }
    }

    /// Run the compiled simulator executable with deterministic terminal settings.
    fn execute<I, S>(args: I) -> Output
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        let executable = Path::new(env!("CARGO_BIN_EXE_nanalogue_sim_bam"));
        assert!(executable.exists(), "simulator executable should exist");
        let executable_bytes = fs::metadata(executable)
            .expect("simulator executable metadata should be readable")
            .len();
        assert!(
            executable_bytes > 0,
            "simulator executable should be nonempty"
        );
        assert!(
            executable_bytes < MAX_EXECUTABLE_BYTES,
            "simulator executable should be smaller than {MAX_EXECUTABLE_BYTES} bytes"
        );

        let owned_args = args
            .into_iter()
            .map(|arg| arg.as_ref().to_owned())
            .collect::<Vec<_>>();
        // The zero-positional, no-flag matrix cell intentionally has no arguments.
        assert!(
            owned_args.len() < MAX_ARGUMENTS,
            "CLI should receive fewer than {MAX_ARGUMENTS} arguments"
        );
        for arg in &owned_args {
            let arg_len = arg.as_encoded_bytes().len();
            assert!(arg_len > 0, "CLI arguments should be nonempty");
            assert!(
                arg_len < MAX_TEXT_BYTES,
                "CLI arguments should be shorter than {MAX_TEXT_BYTES} bytes"
            );
        }

        Command::new(executable)
            .args(&owned_args)
            .env("TERM", "dumb")
            .env("NO_COLOR", "1")
            .env("CLICOLOR", "0")
            .env("CLICOLOR_FORCE", "0")
            .env("COLUMNS", "100")
            .env("LANG", "C")
            .env("LC_ALL", "C")
            .output()
            .expect("nanalogue_sim_bam executable should run")
    }

    /// Assert a successful display action emits one applicable exact output.
    fn assert_display(output: &Output, expected: &[&[u8]], context: &str) {
        assert!(
            !output.stdout.is_empty(),
            "{context}: display stdout should be nonempty"
        );
        assert!(
            output.stdout.len() < MAX_TEXT_BYTES,
            "{context}: display stdout should be shorter than {MAX_TEXT_BYTES} bytes"
        );
        assert!(
            output.stderr.len() < MAX_TEXT_BYTES,
            "{context}: display stderr should be shorter than {MAX_TEXT_BYTES} bytes"
        );
        assert!(
            !expected.is_empty(),
            "{context}: expected output is required"
        );
        assert!(
            expected.len() < MAX_TEXT_BYTES,
            "{context}: expected output list should contain fewer than {MAX_TEXT_BYTES} items"
        );
        for candidate in expected {
            assert!(
                !candidate.is_empty(),
                "{context}: expected output should be nonempty"
            );
            assert!(
                candidate.len() < MAX_TEXT_BYTES,
                "{context}: expected output should be shorter than {MAX_TEXT_BYTES} bytes"
            );
        }
        assert!(!context.is_empty(), "display context should be nonempty");
        assert!(
            context.len() < MAX_TEXT_BYTES,
            "display context should be shorter than {MAX_TEXT_BYTES} bytes"
        );
        assert_eq!(
            output.status.code(),
            Some(0),
            "{context}: stderr={}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            output.stderr.is_empty(),
            "{context}: stderr={}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            expected.contains(&output.stdout.as_slice()),
            "{context}: unexpected stdout={}",
            String::from_utf8_lossy(&output.stdout)
        );
    }

    /// Assert a Clap argument error has the common stream and guidance contract.
    fn assert_cli_error(output: Output, expected_diagnostics: &[&str], context: &str) {
        assert!(
            output.stdout.len() < MAX_TEXT_BYTES,
            "{context}: error stdout should be shorter than {MAX_TEXT_BYTES} bytes"
        );
        assert!(
            !output.stderr.is_empty(),
            "{context}: error stderr should be nonempty"
        );
        assert!(
            output.stderr.len() < MAX_TEXT_BYTES,
            "{context}: error stderr should be shorter than {MAX_TEXT_BYTES} bytes"
        );
        assert!(
            !expected_diagnostics.is_empty(),
            "{context}: expected diagnostics are required"
        );
        assert!(
            expected_diagnostics.len() < MAX_TEXT_BYTES,
            "{context}: expected diagnostic list should contain fewer than {MAX_TEXT_BYTES} items"
        );
        for diagnostic in expected_diagnostics {
            assert!(
                !diagnostic.is_empty(),
                "{context}: expected diagnostic should be nonempty"
            );
            assert!(
                diagnostic.len() < MAX_TEXT_BYTES,
                "{context}: expected diagnostic should be shorter than {MAX_TEXT_BYTES} bytes"
            );
        }
        assert!(!context.is_empty(), "CLI error context should be nonempty");
        assert!(
            context.len() < MAX_TEXT_BYTES,
            "CLI error context should be shorter than {MAX_TEXT_BYTES} bytes"
        );
        assert_eq!(
            output.status.code(),
            Some(2),
            "{context}: stderr={}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            output.stdout.is_empty(),
            "{context}: stdout should be empty"
        );
        let stderr = String::from_utf8(output.stderr).expect("diagnostic should be UTF-8");
        for fragment in expected_diagnostics.iter().copied().chain([
            "Usage: nanalogue_sim_bam",
            "For more information, try '--help'.",
        ]) {
            assert!(stderr.contains(fragment), "{context}: missing {fragment:?}");
        }
    }

    /// Generate one deterministic alphanumeric string that never starts with `-`.
    fn random_token(rng: &mut StdRng) -> String {
        let alphabet = b"ABCDEFGIJKLMNOPQRSTUWXYZabcdefgijklmnopqrstuvwxyz0123456789";
        std::iter::repeat_with(|| {
            let offset = rng.random_range(0..alphabet.len());
            char::from(
                *alphabet
                    .get(offset)
                    .expect("generated alphabet offset should be valid"),
            )
        })
        .take(8)
        .collect()
    }

    /// Visit every positional permutation using iterative Heap's algorithm.
    fn visit_permutations<T>(values: &mut [T], mut visit: impl FnMut(&[T])) {
        assert!(
            !values.is_empty() && values.len() <= MAX_PERMUTATION_ARGUMENTS,
            "permutation input should contain between 1 and {MAX_PERMUTATION_ARGUMENTS} values"
        );
        let mut counters = vec![0; values.len()];
        visit(values);

        let mut index = 0;
        let mut loop_counter = 0usize;
        while let Some(&counter) = counters.get(index) {
            loop_counter = loop_counter
                .checked_add(1)
                .expect("permutation loop counter should not overflow");
            assert!(
                loop_counter <= MAX_PERMUTATION_LOOPS,
                "permutation generation should finish within {MAX_PERMUTATION_LOOPS} loops"
            );
            if counter < index {
                let swap_with = if index.is_multiple_of(2) { 0 } else { counter };
                values.swap(swap_with, index);
                visit(values);
                let current = counters
                    .get_mut(index)
                    .expect("permutation counter index should remain valid");
                *current = current
                    .checked_add(1)
                    .expect("permutation counter should not overflow");
                index = 0;
            } else {
                *counters
                    .get_mut(index)
                    .expect("permutation counter index should remain valid") = 0;
                index = index
                    .checked_add(1)
                    .expect("permutation index should not overflow");
            }
        }
    }

    /// Exercise the complete 4-by-7 matrix and one excess-positional case.
    #[test]
    fn simulator_cli_matrix() {
        let expected = ExpectedOutputs::load();
        for case_index in 0..29 {
            run_case(case_index, &expected);
        }
    }

    /// Run one of the 29 explicit CLI cases through the compiled executable.
    #[expect(
        clippy::integer_division,
        clippy::integer_division_remainder_used,
        clippy::too_many_lines,
        reason = "one cohesive function verifies the matrix and excess-positional case"
    )]
    fn run_case(case_index: usize, expected: &ExpectedOutputs) {
        assert!(
            case_index <= 28,
            "simulator CLI case index must be at most 28"
        );

        if case_index == 28 {
            let mut rng = StdRng::seed_from_u64(CLI_MATRIX_SEED);
            let args = std::iter::repeat_with(|| OsString::from(random_token(&mut rng)))
                .take(4)
                .collect::<Vec<_>>();
            let unexpected = args
                .get(3)
                .expect("four-positional case should have a fourth argument")
                .to_string_lossy()
                .into_owned();
            let context = format!(
                "seed={CLI_MATRIX_SEED:#x}, case={case_index}, argv={}",
                args.iter()
                    .map(|arg| arg.to_string_lossy())
                    .collect::<Vec<_>>()
                    .join(" ")
            );
            let output = execute(&args);
            assert_cli_error(output, &["unexpected argument", &unexpected], &context);
            return;
        }

        let positional_count = case_index / 7;
        let column = case_index % 7;
        let temp = (positional_count == 3).then(|| TempDir::new(&format!("case_{case_index}")));
        let paths = temp.as_ref().map(SimulationPaths::new);
        let mut rng = StdRng::seed_from_u64(CLI_MATRIX_SEED);
        let mut args = if let Some(simulation_paths) = paths.as_ref() {
            vec![
                simulation_paths.input.as_os_str().to_owned(),
                simulation_paths.bam.as_os_str().to_owned(),
                simulation_paths.fasta.as_os_str().to_owned(),
            ]
        } else {
            std::iter::repeat_with(|| OsString::from(random_token(&mut rng)))
                .take(positional_count)
                .collect()
        };
        let unknown_suffix = random_token(&mut rng);
        let (action_args, diagnosed_unknown): (Vec<OsString>, Option<String>) = match column {
            0 => (vec![OsString::from("-h")], None),
            1 => (vec![OsString::from("--help")], None),
            2 => (vec![OsString::from("-V")], None),
            3 => (vec![OsString::from("--version")], None),
            4 => {
                let actions: &[&str] = match positional_count {
                    0 => &["-h", "--version"],
                    1 => &["--version", "--help"],
                    2 => &["-V", "-h"],
                    3 => &["--help", "-V"],
                    _ => unreachable!("matrix has only four positional rows"),
                };
                (actions.iter().map(OsString::from).collect(), None)
            }
            5 if matches!(positional_count, 0 | 2) => {
                let first = unknown_suffix
                    .chars()
                    .next()
                    .expect("unknown suffix should be nonempty");
                (
                    vec![OsString::from(format!("-{unknown_suffix}"))],
                    Some(format!("-{first}")),
                )
            }
            5 => {
                let unknown = format!("--unknown-{unknown_suffix}");
                (vec![OsString::from(&unknown)], Some(unknown))
            }
            6 => (Vec::new(), None),
            _ => unreachable!("case index maps to a column from zero through six"),
        };
        args.extend(action_args);

        if case_index == 27 {
            let context = format!(
                "seed={CLI_MATRIX_SEED:#x}, case={case_index}, argv={}",
                args.iter()
                    .map(|arg| arg.to_string_lossy())
                    .collect::<Vec<_>>()
                    .join(" ")
            );
            let output = execute(&args);
            assert_simulation_output(
                &output,
                paths.as_ref().expect("row three should have output paths"),
                &context,
            );
            return;
        }

        if args.is_empty() {
            assert_eq!(
                case_index, 6,
                "only the zero-positional, no-flag case should have no arguments"
            );
            let context = format!("seed={CLI_MATRIX_SEED:#x}, case={case_index}, argv=");
            let output = execute(&args);
            assert_cli_error(output, &["required arguments were not provided"], &context);
            return;
        }

        let expected_permutation_count = (1..=args.len()).product::<usize>();
        let mut visited_permutations = HashSet::new();
        visit_permutations(&mut args, |permutation| {
            assert!(
                visited_permutations.insert(permutation.to_vec()),
                "case {case_index} generated a duplicate argument permutation"
            );
            let context = format!(
                "seed={CLI_MATRIX_SEED:#x}, case={case_index}, argv={}",
                permutation
                    .iter()
                    .map(|arg| arg.to_string_lossy())
                    .collect::<Vec<_>>()
                    .join(" ")
            );
            let output = execute(permutation);

            match column {
                0 => {
                    assert_display(&output, &[&expected.short_help], &context);
                }
                1 => {
                    assert_display(&output, &[&expected.long_help], &context);
                }
                2 | 3 => {
                    assert_display(&output, &[&expected.version], &context);
                }
                4 => {
                    let applicable_help = if permutation.iter().any(|arg| arg == "-h") {
                        expected.short_help.as_slice()
                    } else {
                        expected.long_help.as_slice()
                    };
                    assert_display(&output, &[applicable_help, &expected.version], &context);
                }
                5 => {
                    let unknown = diagnosed_unknown
                        .as_deref()
                        .expect("unknown column should have a token");
                    assert_cli_error(output, &["unexpected argument", unknown], &context);
                }
                6 if positional_count < 3 => {
                    assert_cli_error(output, &["required arguments were not provided"], &context);
                }
                _ => unreachable!("successful simulation is handled before permutation testing"),
            }

            if let Some(simulation_paths) = paths.as_ref() {
                for path in [
                    &simulation_paths.bam,
                    &simulation_paths.index,
                    &simulation_paths.fasta,
                ] {
                    assert!(
                        !path.exists(),
                        "{context}: {} must not exist",
                        path.display()
                    );
                }
            }
        });
        assert_eq!(
            visited_permutations.len(),
            expected_permutation_count,
            "case {case_index} should exercise every argument permutation"
        );
    }

    /// Validate and decode one successful simulator process output.
    fn assert_simulation_output(output: &Output, paths: &SimulationPaths, context: &str) {
        assert!(
            output.stdout.len() < MAX_TEXT_BYTES,
            "{context}: simulator stdout should be shorter than {MAX_TEXT_BYTES} bytes"
        );
        assert!(
            output.stderr.len() < MAX_TEXT_BYTES,
            "{context}: simulator stderr should be shorter than {MAX_TEXT_BYTES} bytes"
        );
        assert!(
            !EXPECTED_FASTA.is_empty(),
            "expected FASTA should be nonempty"
        );
        assert!(
            EXPECTED_FASTA.len() < MAX_TEXT_BYTES,
            "expected FASTA should be shorter than {MAX_TEXT_BYTES} bytes"
        );
        assert!(!context.is_empty(), "simulation context should be nonempty");
        assert!(
            context.len() < MAX_TEXT_BYTES,
            "simulation context should be shorter than {MAX_TEXT_BYTES} bytes"
        );
        assert_eq!(
            output.status.code(),
            Some(0),
            "{context}: simulator run should succeed: stderr={}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            output.stdout.is_empty(),
            "{context}: successful run should have empty stdout: {}",
            String::from_utf8_lossy(&output.stdout)
        );
        assert!(
            output.stderr.is_empty(),
            "{context}: successful run should have empty stderr: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        for path in [&paths.bam, &paths.index, &paths.fasta] {
            assert!(path.exists(), "{} should exist", path.display());
            assert!(
                fs::metadata(path)
                    .expect("output metadata should be readable")
                    .len()
                    > 0,
                "{} should be nonempty",
                path.display()
            );
        }

        let mut reader = bam::Reader::from_path(&paths.bam).expect("BAM should open");
        assert_eq!(reader.header().target_count(), 1, "expected one reference");
        assert_eq!(reader.header().target_names(), [b"contig_00000".as_slice()]);
        assert_eq!(reader.header().target_len(0), Some(64));
        let records = reader
            .records()
            .collect::<Result<Vec<_>, _>>()
            .expect("BAM records should decode");
        assert_eq!(records.len(), 4, "fixture should produce four records");

        let indexed =
            bam::IndexedReader::from_path(&paths.bam).expect("BAM index should be readable");
        assert_eq!(
            indexed.header().target_count(),
            1,
            "indexed BAM should open"
        );

        let fasta = fs::read(&paths.fasta).expect("FASTA should be readable");
        assert_eq!(fasta, EXPECTED_FASTA, "FASTA should match fixed fixture");
    }
}
