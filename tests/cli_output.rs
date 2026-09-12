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
