#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Process-level tests for the optional BAM viewer command-line interface.

/// Tests that invoke the compiled viewer executable.
#[cfg(test)]
#[cfg(feature = "bam-viewer")]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use std::{
        ffi::OsStr,
        process::{Command, Output},
    };

    /// The first line shared by help output and argument-error usage output.
    const USAGE_LINE: &str = concat!(
        "Usage: nanalogue_bam_viewer <BAM> <CONTIG:START> ",
        "[MOD_TYPE [WINDOW_SIZE individual]]"
    );

    /// Runs the viewer as a child process and captures both output streams.
    fn run_viewer<I, S>(arguments: I) -> Output
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        Command::new(env!("CARGO_BIN_EXE_nanalogue_bam_viewer"))
            .args(arguments)
            .output()
            .expect("viewer executable should run")
    }

    /// Checks the common contract for failures detected during argument parsing.
    fn assert_argument_error<const N: usize>(arguments: [&str; N], expected_error: &str) {
        let output = run_viewer(arguments);
        let stderr = String::from_utf8_lossy(&output.stderr);

        assert_eq!(
            output.status.code(),
            Some(2),
            "argument error should exit with status 2; stderr: {stderr}"
        );
        assert!(
            output.stdout.is_empty(),
            "argument error should not write to stdout"
        );
        assert!(
            stderr.starts_with(&format!("Error: {expected_error}\n")),
            "stderr should begin with the decisive parser error; got: {stderr}"
        );
        assert!(
            stderr.contains(USAGE_LINE),
            "argument error should include usage; got: {stderr}"
        );
    }

    /// Help is a successful request written solely to stdout.
    #[test]
    fn help_prints_usage_and_succeeds() {
        let output = run_viewer(["--help"]);
        let stdout = String::from_utf8_lossy(&output.stdout);

        assert_eq!(
            output.status.code(),
            Some(0),
            "help should exit successfully"
        );
        assert!(output.stderr.is_empty(), "help should not write to stderr");
        assert!(
            stdout.starts_with(USAGE_LINE),
            "help should begin with usage; got: {stdout}"
        );
        assert!(
            stdout.contains("START is zero-based; displayed coordinates are one-based."),
            "help should explain the position convention; got: {stdout}"
        );
    }

    /// The two compulsory positional arguments are diagnosed independently.
    #[test]
    fn missing_bam_or_position_is_an_argument_error() {
        assert_argument_error([], "missing BAM file");
        assert_argument_error(["reads.bam"], "missing CONTIG:START position");
    }

    /// Position syntax and coordinate conversion fail before any BAM is opened.
    #[test]
    fn malformed_positions_are_rejected() {
        assert_argument_error(
            ["reads.bam", "chr1"],
            "position must have the form CONTIG:START",
        );
        assert_argument_error(
            ["reads.bam", "chr1:-7"],
            "START must be a non-negative integer",
        );
    }

    /// Individual mode requires a complete positive-window triplet with its literal keyword.
    #[test]
    fn incomplete_or_bad_individual_mode_is_rejected() {
        assert_argument_error(
            ["reads.bam", "chr1:7", "m", "25"],
            "WINDOW_SIZE and literal 'individual' must be supplied together",
        );
        assert_argument_error(
            ["reads.bam", "chr1:7", "m", "0", "individual"],
            "WINDOW_SIZE must be a positive integer",
        );
        assert_argument_error(
            ["reads.bam", "chr1:7", "m", "25", "table"],
            "expected literal 'individual' after WINDOW_SIZE",
        );
    }

    /// Arguments after a complete individual-mode invocation are not ignored.
    #[test]
    fn trailing_arguments_are_rejected() {
        assert_argument_error(
            ["reads.bam", "chr1:7", "m", "25", "individual", "extra"],
            "unexpected trailing arguments",
        );
    }

    /// Unix permits a non-UTF-8 argument, which must be rejected without lossy parsing.
    #[cfg(unix)]
    #[test]
    fn non_utf8_position_is_rejected() {
        use std::{ffi::OsString, os::unix::ffi::OsStringExt as _};

        let invalid_position = OsString::from_vec(b"chr1:\xff".to_vec());
        let output = run_viewer([OsString::from("reads.bam"), invalid_position]);
        let stderr = String::from_utf8_lossy(&output.stderr);

        assert_eq!(
            output.status.code(),
            Some(2),
            "invalid UTF-8 should be an argument error; stderr: {stderr}"
        );
        assert!(
            output.stdout.is_empty(),
            "invalid UTF-8 should not write to stdout"
        );
        assert!(
            stderr.starts_with("Error: position must be valid UTF-8\n"),
            "stderr should identify the invalid position encoding; got: {stderr}"
        );
        assert!(
            stderr.contains(USAGE_LINE),
            "invalid UTF-8 should include usage; got: {stderr}"
        );
    }

    /// Test-owned `tput` fallback that makes terminal-size discovery deterministic.
    #[cfg(unix)]
    #[derive(Debug)]
    struct FakeTput {
        /// Temporary directory prepended to the child process's search path.
        directory: std::path::PathBuf,
    }

    #[cfg(unix)]
    impl FakeTput {
        /// Creates an executable that reports fixed columns and lines.
        fn new() -> Self {
            use std::os::unix::fs::PermissionsExt as _;

            let directory = std::env::temp_dir()
                .join(format!("nanalogue-viewer-cli-tput-{}", std::process::id()));
            std::fs::create_dir_all(&directory)
                .expect("temporary tput directory should be created");
            let executable = directory.join("tput");
            std::fs::write(
                &executable,
                "#!/bin/sh\ncase \"$1\" in\n  cols) echo 80;;\n  lines) echo 24;;\n  *) exit 1;;\nesac\n",
            )
            .expect("fake tput should be written");
            let mut permissions = std::fs::metadata(&executable)
                .expect("fake tput metadata should be readable")
                .permissions();
            permissions.set_mode(0o755);
            std::fs::set_permissions(executable, permissions)
                .expect("fake tput should be executable");
            Self { directory }
        }

        /// Builds a search path with the fake executable before inherited entries.
        fn child_path(&self) -> std::ffi::OsString {
            let inherited = std::env::var_os("PATH").unwrap_or_default();
            std::env::join_paths(
                std::iter::once(self.directory.clone()).chain(std::env::split_paths(&inherited)),
            )
            .expect("child search path should be valid")
        }
    }

    #[cfg(unix)]
    impl Drop for FakeTput {
        fn drop(&mut self) {
            drop(std::fs::remove_dir_all(&self.directory));
        }
    }

    /// A syntactically valid invocation reaches BAM opening and reports that runtime failure.
    #[cfg(unix)]
    #[test]
    fn missing_bam_after_valid_parse_is_a_runtime_error() {
        let missing_bam = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/this-viewer-input-does-not-exist.bam");
        let fake_tput = FakeTput::new();
        let output = Command::new(env!("CARGO_BIN_EXE_nanalogue_bam_viewer"))
            .args([missing_bam.as_os_str(), OsStr::new("chr7:19")])
            .env("PATH", fake_tput.child_path())
            .env_remove("TERM")
            .output()
            .expect("viewer executable should run");
        let stderr = String::from_utf8_lossy(&output.stderr);

        assert_eq!(
            output.status.code(),
            Some(1),
            "BAM-open failure should exit with status 1; stderr: {stderr}"
        );
        assert!(
            output.stdout.is_empty(),
            "BAM-open failure should not write to stdout"
        );
        assert!(
            stderr.starts_with("Error: rust_htslib error: `file not found:"),
            "stderr should identify the BAM-open failure; got: {stderr}"
        );
        assert!(
            stderr.contains("this-viewer-input-does-not-exist.bam"),
            "stderr should name the missing BAM; got: {stderr}"
        );
        assert!(
            !stderr.contains(USAGE_LINE),
            "runtime failures should not be presented as usage errors; got: {stderr}"
        );
    }
}
