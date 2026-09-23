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

    /// Checks a captured failure detected during argument parsing.
    fn assert_argument_error_output(output: &Output, expected_error: &str) {
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

    /// Checks the common contract for failures detected during argument parsing.
    fn assert_argument_error<const N: usize>(arguments: [&str; N], expected_error: &str) {
        let output = run_viewer(arguments);
        assert_argument_error_output(&output, expected_error);
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

    /// The short help alias has the same successful stream and status contract.
    #[test]
    fn short_help_prints_usage_and_succeeds() {
        let output = run_viewer(["-h"]);
        let stdout = String::from_utf8_lossy(&output.stdout);

        assert_eq!(
            output.status.code(),
            Some(0),
            "short help should exit successfully"
        );
        assert!(
            output.stderr.is_empty(),
            "short help should not write to stderr"
        );
        assert!(
            stdout.starts_with(USAGE_LINE),
            "short help should begin with usage; got: {stdout}"
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
            ["reads.bam", ":7"],
            "position must have the form CONTIG:START",
        );
        assert_argument_error(
            ["reads.bam", "chr1:"],
            "position must have the form CONTIG:START",
        );
        assert_argument_error(
            ["reads.bam", "chr1:-7"],
            "START must be a non-negative integer",
        );
        assert_argument_error(
            ["reads.bam", "chr1:4294967296"],
            "START must be a non-negative integer",
        );
    }

    /// Modification types report their own syntax errors before mode parsing.
    #[test]
    fn malformed_modification_types_are_rejected() {
        assert_argument_error(["reads.bam", "chr1:7", ""], "empty mod type: ``");
        assert_argument_error(["reads.bam", "chr1:7", "@123"], "invalid mod type: `@123`");
        assert_argument_error(
            ["reads.bam", "chr1:7", "1114112"],
            "invalid mod type: `1114112`",
        );
        assert_argument_error(
            ["reads.bam", "chr1:7", "4294967296"],
            "integer parsing error: `number too large to fit in target type`",
        );
    }

    /// Window parsing distinguishes malformed and overflowing values from zero.
    #[test]
    fn malformed_or_overflowing_windows_are_rejected() {
        for window in ["", "seven", "4294967296"] {
            assert_argument_error(
                ["reads.bam", "chr1:7", "m", window, "individual"],
                "WINDOW_SIZE must be a positive integer",
            );
        }
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
        assert_argument_error_output(&output, "position must be valid UTF-8");
    }

    /// Unix-only non-UTF-8 modification types fail before modification parsing.
    #[cfg(unix)]
    #[test]
    fn non_utf8_modification_type_is_rejected() {
        use std::{ffi::OsString, os::unix::ffi::OsStringExt as _};

        let invalid_mod_type = OsString::from_vec(vec![0xff]);
        let output = run_viewer([
            OsString::from("reads.bam"),
            OsString::from("chr1:7"),
            invalid_mod_type,
        ]);
        assert_argument_error_output(&output, "MOD_TYPE must be valid UTF-8");
    }

    /// Unix-only non-UTF-8 windows fail before integer parsing.
    #[cfg(unix)]
    #[test]
    fn non_utf8_window_size_is_rejected() {
        use std::{ffi::OsString, os::unix::ffi::OsStringExt as _};

        let invalid_window = OsString::from_vec(vec![0xff]);
        let output = run_viewer([
            OsString::from("reads.bam"),
            OsString::from("chr1:7"),
            OsString::from("m"),
            invalid_window,
            OsString::from("individual"),
        ]);
        assert_argument_error_output(&output, "WINDOW_SIZE must be valid UTF-8");
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
            use std::sync::atomic::{AtomicU32, Ordering};

            static NEXT_DIRECTORY_ID: AtomicU32 = AtomicU32::new(0);

            let directory = std::env::temp_dir().join(format!(
                "nanalogue-viewer-cli-tput-{}-{}",
                std::process::id(),
                NEXT_DIRECTORY_ID.fetch_add(1, Ordering::Relaxed)
            ));
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

    /// Valid table and individual options both proceed beyond argument parsing.
    #[cfg(unix)]
    #[test]
    fn valid_optional_modes_reach_bam_opening() {
        let missing_bam = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/this-viewer-input-does-not-exist.bam");
        let fake_tput = FakeTput::new();

        for optional_arguments in [&["472232"][..], &["m", "4294967295", "individual"][..]] {
            let output = Command::new(env!("CARGO_BIN_EXE_nanalogue_bam_viewer"))
                .arg(missing_bam.as_os_str())
                .arg("chr7:4294967295")
                .args(optional_arguments)
                .env("PATH", fake_tput.child_path())
                .env_remove("TERM")
                .output()
                .expect("viewer executable should run");
            let stderr = String::from_utf8_lossy(&output.stderr);

            assert_eq!(
                output.status.code(),
                Some(1),
                "valid options should reach BAM opening; stderr: {stderr}"
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
                "runtime error should name the missing BAM; got: {stderr}"
            );
            assert!(
                !stderr.contains(USAGE_LINE),
                "valid options should not produce a usage error; got: {stderr}"
            );
        }
    }

    /// Result captured from a viewer session running in a real Linux pseudo-terminal.
    #[cfg(target_os = "linux")]
    #[derive(Debug)]
    struct PtySession {
        status: std::process::ExitStatus,
        output: Vec<u8>,
        stderr: Vec<u8>,
        input_error: Option<String>,
    }

    /// Runs the viewer under util-linux `script`, with an inner GNU `timeout`
    /// directly supervising the viewer and an outer timeout supervising `script`.
    #[cfg(target_os = "linux")]
    #[expect(
        clippy::too_many_lines,
        reason = "the PTY lifecycle is clearer in one test harness than split across stateful helpers"
    )]
    fn run_viewer_in_pty(bam_name: &str, cols: u16, key_chunks: &[&[u8]]) -> PtySession {
        use std::{
            fs::File,
            io::Write as _,
            process::Stdio,
            sync::atomic::{AtomicU32, Ordering},
            thread,
            time::{Duration, Instant},
        };

        static NEXT_SESSION_ID: AtomicU32 = AtomicU32::new(0);

        let session_id = NEXT_SESSION_ID.fetch_add(1, Ordering::Relaxed);
        let file_prefix = std::env::temp_dir().join(format!(
            "nanalogue-viewer-pty-{}-{session_id}",
            std::process::id()
        ));
        let output_path = file_prefix.with_extension("out");
        let stderr_path = file_prefix.with_extension("err");
        let output_file = File::create(&output_path).expect("PTY output file should be created");
        let stderr_file = File::create(&stderr_path).expect("PTY stderr file should be created");
        let bam = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("examples")
            .join(bam_name);

        // The command is constant: paths travel through environment variables,
        // so spaces or shell metacharacters in the checkout cannot alter it.
        let mut child = Command::new("timeout")
            .args([
                "--kill-after=3s",
                "12s",
                "script",
                "--quiet",
                "--return",
                "--flush",
                "--echo",
                "never",
                "--output-limit",
                "2M",
                "--command",
                concat!(
                    "stty rows 24 cols \"$NANALOGUE_COLS\"; ",
                    "before=$(stty -g) || exit 125; ",
                    "timeout --foreground --kill-after=1s 7s \"$NANALOGUE_VIEWER\" ",
                    "\"$NANALOGUE_BAM\" dummyI:0 m; status=$?; ",
                    "after=$(stty -g) || exit 125; ",
                    "if [ \"$before\" = \"$after\" ]; then ",
                    "printf '\\nNANALOGUE_STTY_RESTORED=yes\\n'; else ",
                    "printf '\\nNANALOGUE_STTY_RESTORED=no\\n'; fi; ",
                    "exit \"$status\""
                ),
                "/dev/null",
            ])
            .env(
                "NANALOGUE_VIEWER",
                env!("CARGO_BIN_EXE_nanalogue_bam_viewer"),
            )
            .env("NANALOGUE_BAM", bam)
            .env("NANALOGUE_COLS", cols.to_string())
            .stdin(Stdio::piped())
            .stdout(Stdio::from(output_file))
            .stderr(Stdio::from(stderr_file))
            .spawn()
            .expect("util-linux script and GNU timeout should be installed");

        let mut input = child.stdin.take().expect("PTY input pipe should exist");
        let mut input_error = None;
        let startup_deadline = Instant::now()
            .checked_add(Duration::from_secs(3))
            .expect("short startup timeout should fit in Instant");
        loop {
            let started = terminal_text(&std::fs::read(&output_path).unwrap_or_default())
                .contains("nanalogue BAM viewer");
            if started {
                break;
            }
            if Instant::now() >= startup_deadline {
                input_error = Some(String::from("initial viewer frame did not render"));
                break;
            }
            thread::sleep(Duration::from_millis(20));
        }
        let mut cancel_controls_baseline = None;
        for keys in key_chunks {
            let prompts_before_goto = (keys == b"g").then(|| {
                terminal_text(&std::fs::read(&output_path).unwrap_or_default())
                    .matches("Go to CONTIG:START: ")
                    .count()
            });
            if let Err(error) = input.write_all(keys).and_then(|()| input.flush()) {
                input_error = Some(format!("viewer keys could not be written: {error}"));
                break;
            }
            if let Some(previous_prompt_frames) = prompts_before_goto {
                let prompt_deadline = Instant::now()
                    .checked_add(Duration::from_secs(2))
                    .expect("short prompt timeout should fit in Instant");
                loop {
                    let text = terminal_text(&std::fs::read(&output_path).unwrap_or_default());
                    if text.matches("Go to CONTIG:START: ").count() > previous_prompt_frames {
                        cancel_controls_baseline = Some(text.matches("h/l 10 bp").count());
                        break;
                    }
                    if Instant::now() >= prompt_deadline {
                        input_error = Some(String::from(
                            "goto key did not render the prompt before cancellation",
                        ));
                        break;
                    }
                    thread::sleep(Duration::from_millis(20));
                }
            }
            if keys == b"\x1b" {
                let Some(previous_control_frames) = cancel_controls_baseline else {
                    input_error = Some(String::from(
                        "goto prompt was not observed before cancellation",
                    ));
                    break;
                };
                let redraw_deadline = Instant::now()
                    .checked_add(Duration::from_secs(2))
                    .expect("short redraw timeout should fit in Instant");
                loop {
                    let controls_after_cancel =
                        terminal_text(&std::fs::read(&output_path).unwrap_or_default())
                            .matches("h/l 10 bp")
                            .count();
                    if controls_after_cancel > previous_control_frames {
                        break;
                    }
                    if Instant::now() >= redraw_deadline {
                        input_error = Some(String::from(
                            "goto cancellation did not redraw the controls before quit",
                        ));
                        break;
                    }
                    thread::sleep(Duration::from_millis(20));
                }
            }
            if input_error.is_some() {
                break;
            }
            thread::sleep(Duration::from_millis(80));
        }

        // Keep stdin open while waiting: script maps pipe EOF to Ctrl-D, which
        // must not provide a second, accidental way for either quit test to pass.
        let outer_deadline = Instant::now()
            .checked_add(Duration::from_secs(17))
            .expect("short outer timeout should fit in Instant");
        let status = loop {
            if let Some(status) = child
                .try_wait()
                .expect("PTY child status should be readable")
            {
                break status;
            }
            if Instant::now() >= outer_deadline {
                drop(child.kill());
                input_error = Some(String::from(
                    "GNU timeout failed to terminate the PTY process group",
                ));
                break child.wait().expect("killed PTY child should be reaped");
            }
            thread::sleep(Duration::from_millis(20));
        };
        drop(input);

        let output = std::fs::read(&output_path).expect("PTY output should be readable");
        let stderr = std::fs::read(&stderr_path).expect("PTY stderr should be readable");
        std::fs::remove_file(output_path).expect("PTY output should be removed");
        std::fs::remove_file(stderr_path).expect("PTY stderr should be removed");
        PtySession {
            status,
            output,
            stderr,
            input_error,
        }
    }

    /// Removes CSI escape sequences from renderer output while retaining the
    /// printable cell contents in draw order for decisive text assertions.
    #[cfg(target_os = "linux")]
    fn terminal_text(output: &[u8]) -> String {
        let mut printable = Vec::with_capacity(output.len());
        let mut bytes = output.iter().copied();
        while let Some(byte) = bytes.next() {
            if byte == 0x1b {
                if bytes.next() == Some(b'[') {
                    for sequence_byte in bytes.by_ref() {
                        if (0x40..=0x7e).contains(&sequence_byte) {
                            break;
                        }
                    }
                }
                continue;
            }
            if byte >= b' ' {
                printable.push(byte);
            }
        }
        String::from_utf8_lossy(&printable).into_owned()
    }

    /// Checks both the rendered frame and `TerminalGuard`'s normal-exit cleanup.
    #[cfg(target_os = "linux")]
    fn assert_clean_pty_exit(session: &PtySession) {
        assert_eq!(
            session.status.code(),
            Some(0),
            "PTY viewer should exit successfully; stderr: {}",
            String::from_utf8_lossy(&session.stderr)
        );
        assert!(
            session.stderr.is_empty(),
            "successful session has no stderr"
        );
        assert_eq!(
            session.input_error, None,
            "all keys and synchronization points should complete"
        );
        assert!(
            session.output.starts_with(b"\x1b[?1049h\x1b[?25l"),
            "viewer should enter the alternate screen and hide the cursor"
        );
        assert!(
            session
                .output
                .windows(8)
                .any(|bytes| bytes == b"\x1b[?2026h"),
            "GhosttyRenderer should use synchronized terminal updates"
        );
        assert!(
            session
                .output
                .windows(14)
                .any(|bytes| bytes == b"\x1b[?25h\x1b[?1049l"),
            "TerminalGuard should show the cursor and leave the alternate screen"
        );
        assert!(
            terminal_text(&session.output).contains("NANALOGUE_STTY_RESTORED=yes"),
            "TerminalGuard should restore the PTY's raw-mode settings"
        );
    }

    /// A queued quit key is consumed only after raw mode starts and exits cleanly.
    #[cfg(target_os = "linux")]
    #[test]
    fn pty_immediate_quit_draws_and_restores_terminal() {
        let session = run_viewer_in_pty("example_1.bam", 100, &[b"q"]);
        assert_clean_pty_exit(&session);

        let text = terminal_text(&session.output);
        assert!(text.contains("nanalogue BAM viewer"));
        assert!(
            text.contains("h/l 81 bp"),
            "PTY width should be fixed at 100"
        );
        assert!(text.contains("g goto  r full IDs  i show ins  q quit"));
    }

    /// Drives display toggles, arrow navigation, and goto cancellation before quitting.
    #[cfg(target_os = "linux")]
    #[test]
    fn pty_interactive_flow_dispatches_keys_and_quits() {
        let session = run_viewer_in_pty(
            "example_3.bam",
            29,
            &[b"r", b"i", b"\x1b[C", b"\x1b[B", b"g", b"\x1b", b"q"],
        );
        assert_clean_pty_exit(&session);

        let text = terminal_text(&session.output);
        let initial_position = text
            .find("dummyI:1-10")
            .expect("initial position should draw");
        let moved_position = text
            .find("dummyI:11-20")
            .expect("Right should navigate one genomic window");
        assert!(
            initial_position < moved_position,
            "navigation should move right"
        );
        assert!(
            text.contains("read001"),
            "a checked-in BAM read should render"
        );
        assert!(
            text.contains("Go to CONTIG:START: "),
            "g should render the goto prompt before Escape cancels it"
        );
        assert!(
            text.matches("nanalogue BAM viewer").count() >= 6,
            "each dispatched key should redraw the viewer"
        );
    }
}
