#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Executable-level tests for CLI output stream contracts.

/// Tests that invoke the compiled CLI executable.
#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use std::process::{Command, Output};

    /// Runs the compiled main executable with the supplied arguments.
    fn run_nanalogue(args: &[&str]) -> Output {
        Command::new(env!("CARGO_BIN_EXE_nanalogue"))
            .args(args)
            .output()
            .expect("nanalogue executable should run")
    }

    /// Runs the `nanalogue` executable with `args`, writes `stdin_bytes` to its
    /// standard input, and returns the captured output.
    fn run_nanalogue_with_stdin<const N: usize>(args: [&str; N], stdin_bytes: &[u8]) -> Output {
        use std::io::Write as _;
        use std::process::Stdio;

        let mut child = Command::new(env!("CARGO_BIN_EXE_nanalogue"))
            .args(args)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .expect("nanalogue executable should spawn");
        child
            .stdin
            .take()
            .expect("stdin pipe should be available")
            .write_all(stdin_bytes)
            .expect("BAM bytes should be written to stdin");
        child
            .wait_with_output()
            .expect("nanalogue executable should run")
    }

    /// Writes the smallest mapped SAM record whose sequence field is `*`.
    /// The returned directory is unique so parallel tests cannot share a fixture.
    fn zero_length_sam_fixture() -> (std::path::PathBuf, std::path::PathBuf) {
        let root = std::env::temp_dir().join(format!(
            "nanalogue_zero_length_cli_{}",
            nanalogue_core::uuid::v4_random()
        ));
        std::fs::create_dir_all(&root).expect("fixture directory should be created");
        let sam = root.join("zero_length.sam");
        std::fs::write(
            &sam,
            concat!(
                "@HD\tVN:1.6\tSO:unsorted\n",
                "@SQ\tSN:ctg1\tLN:100\n",
                "secondary\t256\tctg1\t11\t37\t8M\t*\t0\t0\t*\t*\n"
            ),
        )
        .expect("zero-length SAM fixture should be written");
        (root, sam)
    }

    /// Runtime failures report an error, rather than looking like a successful
    /// empty result or a command-line parsing failure.
    #[test]
    fn missing_input_reports_runtime_failure_on_stderr() {
        let output = Command::new(env!("CARGO_BIN_EXE_nanalogue"))
            .args([
                "read-info",
                concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/examples/example_1.bam/missing.bam"
                ),
            ])
            .output()
            .expect("nanalogue executable should run");

        assert_eq!(output.status.code(), Some(1));
        assert_eq!(output.stdout, Vec::<u8>::new());
        assert!(String::from_utf8_lossy(&output.stderr).starts_with("Error during execution: "));
    }

    /// Commands using shared input filtering must see no records. Covers both
    /// read-table commands, `read-stats`, `read-info`, both window commands,
    /// all six `find-modified-reads` criteria, and `peek`. `peek` skips
    /// zero-length records locally. Headers and zero/empty summaries are
    /// intentionally retained; they do not describe the filtered SAM record.
    #[test]
    #[expect(
        clippy::too_many_lines,
        reason = "one test intentionally documents every CLI command's empty-input contract"
    )]
    fn zero_length_records_are_filtered_by_all_bam_reading_commands() {
        let (root, sam_file) = zero_length_sam_fixture();
        let sam_path = sam_file.to_string_lossy().into_owned();

        for command in ["read-table-show-mods", "read-table-hide-mods"] {
            let output = run_nanalogue(&[command, &sam_path]);
            assert_eq!(
                output.status.code(),
                Some(1),
                "{command} should see no records"
            );
            assert!(
                String::from_utf8_lossy(&output.stderr)
                    .contains("No records found as input for analysis."),
                "{command} must fail as an empty input, not while parsing the record"
            );
            let stdout = String::from_utf8(output.stdout).expect("table output should be UTF-8");
            assert!(
                stdout
                    .lines()
                    .last()
                    .is_some_and(|line| line.starts_with("read_id\t")),
                "the final line must be the table header, not a data row: {stdout}"
            );
            assert!(
                !stdout.contains("secondary"),
                "record-derived output: {stdout}"
            );
        }

        let stats_output = run_nanalogue(&["read-stats", &sam_path]);
        assert!(
            stats_output.status.success(),
            "read-stats failed: {}",
            String::from_utf8_lossy(&stats_output.stderr)
        );
        let stats_text =
            String::from_utf8(stats_output.stdout).expect("statistics should be UTF-8");
        assert!(
            stats_text.starts_with("key\tvalue\n"),
            "statistics must include its header: {stats_text}"
        );
        assert!(
            stats_text.lines().skip(1).all(|line| line.ends_with("\t0")),
            "statistics must be the empty-input zero summary: {stats_text}"
        );
        assert!(
            stats_text.contains("n_secondary_alignments\t0"),
            "record must not contribute to stats: {stats_text}"
        );

        let read_info = run_nanalogue(&["read-info", &sam_path]);
        assert!(
            read_info.status.success(),
            "read-info failed: {}",
            String::from_utf8_lossy(&read_info.stderr)
        );
        let read_info_json = serde_json::from_slice::<serde_json::Value>(&read_info.stdout)
            .expect("read-info output should be JSON");
        assert_eq!(
            read_info_json,
            serde_json::json!([]),
            "read-info must be empty"
        );

        for command in ["window-dens", "window-grad"] {
            let output = run_nanalogue(&[command, "--win", "8", "--step", "1", &sam_path]);
            assert_eq!(
                output.status.code(),
                Some(1),
                "{command} should see no records"
            );
            assert!(
                String::from_utf8_lossy(&output.stderr)
                    .contains("No records found as input for analysis."),
                "{command} must fail as an empty input, not while parsing the record"
            );
            let stdout = String::from_utf8(output.stdout).expect("window output should be UTF-8");
            assert_eq!(
                stdout.lines().count(),
                1,
                "{command} emitted window data: {stdout}"
            );
            assert!(
                stdout.starts_with("#contig\t"),
                "unexpected header: {stdout}"
            );
        }

        for (criterion, args) in [
            ("all-dens-between", vec!["--dens-limits", "0,1"]),
            ("any-dens-above", vec!["--high", "0"]),
            ("any-dens-below", vec!["--low", "1"]),
            (
                "any-dens-below-and-any-dens-above",
                vec!["--low", "1", "--high", "0"],
            ),
            ("dens-range-above", vec!["--min-range", "0"]),
            ("any-abs-grad-above", vec!["--min-grad", "0"]),
        ] {
            let mut command = vec![
                "find-modified-reads",
                criterion,
                "--win",
                "8",
                "--step",
                "1",
                "--tag",
                "m",
            ];
            command.extend(args);
            command.push(&sam_path);
            let output = run_nanalogue(&command);
            assert_eq!(
                output.status.code(),
                Some(1),
                "{criterion} should see no records: {}",
                String::from_utf8_lossy(&output.stderr)
            );
            assert!(
                String::from_utf8_lossy(&output.stderr)
                    .contains("No records found as input for analysis."),
                "{criterion} must fail as an empty input, not while parsing the record"
            );
            assert!(output.stdout.is_empty(), "{criterion} emitted a read ID");
        }

        let peek = run_nanalogue(&["peek", &sam_path]);
        assert!(
            peek.status.success(),
            "peek failed: {}",
            String::from_utf8_lossy(&peek.stderr)
        );
        assert_eq!(
            peek.stdout, b"contigs_and_lengths:\nctg1\t100\n\nmodifications:\nNone\n",
            "peek may show header metadata but must not derive modifications from the record"
        );

        std::fs::remove_dir_all(root).expect("fixture directory should be cleaned up");
    }

    /// Closing the peer before spawning avoids racing a consumer such as head:
    /// even a small output must fail on flush, silently with the SIGPIPE code.
    #[cfg(unix)]
    #[test]
    fn closed_stdout_exits_silently_with_broken_pipe_code() {
        use std::{os::fd::OwnedFd, os::unix::net::UnixStream, process::Stdio};

        let (reader, writer) = UnixStream::pair().expect("socket pair should open");
        drop(reader);
        let output = Command::new(env!("CARGO_BIN_EXE_nanalogue"))
            .args([
                "read-info",
                "--detailed",
                concat!(env!("CARGO_MANIFEST_DIR"), "/examples/example_1.bam"),
            ])
            .stdout(Stdio::from(OwnedFd::from(writer)))
            .output()
            .expect("nanalogue executable should run");

        assert_eq!(output.status.code(), Some(141));
        assert_eq!(output.stderr, Vec::<u8>::new());
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

    #[test]
    fn full_sequence_rejects_region_only_display_options() {
        let bam = concat!(env!("CARGO_MANIFEST_DIR"), "/examples/example_1.bam");
        for args in [
            ["read-table-show-mods", "--seq-full", "--show-mod-z", bam],
            [
                "read-table-show-mods",
                "--seq-full",
                "--show-ins-lowercase",
                bam,
            ],
            [
                "read-table-hide-mods",
                "--seq-full",
                "--show-ins-lowercase",
                bam,
            ],
        ] {
            let output = run_nanalogue(&args);
            assert_eq!(output.status.code(), Some(2), "failed args: {args:?}");
            assert!(output.stdout.is_empty(), "failed args: {args:?}");
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert!(stderr.contains("cannot be used with"), "{stderr}");
            assert!(!stderr.contains("panicked"), "{stderr}");
        }
    }

    #[test]
    fn gradient_commands_reject_single_position_windows_once() {
        let bam = concat!(env!("CARGO_MANIFEST_DIR"), "/examples/example_10.sam");
        for args in [
            vec!["window-grad", "--win", "1", "--step", "1", bam],
            vec![
                "find-modified-reads",
                "any-abs-grad-above",
                "--win",
                "1",
                "--step",
                "1",
                "--tag",
                "N",
                "--min-grad",
                "0.1",
                bam,
            ],
        ] {
            let output = run_nanalogue(&args);
            assert_eq!(output.status.code(), Some(1), "failed args: {args:?}");
            assert!(output.stdout.is_empty(), "failed args: {args:?}");
            let stderr = String::from_utf8_lossy(&output.stderr);
            assert!(
                stderr.contains("gradient calculations require --win to be at least 2"),
                "{stderr}"
            );
            assert_eq!(stderr.lines().count(), 1, "{stderr}");
            assert!(!stderr.contains("Warning: Skipping"), "{stderr}");
        }
    }

    /// Piping raw BAM bytes through stdin exercises the unindexed stdin reader
    /// and must reproduce the checked-in statistics byte-for-byte.
    #[test]
    fn read_stats_from_stdin_matches_checked_in_output() {
        let bam_bytes = std::fs::read(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/examples/example_3.bam"
        ))
        .expect("example_3.bam should be readable");

        let output = run_nanalogue_with_stdin(["read-stats", "-"], &bam_bytes);

        assert_eq!(
            output.status.code(),
            Some(0),
            "stdin read-stats should succeed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            output.stderr.is_empty(),
            "stdin read-stats should not write to stderr"
        );
        let expected = std::fs::read(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/examples/example_3_read_stats"
        ))
        .expect("example_3_read_stats should be readable");
        assert_eq!(
            output.stdout, expected,
            "stdin statistics should match the checked-in output exactly"
        );
    }

    /// Strict sequence consumers reject an explicitly included `SEQ=*` record
    /// before command-specific processing without incorrectly mentioning mod
    /// data. Read-table commands are omitted because they intentionally recover
    /// zero-length records; `peek` has no `--include-zero-len` option.
    #[test]
    fn strict_commands_have_neutral_error_for_included_zero_length_read() {
        let (root, sam_path) = zero_length_sam_fixture();
        let path = sam_path
            .to_str()
            .expect("temporary paths must be valid UTF-8");

        let mut outputs = vec![
            run_nanalogue(&["read-info", "--include-zero-len", path]),
            run_nanalogue(&["read-stats", "--include-zero-len", path]),
        ];

        for command in ["window-dens", "window-grad"] {
            outputs.push(run_nanalogue(&[
                command,
                "--include-zero-len",
                "--win",
                "8",
                "--step",
                "1",
                path,
            ]));
        }

        for (criterion, args) in [
            ("all-dens-between", vec!["--dens-limits", "0,1"]),
            ("any-dens-above", vec!["--high", "0"]),
            ("any-dens-below", vec!["--low", "1"]),
            (
                "any-dens-below-and-any-dens-above",
                vec!["--low", "1", "--high", "0"],
            ),
            ("dens-range-above", vec!["--min-range", "0"]),
            ("any-abs-grad-above", vec!["--min-grad", "0"]),
        ] {
            let mut command = vec![
                "find-modified-reads",
                criterion,
                "--include-zero-len",
                "--win",
                "8",
                "--step",
                "1",
                "--tag",
                "m",
            ];
            command.extend(args);
            command.push(path);
            outputs.push(run_nanalogue(&command));
        }

        std::fs::remove_dir_all(&root).expect("temporary directory must be removable");

        for output in outputs {
            assert_eq!(output.status.code(), Some(1));
            assert!(!String::from_utf8_lossy(&output.stdout).contains("secondary"));
            let stderr = String::from_utf8(output.stderr).expect("error output must be UTF-8");
            assert!(stderr.contains("cannot process record with `SEQ=*`, read_id: secondary"));
            assert!(!stderr.contains("mod data"));
        }
    }

    /// With stdin input there is no index to consult, so `--region` is applied
    /// through the per-record filter; the result must still match the indexed
    /// file path and the region must actually narrow the statistics.
    #[test]
    fn read_stats_from_stdin_with_region_filters_records() {
        let bam_bytes = std::fs::read(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/examples/example_3.bam"
        ))
        .expect("example_3.bam should be readable");

        let stdin_output =
            run_nanalogue_with_stdin(["read-stats", "--region", "dummyI:1-11", "-"], &bam_bytes);

        assert_eq!(
            stdin_output.status.code(),
            Some(0),
            "stdin region read-stats should succeed: {}",
            String::from_utf8_lossy(&stdin_output.stderr)
        );
        assert!(
            stdin_output.stderr.is_empty(),
            "stdin region read-stats should not write to stderr"
        );
        let stdin_stdout = String::from_utf8_lossy(&stdin_output.stdout);
        assert!(
            stdin_stdout.starts_with("key\tvalue\n"),
            "stdin region output should be a statistics table"
        );
        assert!(
            stdin_stdout.contains("n_primary_alignments\t2\n"),
            "only reads passing through the region should be counted, got: {stdin_stdout}"
        );

        let file_output = Command::new(env!("CARGO_BIN_EXE_nanalogue"))
            .args([
                "read-stats",
                "--region",
                "dummyI:1-11",
                concat!(env!("CARGO_MANIFEST_DIR"), "/examples/example_3.bam"),
            ])
            .output()
            .expect("nanalogue executable should run");
        assert_eq!(
            file_output.status.code(),
            Some(0),
            "indexed region read-stats should succeed: {}",
            String::from_utf8_lossy(&file_output.stderr)
        );
        assert_eq!(
            file_output.stdout, stdin_output.stdout,
            "stdin and indexed region filtering should agree"
        );

        let full_file_stats = std::fs::read(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/examples/example_3_read_stats"
        ))
        .expect("example_3_read_stats should be readable");
        assert_ne!(
            stdin_output.stdout, full_file_stats,
            "the region should narrow the statistics compared to the whole file"
        );
    }

    /// A missing JSON config is a runtime failure reported on stderr, not a
    /// command-line parsing failure, and no output files are created.
    #[test]
    fn simulator_missing_json_reports_runtime_failure() {
        let temp_dir = std::env::temp_dir().join(format!(
            "nanalogue_cli_output_sim_missing_{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&temp_dir).expect("temp dir should be created");
        let json_path = temp_dir.join("missing.json");
        let bam_path = temp_dir.join("out.bam");
        let fasta_path = temp_dir.join("out.fa");

        let output = Command::new(env!("CARGO_BIN_EXE_nanalogue_sim_bam"))
            .arg(&json_path)
            .arg(&bam_path)
            .arg(&fasta_path)
            .output()
            .expect("nanalogue_sim_bam executable should run");

        assert_eq!(
            output.status.code(),
            Some(1),
            "a missing config should exit with status 1"
        );
        assert!(
            output.stdout.is_empty(),
            "a failed simulation should not write to stdout"
        );
        assert!(
            String::from_utf8_lossy(&output.stderr).starts_with("Error during execution: "),
            "a missing config should report a runtime failure on stderr, got: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            !bam_path.exists() && !fasta_path.exists(),
            "no outputs should be created when the config cannot be read"
        );

        std::fs::remove_dir_all(&temp_dir).expect("temp dir should be cleaned up");
    }

    /// A minimal valid simulator config creates the requested BAM and FASTA
    /// outputs, with the BAM written as BGZF.
    #[test]
    fn simulator_creates_bam_and_fasta_from_minimal_config() {
        let temp_dir = std::env::temp_dir().join(format!(
            "nanalogue_cli_output_sim_valid_{}",
            std::process::id()
        ));
        std::fs::create_dir_all(&temp_dir).expect("temp dir should be created");
        let config_path = temp_dir.join("config.json");
        std::fs::write(
            &config_path,
            r#"{
                "contigs": {"number": 1, "len_range": [100, 100], "repeated_seq": "ACGT"},
                "reads": [{
                    "number": 5,
                    "mapq_range": [10, 20],
                    "base_qual_range": [20, 30],
                    "len_range": [0.5, 0.5]
                }],
                "seed": 7
            }"#,
        )
        .expect("config should be writable");
        let bam_path = temp_dir.join("out.bam");
        let fasta_path = temp_dir.join("out.fa");

        let output = Command::new(env!("CARGO_BIN_EXE_nanalogue_sim_bam"))
            .arg(&config_path)
            .arg(&bam_path)
            .arg(&fasta_path)
            .output()
            .expect("nanalogue_sim_bam executable should run");

        assert_eq!(
            output.status.code(),
            Some(0),
            "a valid config should succeed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            output.stderr.is_empty(),
            "a successful simulation should not write to stderr"
        );
        let bam = std::fs::read(&bam_path).expect("BAM output should exist");
        let fasta = std::fs::read(&fasta_path).expect("FASTA output should exist");
        assert!(
            bam.starts_with(&[0x1f, 0x8b]),
            "BAM output should be BGZF-compressed"
        );
        assert!(
            fasta.starts_with(b">"),
            "FASTA output should start with a header"
        );
        assert!(
            temp_dir.join("out.bam.bai").exists(),
            "the alignment index should be created next to the BAM output"
        );

        std::fs::remove_dir_all(&temp_dir).expect("temp dir should be cleaned up");
    }

    /// A corrupt (zero-byte) index forces the sequential fallback even though
    /// an index file exists, and region filtering then produces the same
    /// statistics as an unindexed copy.
    #[test]
    fn read_stats_region_falls_back_from_invalid_index() {
        let invalid_index_output = Command::new(env!("CARGO_BIN_EXE_nanalogue"))
            .args([
                "read-stats",
                "--region",
                "dummyI:1-22",
                concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/examples/example_1_copy_invalid_index.bam"
                ),
            ])
            .output()
            .expect("nanalogue executable should run");

        assert_eq!(
            invalid_index_output.status.code(),
            Some(0),
            "the invalid index should fall back to the sequential reader: {}",
            String::from_utf8_lossy(&invalid_index_output.stderr)
        );
        let expected_region_stats = concat!(
            "key\tvalue\n",
            "n_primary_alignments\t1\n",
            "n_secondary_alignments\t0\n",
            "n_supplementary_alignments\t0\n",
            "n_unmapped_reads\t0\n",
            "n_reversed_reads\t0\n",
            "align_len_mean\t8\n",
            "align_len_max\t8\n",
            "align_len_min\t8\n",
            "align_len_median\t8\n",
            "align_len_n50\t8\n",
            "seq_len_mean\t8\n",
            "seq_len_max\t8\n",
            "seq_len_min\t8\n",
            "seq_len_median\t8\n",
            "seq_len_n50\t8\n",
        );
        assert_eq!(
            invalid_index_output.stdout,
            expected_region_stats.as_bytes(),
            "the fallback should still filter by region"
        );
        let warning = "# cannot find index file. region retrieval could be slower.\n";
        assert_eq!(
            invalid_index_output.stderr,
            warning.as_bytes(),
            "the fallback should warn on stderr"
        );

        let no_index_output = Command::new(env!("CARGO_BIN_EXE_nanalogue"))
            .args([
                "read-stats",
                "--region",
                "dummyI:1-22",
                concat!(
                    env!("CARGO_MANIFEST_DIR"),
                    "/examples/example_1_copy_no_index.bam"
                ),
            ])
            .output()
            .expect("nanalogue executable should run");
        assert_eq!(
            no_index_output.stdout, invalid_index_output.stdout,
            "an invalid index and a missing index should agree"
        );
        assert_eq!(
            no_index_output.stderr, invalid_index_output.stderr,
            "an invalid index and a missing index should warn identically"
        );
    }
}
