#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Executable-level tests for CLI output stream contracts.

/// Tests that invoke the compiled CLI executable.
#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use std::process::{Command, Output};

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
        assert!(output.stdout.is_empty());
        assert!(String::from_utf8_lossy(&output.stderr).starts_with("Error during execution: "));
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
        assert!(output.stderr.is_empty());
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
