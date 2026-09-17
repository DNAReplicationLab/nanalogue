//! Executable-level tests for CLI output stream contracts.

/// Tests that invoke the compiled CLI executable.
#[cfg(test)]
mod tests {
    use std::process::Command;

    /// Runtime failures report an error, rather than looking like a successful
    /// empty result or a command-line parsing failure.
    #[test]
    fn missing_input_reports_runtime_failure_on_stderr() {
        let output = Command::new(env!("CARGO_BIN_EXE_nanalogue"))
            .args(["read-info", "examples/example_1.bam/missing.bam"])
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
}
