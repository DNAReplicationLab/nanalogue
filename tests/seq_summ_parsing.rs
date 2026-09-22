#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Integration tests for sequencing-summary TSV parsing reached through the
//! public [`reads_table::run`] entry point.
//!
//! Each malformed summary is written into a fresh temporary directory, the
//! command runs against `examples/example_1.bam`, and the returned error is
//! compared against an independently derived variant and message. Every
//! temporary directory is removed before the test returns.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::constants::reads_table::MAX_SEQ_SUMM_SIZE_PER_LINE;
    use nanalogue_core::constants::shared::MAX_READ_ID_LEN;
    use nanalogue_core::{Error, SeqDisplayOptions, nanalogue_bam_reader, reads_table, uuid};
    use rust_htslib::bam::Read as _;
    use std::fs;
    use std::path::{Path, PathBuf};

    /// Fixture BAM holding three records across three distinct read ids.
    const BAM_PATH: &str = "examples/example_1.bam";

    /// Header accepted by the sequencing-summary parser.
    const HEADER: &str = "read_id\tsequence_length_template\n";

    /// Prefix of every data row that reaches the parser's field extraction.
    const ROW_PREFIX: &str = "read-1\t1234\t";

    /// Creates `seq_summ.tsv` with `contents` in a fresh temporary directory.
    fn write_seq_summ(label: &str, contents: &str) -> (PathBuf, PathBuf) {
        let root =
            std::env::temp_dir().join(format!("nanalogue_seq_summ_{label}_{}", uuid::v4_random()));
        fs::create_dir_all(&root).expect("temporary directory must be creatable");
        let path = root.join("seq_summ.tsv");
        fs::write(&path, contents).expect("sequencing summary must be writable");
        (root, path)
    }

    /// Removes a temporary directory created by [`write_seq_summ`].
    fn remove_temp_dir(root: &Path) {
        fs::remove_dir_all(root).expect("temporary directory must be removable");
    }

    /// Runs `read-table-hide-mods` over [`BAM_PATH`] with the given summary.
    fn run_with_seq_summ(seq_summ_path: &str) -> Result<String, Error> {
        let mut reader = nanalogue_bam_reader(BAM_PATH).expect("fixture BAM must open");
        let mut output = Vec::new();
        reads_table::run(
            &mut output,
            reader.rc_records(),
            None,
            SeqDisplayOptions::No,
            seq_summ_path,
        )?;
        Ok(String::from_utf8(output).expect("reads_table writes UTF-8"))
    }

    /// The stdin marker and a URL are rejected as sequencing-summary inputs.
    #[test]
    fn stdin_marker_and_url_are_rejected_as_non_paths() {
        for (argument, label) in [
            ("-", "stdin marker"),
            ("https://example.com/seq_summ.tsv", "URL"),
        ] {
            let error = run_with_seq_summ(argument)
                .expect_err("a non-path sequencing summary argument must fail");
            let rendered = format!("{error:?}");
            let expected = format!("{argument} does not look like a path");
            assert!(
                matches!(error, Error::InvalidState(message) if message == expected),
                "{label} must be rejected as a non-path, got {rendered}"
            );
        }
    }

    /// A non-numeric `sequence_length_template` names the field and the read.
    #[test]
    fn invalid_sequence_length_names_field_and_read_id() {
        let (root, path) = write_seq_summ("bad_length", &format!("{HEADER}read-1\tnot-a-number\n"));
        let path_str = path.to_str().expect("temporary paths are valid UTF-8");
        let error = run_with_seq_summ(path_str).expect_err("a bad length must fail");
        remove_temp_dir(&root);

        let rendered = format!("{error:?}");
        let expected = format!(
            "sequencing summary tsv parse error in file `{path_str}` for read `read-1`: \
             invalid `sequence_length_template` (invalid digit found in string)"
        );
        assert!(
            matches!(error, Error::InvalidState(message) if message == expected),
            "the parse error must name the field and the read, got {rendered}"
        );
    }

    /// A repeated read id is reported as `InvalidDuplicates` naming the read.
    #[test]
    fn duplicate_read_ids_are_invalid_duplicates() {
        let (root, path) = write_seq_summ(
            "duplicates",
            &format!("{HEADER}read-1\t1234\nread-1\t5678\n"),
        );
        let path_str = path.to_str().expect("temporary paths are valid UTF-8");
        let error = run_with_seq_summ(path_str).expect_err("duplicate read ids must fail");
        remove_temp_dir(&root);

        let rendered = format!("{error:?}");
        let expected = format!("file: {path_str}, read: read-1");
        assert!(
            matches!(error, Error::InvalidDuplicates(message) if message == expected),
            "the second occurrence of a read id must be InvalidDuplicates, got {rendered}"
        );
    }

    /// An empty file and a header-only file fail with distinct messages.
    ///
    /// The coverage plan expected a header-only summary to succeed with no
    /// `bc_len` data, but production deliberately rejects any summary that
    /// yields an empty read map, as also pinned by the
    /// `process_seq_summ_rejects_invalid_preheader_and_empty_input_cases`
    /// unit test in `src/subcommands/reads_table.rs`. This test therefore pins
    /// the actual contract: an empty file reports a missing header, while a
    /// header-only file reports that no reads were found. Omitting the summary
    /// argument is the only way to obtain a run with no `bc_len` data.
    #[test]
    fn empty_and_header_only_files_report_distinct_errors() {
        let (empty_root, empty_path) = write_seq_summ("empty_file", "");
        let empty_str = empty_path
            .to_str()
            .expect("temporary paths are valid UTF-8");
        let empty_error = run_with_seq_summ(empty_str).expect_err("an empty file must fail");
        remove_temp_dir(&empty_root);
        let empty_rendered = format!("{empty_error:?}");
        assert!(
            matches!(empty_error, Error::InvalidState(message) if message == "sequencing summary tsv missing header row"),
            "an empty file must report the missing header, got {empty_rendered}"
        );

        let (header_root, header_path) = write_seq_summ("header_only", HEADER);
        let header_str = header_path
            .to_str()
            .expect("temporary paths are valid UTF-8");
        let header_error = run_with_seq_summ(header_str).expect_err("a header-only file must fail");
        remove_temp_dir(&header_root);
        let header_rendered = format!("{header_error:?}");
        assert!(
            matches!(header_error, Error::InvalidState(message) if message == "sequencing summary TSV did not contain any reads"),
            "a header-only file must report that no reads were found, got {header_rendered}"
        );

        let output = run_with_seq_summ("").expect("a run without a summary must succeed");
        assert!(
            output.starts_with("read_id\talign_length\tsequence_length_template\talignment_type"),
            "a run without a summary must still print the table header, got {output}"
        );
        assert!(
            !output.contains("# seq summ file:"),
            "no summary comment must be printed without a summary argument"
        );
    }

    /// Blank read ids are rejected and the 200-byte boundary is enforced.
    #[test]
    fn blank_and_overlong_read_ids_are_rejected_at_their_boundaries() {
        let (blank_root, blank_path) = write_seq_summ("blank_id", &format!("{HEADER}\t1234\n"));
        let blank_str = blank_path
            .to_str()
            .expect("temporary paths are valid UTF-8");
        let blank_error = run_with_seq_summ(blank_str).expect_err("a blank read id must fail");
        remove_temp_dir(&blank_root);
        let blank_rendered = format!("{blank_error:?}");
        assert!(
            matches!(blank_error, Error::InvalidReadID(message) if message == "read_id is blank"),
            "an empty read id must be InvalidReadID, got {blank_rendered}"
        );

        let max_id = "r".repeat(usize::from(MAX_READ_ID_LEN));
        let (max_root, max_path) = write_seq_summ("max_id", &format!("{HEADER}{max_id}\t1234\n"));
        let max_str = max_path.to_str().expect("temporary paths are valid UTF-8");
        let output = run_with_seq_summ(max_str).expect("a 200-byte read id must be accepted");
        remove_temp_dir(&max_root);
        assert!(
            !output.contains(&max_id),
            "the unmatched maximum-length read id must not appear in the output"
        );

        let long_id = "r".repeat(usize::from(MAX_READ_ID_LEN) + 1);
        let (long_root, long_path) =
            write_seq_summ("long_id", &format!("{HEADER}{long_id}\t1234\n"));
        let long_str = long_path.to_str().expect("temporary paths are valid UTF-8");
        let long_error = run_with_seq_summ(long_str).expect_err("a 201-byte read id must fail");
        remove_temp_dir(&long_root);
        let long_rendered = format!("{long_error:?}");
        let expected = format!("error in setting read_id, length > {MAX_READ_ID_LEN}");
        assert!(
            matches!(long_error, Error::InvalidReadID(message) if message == expected),
            "a 201-byte read id must fail the length check, got {long_rendered}"
        );
    }

    /// A line longer than the 1000-byte cap is rejected; a line at it is not.
    #[test]
    fn overlong_lines_report_the_cap_and_boundary_lines_are_accepted() {
        let line_cap = usize::from(MAX_SEQ_SUMM_SIZE_PER_LINE);
        let padding = "x".repeat(line_cap.saturating_sub(1).saturating_sub(ROW_PREFIX.len()));
        let at_cap = format!("{HEADER}{ROW_PREFIX}{padding}\n");
        let boundary_row = at_cap
            .lines()
            .nth(1)
            .expect("the boundary file must have a data row");
        assert_eq!(
            boundary_row.len(),
            line_cap.saturating_sub(1),
            "the boundary row plus its newline must consume exactly the cap"
        );

        let (at_root, at_path) = write_seq_summ("at_cap", &at_cap);
        let at_str = at_path.to_str().expect("temporary paths are valid UTF-8");
        let output = run_with_seq_summ(at_str).expect("a line at the cap must be accepted");
        remove_temp_dir(&at_root);
        assert!(
            output.contains("read_id\talign_length\tsequence_length_template\talignment_type"),
            "an accepted boundary row must still produce the table header, got {output}"
        );

        let over_cap = format!("{HEADER}{ROW_PREFIX}{}\n", "x".repeat(line_cap));
        let (over_root, over_path) = write_seq_summ("over_cap", &over_cap);
        let over_str = over_path.to_str().expect("temporary paths are valid UTF-8");
        let error = run_with_seq_summ(over_str).expect_err("an overlong line must fail");
        remove_temp_dir(&over_root);

        let rendered = format!("{error:?}");
        let expected = format!("line is too long (>{line_cap} bytes)");
        assert!(
            matches!(error, Error::InvalidState(message) if message == expected),
            "an overlong line must report the cap, got {rendered}"
        );
    }

    /// A `\r` that is not followed by `\n` is rejected, including at EOF.
    #[test]
    fn carriage_returns_not_followed_by_newline_are_rejected() {
        let cases = [
            (
                "lone_cr",
                format!("{HEADER}read-1\t12\r34\n"),
                "a `\\r` inside a line",
            ),
            (
                "eof_cr",
                format!("{HEADER}read-1\t1234\r"),
                "a `\\r` at end of file",
            ),
        ];
        for (label, contents, description) in cases {
            let (root, path) = write_seq_summ(label, &contents);
            let path_str = path.to_str().expect("temporary paths are valid UTF-8");
            let error = run_with_seq_summ(path_str).expect_err("a dangling `\\r` must fail");
            remove_temp_dir(&root);
            let rendered = format!("{error:?}");
            assert!(
                matches!(error, Error::InvalidState(message) if message == "\\r must be followed by \\n"),
                "{description} must be rejected, got {rendered}"
            );
        }
    }
}
