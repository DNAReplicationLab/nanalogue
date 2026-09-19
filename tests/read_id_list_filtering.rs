#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Integration tests for the read-ID-list loading path reached through the
//! public [`commands::run`] and [`BamRcRecords::new`] APIs.
//!
//! Every list file is written into a unique temporary directory, and the
//! assertions compare exact read-ID sets or per-record counts instead of
//! checking only that a command succeeded.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use clap::Parser as _;
    use nanalogue_core::constants::shared::MAX_READ_ID_LEN;
    use nanalogue_core::{
        BamRcRecords, Error, InputBamBuilder, InputMods, OptionalTag, PathOrURLOrStdin, commands,
        nanalogue_bam_reader, reads_table::sort_output_lines, uuid,
    };
    use std::collections::HashSet;
    use std::fs;
    use std::io;
    use std::path::{Path, PathBuf};

    /// Fixture with four records: three distinct read IDs, one of which
    /// ([`DUPLICATED_READ_ID`]) occurs on two alignments.
    const BAM_PATH: &str = "examples/example_1.bam";

    /// Read ID of the first alignment in [`BAM_PATH`]; occurs once.
    const SINGLE_READ_ID: &str = "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575";

    /// Read ID of the reverse-strand alignment in [`BAM_PATH`]; occurs once.
    const OTHER_SINGLE_READ_ID: &str = "fffffff1-10d2-49cb-8ca3-e8d48979001b";

    /// Read ID shared by a mapped and an unmapped alignment in [`BAM_PATH`].
    const DUPLICATED_READ_ID: &str = "a4f36092-b4d5-47a9-813e-c22c3b477a0c";

    /// Writes `contents` to `read_ids.txt` inside a fresh temporary directory.
    fn write_read_id_list(label: &str, contents: &[u8]) -> (PathBuf, PathBuf) {
        let root = std::env::temp_dir().join(format!(
            "nanalogue_read_id_list_{label}_{}",
            uuid::v4_random()
        ));
        fs::create_dir_all(&root).expect("temporary directory must be creatable");
        let list = root.join("read_ids.txt");
        fs::write(&list, contents).expect("read id list must be writable");
        (root, list)
    }

    /// Deletes the temporary directory that holds a list created above.
    fn remove_temp_dir(root: &Path) {
        fs::remove_dir_all(root).expect("temporary directory must be removable");
    }

    /// Runs `subcommand` through the public CLI with the given read-id list.
    fn run_cli(subcommand: &str, list: &Path) -> Result<String, Error> {
        let cli = commands::Cli::parse_from([
            "",
            subcommand,
            BAM_PATH,
            "--read-id-list",
            list.to_str().expect("temporary paths are valid UTF-8"),
        ]);
        let mut output = Vec::new();
        commands::run(cli, &mut output)?;
        Ok(String::from_utf8(output).expect("subcommand output is UTF-8"))
    }

    /// Runs `read-table-hide-mods` and returns the expected failure.
    fn hide_mods_failure(list: &Path) -> Error {
        run_cli("read-table-hide-mods", list).expect_err("the run must fail")
    }

    /// Read IDs in `read-table` output, sorted for order-independent checks.
    fn table_read_ids(output: &str) -> Vec<String> {
        sort_output_lines(output)
            .into_iter()
            .filter(|line| !line.starts_with('#') && !line.starts_with("read_id"))
            .map(|line| line.split('\t').next().unwrap_or_default().to_owned())
            .collect()
    }

    /// Read IDs reported by `read-info`, one entry per BAM record.
    fn info_read_ids(output: &str) -> Vec<String> {
        let entries: Vec<serde_json::Value> =
            serde_json::from_str(output).expect("read-info emits a JSON array");
        entries
            .iter()
            .map(|entry| {
                entry
                    .get("read_id")
                    .and_then(serde_json::Value::as_str)
                    .expect("every read-info entry carries a read id")
                    .to_owned()
            })
            .collect()
    }

    /// Two LF-separated IDs select exactly those reads from the fixture.
    #[test]
    fn lf_delimited_list_selects_exactly_the_requested_reads() {
        let (root, list) = write_read_id_list(
            "lf",
            format!("{SINGLE_READ_ID}\n{OTHER_SINGLE_READ_ID}\n").as_bytes(),
        );

        let output =
            run_cli("read-table-hide-mods", &list).expect("an LF-delimited list must load");
        let ids = table_read_ids(&output);

        assert_eq!(
            ids,
            vec![SINGLE_READ_ID.to_owned(), OTHER_SINGLE_READ_ID.to_owned()],
            "the read table must contain exactly the two listed read ids"
        );
        assert!(
            !output.contains(DUPLICATED_READ_ID),
            "the unlisted read id must be filtered out"
        );
        remove_temp_dir(&root);
    }

    /// CRLF line endings select the same reads as LF endings.
    #[test]
    fn crlf_delimited_list_matches_lf_output() {
        let (lf_root, lf_list) = write_read_id_list(
            "lf_compare",
            format!("{SINGLE_READ_ID}\n{OTHER_SINGLE_READ_ID}\n").as_bytes(),
        );
        let (crlf_root, crlf_list) = write_read_id_list(
            "crlf_compare",
            format!("{SINGLE_READ_ID}\r\n{OTHER_SINGLE_READ_ID}\r\n").as_bytes(),
        );

        let lf_ids =
            table_read_ids(&run_cli("read-table-hide-mods", &lf_list).expect("LF list loads"));
        let crlf_ids =
            table_read_ids(&run_cli("read-table-hide-mods", &crlf_list).expect("CRLF list loads"));

        assert_eq!(
            crlf_ids, lf_ids,
            "line endings must not change the selected reads"
        );
        assert_eq!(
            crlf_ids,
            vec![SINGLE_READ_ID.to_owned(), OTHER_SINGLE_READ_ID.to_owned()],
            "CRLF must not leave a trailing `\\r` on any read id"
        );
        remove_temp_dir(&lf_root);
        remove_temp_dir(&crlf_root);
    }

    /// A repeated ID in the list must not duplicate the filtered records.
    #[test]
    fn duplicate_list_entries_do_not_duplicate_records() {
        let (root, list) = write_read_id_list(
            "duplicates",
            format!("{DUPLICATED_READ_ID}\n{SINGLE_READ_ID}\n{DUPLICATED_READ_ID}\n").as_bytes(),
        );
        let (unique_root, unique_list) = write_read_id_list(
            "deduplicated",
            format!("{DUPLICATED_READ_ID}\n{SINGLE_READ_ID}\n").as_bytes(),
        );

        let duplicated_ids =
            info_read_ids(&run_cli("read-info", &list).expect("duplicate ids must load"));
        let unique_ids =
            info_read_ids(&run_cli("read-info", &unique_list).expect("unique ids must load"));

        assert_eq!(
            duplicated_ids, unique_ids,
            "a repeated id must not change the output"
        );
        assert_eq!(
            duplicated_ids.len(),
            3,
            "two alignments of the shared id plus one single-alignment read"
        );
        assert_eq!(
            duplicated_ids
                .iter()
                .filter(|id| id.as_str() == DUPLICATED_READ_ID)
                .count(),
            2,
            "the shared read id must appear once per alignment"
        );
        assert_eq!(
            duplicated_ids
                .iter()
                .filter(|id| id.as_str() == SINGLE_READ_ID)
                .count(),
            1,
            "the single-alignment read must appear exactly once"
        );
        remove_temp_dir(&root);
        remove_temp_dir(&unique_root);
    }

    /// A blank line aborts both read-table subcommands with `InvalidState`.
    #[test]
    fn blank_line_is_rejected() {
        let (root, list) = write_read_id_list(
            "blank",
            format!("{SINGLE_READ_ID}\n\n{OTHER_SINGLE_READ_ID}\n").as_bytes(),
        );

        for (subcommand, result) in [
            (
                "read-table-hide-mods",
                run_cli("read-table-hide-mods", &list),
            ),
            (
                "read-table-show-mods",
                run_cli("read-table-show-mods", &list),
            ),
        ] {
            let error = match result {
                Err(error) => error,
                Ok(output) => unreachable!("{subcommand} unexpectedly produced {output}"),
            };
            let rendered = format!("{error:?}");
            assert!(
                matches!(error, Error::InvalidState(message) if message == "blank line found in read id file!"),
                "{subcommand} must report the documented InvalidState, got {rendered}"
            );
        }
        remove_temp_dir(&root);
    }

    /// The 200-byte maximum loads; a 201-byte ID fails the length check.
    #[test]
    fn overlong_read_ids_are_rejected() {
        let max_length_id = "r".repeat(usize::from(MAX_READ_ID_LEN));
        let (max_root, max_list) = write_read_id_list(
            "max_length",
            format!("{max_length_id}\n{SINGLE_READ_ID}\n").as_bytes(),
        );

        let output = run_cli("read-table-hide-mods", &max_list)
            .expect("an id of exactly 200 bytes must load");
        assert_eq!(
            table_read_ids(&output),
            vec![SINGLE_READ_ID.to_owned()],
            "the maximum-length id loads but matches no record in the fixture"
        );
        assert!(
            !output.contains(&max_length_id),
            "an unmatched read id must not appear in the output"
        );
        remove_temp_dir(&max_root);

        let too_long_id = "r".repeat(usize::from(MAX_READ_ID_LEN) + 1);
        let (long_root, long_list) =
            write_read_id_list("too_long", format!("{too_long_id}\n").as_bytes());
        let error = hide_mods_failure(&long_list);
        let rendered = format!("{error:?}");
        assert!(
            matches!(error, Error::InvalidReadID(message) if message == "error in setting read_id, length > 200"),
            "a 201-byte id must fail the length check, got {rendered}"
        );
        remove_temp_dir(&long_root);
    }

    /// A byte outside the read-id alphabet is rejected.
    #[test]
    fn invalid_alphabet_byte_is_rejected() {
        // A space passes the line reader but is forbidden by the identifier check.
        let (root, list) = write_read_id_list("space", b"read id\n");

        let error = hide_mods_failure(&list);
        let rendered = format!("{error:?}");
        assert!(
            matches!(error, Error::InvalidReadID(message) if message == "read_id contains forbidden characters"),
            "a space must fail the identifier alphabet check, got {rendered}"
        );
        remove_temp_dir(&root);
    }

    /// A missing list file surfaces as `InputOutputError` with `NotFound`.
    #[test]
    fn nonexistent_list_path_is_an_io_error() {
        let (root, list) = write_read_id_list("missing", format!("{SINGLE_READ_ID}\n").as_bytes());
        fs::remove_file(&list).expect("the list file must be removable");
        assert!(!list.exists(), "the list must be absent before the run");

        let error = hide_mods_failure(&list);
        let Error::InputOutputError(source) = error else {
            unreachable!("expected InputOutputError, got {error:?}")
        };

        assert_eq!(
            source.kind(),
            io::ErrorKind::NotFound,
            "the underlying io error must be preserved"
        );
        remove_temp_dir(&root);
    }

    /// Setting both `read_id_set` and `read_id_list` is rejected before loading.
    #[test]
    fn read_id_set_and_read_id_list_conflict_is_rejected() {
        let (root, list) = write_read_id_list("conflict", format!("{SINGLE_READ_ID}\n").as_bytes());
        let mut bam_opts = InputBamBuilder::default()
            .bam_path(PathOrURLOrStdin::Path(BAM_PATH.into()))
            .read_id_list(list.to_str().expect("temporary paths are valid UTF-8"))
            .build()
            .expect("a read-id list alone is a valid option set");

        let mut read_id_set = HashSet::new();
        let _: bool = read_id_set.insert(SINGLE_READ_ID.to_owned());
        bam_opts.read_id_set = Some(read_id_set);

        let mut reader = nanalogue_bam_reader(BAM_PATH).expect("the fixture must open");
        let error = BamRcRecords::new(
            &mut reader,
            &mut bam_opts,
            &mut InputMods::<OptionalTag>::default(),
        )
        .expect_err("setting both inputs must be rejected");
        let rendered = format!("{error:?}");

        assert!(
            matches!(error, Error::InvalidState(message) if message == "cannot set both `read_id_set` and `read_id_list` in `bam_opts` in `BamRcRecords`"),
            "the conflicting inputs must be InvalidState, got {rendered}"
        );
        remove_temp_dir(&root);
    }

    /// Leading comment lines are documented and must be ignored.
    #[test]
    fn leading_comments_are_ignored() {
        let (comment_root, comment_list) = write_read_id_list(
            "comments",
            format!(
                "# first comment\n#second comment\r\n{SINGLE_READ_ID}\n{OTHER_SINGLE_READ_ID}\r\n"
            )
            .as_bytes(),
        );
        let (plain_root, plain_list) = write_read_id_list(
            "no_comments",
            format!("{SINGLE_READ_ID}\n{OTHER_SINGLE_READ_ID}\n").as_bytes(),
        );

        let comment_ids = table_read_ids(
            &run_cli("read-table-hide-mods", &comment_list).expect("comments must be ignored"),
        );
        let plain_ids =
            table_read_ids(&run_cli("read-table-hide-mods", &plain_list).expect("plain list"));

        assert_eq!(
            comment_ids, plain_ids,
            "comment lines must not change the selected reads"
        );
        assert_eq!(
            comment_ids,
            vec![SINGLE_READ_ID.to_owned(), OTHER_SINGLE_READ_ID.to_owned()],
            "the two read ids after the comments must still filter"
        );
        remove_temp_dir(&comment_root);
        remove_temp_dir(&plain_root);
    }
}
