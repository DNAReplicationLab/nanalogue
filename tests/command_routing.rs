#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Integration coverage for public command routing and sequence-display options.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use clap::Parser as _;
    use nanalogue_core::{
        Error, InputBamBuilder, PathOrURLOrStdin,
        commands::{self, Cli, CliBuilder, Commands},
    };
    use rust_htslib::errors::Error as HtslibError;
    use url::Url;

    /// Runs a parsed command and returns its UTF-8 output.
    fn run_cli<const N: usize>(args: [&str; N]) -> String {
        let cli = Cli::parse_from(args);
        let mut output = Vec::new();
        commands::run(cli, &mut output).expect("command should succeed");
        String::from_utf8(output).expect("command output should be UTF-8")
    }

    /// Checks command output byte-for-byte against a single-row checked-in table.
    fn assert_matches_example(actual: &str, expected_path: &str) {
        let expected = std::fs::read_to_string(expected_path)
            .expect("checked-in command output should be readable");
        assert_eq!(actual, expected, "output differs from {expected_path}");
    }

    /// A programmatically constructed unsupported URL reaches the ordinary URL
    /// arm and must preserve the shared URL-policy error without output.
    #[test]
    fn run_routes_unindexed_url_errors() {
        let url = Url::parse("file:///not/a/network/input.bam").expect("URL should parse");
        let bam = InputBamBuilder::default()
            .bam_path(PathOrURLOrStdin::URL(url))
            .build()
            .expect("input options should build");
        let cli = CliBuilder::default()
            .command(Commands::ReadStats { bam })
            .build()
            .expect("CLI should build");
        let mut output = Vec::new();

        let error =
            commands::run(cli, &mut output).expect_err("unsupported URL scheme should be rejected");

        assert!(output.is_empty(), "failed routing must not emit a table");
        assert!(
            matches!(error, Error::InvalidState(message)
                if message == "URL scheme `file` is not in the allow-list (http, https, ftp)"),
            "routing should return the URL reader's precise policy error"
        );
    }

    /// Adding a region reaches the indexed URL arm and its shared policy failure
    /// must likewise be returned before any command output is written.
    #[test]
    fn run_routes_indexed_url_errors() {
        let url = Url::parse("ssh://example.invalid/input.bam").expect("URL should parse");
        let bam = InputBamBuilder::default()
            .bam_path(PathOrURLOrStdin::URL(url))
            .region("chr1:1-10".to_owned())
            .build()
            .expect("region input options should build");
        let cli = CliBuilder::default()
            .command(Commands::ReadStats { bam })
            .build()
            .expect("CLI should build");
        let mut output = Vec::new();

        let error = commands::run(cli, &mut output)
            .expect_err("unsupported indexed URL scheme should be rejected");

        assert!(output.is_empty(), "failed routing must not emit a table");
        assert!(
            matches!(error, Error::InvalidState(message)
                if message == "URL scheme `ssh` is not in the allow-list (http, https, ftp)"),
            "indexed routing should return the URL reader's precise policy error"
        );
    }

    /// A fetch failure from an indexed BAM must be propagated rather than
    /// falling back as though only its index were absent.
    #[test]
    fn run_propagates_indexed_path_open_errors() {
        let cli = Cli::parse_from([
            "nanalogue",
            "read-stats",
            "--region",
            "absent:1-10",
            "./examples/example_1.bam",
        ]);
        let mut output = Vec::new();

        let error = commands::run(cli, &mut output)
            .expect_err("fetching an absent contig from an indexed BAM should fail");

        assert!(
            output.is_empty(),
            "an indexed fetch failure must not emit statistics"
        );
        assert!(
            matches!(&error, Error::RustHtslibError(source)
                if matches!(source.as_ref(), HtslibError::Fetch)),
            "the indexed reader's fetch error should be propagated, got {error:?}"
        );
    }

    /// A valid BAM without an index must fall back to sequential reading while
    /// retaining the requested region filter and exact statistics semantics.
    #[test]
    fn run_falls_back_for_unindexed_regional_path() {
        let actual = run_cli([
            "nanalogue",
            "read-stats",
            "--region",
            "dummyI:1-22",
            "./examples/example_1_copy_no_index.bam",
        ]);
        let expected = concat!(
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

        assert_eq!(actual, expected, "fallback must preserve region filtering");
    }

    /// Full-sequence routing forwards the quality flag for the show-mods
    /// command, producing both complete sequence and quality columns.
    #[test]
    fn show_mods_routes_full_sequence_with_qualities() {
        let actual = run_cli([
            "nanalogue",
            "read-table-show-mods",
            "--seq-full",
            "--show-base-qual",
            "./examples/example_5_valid_basequal.sam",
        ]);

        assert_matches_example(
            &actual,
            "./examples/example_5_valid_basequal_read_table_show_mods",
        );
        assert!(
            actual.contains("\tsequence\tqualities\n"),
            "both requested columns should be present"
        );
        assert!(
            actual.contains("\tTCGTTTCT\t32.0.69.80.79.81.29.30\n"),
            "the full sequence and corresponding qualities should be emitted"
        );
    }

    /// Region-sequence routing independently forwards the quality flag for the
    /// show-mods command and limits both values to the requested coordinates.
    #[test]
    fn show_mods_routes_region_sequence_with_qualities() {
        let actual = run_cli([
            "nanalogue",
            "read-table-show-mods",
            "--seq-region",
            "dummyI:10-12",
            "--show-base-qual",
            "./examples/example_5_valid_basequal.sam",
        ]);

        assert_matches_example(
            &actual,
            "./examples/example_5_valid_basequal_read_table_show_mods_subset",
        );
        assert!(
            actual.contains("\tCG\t0.69\n"),
            "sequence and qualities should cover only the two-base region"
        );
        assert!(
            !actual.contains("\tTCGTTTCT\t"),
            "region routing must not emit the full sequence"
        );
    }

    /// Region routing forwards insertion casing for the hide-mods command;
    /// inserted bases become lowercase while reference gaps remain dots.
    #[test]
    fn hide_mods_routes_region_with_lowercase_insertions() {
        let actual = run_cli([
            "nanalogue",
            "read-table-hide-mods",
            "--seq-region",
            "dummyI",
            "--show-ins-lowercase",
            "./examples/example_7.sam",
        ]);

        assert_matches_example(
            &actual,
            "./examples/example_7_table_hide_mods_ins_lowercase",
        );
        assert!(
            actual.contains("\tT..aTTTGT\n"),
            "the insertion should be lowercase and deletions should remain dots"
        );
        assert!(
            !actual.contains("\tmod_count"),
            "hide-mods routing must not add modification columns"
        );
    }
}
