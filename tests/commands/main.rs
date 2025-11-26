//! Integration tests for [`commands::run`]

/// Integration tests for [`commands::run`]
#[cfg(test)]
mod tests {

    use clap::Parser as _;
    use nanalogue_core::{commands, reads_table::sort_output_lines};

    #[test]
    fn run_commands_read_table_show_mods_with_seq_summ() {
        // Construct CLI
        let cli = commands::Cli::parse_from([
            "",
            "read-table-show-mods",
            "./examples/example_1.bam",
            "./examples/example_1_sequencing_summary",
        ]);

        // Capture output to buffer (same as existing test)
        let mut output = Vec::new();

        // Run the command
        commands::run(cli, &mut output).expect("no error");

        // Compare with expected output (same as existing test)
        let actual_output = String::from_utf8(output).expect("Invalid UTF-8");
        let expected_output =
            std::fs::read_to_string("./examples/example_1_read_table_show_mods_seq_summ")
                .expect("Failed to read expected output file");

        let actual_lines = sort_output_lines(&actual_output);
        let expected_lines = sort_output_lines(&expected_output);

        assert_eq!(actual_lines, expected_lines);
    }
}
