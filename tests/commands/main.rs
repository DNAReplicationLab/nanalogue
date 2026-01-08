//! Integration tests for [`commands::run`]

/// Integration tests for [`commands::run`]
#[cfg(test)]
mod tests {

    use clap::Parser as _;
    use nanalogue_core::{commands, reads_table::sort_output_lines};

    /// Helper function that runs a CLI command and compares its sorted output with a file
    fn assert_command_output_matches_file(cli: commands::Cli, expected_output_path: &str) {
        // Capture output to buffer
        let mut output = Vec::new();

        // Run the command
        commands::run(cli, &mut output).expect("no error");

        // Compare with expected output
        let actual_output = String::from_utf8(output).expect("Invalid UTF-8");
        let expected_output = std::fs::read_to_string(expected_output_path)
            .expect("Failed to read expected output file");

        let actual_lines = sort_output_lines(&actual_output);
        let expected_lines = sort_output_lines(&expected_output);

        assert_eq!(actual_lines, expected_lines);
    }

    #[test]
    fn run_commands_read_table_show_mods_with_seq_summ() {
        let cli = commands::Cli::parse_from([
            "",
            "read-table-show-mods",
            "./examples/example_1.bam",
            "./examples/example_1_sequencing_summary",
        ]);

        assert_command_output_matches_file(
            cli,
            "./examples/example_1_read_table_show_mods_seq_summ",
        );
    }

    #[test]
    fn run_commands_read_table_hide_mods() {
        let cli =
            commands::Cli::parse_from(["", "read-table-hide-mods", "./examples/example_1.bam"]);

        assert_command_output_matches_file(cli, "./examples/example_1_read_table_hide_mods");
    }

    #[test]
    fn run_commands_read_info_detailed_print() {
        let cli =
            commands::Cli::parse_from(["", "read-info", "--detailed", "./examples/example_1.bam"]);

        assert_command_output_matches_file(cli, "./examples/example_1_detailed_print");
    }

    #[test]
    fn run_commands_read_info_detailed_pretty_print() {
        let cli = commands::Cli::parse_from([
            "",
            "read-info",
            "--detailed-pretty",
            "./examples/example_1.bam",
        ]);

        assert_command_output_matches_file(cli, "./examples/example_1_detailed_pretty_print");
    }

    #[test]
    fn run_commands_window_dens_win_2_step_1() {
        let cli = commands::Cli::parse_from([
            "",
            "window-dens",
            "--win",
            "2",
            "--step",
            "1",
            "./examples/example_1.bam",
        ]);

        assert_command_output_matches_file(cli, "./examples/example_1_window_reads");
    }

    #[test]
    fn run_commands_peek() {
        let cli = commands::Cli::parse_from(["", "peek", "./examples/example_1.bam"]);

        assert_command_output_matches_file(cli, "./examples/example_1_peek");
    }
}

/// Tests for BAM input option retrieval from `Commands` and `FindModReadsSubcommands`
#[cfg(test)]
mod bam_input_option_retrieval_from_commands {
    use nanalogue_core::{
        F32Bw0and1, InputBam, InputBamBuilder, InputMods, InputModsBuilder, InputWindowing,
        OptionalTag, OrdPair, PathOrURLOrStdin, RequiredTag,
        commands::{Commands, FindModReadsSubcommands},
    };
    use std::str::FromStr as _;

    /// Helper function to create a test `InputBam` with distinctive values
    fn create_test_input_bam() -> InputBam {
        InputBamBuilder::default()
            .bam_path(PathOrURLOrStdin::Path("/test/path/to/test_bam.bam".into()))
            .min_seq_len(12345u64)
            .build()
            .expect("Failed to build test InputBam")
    }

    /// Helper function to verify two `InputBam` instances match on key fields
    fn assert_bam_matches(actual: &InputBam, expected: &InputBam) {
        assert_eq!(
            actual.bam_path.to_string(),
            expected.bam_path.to_string(),
            "bam_path mismatch"
        );
        assert_eq!(
            actual.min_seq_len, expected.min_seq_len,
            "min_seq_len mismatch"
        );
    }

    #[test]
    fn find_mod_reads_all_dens_between_bam_retrieval() {
        let test_bam = create_test_input_bam();
        let command = FindModReadsSubcommands::AllDensBetween {
            bam: test_bam.clone(),
            win: InputWindowing::default(),
            mods: InputModsBuilder::<RequiredTag>::default()
                .tag(RequiredTag::from_str("m").unwrap())
                .build()
                .unwrap(),
            dens_limits: OrdPair::new(F32Bw0and1::new(0.2).unwrap(), F32Bw0and1::new(0.8).unwrap())
                .unwrap(),
        };

        let retrieved_bam = command.bam();
        assert_bam_matches(&retrieved_bam, &test_bam);
    }

    #[test]
    fn find_mod_reads_any_dens_above_bam_retrieval() {
        let test_bam = create_test_input_bam();
        let command = FindModReadsSubcommands::AnyDensAbove {
            bam: test_bam.clone(),
            win: InputWindowing::default(),
            mods: InputModsBuilder::<RequiredTag>::default()
                .tag(RequiredTag::from_str("m").unwrap())
                .build()
                .unwrap(),
            high: F32Bw0and1::new(0.7).unwrap(),
        };

        let retrieved_bam = command.bam();
        assert_bam_matches(&retrieved_bam, &test_bam);
    }

    #[test]
    fn find_mod_reads_any_dens_below_bam_retrieval() {
        let test_bam = create_test_input_bam();
        let command = FindModReadsSubcommands::AnyDensBelow {
            bam: test_bam.clone(),
            win: InputWindowing::default(),
            mods: InputModsBuilder::<RequiredTag>::default()
                .tag(RequiredTag::from_str("m").unwrap())
                .build()
                .unwrap(),
            low: F32Bw0and1::new(0.3).unwrap(),
        };

        let retrieved_bam = command.bam();
        assert_bam_matches(&retrieved_bam, &test_bam);
    }

    #[test]
    fn find_mod_reads_any_dens_below_and_any_dens_above_bam_retrieval() {
        let test_bam = create_test_input_bam();
        let command = FindModReadsSubcommands::AnyDensBelowAndAnyDensAbove {
            bam: test_bam.clone(),
            win: InputWindowing::default(),
            mods: InputModsBuilder::<RequiredTag>::default()
                .tag(RequiredTag::from_str("m").unwrap())
                .build()
                .unwrap(),
            low: F32Bw0and1::new(0.2).unwrap(),
            high: F32Bw0and1::new(0.8).unwrap(),
        };

        let retrieved_bam = command.bam();
        assert_bam_matches(&retrieved_bam, &test_bam);
    }

    #[test]
    fn find_mod_reads_dens_range_above_bam_retrieval() {
        let test_bam = create_test_input_bam();
        let command = FindModReadsSubcommands::DensRangeAbove {
            bam: test_bam.clone(),
            win: InputWindowing::default(),
            mods: InputModsBuilder::<RequiredTag>::default()
                .tag(RequiredTag::from_str("m").unwrap())
                .build()
                .unwrap(),
            min_range: F32Bw0and1::new(0.5).unwrap(),
        };

        let retrieved_bam = command.bam();
        assert_bam_matches(&retrieved_bam, &test_bam);
    }

    #[test]
    fn find_mod_reads_any_abs_grad_above_bam_retrieval() {
        let test_bam = create_test_input_bam();
        let command = FindModReadsSubcommands::AnyAbsGradAbove {
            bam: test_bam.clone(),
            win: InputWindowing::default(),
            mods: InputModsBuilder::<RequiredTag>::default()
                .tag(RequiredTag::from_str("m").unwrap())
                .build()
                .unwrap(),
            min_grad: F32Bw0and1::new(0.1).unwrap(),
        };

        let retrieved_bam = command.bam();
        assert_bam_matches(&retrieved_bam, &test_bam);
    }

    #[test]
    fn commands_read_table_show_mods_bam_retrieval() {
        let test_bam = create_test_input_bam();
        let command = Commands::ReadTableShowMods {
            bam: test_bam.clone(),
            mods: InputMods::<OptionalTag>::default(),
            seq_region: None,
            seq_full: false,
            show_base_qual: false,
            show_ins_lowercase: false,
            show_mod_z: false,
            seq_summ_file: String::new(),
        };

        let retrieved_bam = command.bam();
        assert_bam_matches(&retrieved_bam, &test_bam);
    }

    #[test]
    fn commands_read_table_hide_mods_bam_retrieval() {
        let test_bam = create_test_input_bam();
        let command = Commands::ReadTableHideMods {
            bam: test_bam.clone(),
            seq_region: None,
            seq_full: false,
            show_base_qual: false,
            show_ins_lowercase: false,
            seq_summ_file: String::new(),
        };

        let retrieved_bam = command.bam();
        assert_bam_matches(&retrieved_bam, &test_bam);
    }

    #[test]
    fn commands_read_stats_bam_retrieval() {
        let test_bam = create_test_input_bam();
        let command = Commands::ReadStats {
            bam: test_bam.clone(),
        };

        let retrieved_bam = command.bam();
        assert_bam_matches(&retrieved_bam, &test_bam);
    }

    #[test]
    fn commands_read_info_bam_retrieval() {
        let test_bam = create_test_input_bam();
        let command = Commands::ReadInfo {
            bam: test_bam.clone(),
            mods: InputMods::<OptionalTag>::default(),
            detailed: false,
            detailed_pretty: false,
        };

        let retrieved_bam = command.bam();
        assert_bam_matches(&retrieved_bam, &test_bam);
    }

    #[test]
    fn commands_find_modified_reads_bam_retrieval() {
        let test_bam = create_test_input_bam();
        let subcommand = FindModReadsSubcommands::DensRangeAbove {
            bam: test_bam.clone(),
            win: InputWindowing::default(),
            mods: InputModsBuilder::<RequiredTag>::default()
                .tag(RequiredTag::from_str("m").unwrap())
                .build()
                .unwrap(),
            min_range: F32Bw0and1::new(0.5).unwrap(),
        };

        let command = Commands::FindModifiedReads {
            command: subcommand,
        };

        let retrieved_bam = command.bam();
        assert_bam_matches(&retrieved_bam, &test_bam);
    }

    #[test]
    fn commands_window_dens_bam_retrieval() {
        let test_bam = create_test_input_bam();
        let command = Commands::WindowDens {
            bam: test_bam.clone(),
            win: InputWindowing::default(),
            mods: InputMods::<OptionalTag>::default(),
        };

        let retrieved_bam = command.bam();
        assert_bam_matches(&retrieved_bam, &test_bam);
    }

    #[test]
    fn commands_window_grad_bam_retrieval() {
        let test_bam = create_test_input_bam();
        let command = Commands::WindowGrad {
            bam: test_bam.clone(),
            win: InputWindowing::default(),
            mods: InputMods::<OptionalTag>::default(),
        };

        let retrieved_bam = command.bam();
        assert_bam_matches(&retrieved_bam, &test_bam);
    }
}
