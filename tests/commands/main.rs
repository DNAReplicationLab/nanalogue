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

    #[test]
    fn run_commands_read_stats_example_3() {
        let cli = commands::Cli::parse_from(["", "read-stats", "./examples/example_3.bam"]);

        assert_command_output_matches_file(cli, "./examples/example_3_read_stats");
    }

    #[test]
    fn run_commands_read_info() {
        let cli = commands::Cli::parse_from(["", "read-info", "./examples/example_1.bam"]);

        assert_command_output_matches_file(cli, "./examples/example_1_read_info");
    }

    #[test]
    fn run_commands_window_grad_example_10_win_10_step_1() {
        let cli = commands::Cli::parse_from([
            "",
            "window-grad",
            "--win",
            "10",
            "--step",
            "1",
            "./examples/example_10.sam",
        ]);

        assert_command_output_matches_file(cli, "./examples/example_10_win_grad_win_10_step_1");
    }

    #[test]
    fn run_commands_window_grad_example_10_win_20_step_2() {
        let cli = commands::Cli::parse_from([
            "",
            "window-grad",
            "--win",
            "20",
            "--step",
            "2",
            "./examples/example_10.sam",
        ]);

        assert_command_output_matches_file(cli, "./examples/example_10_win_grad_win_20_step_2");
    }

    #[test]
    fn run_commands_window_grad_example_11_win_10_step_1() {
        let cli = commands::Cli::parse_from([
            "",
            "window-grad",
            "--win",
            "10",
            "--step",
            "1",
            "./examples/example_11.sam",
        ]);

        assert_command_output_matches_file(cli, "./examples/example_11_win_grad_win_10_step_1");
    }

    #[test]
    fn run_commands_window_grad_example_11_win_20_step_2() {
        let cli = commands::Cli::parse_from([
            "",
            "window-grad",
            "--win",
            "20",
            "--step",
            "2",
            "./examples/example_11.sam",
        ]);

        assert_command_output_matches_file(cli, "./examples/example_11_win_grad_win_20_step_2");
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

/// Integration tests for `find-modified-reads` subcommands using simulated BAM files.
/// Each test creates a BAM file with two read groups:
/// - Group 0 (1000 reads): modification pattern that does NOT match the filter
/// - Group 1 (1000 reads): modification pattern that DOES match the filter
///   Read names are formatted as `{group_number}.{uuid}`, so we verify that all output
///   starts with "1." and has at least ~1000 lines.
#[cfg(test)]
mod find_modified_reads_tests {
    use clap::Parser as _;
    use nanalogue_core::commands;
    use nanalogue_core::simulate_mod_bam::{SimulationConfig, TempBamSimulation};
    use nanalogue_core::{F32Bw0and1, OrdPair};

    /// Creates a probability interval [x1, x2] that will produce the desired modification density.
    ///
    /// Since probabilities >= 0.5 are considered "modified" and < 0.5 are "unmodified",
    /// this function returns an interval such that when probabilities are uniformly
    /// distributed within it, the expected fraction of modified bases equals the input density.
    ///
    /// The interval [d/2, 0.5 + d/2] has width 0.5 and the fraction >= 0.5 is exactly d.
    fn make_interval_with_desired_mod_density(desired_density: F32Bw0and1) -> OrdPair<F32Bw0and1> {
        let d = desired_density.val();
        (d / 2.0, 0.5 + d / 2.0).try_into().expect("valid interval")
    }

    /// Configuration for a single read group's modification pattern.
    #[derive(Clone, Copy)]
    enum GroupModConfig {
        /// Single window with uniform modification density.
        SingleWindow(F32Bw0and1),
        /// Two windows with different modification densities (for gradient/range tests).
        DualWindow(OrdPair<F32Bw0and1>),
    }

    /// Creates a BAM simulation with two read groups having the specified modification patterns.
    ///
    /// Group 0 (1000 reads) will NOT match the filter criteria.
    /// Group 1 (1000 reads) will match the filter criteria.
    fn create_two_group_simulation(
        group0: GroupModConfig,
        group1: GroupModConfig,
    ) -> TempBamSimulation {
        fn format_mod_spec(
            config: &GroupModConfig,
            make_interval: fn(F32Bw0and1) -> OrdPair<F32Bw0and1>,
        ) -> (String, String) {
            match *config {
                GroupModConfig::SingleWindow(density) => {
                    let interval = make_interval(density);
                    ("[100]".to_string(), format!("[[{interval}]]"))
                }
                GroupModConfig::DualWindow(densities) => {
                    let low = make_interval(densities.low());
                    let high = make_interval(densities.high());
                    ("[100, 100]".to_string(), format!("[[{low}], [{high}]]"))
                }
            }
        }

        let (win0, range0) = format_mod_spec(&group0, make_interval_with_desired_mod_density);
        let (win1, range1) = format_mod_spec(&group1, make_interval_with_desired_mod_density);

        let config_json = format!(
            r#"{{
                "contigs": {{ "number": 1, "len_range": [2000, 2000], "repeated_seq": "ACGTACGT" }},
                "reads": [
                    {{ "number": 1000, "len_range": [1.0, 1.0],
                       "mods": [{{ "base": "C", "mod_code": "m",
                                   "win": {win0}, "mod_range": {range0} }}] }},
                    {{ "number": 1000, "len_range": [1.0, 1.0],
                       "mods": [{{ "base": "C", "mod_code": "m",
                                   "win": {win1}, "mod_range": {range1} }}] }}
                ]
            }}"#
        );

        let config: SimulationConfig = serde_json::from_str(&config_json).unwrap();
        TempBamSimulation::new(config).unwrap()
    }

    /// Builds a `find-modified-reads` CLI command with the given parameters.
    fn build_find_modified_reads_cli(
        subcommand: &str,
        win_size: u8,
        step_size: u8,
        extras: &[&str],
        bam_path: &str,
    ) -> commands::Cli {
        let win = win_size.to_string();
        let step = step_size.to_string();

        let mut args: Vec<&str> = vec![
            "",
            "find-modified-reads",
            subcommand,
            "--win",
            &win,
            "--step",
            &step,
            "--tag",
            "m",
        ];
        args.extend(extras.iter().copied());
        args.push(bam_path);
        commands::Cli::parse_from(args)
    }

    /// Helper function that runs a find-modified-reads command and verifies:
    /// 1. All output lines start with "1." (from read group 1)
    /// 2. There are at least (1000 - tolerance%) lines (allowing for statistical variance)
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "converting f32 line count to usize is safe for test assertions"
    )]
    fn assert_find_mod_reads_selects_group_1(cli: commands::Cli, tolerance_fraction: F32Bw0and1) {
        let mut output = Vec::new();
        commands::run(cli, &mut output).expect("command should succeed");

        let output_str = String::from_utf8(output).expect("Invalid UTF-8");
        let lines: Vec<&str> = output_str.lines().collect();

        let min_expected = (1000.0 * (1.0 - tolerance_fraction.val())).floor() as usize;
        assert!(
            lines.len() >= min_expected,
            "Expected at least {min_expected} lines (all from group 1), got {}",
            lines.len()
        );

        for (i, line) in lines.iter().enumerate() {
            assert!(
                line.starts_with("1."),
                "Line {i} should start with '1.' but was: {line}"
            );
        }
    }

    /// Test `all-dens-between`: finds reads where ALL windowed densities are within [low, high].
    /// - Group 0: All densities at 0.9 (outside range 0.3-0.7)
    /// - Group 1: All densities at 0.5 (within range 0.3-0.7)
    #[test]
    fn find_modified_reads_all_dens_between() {
        let sim = create_two_group_simulation(
            GroupModConfig::SingleWindow(0.9.try_into().unwrap()),
            GroupModConfig::SingleWindow(0.5.try_into().unwrap()),
        );

        let cli = build_find_modified_reads_cli(
            "all-dens-between",
            100,
            100,
            &["--dens-limits", "0.3,0.7"],
            sim.bam_path(),
        );

        assert_find_mod_reads_selects_group_1(cli, 0.05.try_into().unwrap());
    }

    /// Test `any-dens-above`: finds reads where at least one window has density >= high.
    /// - Group 0: All densities at 0.3 (none >= 0.7)
    /// - Group 1: All densities at 0.9 (all >= 0.7)
    #[test]
    fn find_modified_reads_any_dens_above() {
        let sim = create_two_group_simulation(
            GroupModConfig::SingleWindow(0.3.try_into().unwrap()),
            GroupModConfig::SingleWindow(0.9.try_into().unwrap()),
        );

        let cli = build_find_modified_reads_cli(
            "any-dens-above",
            100,
            100,
            &["--high", "0.7"],
            sim.bam_path(),
        );

        assert_find_mod_reads_selects_group_1(cli, 0.05.try_into().unwrap());
    }

    /// Test `any-dens-below`: finds reads where at least one window has density <= low.
    /// - Group 0: All densities at 0.7 (none <= 0.3)
    /// - Group 1: All densities at 0.1 (all <= 0.3)
    #[test]
    fn find_modified_reads_any_dens_below() {
        let sim = create_two_group_simulation(
            GroupModConfig::SingleWindow(0.7.try_into().unwrap()),
            GroupModConfig::SingleWindow(0.1.try_into().unwrap()),
        );

        let cli = build_find_modified_reads_cli(
            "any-dens-below",
            100,
            100,
            &["--low", "0.3"],
            sim.bam_path(),
        );

        assert_find_mod_reads_selects_group_1(cli, 0.05.try_into().unwrap());
    }

    /// Test `any-dens-below-and-any-dens-above`: finds reads where at least one window
    /// is <= low AND at least one window is >= high.
    /// - Group 0: All densities at 0.5 (no extremes)
    /// - Group 1: Alternating pattern with densities at 0.1 and 0.9 (both extremes present)
    #[test]
    fn find_modified_reads_any_dens_below_and_any_dens_above() {
        let sim = create_two_group_simulation(
            GroupModConfig::SingleWindow(0.5.try_into().unwrap()),
            GroupModConfig::DualWindow((0.1, 0.9).try_into().unwrap()),
        );

        let cli = build_find_modified_reads_cli(
            "any-dens-below-and-any-dens-above",
            100,
            100,
            &["--low", "0.3", "--high", "0.7"],
            sim.bam_path(),
        );

        assert_find_mod_reads_selects_group_1(cli, 0.05.try_into().unwrap());
    }

    /// Test `dens-range-above`: finds reads where max(densities) - min(densities) >= threshold.
    /// - Group 0: All densities at 0.5 (range = 0)
    /// - Group 1: Alternating pattern with densities at 0.1 and 0.9 (range = 0.8)
    #[test]
    fn find_modified_reads_dens_range_above() {
        let sim = create_two_group_simulation(
            GroupModConfig::SingleWindow(0.5.try_into().unwrap()),
            GroupModConfig::DualWindow((0.1, 0.9).try_into().unwrap()),
        );

        let cli = build_find_modified_reads_cli(
            "dens-range-above",
            100,
            100,
            &["--min-range", "0.5"],
            sim.bam_path(),
        );

        assert_find_mod_reads_selects_group_1(cli, 0.05.try_into().unwrap());
    }

    /// Test `any-abs-grad-above`: finds reads where |gradient| >= threshold.
    /// - Group 0: Uniform densities at 0.5 (gradient = 0)
    /// - Group 1: Sharp transition from 0.1 to 0.9 (high gradient)
    #[test]
    fn find_modified_reads_any_abs_grad_above() {
        let sim = create_two_group_simulation(
            GroupModConfig::SingleWindow(0.5.try_into().unwrap()),
            GroupModConfig::DualWindow((0.1, 0.9).try_into().unwrap()),
        );

        let cli = build_find_modified_reads_cli(
            "any-abs-grad-above",
            100,
            10,
            &["--min-grad", "0.01"],
            sim.bam_path(),
        );

        assert_find_mod_reads_selects_group_1(cli, 0.05.try_into().unwrap());
    }
}
