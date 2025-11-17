//! # Nanalogue Simulate BAM
//!
//! Companion tool to nanalogue which creates artificial BAM or mod BAM files
//! for developers wishing to test BAM parsing or BAM modification data parsing.
use clap::Parser;
use nanalogue_core::{Error, SimulationConfig, simulate_mod_bam};

/// Main command line parsing struct that gets paths to files to be created.
#[derive(Parser, Debug)]
#[command(author, version,
    about = "Simulate BAM with or without modifications. Aimed at developers who wish to test their BAM parsers",
    long_about = None)]
struct Cli {
    /// Input JSON file path
    json: String,
    /// Output mod BAM file path; if pre-existing, the file will be overwritten.
    bam: String,
    /// Output fasta file path; if pre-existing, the file will be overwritten.
    fasta: String,
}

/// Main function, run the program. All business logic handled by `run`
///
/// This separation of function between `main` and `run` is so that we
/// can test the functionality of `run` easily through our code without
/// actually running the program on the command line like an external user.
fn main() {
    // Parse command line options
    let cli = Cli::parse();

    // call the run function and get the result
    match run(&cli) {
        Ok(()) => {}
        Err(e) => {
            eprintln!("Error during execution: {e}");
            std::process::exit(1);
        }
    }
}

/// Simple wrapper around `simulate_mod_bam`.
///
/// # Errors
/// Returns errors from simulating BAM files
fn run(cli: &Cli) -> Result<(), Error> {
    let json_str = std::fs::read_to_string(&cli.json)?;
    let config: SimulationConfig = serde_json::from_str(&json_str)?;
    simulate_mod_bam::run(config, &cli.bam, &cli.fasta)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::env;
    use std::fs;
    use uuid::Uuid;

    /// Test that the run function successfully creates BAM, BAI, and FASTA files from valid JSON config
    #[test]
    fn run_creates_output_files() {
        // Create temporary directory for test files
        let temp_path = env::temp_dir().join(format!("nanalogue_test_{}", Uuid::new_v4()));
        fs::create_dir_all(&temp_path).expect("failed to create temp dir");

        // Create JSON configuration file
        let json_config = r#"{
            "contigs": {
                "number": 2,
                "len_range": [100, 200],
                "repeated_seq": "ACGTACGT"
            },
            "reads": [
                {
                    "number": 10,
                    "mapq_range": [10, 30],
                    "base_qual_range": [20, 40],
                    "len_range": [0.1, 0.9],
                    "barcode": "ACGTAA",
                    "mods": [{
                        "base": "C",
                        "is_strand_plus": true,
                        "mod_code": "m",
                        "win": [5, 3],
                        "mod_range": [[0.3, 0.7], [0.1, 0.5]]
                    }]
                }
            ]
        }"#;

        let json_path = temp_path.join("config.json");
        fs::write(&json_path, json_config).expect("failed to write JSON config");

        // Create paths for output files
        let bam_path = temp_path.join("output.bam");
        let fasta_path = temp_path.join("output.fasta");

        // Create Cli struct with our test paths
        let cli = Cli {
            json: json_path.to_string_lossy().to_string(),
            bam: bam_path.to_string_lossy().to_string(),
            fasta: fasta_path.to_string_lossy().to_string(),
        };

        // Run the function
        run(&cli).expect("run should succeed with valid config");

        // Verify output files were created
        assert!(bam_path.exists(), "BAM file should be created");
        assert!(fasta_path.exists(), "FASTA file should be created");

        // Verify BAI index file was created
        let bai_path = temp_path.join("output.bam.bai");
        assert!(bai_path.exists(), "BAI index file should be created");

        // Verify files are not empty
        let bam_meta = fs::metadata(&bam_path).expect("failed to get BAM metadata");
        let fasta_meta = fs::metadata(&fasta_path).expect("failed to get FASTA metadata");
        let index_meta = fs::metadata(&bai_path).expect("failed to get BAI metadata");

        assert!(bam_meta.len() > 0, "BAM file should not be empty");
        assert!(fasta_meta.len() > 0, "FASTA file should not be empty");
        assert!(index_meta.len() > 0, "BAI index file should not be empty");

        // Clean up
        fs::remove_dir_all(&temp_path).expect("failed to clean up temp dir");
    }

    /// Test that the run function returns an error with invalid JSON
    #[test]
    fn run_fails_with_invalid_json() {
        let temp_path = env::temp_dir().join(format!("nanalogue_test_{}", Uuid::new_v4()));
        fs::create_dir_all(&temp_path).expect("failed to create temp dir");

        // Create invalid JSON file
        let json_path = temp_path.join("invalid.json");
        fs::write(&json_path, "{ invalid json }").expect("failed to write invalid JSON");

        let bam_path = temp_path.join("output.bam");
        let fasta_path = temp_path.join("output.fasta");

        let cli = Cli {
            json: json_path.to_string_lossy().to_string(),
            bam: bam_path.to_string_lossy().to_string(),
            fasta: fasta_path.to_string_lossy().to_string(),
        };

        // Run should fail with invalid JSON
        let result = run(&cli);
        assert!(result.is_err(), "run should fail with invalid JSON");

        // Clean up
        fs::remove_dir_all(&temp_path).expect("failed to clean up temp dir");
    }

    /// Test that the run function returns an error when JSON file doesn't exist
    #[test]
    fn run_fails_with_missing_json_file() {
        let temp_path = env::temp_dir().join(format!("nanalogue_test_{}", Uuid::new_v4()));
        fs::create_dir_all(&temp_path).expect("failed to create temp dir");

        let cli = Cli {
            json: temp_path
                .join("nonexistent.json")
                .to_string_lossy()
                .to_string(),
            bam: temp_path.join("output.bam").to_string_lossy().to_string(),
            fasta: temp_path.join("output.fasta").to_string_lossy().to_string(),
        };

        let result = run(&cli);
        assert!(
            result.is_err(),
            "run should fail when JSON file doesn't exist"
        );

        // Clean up
        fs::remove_dir_all(&temp_path).expect("failed to clean up temp dir");
    }
}
