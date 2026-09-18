#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Simulation configuration validation through raw JSON, bypassing the
//! builder so runtime validation at the public `run` boundary is exercised.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::simulate_mod_bam::{ModConfig, PerfectSeqMatchToNot, run};
    use nanalogue_core::{Error, F32Bw0and1, ReadState, SimulationConfig, uuid};
    use rand::{SeedableRng as _, rngs::StdRng};
    use rust_htslib::bam;
    use rust_htslib::bam::Read as _;
    use std::path::{Path, PathBuf};

    /// Valid one-contig, one-read-group fixture with one modification config.
    /// Invalid cases replace exactly one of `win`, `mod_range`, or `drop`.
    const BASE_MOD_CONFIG: &str = r#"{
        "contigs": { "number": 1, "len_range": [40, 40] },
        "reads": [{ "number": 2, "mapq_range": [10, 20], "base_qual_range": [10, 20],
            "len_range": [0.5, 0.5],
            "mods": [{ "base": "C", "is_strand_plus": true, "mod_code": "m",
                "win": [3], "mod_range": [[0.5, 0.5]], "drop": [] }] }]
    }"#;

    /// Valid configuration without modifications, used for output-path checks.
    const PLAIN_CONFIG: &str = r#"{
        "contigs": { "number": 1, "len_range": [40, 40] },
        "reads": [{ "number": 2, "mapq_range": [10, 20], "base_qual_range": [10, 20],
            "len_range": [0.5, 0.5] }]
    }"#;

    /// Returns [`BASE_MOD_CONFIG`] with `from` replaced by `to`.
    fn with_field(from: &str, to: &str) -> String {
        BASE_MOD_CONFIG.replace(from, to)
    }

    /// Removes its unique directory when dropped so failed assertions do not
    /// leave simulated outputs behind.
    struct TempDir {
        /// Unique directory created for one test.
        path: PathBuf,
    }

    impl TempDir {
        /// Creates a fresh uniquely named directory under the system temp dir.
        fn new(tag: &str) -> Self {
            let path =
                std::env::temp_dir().join(format!("nanalogue_sim_{tag}_{}", uuid::v4_random()));
            std::fs::create_dir_all(&path).expect("test temp directory must be creatable");
            Self { path }
        }
    }

    impl Drop for TempDir {
        fn drop(&mut self) {
            drop(std::fs::remove_dir_all(&self.path));
        }
    }

    /// Deserializes a fixture that is structurally valid JSON, asserting that
    /// deserialization succeeds so a serde failure cannot be mistaken for the
    /// runtime validation error under test. Returns the config and its first
    /// modification config.
    fn deserialize_mod_config(json: &str) -> (SimulationConfig, ModConfig) {
        let config: SimulationConfig = serde_json::from_str(json)
            .expect("fixture must deserialize; this case targets runtime validation, not serde");
        let mod_config = config
            .reads
            .first()
            .and_then(|read| read.mods.first())
            .cloned()
            .expect("fixture must contain one modification config");
        (config, mod_config)
    }

    /// Asserts that `error` is an `InvalidState` carrying exactly `expected`.
    fn assert_invalid_state(error: &Error, expected: &str) {
        assert!(
            matches!(error, Error::InvalidState(message) if message == expected),
            "expected InvalidState({expected:?}), got {error:?}"
        );
    }

    /// Runs `config` inside `dir`, expecting rejection, and asserts that no
    /// alignment or FASTA output was created before the error.
    fn run_expecting_rejection(config: SimulationConfig, dir: &Path) -> Error {
        let bam_path = dir.join("simulation.bam");
        let fasta_path = dir.join("reference.fa");
        let error = run(config, &bam_path, &fasta_path)
            .expect_err("invalid simulation configuration must be rejected");
        assert!(
            !bam_path.exists() && !fasta_path.exists(),
            "validation must fail before writing any output"
        );
        error
    }

    /// An empty `win` schedule deserializes, but `run` rejects it at the
    /// `SimulationConfig::validate` boundary. With both schedules empty, that
    /// runtime win check still runs before read generation, where an empty
    /// `mod_range` would otherwise be checked.
    #[test]
    fn json_empty_win_fails_at_runtime_validation() {
        let dir = TempDir::new("empty_win");
        let json = with_field(r#""win": [3]"#, r#""win": []"#);
        let (config, mod_config) = deserialize_mod_config(&json);
        assert!(
            mod_config.win.is_empty(),
            "the empty win must survive deserialization to reach runtime validation"
        );
        let error = run_expecting_rejection(config, &dir.path);
        assert_invalid_state(&error, "win must contain at least one window value");

        let both_empty = json.replace(r#""mod_range": [[0.5, 0.5]]"#, r#""mod_range": []"#);
        let (precedence_config, precedence_mod) = deserialize_mod_config(&both_empty);
        assert!(
            precedence_mod.win.is_empty() && precedence_mod.mod_range.is_empty(),
            "both schedules must be empty after deserialization"
        );
        let precedence_error = run_expecting_rejection(precedence_config, &dir.path);
        assert_invalid_state(
            &precedence_error,
            "win must contain at least one window value",
        );
    }

    /// `drop: [5]` with `win: [3]` deserializes, but every drop value must be
    /// at most the smallest window, so `run` rejects it at the runtime
    /// validation boundary with both values in the message.
    #[test]
    fn json_drop_exceeding_min_win_is_rejected_at_runtime_validation() {
        let dir = TempDir::new("drop_too_large");
        let json = with_field(r#""drop": []"#, r#""drop": [5]"#);
        let (config, mod_config) = deserialize_mod_config(&json);
        assert_eq!(
            mod_config.drop,
            vec![5],
            "the drop must survive deserialization"
        );
        assert_eq!(
            mod_config.win,
            vec![std::num::NonZeroU32::new(3).expect("3 is nonzero")],
            "the fixture window must be smaller than the drop value"
        );

        let error = run_expecting_rejection(config, &dir.path);
        assert_invalid_state(
            &error,
            "drop value 5 exceeds minimum window size 3; \
             no drop value may exceed the smallest window",
        );
    }

    /// An empty `mod_range` passes `SimulationConfig::validate` but is rejected
    /// while reads are generated, before any output is written.
    #[test]
    fn json_empty_mod_range_is_rejected_during_read_generation() {
        let dir = TempDir::new("empty_mod_range");
        let json = with_field(r#""mod_range": [[0.5, 0.5]]"#, r#""mod_range": []"#);
        let (config, mod_config) = deserialize_mod_config(&json);
        assert!(
            mod_config.mod_range.is_empty() && !mod_config.win.is_empty(),
            "only the modification range may be empty in this fixture"
        );

        let error = run_expecting_rejection(config, &dir.path);
        assert_invalid_state(
            &error,
            "modification config 0 has an empty mod_range schedule",
        );
    }

    /// `number` fields are `NonZeroU32`, so JSON asking for zero contigs or
    /// zero reads fails during deserialization rather than at runtime.
    #[test]
    fn json_zero_number_fails_deserialization_as_nonzero_u32() {
        let cases = [
            (
                r#"{"contigs":{"number":0,"len_range":[20,20]},"reads":[]}"#,
                "contig count",
            ),
            (
                r#"{"contigs":{"number":1,"len_range":[20,20]},"reads":[{"number":0}]}"#,
                "read count",
            ),
        ];
        for (json, context) in cases {
            let error = serde_json::from_str::<SimulationConfig>(json)
                .expect_err("NonZeroU32 fields must reject zero at deserialization");
            let text = error.to_string();
            assert!(
                text.contains("invalid value: integer `0`") && text.contains("a nonzero u32"),
                "{context} must report a rejected zero for NonZeroU32, got: {text}"
            );
        }
    }

    /// A deserialized configuration always has at least one contig because
    /// `contigs.number` is a `NonZeroU32`, so the generator's
    /// `UnavailableData("no contigs found")` branch (which needs an empty
    /// contig slice) is unreachable from JSON. Exposing it would require
    /// calling the crate-private generator with an empty slice or weakening
    /// `NonZeroU32`, so this asserts the deserialization guard and the fact
    /// that an empty `reads` array still runs instead of fabricating a change.
    #[test]
    fn json_cannot_reach_the_no_contigs_runtime_state() {
        let error = serde_json::from_str::<SimulationConfig>(
            r#"{"contigs":{"number":0,"len_range":[20,20]},"reads":[]}"#,
        )
        .expect_err("zero contigs cannot be expressed through JSON");
        assert!(
            error.to_string().contains("a nonzero u32"),
            "zero contigs must fail at the NonZeroU32 boundary, got: {error}"
        );

        // Empty read groups are valid and produce an empty alignment; they are
        // not the "no contigs found" path.
        let dir = TempDir::new("no_contigs");
        let empty_reads: SimulationConfig =
            serde_json::from_str(r#"{"contigs":{"number":1,"len_range":[20,20]},"reads":[]}"#)
                .expect("empty reads must deserialize");
        let bam_path = dir.path.join("simulation.bam");
        run(empty_reads, &bam_path, &dir.path.join("reference.fa"))
            .expect("empty read groups must still write a valid alignment");
        assert!(
            bam_path.exists(),
            "an empty-reads run must still create the alignment output"
        );
    }

    /// A configuration whose `drop` equals `win` drops every base of interest:
    /// ten bases with a window of two drop two per window, no window modifies a
    /// base, and the final window consumes the last base. The same JSON runs
    /// with and without a seed, exercising both `run` RNG paths. With the
    /// implicit no-suffix style, the all-dropped group is emitted as an empty
    /// MM group and produces no ML values.
    #[test]
    fn json_drop_equal_to_win_drops_every_base_in_both_rng_paths() {
        for seed_field in ["", r#", "seed": 7"#] {
            let dir = TempDir::new("all_dropped");
            let json = format!(
                r#"{{
                    "contigs": {{ "number": 1, "len_range": [10, 10],
                        "repeated_seq": "ACGTACGTAC" }},
                    "reads": [{{ "number": 1, "len_range": [1.0, 1.0],
                        "mods": [{{ "base": "N", "is_strand_plus": false,
                            "mod_code": "m", "mm_suffix": "none",
                            "win": [2], "drop": [2],
                            "mod_range": [[0.5, 0.5]] }}] }}]{seed_field}
                }}"#
            );
            let config: SimulationConfig =
                serde_json::from_str(&json).expect("all-dropped config must deserialize");
            let bam_path = dir.path.join("simulation.bam");
            run(config, &bam_path, &dir.path.join("reference.fa"))
                .expect("all-dropped config must run");
            let mut reader = bam::Reader::from_path(&bam_path).expect("output must be valid BAM");
            let record = reader
                .records()
                .next()
                .expect("one read was requested")
                .expect("the read must decode");
            assert!(
                matches!(record.aux(b"MM"), Ok(bam::record::Aux::String("N-m;"))),
                "an all-dropped group with no suffix must emit an empty MM group"
            );
            drop(
                record
                    .aux(b"ML")
                    .expect_err("all-dropped bases must emit no ML values"),
            );
        }
    }

    /// Public simulator entry points that the JSON-only cases above do not
    /// reach: `run` with `String` paths (both RNG paths), the temporary
    /// simulation accessors, the builder-only drop validation, and the
    /// insertion setter.
    #[test]
    fn public_simulator_api_smoke_covers_string_paths_and_accessors() {
        use nanalogue_core::DNARestrictive;
        use nanalogue_core::simulate_mod_bam::{
            AlignmentFormat, ModConfigBuilder, TempBamSimulation,
        };

        let dir = TempDir::new("api_smoke");
        for seed_field in ["", r#", "seed": 7"#] {
            let json = [PLAIN_CONFIG.trim_end_matches('}'), seed_field, "}"].concat();
            let config: SimulationConfig =
                serde_json::from_str(&json).expect("fixture with optional seed must deserialize");
            let bam_path = dir
                .path
                .join("simulation.bam")
                .to_string_lossy()
                .to_string();
            let fasta_path = dir.path.join("reference.fa").to_string_lossy().to_string();
            run(config, &bam_path, &fasta_path).expect("String paths must be accepted");
        }

        let config: SimulationConfig =
            serde_json::from_str(PLAIN_CONFIG).expect("plain fixture must deserialize");
        let simulation = TempBamSimulation::new(config, AlignmentFormat::Bam)
            .expect("temporary simulation must build");
        assert!(
            matches!(simulation.format(), AlignmentFormat::Bam),
            "the remembered format must be BAM"
        );
        assert!(
            simulation.fasta_path().ends_with("simulation.fa"),
            "the temporary simulation must remember its FASTA path"
        );

        let mod_config = ModConfigBuilder::default()
            .base('C')
            .mod_code("m".into())
            .win(vec![2])
            .drop(vec![2])
            .mod_range(vec![(0.5, 0.5)])
            .build()
            .expect("drop equal to the smallest window must build");
        assert_eq!(mod_config.drop, vec![2], "the builder must keep the drop");

        let insert = "TT".parse::<DNARestrictive>().expect("valid DNA");
        let mut thread_rng = rand::rng();
        let (sequence, _cigar) = PerfectSeqMatchToNot::seq(b"ACGTACGT".to_vec())
            .expect("non-empty sequence")
            .insert_middle(insert)
            .build(ReadState::PrimaryFwd, &mut thread_rng)
            .expect("middle insertion must build");
        assert_eq!(
            sequence, b"ACGTTTACGT",
            "the inserted bases must appear at the middle"
        );

        // The empty-sequence guard is the only safeguard that runs before any
        // random state is consumed.
        assert!(
            matches!(
                PerfectSeqMatchToNot::seq(Vec::new()),
                Err(Error::InvalidState(_))
            ),
            "an empty sequence must be rejected"
        );

        // A FASTA output path that is an existing directory cannot be created;
        // the failure must propagate instead of being silently ignored.
        let blocked_fasta = dir.path.join("blocked.fa");
        std::fs::create_dir_all(&blocked_fasta).expect("test directory must be creatable");
        let blocked_config: SimulationConfig =
            serde_json::from_str(PLAIN_CONFIG).expect("plain fixture must deserialize");
        let blocked_error = run(
            blocked_config,
            &dir.path.join("blocked.bam"),
            &blocked_fasta,
        )
        .expect_err("an unusable FASTA path must fail the run");
        assert!(
            !matches!(blocked_error, Error::InvalidState(_)),
            "the failure must be an output-write error, got {blocked_error:?}"
        );

        // The CRAM branch writes its FASTA before indexing, so an unusable
        // FASTA path must propagate the write error on that branch as well.
        let cram_config: SimulationConfig =
            serde_json::from_str(PLAIN_CONFIG).expect("plain fixture must deserialize");
        let cram_error = run(cram_config, &dir.path.join("blocked.cram"), &blocked_fasta)
            .expect_err("an unusable FASTA path must fail the CRAM run");
        assert!(
            !matches!(cram_error, Error::InvalidState(_)),
            "the CRAM failure must be an output-write error, got {cram_error:?}"
        );

        // If the FASTA index path is unusable, the CRAM branch fails while
        // building the index after the FASTA itself was written.
        std::fs::create_dir_all(dir.path.join("blocked_fai.fa.fai"))
            .expect("test directory must be creatable");
        let index_config: SimulationConfig =
            serde_json::from_str(PLAIN_CONFIG).expect("plain fixture must deserialize");
        let index_error = run(
            index_config,
            &dir.path.join("indexed.cram"),
            &dir.path.join("blocked_fai.fa"),
        )
        .expect_err("an unusable FASTA index path must fail the CRAM run");
        assert!(
            !matches!(index_error, Error::InvalidState(_)),
            "the FASTA index failure must be an output-write error, got {index_error:?}"
        );
    }

    /// A `.sam` alignment path is rejected by extension, while an uppercase
    /// `.BAM` path is accepted case-insensitively and produces a real BAM.
    #[test]
    fn json_config_rejects_sam_and_accepts_uppercase_bam_output() {
        let dir = TempDir::new("output_extension");
        let sam_config: SimulationConfig =
            serde_json::from_str(PLAIN_CONFIG).expect("plain fixture must deserialize");
        let sam_path = dir.path.join("simulation.sam");
        let error = run(sam_config, &sam_path, &dir.path.join("reference.fa"))
            .expect_err(".sam output must be rejected");
        assert_invalid_state(&error, "alignment output path must end in .bam or .cram");
        assert!(
            !sam_path.exists(),
            "the rejected .sam path must not be created"
        );

        let bam_config: SimulationConfig =
            serde_json::from_str(PLAIN_CONFIG).expect("plain fixture must deserialize");
        let uppercase_bam = dir.path.join("simulation.BAM");
        run(bam_config, &uppercase_bam, &dir.path.join("reference.fa"))
            .expect("uppercase .BAM must be accepted case-insensitively");
        assert!(
            uppercase_bam.exists() && dir.path.join("simulation.BAM.bai").exists(),
            "the uppercase .BAM output and its index sidecar must be written"
        );
        let mut reader = bam::Reader::from_path(&uppercase_bam).expect("output must be valid BAM");
        assert_eq!(
            reader.header().target_count(),
            1,
            "one contig was requested"
        );
        assert_eq!(reader.records().count(), 2, "two reads were requested");
    }

    /// A mismatch pass that covers every position still preserves an `N` base,
    /// because `N` has no alternative canonical base. The input places `N`
    /// off-centre and the fraction covers all six positions, so exactly the
    /// five non-`N` positions change.
    #[test]
    fn mismatch_preserves_an_n_base_covered_by_the_range() {
        let input = b"ACNGTA".to_vec();
        let mismatch = F32Bw0and1::try_from(1.0).expect("1.0 is a valid fraction");
        let mut rng = StdRng::seed_from_u64(7);
        let (mutated, mutated_cigar) = PerfectSeqMatchToNot::seq(input.clone())
            .expect("non-empty sequence")
            .mismatch(mismatch)
            .build(ReadState::PrimaryFwd, &mut rng)
            .expect("mismatches must not consume the whole read");

        assert_eq!(
            mutated_cigar.expect("mapped read has a CIGAR").to_string(),
            "6M",
            "mismatches keep match operations"
        );
        assert_eq!(
            mutated.get(2).copied(),
            Some(b'N'),
            "the N at index 2 must survive a pass that covers every position"
        );
        let changed: Vec<usize> = mutated
            .iter()
            .zip(input.iter())
            .enumerate()
            .filter_map(|(index, (actual, original))| (actual != original).then_some(index))
            .collect();
        assert_eq!(
            changed,
            vec![0, 1, 3, 4, 5],
            "only non-N positions may change"
        );
        for (index, &base) in mutated.iter().enumerate() {
            if index != 2 {
                assert!(
                    b"ACGT".contains(&base),
                    "position {index} must be mutated to a canonical base, got {base}"
                );
            }
        }

        // A single-N read is the smallest exact check: its only base cannot be
        // mutated, so the resulting sequence is unchanged byte for byte. This
        // uses the thread RNG to cover the production `rand::rng()` code path,
        // which the seeded case above does not.
        let mut thread_rng = rand::rng();
        let (single_n, single_n_cigar) = PerfectSeqMatchToNot::seq(b"N".to_vec())
            .expect("non-empty sequence")
            .mismatch(mismatch)
            .build(ReadState::PrimaryFwd, &mut thread_rng)
            .expect("a single-N read is valid");
        assert_eq!(single_n, b"N", "the only base is N and must be preserved");
        assert_eq!(
            single_n_cigar.expect("mapped read has a CIGAR").to_string(),
            "1M",
            "the preserved base is still a match"
        );
    }
}
