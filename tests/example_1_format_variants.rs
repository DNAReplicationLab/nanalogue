#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Golden alignment-format variants of `examples/example_1.bam`.
//!
//! `examples/example_1_formats/` holds the same four records re-encoded as
//! BAM, SAM, and CRAM with several indices, CRAM versions, reference modes,
//! codecs, slice layouts, and interleaved per-record encodings. The files are
//! committed goldens made by `generate.sh` in that directory; tests only read
//! them. Every variant must behave exactly like `examples/example_1.bam`,
//! apart from the documented `HTSlib` limitations pinned below.
//!
//! External-reference CRAMs name the reference through the relative header
//! tag `UR:examples/contigs.fa`. `HTSlib` tries a header M5 lookup through
//! `REF_PATH` (by default, the EBI server) before `UR`, so tests that decode
//! them run the executable with `REF_PATH` pointed away from the network and
//! `REF_CACHE` empty, so `HTSlib` neither downloads nor caches references. In-process tests only open variants that need no reference.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::{
        Error, nanalogue_bam_reader, nanalogue_indexed_bam_reader, reads_table::sort_output_lines,
        uuid,
    };
    use rust_htslib::bam::record::Aux;
    use rust_htslib::bam::{FetchDefinition, Read as _};
    use rust_htslib::errors::Error as HtslibError;
    use std::collections::BTreeSet;
    use std::fs;
    use std::path::{Path, PathBuf};
    use std::process::{Command, Output};

    /// Directory holding the golden variants, relative to the repository root.
    const VARIANT_DIR: &str = "examples/example_1_formats";

    /// The original file that every variant re-encodes.
    const ORIGINAL: &str = "examples/example_1.bam";

    /// The reference FASTA named by external-reference CRAM headers.
    const REFERENCE: &str = "examples/contigs.fa";

    /// A nonexistent `REF_PATH`, which disables `HTSlib`'s reference server
    /// lookups. `HTSlib` only reads from `REF_PATH`.
    const NO_REFERENCE_SEARCH: &str = "/nanalogue/example_1_formats/does/not/exist";

    /// Message printed to stderr when a region query falls back to a full scan.
    const NO_INDEX_WARNING: &str = "cannot find index file";

    /// Contig names and their header M5 checksums in `examples/contigs.fa`.
    const CONTIG_MD5S: [(&str, &str); 3] = [
        ("dummyI", "38365e17c1a365da4a246503ae46a622"),
        ("dummyII", "f76089038189b6a9c29ae23c6feb11bb"),
        ("dummyIII", "1d01faae0234a2ebdf640e80ecad866d"),
    ];

    /// First bytes of a BGZF file (BAM, bgzipped SAM, CSI).
    const BGZF_MAGIC: &[u8] = b"\x1f\x8b\x08\x04";

    /// Alignment index formats present next to a variant.
    #[derive(Clone, Copy)]
    enum Index {
        /// BAM index, `<file>.bai`.
        Bai,
        /// Coordinate-sorted index, `<file>.csi`.
        Csi,
        /// CRAM index, `<file>.crai`.
        Crai,
    }

    impl Index {
        /// Suffix appended to the alignment file name.
        const fn suffix(self) -> &'static str {
            match self {
                Self::Bai => ".bai",
                Self::Csi => ".csi",
                Self::Crai => ".crai",
            }
        }

        /// Leading bytes of a valid index of this kind.
        const fn magic(self) -> &'static [u8] {
            match self {
                Self::Bai => b"BAI\x01",
                Self::Csi => BGZF_MAGIC,
                // Plain gzip (no BGZF extra field) holding text.
                Self::Crai => b"\x1f\x8b\x08\x00",
            }
        }
    }

    /// How a variant stores the bases that match the reference.
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Reference {
        /// Decodes without a FASTA: BAM, SAM, and CRAM containers that embed
        /// a reference or store bases without one.
        SelfContained,
        /// Some CRAM slices need `examples/contigs.fa`. Each such slice
        /// records an MD5 of its reference span, so a wrong FASTA is rejected.
        External,
        /// External-reference CRAM made only of multi-reference slices. The
        /// CRAM specification zeroes their reference MD5, so `HTSlib` cannot
        /// detect a wrong FASTA and decodes the wrong bases.
        ExternalUnchecked,
    }

    /// One golden file and the properties that its tests rely on.
    struct Variant {
        /// File name inside [`VARIANT_DIR`].
        file: &'static str,
        /// Leading bytes identifying the container format and version.
        magic: &'static [u8],
        /// Index stored next to the file, if any.
        index: Option<Index>,
        /// Whether decoding needs the external reference.
        reference: Reference,
        /// False where `HTSlib` replaces absent base qualities (`*`) with
        /// stored values. All `example_1` records lack base qualities.
        keeps_absent_base_quals: bool,
    }

    impl Variant {
        /// A variant that decodes to records identical to the original.
        const fn new(
            file: &'static str,
            magic: &'static [u8],
            index: Option<Index>,
            reference: Reference,
        ) -> Self {
            Self {
                file,
                magic,
                index,
                reference,
                keeps_absent_base_quals: true,
            }
        }

        /// Marks a variant whose decoder does not restore absent qualities.
        const fn losing_absent_base_quals(self) -> Self {
            Self {
                keeps_absent_base_quals: false,
                ..self
            }
        }

        /// Path relative to the repository root.
        fn path(&self) -> String {
            format!("{VARIANT_DIR}/{}", self.file)
        }

        /// Whether this variant can be decoded without `examples/contigs.fa`.
        fn is_self_contained(&self) -> bool {
            self.reference == Reference::SelfContained
        }
    }

    /// CRAM magic followed by major and minor version bytes.
    const CRAM_2_0: &[u8] = b"CRAM\x02\x00";
    /// See [`CRAM_2_0`].
    const CRAM_2_1: &[u8] = b"CRAM\x02\x01";
    /// See [`CRAM_2_0`].
    const CRAM_3_0: &[u8] = b"CRAM\x03\x00";
    /// See [`CRAM_2_0`].
    const CRAM_3_1: &[u8] = b"CRAM\x03\x01";

    /// Every golden variant, described in `generate.sh`.
    const VARIANTS: [Variant; 29] = {
        use Index::{Bai, Crai, Csi};
        use Reference::{External, ExternalUnchecked, SelfContained};
        [
            Variant::new("bam_csi.bam", BGZF_MAGIC, Some(Csi), SelfContained),
            Variant::new("bam_uncompressed.bam", BGZF_MAGIC, Some(Bai), SelfContained),
            Variant::new(
                "bam_interleaved_levels.bam",
                BGZF_MAGIC,
                Some(Bai),
                SelfContained,
            ),
            Variant::new("sam_plain.sam", b"@HD\t", None, SelfContained),
            Variant::new("sam_bgzf.sam.gz", BGZF_MAGIC, Some(Csi), SelfContained),
            Variant::new(
                "cram_v2.0_external_ref.cram",
                CRAM_2_0,
                Some(Crai),
                External,
            ),
            Variant::new(
                "cram_v2.0_embedded_ref.cram",
                CRAM_2_0,
                Some(Crai),
                SelfContained,
            ),
            Variant::new("cram_v2.0_no_ref.cram", CRAM_2_0, Some(Crai), SelfContained)
                .losing_absent_base_quals(),
            Variant::new(
                "cram_v2.1_external_ref.cram",
                CRAM_2_1,
                Some(Crai),
                External,
            ),
            Variant::new(
                "cram_v2.1_embedded_ref.cram",
                CRAM_2_1,
                Some(Crai),
                SelfContained,
            ),
            Variant::new("cram_v2.1_no_ref.cram", CRAM_2_1, Some(Crai), SelfContained)
                .losing_absent_base_quals(),
            Variant::new("cram_v2.1_interleaved.cram", CRAM_2_1, Some(Crai), External)
                .losing_absent_base_quals(),
            Variant::new(
                "cram_v3.0_external_ref.cram",
                CRAM_3_0,
                Some(Crai),
                External,
            ),
            Variant::new(
                "cram_v3.0_embedded_ref.cram",
                CRAM_3_0,
                Some(Crai),
                SelfContained,
            ),
            Variant::new("cram_v3.0_no_ref.cram", CRAM_3_0, Some(Crai), SelfContained),
            Variant::new("cram_v3.0_interleaved.cram", CRAM_3_0, Some(Crai), External),
            Variant::new(
                "cram_v3.1_external_ref.cram",
                CRAM_3_1,
                Some(Crai),
                External,
            ),
            Variant::new(
                "cram_v3.1_embedded_ref.cram",
                CRAM_3_1,
                Some(Crai),
                SelfContained,
            ),
            Variant::new("cram_v3.1_no_ref.cram", CRAM_3_1, Some(Crai), SelfContained),
            Variant::new("cram_v3.1_interleaved.cram", CRAM_3_1, Some(Crai), External),
            Variant::new(
                "cram_v3.1_consensus_ref.cram",
                CRAM_3_1,
                Some(Crai),
                SelfContained,
            ),
            Variant::new(
                "cram_v3.1_profile_fast.cram",
                CRAM_3_1,
                Some(Crai),
                External,
            ),
            Variant::new(
                "cram_v3.1_profile_archive.cram",
                CRAM_3_1,
                Some(Crai),
                External,
            ),
            Variant::new("cram_v3.1_gzip_only.cram", CRAM_3_1, Some(Crai), External),
            Variant::new(
                "cram_v3.1_uncompressed.cram",
                CRAM_3_1,
                Some(Crai),
                External,
            ),
            Variant::new(
                "cram_v3.1_multi_ref_slice.cram",
                CRAM_3_1,
                Some(Crai),
                ExternalUnchecked,
            ),
            Variant::new(
                "cram_v3.1_slice_per_record.cram",
                CRAM_3_1,
                Some(Crai),
                ExternalUnchecked,
            ),
            Variant::new(
                "cram_v3.1_container_per_record.cram",
                CRAM_3_1,
                Some(Crai),
                External,
            ),
            Variant::new("cram_v3.1_unsorted.cram", CRAM_3_1, None, External),
        ]
    };

    /// A CLI invocation, without the input path, compared across variants.
    struct CliCase {
        /// Arguments placed before the input path.
        args: &'static [&'static str],
        /// Whether the output reports base qualities.
        shows_base_quals: bool,
    }

    /// Subcommands covering every output path that reads record contents.
    const CLI_CASES: [CliCase; 13] = [
        CliCase {
            args: &["read-table-show-mods"],
            shows_base_quals: false,
        },
        CliCase {
            args: &["read-table-show-mods", "--seq-full"],
            shows_base_quals: false,
        },
        CliCase {
            args: &["read-table-show-mods", "--seq-full", "--show-base-qual"],
            shows_base_quals: true,
        },
        CliCase {
            args: &["read-table-show-mods", "--mod-region", "dummyIII:30-50"],
            shows_base_quals: false,
        },
        CliCase {
            args: &["read-table-hide-mods"],
            shows_base_quals: false,
        },
        CliCase {
            args: &["read-info"],
            shows_base_quals: false,
        },
        CliCase {
            args: &["read-info", "--detailed"],
            shows_base_quals: false,
        },
        CliCase {
            args: &["read-stats"],
            shows_base_quals: false,
        },
        CliCase {
            args: &["window-dens", "--win", "2", "--step", "1"],
            shows_base_quals: true,
        },
        CliCase {
            args: &["window-grad", "--win", "4", "--step", "1"],
            shows_base_quals: true,
        },
        CliCase {
            args: &[
                "find-modified-reads",
                "any-dens-above",
                "--high",
                "0.5",
                "--tag",
                "T",
                "--win",
                "2",
                "--step",
                "1",
            ],
            shows_base_quals: false,
        },
        CliCase {
            args: &["peek"],
            shows_base_quals: false,
        },
        CliCase {
            args: &["read-info", "--detailed", "--mapq-filter", "1"],
            shows_base_quals: false,
        },
    ];

    /// Regions (samtools notation) covering whole contigs, partial overlaps,
    /// a gap without reads, and the reverse-strand read.
    const REGIONS: [&str; 5] = [
        "dummyI",
        "dummyII:5-10",
        "dummyIII:1-10",
        "dummyIII:20-30",
        "dummyIII:70-76",
    ];

    /// A uniquely named temporary directory removed on drop.
    struct TempDir(PathBuf);

    impl TempDir {
        /// Creates an empty directory.
        fn new() -> Self {
            let path = std::env::temp_dir()
                .join(format!("nanalogue_example_1_formats_{}", uuid::v4_random()));
            fs::create_dir_all(&path).expect("temporary directory should be creatable");
            Self(path)
        }

        /// Copies a variant and its index into this directory, returning its
        /// file name, which resolves relative to this directory.
        fn copy_variant(&self, variant: &Variant) -> &'static str {
            let _: u64 = fs::copy(variant.path(), self.0.join(variant.file))
                .expect("variant should be copyable");
            if let Some(index) = variant.index {
                let name = format!("{}{}", variant.file, index.suffix());
                let _: u64 = fs::copy(format!("{VARIANT_DIR}/{name}"), self.0.join(name))
                    .expect("index should be copyable");
            }
            variant.file
        }
    }

    impl Drop for TempDir {
        fn drop(&mut self) {
            let removal = fs::remove_dir_all(&self.0);
            // A second panic during unwinding would abort and hide the failure.
            if !std::thread::panicking() {
                removal.expect("temporary directory should be removable");
            }
        }
    }

    /// Runs the executable in `dir` with `HTSlib`'s reference search limited to
    /// `ref_path`, so no test consults a reference server. An empty
    /// `REF_CACHE` disables the local reference cache, which `HTSlib` would
    /// otherwise create and fill (changing permissions along the way) whenever
    /// a reference is found through `REF_PATH`.
    fn run_nanalogue(args: &[&str], dir: &Path, ref_path: &str) -> Output {
        Command::new(env!("CARGO_BIN_EXE_nanalogue"))
            .args(args)
            .current_dir(dir)
            .env("REF_PATH", ref_path)
            .env("REF_CACHE", "")
            .output()
            .expect("nanalogue executable should run")
    }

    /// Runs the executable from the repository root, where the relative `UR`
    /// header tag of external-reference CRAMs resolves.
    fn run_in_repo(args: &[&str]) -> Output {
        run_nanalogue(
            args,
            Path::new(env!("CARGO_MANIFEST_DIR")),
            NO_REFERENCE_SEARCH,
        )
    }

    /// Order-insensitive form of a successful run's stdout. JSON arrays are
    /// compared record by record; other outputs line by line.
    fn normalized_stdout(output: &Output, context: &str) -> Vec<String> {
        assert!(
            output.status.success(),
            "{context} should succeed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let stdout = String::from_utf8(output.stdout.clone()).expect("stdout should be UTF-8");
        if stdout.trim_start().starts_with('[') {
            let values: Vec<serde_json::Value> =
                serde_json::from_str(&stdout).expect("JSON output should parse");
            let mut records: Vec<String> = values.iter().map(ToString::to_string).collect();
            records.sort();
            records
        } else {
            sort_output_lines(&stdout)
        }
    }

    /// Normalized stdout of a successful run, or the error lines of a failed
    /// one without `HTSlib` diagnostics and nanalogue warnings, which differ
    /// between formats.
    fn outcome(output: &Output, context: &str) -> Result<Vec<String>, Vec<String>> {
        if output.status.success() {
            Ok(normalized_stdout(output, context))
        } else {
            Err(String::from_utf8_lossy(&output.stderr)
                .lines()
                .filter(|line| !line.starts_with('[') && !line.starts_with('#'))
                .map(str::to_owned)
                .collect())
        }
    }

    /// Drops the last tab-separated column, which holds base qualities in
    /// every output that prints them.
    fn without_last_column(lines: &[String]) -> Vec<String> {
        lines
            .iter()
            .map(|line| line.rsplit_once('\t').map_or(line.as_str(), |pair| pair.0))
            .map(str::to_owned)
            .collect()
    }

    /// Asserts a run failed with an ordinary error rather than a panic or
    /// partial output.
    fn assert_clean_failure(output: &Output, context: &str) {
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(
            output.status.code(),
            Some(1),
            "{context} should exit with an error: {stderr}"
        );
        assert!(output.stdout.is_empty(), "{context} should not emit output");
        assert!(
            stderr.contains("rust_htslib error"),
            "{context} should report the decoding error: {stderr}"
        );
        assert!(!stderr.contains("panicked"), "{context} must not panic");
    }

    /// Record fields that every variant must reproduce.
    #[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
    struct RecordSummary {
        /// Read name.
        name: Vec<u8>,
        /// SAM flags.
        flags: u16,
        /// Reference id.
        tid: i32,
        /// Zero-based position.
        pos: i64,
        /// Mapping quality.
        mapq: u8,
        /// CIGAR string.
        cigar: String,
        /// Read bases.
        seq: Vec<u8>,
        /// Base qualities, omitted for variants that do not keep them.
        qual: Option<Vec<u8>>,
        /// `MM` tag.
        mm: String,
        /// `ML` tag.
        ml: Vec<u8>,
    }

    /// Reads every record of a reader into sorted summaries.
    fn summaries<R: rust_htslib::bam::Read>(mut reader: R, keep_qual: bool) -> Vec<RecordSummary> {
        let mut summaries: Vec<RecordSummary> = reader
            .records()
            .map(|result| {
                let record = result.expect("record should decode");
                let Ok(Aux::String(mm)) = record.aux(b"MM") else {
                    unreachable!("MM should be a string tag")
                };
                let Ok(Aux::ArrayU8(ml)) = record.aux(b"ML") else {
                    unreachable!("ML should be a u8 array tag")
                };
                RecordSummary {
                    name: record.qname().to_vec(),
                    flags: record.flags(),
                    tid: record.tid(),
                    pos: record.pos(),
                    mapq: record.mapq(),
                    cigar: record.cigar().to_string(),
                    seq: record.seq().as_bytes(),
                    qual: keep_qual.then(|| record.qual().to_vec()),
                    mm: mm.to_owned(),
                    ml: ml.iter().collect(),
                }
            })
            .collect();
        summaries.sort();
        summaries
    }

    /// Read names returned by an indexed fetch, sorted.
    fn fetched_names(path: &str, fetch: FetchDefinition<'_>) -> Result<Vec<Vec<u8>>, Error> {
        let mut reader = nanalogue_indexed_bam_reader(path, fetch)?;
        let mut names = reader
            .records()
            .map(|result| Ok(result?.qname().to_vec()))
            .collect::<Result<Vec<_>, Error>>()?;
        names.sort();
        Ok(names)
    }

    /// Fetches that exercise each index on mapped and unmapped reads.
    fn fetches() -> [FetchDefinition<'static>; 7] {
        [
            FetchDefinition::All,
            FetchDefinition::String(b"dummyI"),
            FetchDefinition::CompleteTid(1),
            FetchDefinition::RegionString(b"dummyIII", 10, 20),
            FetchDefinition::RegionString(b"dummyIII", 20, 30),
            FetchDefinition::RegionString(b"dummyII", 30, 60),
            FetchDefinition::Unmapped,
        ]
    }

    /// Writes `examples/contigs.fa` with every `C` replaced by `A`, keeping
    /// contig names and lengths.
    fn write_wrong_reference(dir: &Path) {
        let reference = fs::read_to_string(REFERENCE).expect("reference should be readable");
        let wrong: String = reference
            .lines()
            .map(|line| {
                if line.starts_with('>') {
                    format!("{line}\n")
                } else {
                    format!("{}\n", line.replace('C', "A"))
                }
            })
            .collect();
        assert_ne!(wrong, reference, "the wrong reference must differ");
        fs::create_dir_all(dir.join("examples")).expect("examples directory should be creatable");
        fs::write(dir.join(REFERENCE), wrong).expect("wrong reference should be writable");
    }

    /// The golden directory holds exactly the listed variants and indices,
    /// each with the expected container format, version, and index format.
    #[test]
    fn golden_directory_matches_variant_table() {
        let mut expected = BTreeSet::from(["generate.sh".to_owned()]);
        for variant in &VARIANTS {
            let contents = fs::read(variant.path()).expect("variant should be readable");
            assert!(
                contents.starts_with(variant.magic),
                "{} has an unexpected format or version",
                variant.file
            );
            assert!(
                expected.insert(variant.file.to_owned()),
                "duplicate variant"
            );
            if let Some(index) = variant.index {
                let name = format!("{}{}", variant.file, index.suffix());
                let index_contents =
                    fs::read(format!("{VARIANT_DIR}/{name}")).expect("index should be readable");
                assert!(
                    index_contents.starts_with(index.magic()),
                    "{name} has an unexpected format"
                );
                assert!(expected.insert(name), "duplicate index");
            }
        }
        let actual: BTreeSet<String> = fs::read_dir(VARIANT_DIR)
            .expect("variant directory should be readable")
            .map(|entry| {
                entry
                    .expect("directory entry should be readable")
                    .file_name()
                    .into_string()
                    .expect("file names should be UTF-8")
            })
            .collect();
        assert_eq!(actual, expected, "golden files and VARIANTS must agree");
    }

    /// Every CLI subcommand prints the same output for every variant as for
    /// the original BAM. For variants that do not keep absent base qualities,
    /// the last column, which holds qualities, is excluded from outputs that
    /// print them; see [`cram_v2_no_ref_replaces_absent_base_quals`].
    #[test]
    fn cli_output_matches_original_for_every_variant() {
        for case in &CLI_CASES {
            let original = run_in_repo(&[case.args, &[ORIGINAL]].concat());
            let expected = normalized_stdout(&original, ORIGINAL);
            assert!(!expected.is_empty(), "the original should produce output");
            for variant in &VARIANTS {
                let path = variant.path();
                let context = format!("{} on {path}", case.args.join(" "));
                let output = run_in_repo(&[case.args, &[path.as_str()]].concat());
                let actual = normalized_stdout(&output, &context);
                if case.shows_base_quals && !variant.keeps_absent_base_quals {
                    assert_eq!(
                        without_last_column(&actual),
                        without_last_column(&expected),
                        "{context} differs from the original outside base qualities"
                    );
                } else {
                    assert_eq!(actual, expected, "{context} differs from the original");
                }
            }
        }
    }

    /// Region queries match the original. Indexed variants must use their
    /// BAI, CSI, or CRAI; unindexed ones fall back to a full scan and warn.
    #[test]
    fn cli_region_queries_match_original_and_use_indices() {
        for region in REGIONS {
            for full_region in [false, true] {
                let mut args = vec!["read-table-show-mods", "--seq-full", "--region", region];
                if full_region {
                    args.push("--full-region");
                }
                let original = run_in_repo(&[args.as_slice(), &[ORIGINAL]].concat());
                let expected = outcome(&original, ORIGINAL);
                for variant in &VARIANTS {
                    let path = variant.path();
                    let context = format!("{} on {path}", args.join(" "));
                    let output = run_in_repo(&[args.as_slice(), &[path.as_str()]].concat());
                    assert_eq!(
                        outcome(&output, &context),
                        expected,
                        "{context} differs from the original"
                    );
                    let warned = String::from_utf8_lossy(&output.stderr).contains(NO_INDEX_WARNING);
                    assert_eq!(
                        warned,
                        variant.index.is_none(),
                        "{context}: index use does not match the variant's index"
                    );
                }
            }
        }
    }

    /// Variants that need no reference decode in process to the original
    /// records, field by field.
    #[test]
    fn self_contained_variants_decode_to_original_records() {
        let original = summaries(
            nanalogue_bam_reader(ORIGINAL).expect("original should open"),
            true,
        );
        assert_eq!(original.len(), 4, "the original has four records");
        for variant in VARIANTS
            .iter()
            .filter(|variant| variant.is_self_contained())
        {
            let reader = nanalogue_bam_reader(&variant.path()).expect("variant should open");
            let keep_qual = variant.keeps_absent_base_quals;
            let expected = summaries(
                nanalogue_bam_reader(ORIGINAL).expect("original should open"),
                keep_qual,
            );
            assert_eq!(
                summaries(reader, keep_qual),
                expected,
                "{} records differ from the original",
                variant.file
            );
        }
    }

    /// CRAM 2.x reference-free encoding does not restore absent base
    /// qualities: `HTSlib` decodes a quality of 30 at the first base and 255
    /// elsewhere instead of 255 everywhere. Everything else round-trips, as
    /// checked by the other tests. This pins the `HTSlib` 1.19 behaviour.
    ///
    /// NOTE: this test deliberately asserts incorrect behaviour. The correct
    /// result would be 255 at every base. If `HTSlib` fixes this, the test
    /// fails; then drop `losing_absent_base_quals` from those variants.
    #[test]
    fn cram_v2_no_ref_replaces_absent_base_quals() {
        let lossy: Vec<&Variant> = VARIANTS
            .iter()
            .filter(|variant| !variant.keeps_absent_base_quals)
            .collect();
        assert_eq!(lossy.len(), 3, "CRAM 2.0, 2.1 and 2.1 interleaved");
        for variant in lossy.iter().filter(|variant| variant.is_self_contained()) {
            let reader = nanalogue_bam_reader(&variant.path()).expect("variant should open");
            for summary in summaries(reader, true) {
                let qual = summary.qual.expect("qualities were requested");
                // Incorrect on purpose: mapped reads should also start at 255.
                let expected_first = if summary.flags & 4 == 0 { 30 } else { 255 };
                assert_eq!(
                    qual.first(),
                    Some(&expected_first),
                    "{}: unexpected first quality",
                    variant.file
                );
                assert!(
                    qual.iter().skip(1).all(|&value| value == 255),
                    "{}: later qualities stay absent",
                    variant.file
                );
            }
        }

        // Printed qualities are the only CLI difference.
        let args = ["read-table-show-mods", "--seq-full", "--show-base-qual"];
        let original = run_in_repo(&[args.as_slice(), &[ORIGINAL]].concat());
        let expected = normalized_stdout(&original, ORIGINAL);
        for variant in lossy {
            let path = variant.path();
            let output = run_in_repo(&[args.as_slice(), &[path.as_str()]].concat());
            let actual = normalized_stdout(&output, &path);
            // Incorrect on purpose: the qualities should match the original.
            assert_ne!(actual, expected, "{path}: qualities should differ");
            assert_eq!(
                without_last_column(&actual),
                without_last_column(&expected),
                "{path}: only the qualities column may differ"
            );
        }
    }

    /// Indexed self-contained variants fetch the same reads as the original
    /// through the public indexed reader; unindexed ones report a missing
    /// index.
    #[test]
    fn self_contained_variants_fetch_like_original() {
        for variant in VARIANTS
            .iter()
            .filter(|variant| variant.is_self_contained())
        {
            let path = variant.path();
            for (original_fetch, fetch) in fetches().into_iter().zip(fetches()) {
                let expected = fetched_names(ORIGINAL, original_fetch).expect("original fetches");
                match (variant.index, fetched_names(&path, fetch)) {
                    (Some(_), Ok(names)) => assert_eq!(
                        names, expected,
                        "{path}: fetched reads differ from the original"
                    ),
                    (None, Err(Error::RustHtslibError(error))) => assert!(
                        matches!(*error, HtslibError::BamInvalidIndex { .. }),
                        "{path}: expected a missing-index error"
                    ),
                    (Some(_), Err(error)) => {
                        unreachable!("{path}: indexed fetch failed: {error}")
                    }
                    (None, Ok(_)) => unreachable!("{path}: unindexed fetch should fail"),
                    (None, Err(error)) => {
                        unreachable!("{path}: expected a missing-index error, got {error}")
                    }
                }
            }
        }
    }

    /// Outside the repository root, the relative `UR` tag no longer resolves:
    /// self-contained variants still decode, and every external-reference
    /// CRAM fails cleanly instead of guessing bases.
    #[test]
    fn only_self_contained_variants_decode_without_reference() {
        let args = ["read-table-show-mods", "--seq-full"];
        let original = run_in_repo(&[args.as_slice(), &[ORIGINAL]].concat());
        let expected = normalized_stdout(&original, ORIGINAL);
        for variant in &VARIANTS {
            let dir = TempDir::new();
            let file = dir.copy_variant(variant);
            let output = run_nanalogue(
                &[args.as_slice(), &[file]].concat(),
                &dir.0,
                NO_REFERENCE_SEARCH,
            );
            let context = format!("{file} without a reference");
            if variant.is_self_contained() {
                assert_eq!(
                    normalized_stdout(&output, &context),
                    expected,
                    "{context} differs from the original"
                );
            } else {
                assert_clean_failure(&output, &context);
            }
        }
    }

    /// External-reference CRAMs also resolve their reference through an M5
    /// cache named by `REF_PATH`, as used for large reference collections.
    #[test]
    fn external_reference_resolves_through_ref_path_md5_cache() {
        let reference = fs::read_to_string(REFERENCE).expect("reference should be readable");
        let dir = TempDir::new();
        let cache = dir.0.join("md5_cache");
        fs::create_dir_all(&cache).expect("cache directory should be creatable");
        let mut sequences: Vec<(&str, String)> = Vec::new();
        for line in reference.lines() {
            if let Some(name) = line.strip_prefix('>') {
                sequences.push((name, String::new()));
            } else {
                let last = sequences
                    .last_mut()
                    .expect("reference starts with a header");
                last.1.push_str(line);
            }
        }
        assert_eq!(sequences.len(), CONTIG_MD5S.len(), "three contigs");
        for ((name, sequence), (expected_name, md5)) in sequences.into_iter().zip(CONTIG_MD5S) {
            assert_eq!(name, expected_name, "contig order");
            fs::write(cache.join(md5), sequence).expect("cache entry should be writable");
        }
        let ref_path = format!("{}/%s", cache.display());

        let args = ["read-table-show-mods", "--seq-full"];
        let original = run_in_repo(&[args.as_slice(), &[ORIGINAL]].concat());
        let expected = normalized_stdout(&original, ORIGINAL);
        for variant in VARIANTS
            .iter()
            .filter(|variant| !variant.is_self_contained())
        {
            let file = dir.copy_variant(variant);
            let output = run_nanalogue(&[args.as_slice(), &[file]].concat(), &dir.0, &ref_path);
            let context = format!("{file} with an M5 cache");
            assert_eq!(
                normalized_stdout(&output, &context),
                expected,
                "{context} differs from the original"
            );
        }
    }

    /// A wrong FASTA at the `UR` path, with the right names and lengths, is
    /// rejected wherever a slice records its reference MD5. In multi-reference
    /// slices the CRAM specification leaves that MD5 empty, so `HTSlib` decodes
    /// wrong bases without an error. This pins that limitation: if it starts
    /// failing because the wrong reference is detected, move those variants
    /// to [`Reference::External`].
    ///
    /// NOTE: this test deliberately asserts incorrect behaviour. The correct
    /// result would be a clean failure for every external-reference variant.
    #[test]
    fn wrong_reference_is_rejected_except_in_multi_reference_slices() {
        let args = ["read-table-show-mods", "--seq-full"];
        let original = run_in_repo(&[args.as_slice(), &[ORIGINAL]].concat());
        let expected = normalized_stdout(&original, ORIGINAL);
        for variant in VARIANTS
            .iter()
            .filter(|variant| !variant.is_self_contained())
        {
            let dir = TempDir::new();
            write_wrong_reference(&dir.0);
            let file = dir.copy_variant(variant);
            let output = run_nanalogue(
                &[args.as_slice(), &[file]].concat(),
                &dir.0,
                NO_REFERENCE_SEARCH,
            );
            let context = format!("{file} with a wrong reference");
            if variant.reference == Reference::ExternalUnchecked {
                // Incorrect on purpose: these should fail cleanly as well.
                let actual = normalized_stdout(&output, &context);
                assert_eq!(actual.len(), expected.len(), "{context}: same reads");
                assert_ne!(
                    actual, expected,
                    "{context}: bases come from the wrong FASTA"
                );
            } else {
                assert_clean_failure(&output, &context);
            }
        }
    }
}
