#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Public local-file reader, writer, index, fetch, and output-path contracts.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::{
        DNARestrictive, Error, nanalogue_bam_reader, nanalogue_indexed_bam_reader, uuid,
        write_bam_denovo, write_cram_denovo, write_fasta,
    };
    use rust_htslib::bam::record::{Cigar, CigarString};
    use rust_htslib::bam::{FetchDefinition, Read as _, Record};
    use rust_htslib::errors::Error as HtslibError;
    use rust_htslib::faidx;
    use std::fs;
    use std::num::NonZeroU32;
    use std::path::{Path, PathBuf};
    use std::str::FromStr as _;

    /// A uniquely named temporary directory removed even when a test returns early.
    struct TempDir(PathBuf);

    impl TempDir {
        fn new(label: &str) -> Self {
            let path = std::env::temp_dir()
                .join(format!("nanalogue_file_io_{label}_{}", uuid::v4_random()));
            fs::create_dir_all(&path).expect("temporary directory should be creatable");
            Self(path)
        }

        fn join(&self, name: &str) -> PathBuf {
            self.0.join(name)
        }
    }

    impl Drop for TempDir {
        fn drop(&mut self) {
            fs::remove_dir_all(&self.0).expect("temporary directory should be removable");
        }
    }

    /// Builds a four-base mapped record with an asymmetric identity and position.
    fn mapped_record(read_id: &str, tid: i32, pos: i64, reverse: bool) -> Record {
        let mut record = Record::new();
        record.set(
            read_id.as_bytes(),
            Some(&CigarString(vec![Cigar::Match(4)])),
            b"ACGT",
            &[31, 32, 33, 34],
        );
        record.set_flags(0);
        record.set_tid(tid);
        record.set_pos(pos);
        record.set_mapq(42);
        if reverse {
            record.set_reverse();
        }
        record
    }

    /// Builds an unmapped record which must sort after every mapped record.
    fn unmapped_record() -> Record {
        let mut record = Record::new();
        record.set(b"unmapped", None, b"TGCA", &[21, 22, 23, 24]);
        record.set_unmapped();
        record.set_tid(-1);
        record.set_pos(-1);
        record
    }

    /// Returns read names after a public indexed-reader fetch.
    fn fetched_names(path: &Path, fetch: FetchDefinition<'_>) -> Result<Vec<String>, Error> {
        let mut reader = nanalogue_indexed_bam_reader(path, fetch)?;
        reader
            .records()
            .map(|result| {
                let record = result?;
                Ok(String::from_utf8(record.qname().to_vec())
                    .expect("test read names are valid UTF-8"))
            })
            .collect()
    }

    /// One written BAM exercises numeric, named, interval, all, and unmapped
    /// fetches. Boundary intervals distinguish overlap from mere adjacency.
    #[test]
    fn written_bam_supports_all_local_fetch_forms() -> Result<(), Error> {
        let temp = TempDir::new("fetches");
        let bam_path = temp.join("asymmetric.bam");
        let records = vec![
            mapped_record("same_pos_forward", 0, 5, false),
            mapped_record("same_pos_reverse", 0, 5, true),
            mapped_record("later_chr_a", 0, 30, false),
            mapped_record("chr_b", 1, 2, false),
            unmapped_record(),
        ];

        write_bam_denovo(
            records,
            [("chrA".to_owned(), 50), ("chrB".to_owned(), 17)],
            ["flowcell-7".to_owned()],
            ["asymmetric integration fixture".to_owned()],
            &bam_path,
        )?;

        let index_path = temp.join("asymmetric.bam.bai");
        assert!(
            index_path.is_file(),
            "the BAI must be written beside the BAM"
        );
        assert!(
            fs::metadata(&index_path)?.len() > 0,
            "the BAI must be nonempty"
        );

        let reader = nanalogue_bam_reader(&bam_path)?;
        assert_eq!(reader.header().target_names(), [b"chrA", b"chrB"]);
        let header = std::str::from_utf8(reader.header().as_bytes())
            .expect("generated header must be UTF-8");
        assert!(header.contains("@RG\tID:flowcell-7\tPL:ONT"));
        assert!(header.contains("@CO\tasymmetric integration fixture"));

        assert_eq!(
            fetched_names(&bam_path, FetchDefinition::All)?,
            [
                "same_pos_forward",
                "same_pos_reverse",
                "later_chr_a",
                "chr_b",
                "unmapped",
            ]
        );
        assert_eq!(
            fetched_names(&bam_path, FetchDefinition::CompleteTid(1))?,
            ["chr_b"]
        );
        assert_eq!(
            fetched_names(&bam_path, FetchDefinition::String(b"chrA"))?,
            ["same_pos_forward", "same_pos_reverse", "later_chr_a"]
        );
        assert_eq!(
            fetched_names(&bam_path, FetchDefinition::Region(0, 6, 7))?,
            ["same_pos_forward", "same_pos_reverse"]
        );
        assert!(
            fetched_names(&bam_path, FetchDefinition::RegionString(b"chrA", 9, 30))?.is_empty(),
            "a read ending at 9 and one starting at 30 are adjacent, not overlapping"
        );
        assert_eq!(
            fetched_names(&bam_path, FetchDefinition::Unmapped)?,
            ["unmapped"]
        );
        Ok(())
    }

    /// An invalid fetch and absent local inputs preserve their public error
    /// categories without creating files in the temporary directory.
    #[test]
    fn local_reader_errors_do_not_create_inputs() {
        let temp = TempDir::new("reader_errors");
        let missing = temp.join("missing.bam");
        let ordinary_error = nanalogue_bam_reader(&missing).expect_err("missing BAM must fail");
        assert!(matches!(
            ordinary_error,
            Error::RustHtslibError(source)
                if matches!(source.as_ref(), HtslibError::FileNotFound { path } if path == &missing)
        ));
        assert!(!missing.exists());

        let indexed_error = nanalogue_indexed_bam_reader(&missing, FetchDefinition::All)
            .expect_err("missing indexed BAM must fail");
        assert!(matches!(
            indexed_error,
            Error::RustHtslibError(source)
                if matches!(source.as_ref(), HtslibError::FileNotFound { path } if path == &missing)
        ));
        assert!(!missing.exists());

        let existing = Path::new("examples/example_1.bam");
        let fetch_error = nanalogue_indexed_bam_reader(existing, FetchDefinition::CompleteTid(99))
            .expect_err("an out-of-range numeric target must fail");
        assert!(matches!(
            fetch_error,
            Error::RustHtslibError(source) if matches!(source.as_ref(), HtslibError::Fetch)
        ));
        assert_eq!(
            fs::read_dir(&temp.0)
                .expect("temporary directory remains readable")
                .count(),
            0
        );
    }

    /// Index creation happens only after the BAM stream closes. Blocking the
    /// conventional BAI path therefore reports an index error while leaving a
    /// complete, sequentially readable BAM and the pre-existing path intact.
    #[test]
    fn blocked_bai_path_preserves_completed_bam() -> Result<(), Error> {
        let temp = TempDir::new("blocked_bai");
        let bam_path = temp.join("output.bam");
        let bai_path = temp.join("output.bam.bai");
        fs::create_dir_all(&bai_path)?;

        let error = write_bam_denovo(
            vec![mapped_record("persisted", 0, 11, false)],
            vec![("chrA".to_owned(), 30)],
            Vec::<String>::new(),
            Vec::<String>::new(),
            &bam_path,
        )
        .expect_err("a directory at the BAI path must block index creation");
        assert!(matches!(
            error,
            Error::RustHtslibError(source)
                if matches!(source.as_ref(), HtslibError::BamWriteIndex)
        ));
        assert!(bai_path.is_dir(), "the blocker must not be replaced");

        let mut reader = nanalogue_bam_reader(&bam_path)?;
        let records = reader.records().collect::<Result<Vec<_>, _>>()?;
        assert_eq!(records.len(), 1);
        let record = records.first().expect("one record was read");
        assert_eq!(record.qname(), b"persisted");
        assert_eq!(record.tid(), 0);
        assert_eq!(record.pos(), 11);
        Ok(())
    }

    /// The public CRAM writer rejects a thread count beyond `HTSlib`'s signed
    /// range before creating output, while preserving its reference.
    #[test]
    fn cram_thread_count_must_fit_htslib() -> Result<(), Error> {
        let temp = TempDir::new("cram_threads");
        let reference = temp.join("reference.fa");
        write_fasta(
            [(
                "chrA".to_owned(),
                DNARestrictive::from_str("ACGTACGT").expect("valid DNA"),
            )],
            &reference,
        )?;
        faidx::build(&reference).expect("reference index should be written");
        let original_reference = fs::read(&reference)?;
        let cram_path = temp.join("output.cram");
        let too_many_threads = NonZeroU32::new(
            u32::try_from(i32::MAX)
                .expect("i32::MAX fits u32")
                .checked_add(1)
                .expect("i32::MAX + 1 fits u32"),
        )
        .expect("the out-of-range count remains nonzero");

        let error = write_cram_denovo(
            Vec::<Record>::new(),
            vec![("chrA".to_owned(), 8)],
            Vec::<String>::new(),
            Vec::<String>::new(),
            &cram_path,
            &reference,
            too_many_threads,
        )
        .expect_err("HTSlib's thread count must fit i32");
        assert!(matches!(error, Error::IntConversionError(_)));
        assert!(
            !cram_path.exists(),
            "invalid thread counts must be rejected before opening the CRAM"
        );
        assert!(!temp.join("output.cram.crai").exists());
        assert_eq!(fs::read(reference)?, original_reference);
        Ok(())
    }

    /// FASTA output truncates an existing file, keeps sequence order exactly,
    /// and rejects a directory output without modifying that directory.
    #[test]
    fn fasta_output_replaces_files_but_not_directories() -> Result<(), Error> {
        let temp = TempDir::new("fasta_paths");
        let fasta = temp.join("reference.fa");
        fs::write(&fasta, b"stale bytes that must disappear")?;
        write_fasta(
            [
                (
                    "short".to_owned(),
                    DNARestrictive::from_str("ACG").expect("valid DNA"),
                ),
                (
                    "asymmetric".to_owned(),
                    DNARestrictive::from_str("TTGCA").expect("valid DNA"),
                ),
            ],
            &fasta,
        )?;
        assert_eq!(fs::read(&fasta)?, b">short\nACG\n>asymmetric\nTTGCA\n");

        let blocked = temp.join("blocked.fa");
        fs::create_dir_all(&blocked)?;
        let error = write_fasta(
            [(
                "unused".to_owned(),
                DNARestrictive::from_str("A").expect("valid DNA"),
            )],
            &blocked,
        )
        .expect_err("a directory cannot be replaced with FASTA data");
        assert!(matches!(error, Error::InputOutputError(_)));
        assert!(blocked.is_dir());
        assert_eq!(fs::read_dir(blocked)?.count(), 0);
        Ok(())
    }
}
