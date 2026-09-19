#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Offline coverage for BAM readers opened through absolute `file:` URLs.
//!
//! The fixture and its adjacent BAI are repository data, so these tests cover
//! `HTSlib`'s URL path without a server, network access, or TLS configuration.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::{
        Error, nanalogue_bam_reader_from_url, nanalogue_indexed_bam_reader_from_url,
    };
    use rust_htslib::bam::{FetchDefinition, Read as _};
    use rust_htslib::errors::Error as HtslibError;
    use std::path::Path;
    use url::Url;

    /// `example_1.bam` holds four alignments belonging to three read IDs.
    const EXAMPLE_1_IDS: [&str; 4] = [
        "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575",
        "a4f36092-b4d5-47a9-813e-c22c3b477a0c",
        "a4f36092-b4d5-47a9-813e-c22c3b477a0c",
        "fffffff1-10d2-49cb-8ca3-e8d48979001b",
    ];

    /// The `dummyIII:20-30` half-open interval contains this one alignment.
    const NARROW_INTERVAL_RECORD: (&str, i64) = ("a4f36092-b4d5-47a9-813e-c22c3b477a0c", 23);

    /// Resolves a repository fixture to an absolute file URL as `HTSlib` sees it.
    fn fixture_url(name: &str) -> Url {
        let path = Path::new("examples")
            .join(name)
            .canonicalize()
            .expect("checked-in fixture must exist");
        Url::from_file_path(path).expect("absolute fixture path must form a file URL")
    }

    /// Produces a unique, absent absolute path without creating a resource.
    fn missing_url() -> Url {
        let path = std::env::temp_dir().join(format!(
            "nanalogue_missing_url_reader_fixture_{}.bam",
            nanalogue_core::uuid::v4_random()
        ));
        assert!(
            !path.exists(),
            "the generated missing fixture path is unused"
        );
        Url::from_file_path(path).expect("temporary directory paths are absolute")
    }

    /// Reads names through the ordinary reader and sorts them independently of
    /// fixture order, preserving the intentional duplicated alignment ID.
    fn ordinary_ids(url: &Url) -> Result<Vec<String>, Error> {
        let mut reader = nanalogue_bam_reader_from_url(url)?;
        let mut ids: Vec<String> = reader
            .records()
            .map(|result| result.map(|record| String::from_utf8_lossy(record.qname()).into_owned()))
            .collect::<Result<_, _>>()?;
        ids.sort();
        Ok(ids)
    }

    /// Reads ID and zero-based alignment start pairs from an indexed fetch.
    fn fetched_records(
        url: &Url,
        definition: FetchDefinition<'_>,
    ) -> Result<Vec<(String, i64)>, Error> {
        let mut reader = nanalogue_indexed_bam_reader_from_url(url, definition)?;
        reader
            .records()
            .map(|result| {
                result.map(|record| {
                    (
                        String::from_utf8_lossy(record.qname()).into_owned(),
                        record.pos(),
                    )
                })
            })
            .collect::<Result<_, _>>()
            .map_err(Error::from)
    }

    /// A file URL to the indexed fixture opens with the ordinary reader and
    /// yields exactly the documented four-alignment ID multiset.
    #[test]
    fn ordinary_file_url_reads_the_documented_four_alignment_ids() -> Result<(), Error> {
        let mut expected = EXAMPLE_1_IDS.map(str::to_owned).to_vec();
        expected.sort();

        assert_eq!(
            ordinary_ids(&fixture_url("example_1.bam"))?,
            expected,
            "the file URL reader must preserve all four fixture alignments"
        );
        Ok(())
    }

    /// A file URL also locates the adjacent BAI: fetching all records and a
    /// narrow half-open interval produces independently enumerated records.
    #[test]
    fn indexed_file_url_fetches_all_and_a_narrow_interval() -> Result<(), Error> {
        let url = fixture_url("example_1.bam");
        let mut expected_all = EXAMPLE_1_IDS.map(str::to_owned).to_vec();
        expected_all.sort();
        let mut all: Vec<String> = fetched_records(&url, FetchDefinition::All)?
            .into_iter()
            .map(|(id, _position)| id)
            .collect();
        all.sort();
        assert_eq!(
            all, expected_all,
            "FetchDefinition::All returns every alignment"
        );

        assert_eq!(
            fetched_records(&url, FetchDefinition::RegionString(b"dummyIII", 20, 30),)?,
            vec![(
                NARROW_INTERVAL_RECORD.0.to_owned(),
                NARROW_INTERVAL_RECORD.1
            )],
            "dummyIII:20-30 includes its one overlapping alignment at position 23"
        );
        // The fixture alignment spans [23, 71). These queries distinguish
        // overlap-based interval fetching from merely selecting the contig.
        assert!(
            fetched_records(&url, FetchDefinition::RegionString(b"dummyIII", 20, 23))?.is_empty(),
            "an interval ending at the alignment start must exclude it"
        );
        assert_eq!(
            fetched_records(&url, FetchDefinition::RegionString(b"dummyIII", 70, 71))?,
            vec![(
                NARROW_INTERVAL_RECORD.0.to_owned(),
                NARROW_INTERVAL_RECORD.1
            )],
            "the last aligned base overlaps even though the read starts before the query"
        );
        assert!(
            fetched_records(&url, FetchDefinition::RegionString(b"dummyIII", 71, 72))?.is_empty(),
            "an interval starting at the exclusive alignment end must exclude it"
        );
        Ok(())
    }

    /// URL opening errors retain their exact `HTSlib` variants and URL targets.
    #[test]
    fn url_readers_wrap_missing_and_unindexed_file_errors() {
        let missing = missing_url();
        let missing_error = nanalogue_bam_reader_from_url(&missing)
            .expect_err("a missing file URL cannot open as BAM");
        assert!(
            matches!(&missing_error, Error::RustHtslibError(error)
                if matches!(error.as_ref(), HtslibError::BamOpen { target } if target == missing.as_str())),
            "missing file URL must retain BamOpen and its URL target, got {missing_error:?}"
        );

        let unindexed = fixture_url("example_1_copy_no_index.bam");
        let index_error = nanalogue_indexed_bam_reader_from_url(&unindexed, FetchDefinition::All)
            .expect_err("an unindexed file URL cannot create an indexed reader");
        assert!(
            matches!(&index_error, Error::RustHtslibError(error)
                if matches!(error.as_ref(), HtslibError::BamInvalidIndex { target } if target == unindexed.as_str())),
            "missing BAI must retain BamInvalidIndex and its URL target, got {index_error:?}"
        );
    }

    /// Both an out-of-range numeric target and an unknown contig fail during
    /// fetch with `HTSlib`'s precise `Fetch` error, not during URL opening.
    #[test]
    fn indexed_file_url_wraps_invalid_and_unknown_fetches() {
        let url = fixture_url("example_1.bam");
        for definition in [
            FetchDefinition::CompleteTid(99),
            FetchDefinition::String(b"not_a_fixture_contig"),
        ] {
            let error = nanalogue_indexed_bam_reader_from_url(&url, definition)
                .expect_err("an invalid fetch definition must be rejected");
            assert!(
                matches!(&error, Error::RustHtslibError(inner)
                    if matches!(inner.as_ref(), HtslibError::Fetch)),
                "invalid and unknown contig fetches must retain Fetch, got {error:?}"
            );
        }
    }
}
