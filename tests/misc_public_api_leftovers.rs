#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Focused public-API tests for the last uncovered behaviour around the
//! `bedrs` compatibility types, the `u32` contig-length boundary of
//! [`GenomicRegion::try_to_bed3`], `sample_fraction` extremes, zero-length
//! sequences during coordinate retrieval, and the identifier-validation
//! message rewrites shared by `peek::run` and `CurrRead::set_read_state_and_id`.
//!
//! Expectations are derived independently from public constants, the fixture
//! contents, and the documented SAM/BAM formats rather than from the
//! implementation under test.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use clap::Parser as _;
    use nanalogue_core::bedrs::{Bed3, Coordinates as _, Strand, StrandedBed3};
    use nanalogue_core::constants::shared::{MAX_CONTIG_NAME_LENGTH, MAX_READ_ID_LEN};
    use nanalogue_core::{
        BamPreFilt as _, CurrRead, Error, F32Bw0and1, GenomicBed3, GenomicRegion, InputBamBuilder,
        commands, nanalogue_bam_reader, peek,
    };
    use rust_htslib::bam::header::HeaderRecord;
    use rust_htslib::bam::record::{Cigar, CigarString};
    use rust_htslib::bam::{Header, HeaderView, Read as _, Record};
    use std::str::FromStr as _;
    use std::sync::Arc;

    /// Fixture holding four alignments across three distinct read IDs.
    const BAM_PATH: &str = "examples/example_1.bam";

    /// The read IDs of `example_1.bam` in sorted order; the first two entries
    /// come from two different alignments that share one ID.
    const EXAMPLE_1_READ_IDS: [&str; 4] = [
        "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575",
        "a4f36092-b4d5-47a9-813e-c22c3b477a0c",
        "a4f36092-b4d5-47a9-813e-c22c3b477a0c",
        "fffffff1-10d2-49cb-8ca3-e8d48979001b",
    ];

    /// The read IDs of `example_1.bam`, sorted so comparisons ignore BAM order.
    fn expected_ids() -> Vec<String> {
        EXAMPLE_1_READ_IDS.map(str::to_owned).to_vec()
    }

    /// Builds a header with a single `chr1` contig of the requested length.
    fn header_with_contig_length(length: u64) -> HeaderView {
        let text = format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{length}\n");
        HeaderView::from_bytes(text.as_bytes())
    }

    /// Converts `region` with the given header, expecting both steps to succeed.
    fn region_to_bed3(header: &HeaderView, region: &str) -> GenomicBed3 {
        GenomicRegion::from_str(region)
            .expect("region parses")
            .try_to_bed3(header)
            .expect("region converts")
    }

    /// Applies `pre_filt` with a `sample_fraction` built through the public
    /// `InputBamBuilder`, returning the surviving read IDs sorted by name.
    fn sampled_ids(fraction: F32Bw0and1) -> Result<Vec<String>, Error> {
        let bam = InputBamBuilder::default()
            .sample_fraction(fraction)
            .build()?;
        let mut reader = nanalogue_bam_reader(BAM_PATH)?;
        let mut ids = Vec::new();
        for record in reader.records() {
            let rec = record?;
            if rec.pre_filt(&bam) {
                ids.push(String::from_utf8_lossy(rec.qname()).into_owned());
            }
        }
        ids.sort();
        Ok(ids)
    }

    /// Runs `read-info` through the public CLI at the given `--sample-fraction`
    /// and returns the reported read IDs, sorted by name.
    fn cli_read_info_ids(fraction: &str) -> Result<Vec<String>, Error> {
        let cli =
            commands::Cli::parse_from(["", "read-info", BAM_PATH, "--sample-fraction", fraction]);
        let mut output = Vec::new();
        commands::run(cli, &mut output)?;
        let entries: Vec<serde_json::Value> =
            serde_json::from_slice(&output).expect("read-info must emit a JSON array");
        let mut ids: Vec<String> = entries
            .iter()
            .map(|entry| {
                entry
                    .get("read_id")
                    .and_then(serde_json::Value::as_str)
                    .expect("every read-info entry must carry a read_id")
                    .to_owned()
            })
            .collect();
        ids.sort();
        Ok(ids)
    }

    /// `Bed3::empty` is exactly the derived default and reports no strand,
    /// while its stranded twin reports the strand it was built with.
    #[test]
    fn empty_bed3_is_default_and_unstranded() {
        let empty: Bed3<String, u32> = Bed3::empty();

        assert_eq!(
            empty,
            Bed3::<String, u32>::default(),
            "empty must be default"
        );
        // Derived independently from the component defaults.
        assert_eq!(
            empty,
            Bed3::new(String::new(), 0u32, 0u32),
            "empty components"
        );
        assert_eq!(empty.chr(), &String::new(), "empty chromosome is blank");
        assert_eq!(empty.start(), 0, "empty start is zero");
        assert_eq!(empty.end(), 0, "empty end is zero");

        let stranded: StrandedBed3<String, u32> =
            StrandedBed3::new(String::from("chr1"), 10, 20, Strand::Reverse);
        assert_eq!(
            empty.strand(),
            None,
            "an unstranded Bed3 must not invent a strand"
        );
        assert_eq!(
            stranded.strand(),
            Some(Strand::Reverse),
            "strand must survive"
        );
        assert_ne!(
            empty.strand(),
            stranded.strand(),
            "the strand is the only difference"
        );
    }

    /// A contig of exactly `u32::MAX` bases is the largest supported length:
    /// the whole contig converts to the exact interval `[0, u32::MAX)`, while
    /// one base less converts to `[0, u32::MAX - 1)`.
    ///
    /// A header declaring `LN = u32::MAX + 1` is clamped by `HTSlib` before
    /// nanalogue sees it: reference lengths are stored in a `u32`, so
    /// `target_len` reports `u32::MAX`. Consequently the
    /// `contig_len > u32::MAX` guard inside `try_to_bed3` (which would return
    /// `Error::InvalidSeqLength`) is unreachable through any header produced
    /// by the public HTSlib-backed API, and the conversion succeeds with the
    /// clamped interval instead of failing.
    #[test]
    fn region_conversion_handles_the_u32_length_boundary() {
        let max_header = header_with_contig_length(u64::from(u32::MAX));
        assert_eq!(
            region_to_bed3(&max_header, "chr1"),
            GenomicBed3::new(0, 0, u32::MAX),
            "the whole contig is the exact BED interval"
        );
        assert_eq!(
            region_to_bed3(&max_header, "chr1:0-"),
            GenomicBed3::new(0, 0, u32::MAX),
            "an open end is the contig length"
        );

        let below_header = header_with_contig_length(u64::from(u32::MAX) - 1);
        assert_eq!(
            region_to_bed3(&below_header, "chr1"),
            GenomicBed3::new(0, 0, u32::MAX - 1),
            "the end is the contig length, not a constant"
        );

        let over_max_header = header_with_contig_length(u64::from(u32::MAX) + 1);
        assert_eq!(
            over_max_header.target_len(0),
            Some(u64::from(u32::MAX)),
            "HTSlib clamps an oversized reference length to u32::MAX"
        );
        assert_eq!(
            region_to_bed3(&over_max_header, "chr1"),
            GenomicBed3::new(0, 0, u32::MAX),
            "the conversion sees the clamped length"
        );
    }

    /// `sample_fraction` extremes retain no records at zero and every record
    /// (with exact IDs) at one, both through the builder and the CLI.
    #[test]
    fn sample_fraction_extremes_retain_none_and_all_records() -> Result<(), Error> {
        assert_eq!(
            sampled_ids(F32Bw0and1::zero())?,
            Vec::<String>::new(),
            "fraction 0 must retain no records"
        );
        assert_eq!(
            sampled_ids(F32Bw0and1::one())?,
            expected_ids(),
            "fraction 1 must retain every record"
        );

        assert_eq!(
            cli_read_info_ids("0")?,
            Vec::<String>::new(),
            "--sample-fraction 0 must retain no records"
        );
        assert_eq!(
            cli_read_info_ids("1")?,
            expected_ids(),
            "--sample-fraction 1 must retain every record"
        );
        Ok(())
    }

    /// A mapped record with an empty `SEQ` field reaches coordinate retrieval
    /// through the lenient constructor and is then rejected with the dedicated
    /// zero-length `InvalidState`; the strict constructor rejects it earlier
    /// with `ZeroSeqLen`.
    #[test]
    fn empty_sequence_is_rejected_by_coordinate_retrieval() -> Result<(), Error> {
        let mut header = Header::new();
        let _header = header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 100),
        );
        let mut record = Record::new();
        record.set_header(Arc::new(HeaderView::from_header(&header)));
        record.set(
            b"empty-seq-read",
            Some(&CigarString(vec![Cigar::Match(10)])),
            b"",
            &[],
        );
        record.set_flags(0);
        record.set_tid(0);
        record.set_pos(10);

        assert_eq!(record.seq_len(), 0, "the record carries no bases");
        let strict_error = CurrRead::default()
            .try_from_only_alignment(&record)
            .expect_err("the strict constructor must refuse an empty sequence");
        assert!(
            matches!(strict_error, Error::ZeroSeqLen(_)),
            "expected ZeroSeqLen, got {strict_error:?}"
        );

        let read = CurrRead::default().try_from_only_alignment_zero_seq_len(&record)?;
        assert_eq!(read.seq_len()?, 0, "the lenient constructor records zero");
        let region = GenomicBed3::new(0, 10, 20);
        let error = read
            .seq_coords_from_ref_coords(&record, &region)
            .expect_err("coordinate retrieval must refuse an empty sequence");
        assert!(
            matches!(&error, Error::InvalidState(message)
                if message == "zero-len sequences cannot be used here even with valid CIGAR strings"),
            "expected the zero-length InvalidState, got {error:?}"
        );
        Ok(())
    }

    /// Overlong contig and read-ID values exercise both identifier-validator
    /// message rewrites: `peek::run` rewrites the contig failure into
    /// `InvalidContig` and `set_read_state_and_id` rewrites the read-ID failure
    /// into `InvalidReadID`, each naming the identifier that triggered it.
    #[test]
    fn overlong_identifiers_report_identifier_specific_errors() {
        let overlong_contig = "a".repeat(usize::from(MAX_CONTIG_NAME_LENGTH) + 1);
        let header_text =
            format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:{overlong_contig}\tLN:100\n");
        let header = HeaderView::from_bytes(header_text.as_bytes());
        let contig_error = peek::run(&mut Vec::new(), &header, std::iter::empty())
            .expect_err("an overlong contig name must stop peek");
        let contig_rendered = format!("{contig_error}");
        let Error::InvalidContig(contig_message) = contig_error else {
            unreachable!("expected InvalidContig, got {contig_error:?}")
        };
        assert_eq!(
            contig_message,
            format!("error in setting contig, length > {MAX_CONTIG_NAME_LENGTH}"),
            "the contig rewrite must name the contig"
        );
        assert!(
            !contig_message.contains("read_id"),
            "the contig rewrite must not name a read id"
        );

        let overlong_read_id = "r".repeat(usize::from(MAX_READ_ID_LEN) + 1);
        let mut record = Record::new();
        record.set_qname(overlong_read_id.as_bytes());
        let read_error = CurrRead::default()
            .set_read_state_and_id(&record)
            .expect_err("an overlong read id must be rejected");
        let read_rendered = format!("{read_error}");
        let Error::InvalidReadID(read_message) = read_error else {
            unreachable!("expected InvalidReadID, got {read_error:?}")
        };
        assert_eq!(
            read_message,
            format!("error in setting read_id, length > {MAX_READ_ID_LEN}"),
            "the read-id rewrite must name the read id"
        );
        assert!(
            !read_message.contains("contig"),
            "the read-id rewrite must not name a contig"
        );

        assert!(
            contig_rendered.starts_with("invalid contig: "),
            "the contig Display rewrite is identifier-specific"
        );
        assert!(
            read_rendered.starts_with("invalid read id: "),
            "the read-id Display rewrite is identifier-specific"
        );
    }
}
