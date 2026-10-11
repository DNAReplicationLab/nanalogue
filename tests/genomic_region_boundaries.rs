#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Boundary tests for genomic region parsing and header-aware conversion.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::constants::shared::MAX_CONTIG_LEN;
    use nanalogue_core::{Error, GenomicBed3, GenomicRegion};
    use rust_htslib::bam;

    #[test]
    fn distinguishes_the_region_byte_limit_from_coordinate_validation() {
        // Leading zeroes preserve the numeric value while reaching the exact
        // 255-byte limit without exceeding the separate contig-name limit.
        let at_limit = format!("chr1:{}0-1", "0".repeat(247));
        let over_limit = format!("{at_limit} ");
        assert_eq!(at_limit.len(), 255);
        assert_eq!(over_limit.len(), 256);
        let parsed = at_limit.parse::<GenomicRegion>().expect("255 bytes fit");
        assert_eq!(parsed.contig(), "chr1");
        assert_eq!(parsed.start_end(), Some((0, 1)));
        let error = over_limit.parse::<GenomicRegion>().unwrap_err();
        assert!(
            matches!(&error, Error::InvalidState(message)
                if message == "genomic region is too long!"),
            "256-byte input must fail at the length guard: {error}"
        );
        // The extra space would be trimmed by interval parsing, so rejection
        // is due to length, not a malformed number or an invalid contig.
        let short = "chr1:0-1 ".parse::<GenomicRegion>().unwrap();
        assert_eq!(short.start_end(), Some((0, 1)));
    }

    #[test]
    fn counts_utf8_bytes_before_trimming_interval_whitespace() {
        // EM SPACE takes three bytes but is valid whitespace for the interval
        // parser. These inputs have the same number of Unicode characters.
        let at_limit = format!("chr1:{} 0-1", "\u{2003}".repeat(82));
        let over_limit = format!("chr1:{}\u{2002}0-1", "\u{2003}".repeat(82));
        assert_eq!(at_limit.len(), 255);
        assert_eq!(over_limit.len(), 257);
        assert_eq!(at_limit.chars().count(), over_limit.chars().count());
        assert_eq!(
            at_limit.parse::<GenomicRegion>().unwrap().start_end(),
            Some((0, 1))
        );
        let error = over_limit.parse::<GenomicRegion>().unwrap_err();
        assert!(
            matches!(&error, Error::InvalidState(message)
                if message == "genomic region is too long!"),
            "the limit is bytes, not Unicode characters: {error}"
        );
    }

    #[test]
    fn handles_the_maximum_length_exposed_by_the_header() {
        // Header-only fixtures exercise large coordinates without allocating
        // a multi-gigabyte sequence or depending on an external reference.
        let header = bam::HeaderView::from_bytes(
            b"@SQ\tSN:limit\tLN:4294967294\n@SQ\tSN:too_long\tLN:4294967295\n@SQ\tSN:clamped\tLN:4294967296\n",
        );
        assert_eq!(header.target_len(0), Some(u64::from(MAX_CONTIG_LEN)));
        // The pinned HTSlib clamps larger declared lengths to u32::MAX.
        assert_eq!(header.target_len(1), Some(0xFFFF_FFFF));
        assert_eq!(header.target_len(2), Some(0xFFFF_FFFF));
        for (input, start) in [
            ("limit", 0),
            ("limit:4294967293-", 0xFFFF_FFFD),
            ("limit:4294967293-4294967294", 0xFFFF_FFFD),
        ] {
            let bed = input
                .parse::<GenomicRegion>()
                .unwrap()
                .try_to_bed3(&header)
                .expect("MAX_CONTIG_LEN is a supported contig length");
            assert_eq!(bed, GenomicBed3::new(0, start, MAX_CONTIG_LEN).unwrap());
        }
        for contig in ["too_long", "clamped"] {
            for input in [
                contig.to_owned(),
                format!("{contig}:0-1"),
                format!("{contig}:1-"),
            ] {
                let error = input
                    .parse::<GenomicRegion>()
                    .unwrap()
                    .try_to_bed3(&header)
                    .unwrap_err();
                assert!(
                    matches!(error, Error::InvalidSeqLength(_)),
                    "reject oversized contigs even for short bounded regions"
                );
            }
        }
    }

    #[test]
    fn checks_closed_end_but_resolves_open_end_to_header_length() {
        // Use the second contig to catch accidentally hard-coded numeric IDs.
        let header = bam::HeaderView::from_bytes(b"@SQ\tSN:other\tLN:19\n@SQ\tSN:target\tLN:37\n");
        for (input, end) in [
            ("target:11-36", 36),
            ("target:11-37", 37),
            ("target:11-", 37),
        ] {
            let bed = input
                .parse::<GenomicRegion>()
                .unwrap()
                .try_to_bed3(&header)
                .expect("end is inside the contig or an open-end sentinel");
            assert_eq!(bed, GenomicBed3::new(1, 11, end).unwrap(), "input: {input}");
        }
        for (input, invalid_pos) in [
            ("target:11-38", 38),
            ("target:37-38", 37),
            ("target:37-", 37),
        ] {
            let error = input
                .parse::<GenomicRegion>()
                .unwrap()
                .try_to_bed3(&header)
                .unwrap_err();
            assert!(
                matches!(&error, Error::InvalidRegion { region, pos, contig_length }
                    if region == input && *pos == invalid_pos && *contig_length == 37),
                "report the offending bound; start takes precedence: {error}"
            );
        }
    }

    #[test]
    fn reserves_maximum_coordinate_for_omitted_ends() {
        for input in ["chr16:0-4294967295", "chr16:0- 4294967295 \t"] {
            let error = input.parse::<GenomicRegion>().unwrap_err();
            assert!(
                matches!(error, Error::OrdPairConversion(_)),
                "a numeric end above MAX_CONTIG_LEN must fail during parsing"
            );
        }
        let header = bam::HeaderView::from_bytes(b"@SQ\tSN:chr16\tLN:90338345\n");
        let open = "chr16:11- \t".parse::<GenomicRegion>().unwrap();
        assert_eq!(open.start_end(), Some((11, u32::MAX)));
        assert_eq!(
            serde_json::to_value(&open).unwrap(),
            serde_json::json!(["chr16", {"low": 11, "high": u32::MAX}])
        );
        for region in [
            open,
            GenomicRegion::try_from(("chr16".to_owned(), (11, u32::MAX))).unwrap(),
            serde_json::from_str::<GenomicRegion>(r#"["chr16", {"low": 11, "high": 4294967295}]"#)
                .unwrap(),
        ] {
            assert!(
                matches!(bam::FetchDefinition::from(&region),
                    bam::FetchDefinition::RegionString(name, 11, i64::MAX)
                    if name == b"chr16"),
                "the sentinel must remain an open fetch end"
            );
            assert_eq!(
                region.try_to_bed3(&header).unwrap(),
                GenomicBed3::new(0, 11, 90_338_345).unwrap()
            );
        }
    }
}
