#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Bounds at the BAM-to-`CurrRead` alignment boundary.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::{CurrRead, Error, ReadState, constants::shared::MAX_CONTIGS};
    use rust_htslib::bam::{
        Record,
        ext::BamRecordExtensions as _,
        record::{Cigar, CigarString},
    };

    /// Ten query bases span fifteen reference bases because of an insertion
    /// and a deletion. No header is needed for the numeric alignment setters.
    fn mapped_record() -> Record {
        let mut record = Record::new();
        record.set(
            b"boundary-read",
            Some(&CigarString(vec![
                Cigar::Match(3),
                Cigar::Ins(2),
                Cigar::Del(7),
                Cigar::Match(5),
            ])),
            b"ACGTACGTAC",
            &[30; 10],
        );
        record.set_flags(0);
        record.set_tid(2);
        record.set_pos(17);
        record.set_mapq(41);
        record
    }

    #[test]
    fn alignment_end_limit_is_inclusive_and_uses_reference_span() {
        let mut record = mapped_record();
        for start in [0, 17, i64::from(u32::MAX) - 15] {
            record.set_pos(start);
            let read = CurrRead::default()
                .set_read_state_and_id(&record)
                .unwrap()
                .set_align_len(&record)
                .unwrap()
                .set_seq_len(&record)
                .unwrap();
            assert_eq!(read.align_len().unwrap(), 15);
            assert_eq!(read.seq_len().unwrap(), 10);
        }

        // The start still fits u32, but the reference end exceeds it by one.
        // Using query length instead of reference span would wrongly pass.
        record.set_pos(i64::from(u32::MAX) - 14);
        assert_eq!(record.reference_end(), i64::from(u32::MAX) + 1);
        let error = CurrRead::default()
            .set_read_state_and_id(&record)
            .unwrap()
            .set_align_len(&record)
            .unwrap_err();
        assert!(matches!(error, Error::InvalidAlignLength(_)));

        record.set_pos(-1);
        assert_eq!(record.reference_end(), 14);
        let negative_start_error = CurrRead::default()
            .set_read_state_and_id(&record)
            .unwrap()
            .set_align_len(&record)
            .unwrap_err();
        assert!(matches!(negative_start_error, Error::InvalidAlignLength(_)));
    }

    #[test]
    fn numeric_contig_ids_exclude_negative_and_limit_values() {
        let mut record = mapped_record();
        let limit = i32::try_from(MAX_CONTIGS).unwrap();
        for tid in [0, limit - 1] {
            record.set_tid(tid);
            let read = CurrRead::default()
                .set_read_state_and_id(&record)
                .unwrap()
                .set_contig_id_and_start(&record)
                .unwrap();
            assert_eq!(read.contig_id_and_start().unwrap(), (tid, 17));
        }
        for (tid, expected) in [
            (-1, "contig id < 0, seems malformed!".to_owned()),
            (
                limit,
                format!("cannot process contigs more than {MAX_CONTIGS}"),
            ),
        ] {
            record.set_tid(tid);
            let error = CurrRead::default()
                .set_read_state_and_id(&record)
                .unwrap()
                .set_contig_id_and_start(&record)
                .unwrap_err();
            assert!(
                matches!(error, Error::InvalidState(message) if message == expected),
                "invalid contig ID {tid} must fail before it is stored"
            );
        }
    }

    #[test]
    fn omitted_unmapped_sequence_requires_explicit_opt_in() {
        let mut record = mapped_record();
        record.set(b"boundary-read", None, b"", &[]);
        record.set_unmapped();
        record.set_tid(-1);
        record.set_pos(-1);
        assert_eq!(record.seq_len(), 0);

        let error = CurrRead::default()
            .try_from_only_alignment(&record)
            .unwrap_err();
        assert!(matches!(error, Error::ZeroSeqLen(_)));

        let read = CurrRead::default()
            .try_from_only_alignment_zero_seq_len(&record)
            .unwrap();
        assert_eq!(read.seq_len().unwrap(), 0);
        assert_eq!(read.read_state(), ReadState::Unmapped);
        assert_eq!(read.read_id(), "boundary-read");
        assert_eq!(read.mapq(), 41);
        // Opting in to an omitted sequence must not invent an alignment from
        // the sentinel coordinates or turn missing properties into zeroes.
        assert!(matches!(read.align_len().unwrap_err(), Error::Unmapped(_)));
        assert!(matches!(
            read.contig_id_and_start().unwrap_err(),
            Error::Unmapped(_)
        ));
    }
}
