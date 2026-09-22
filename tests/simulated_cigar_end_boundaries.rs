#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Deletion endpoints must retain an aligned base at each end of a simulated read.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::simulate_mod_bam::PerfectSeqMatchToNot;
    use nanalogue_core::{DNARestrictive, Error, OrdPair, ReadState};
    use rand::{SeedableRng as _, rngs::StdRng};

    #[test]
    fn rejects_terminal_deletions_even_when_barcodes_hide_the_end() {
        let mut rng = StdRng::seed_from_u64(6);
        for state in [
            ReadState::PrimaryFwd,
            ReadState::PrimaryRev,
            ReadState::Unmapped,
        ] {
            for barcode in [None, Some("AG".parse::<DNARestrictive>().unwrap())] {
                for bounds in [(0.0, 0.25), (0.5, 1.0)] {
                    let mut builder = PerfectSeqMatchToNot::seq(b"ACGTACGT".to_vec())
                        .unwrap()
                        .delete(OrdPair::try_from(bounds).unwrap());
                    if let Some(value) = barcode.clone() {
                        builder = builder.barcode(value);
                    }
                    let error = builder.build(state, &mut rng).unwrap_err();
                    assert!(
                        matches!(error, Error::SimulateDNASeqCIGAREndProblem(_)),
                        "terminal deletion must fail, including before unmapped CIGAR removal"
                    );
                }
            }
        }
    }

    #[test]
    fn permits_one_surviving_base_after_an_internal_deletion() {
        let mut rng = StdRng::seed_from_u64(6);
        // Delete positions [4,7), leaving the final T as an aligned base.
        // Reverse state changes the barcode orientation, not the input sequence.
        for (state, expected_seq) in [
            (ReadState::PrimaryFwd, b"AGACGTTCT".as_slice()),
            (ReadState::PrimaryRev, b"TCACGTTGA".as_slice()),
        ] {
            let (sequence, cigar) = PerfectSeqMatchToNot::seq(b"ACGTACGT".to_vec())
                .unwrap()
                .delete(OrdPair::try_from((0.5, 0.875)).unwrap())
                .barcode("AG".parse().unwrap())
                .build(state, &mut rng)
                .unwrap();
            assert_eq!(sequence, expected_seq);
            assert_eq!(cigar.unwrap().to_string(), "2S4M3D1M2S");
        }

        let (sequence, cigar) = PerfectSeqMatchToNot::seq(b"ACGTACGT".to_vec())
            .unwrap()
            .delete(OrdPair::try_from((0.5, 0.875)).unwrap())
            .build(ReadState::Unmapped, &mut rng)
            .unwrap();
        assert_eq!(sequence, b"ACGTT");
        assert!(cigar.is_none(), "valid unmapped output has no CIGAR");
    }

    #[test]
    fn rounds_the_deletion_endpoint_before_checking_the_last_base() {
        let mut rng = StdRng::seed_from_u64(6);
        // With eight bases, 0.9375 maps to 7.5 and rounds to 8. Checking
        // only whether the fraction is below 1 would accept an invalid end.
        for (end, expected_sequence, expected_cigar) in [
            (0.875, b"ACGTT".as_slice(), "4M3D1M"),
            (0.9374, b"ACGTT".as_slice(), "4M3D1M"),
        ] {
            let (sequence, cigar) = PerfectSeqMatchToNot::seq(b"ACGTACGT".to_vec())
                .unwrap()
                .delete(OrdPair::try_from((0.5, end)).unwrap())
                .build(ReadState::PrimaryFwd, &mut rng)
                .unwrap();
            assert_eq!(sequence, expected_sequence);
            assert_eq!(cigar.unwrap().to_string(), expected_cigar);
        }
        for end in [0.9375, 1.0] {
            let error = PerfectSeqMatchToNot::seq(b"ACGTACGT".to_vec())
                .unwrap()
                .delete(OrdPair::try_from((0.5, end)).unwrap())
                .build(ReadState::PrimaryFwd, &mut rng)
                .unwrap_err();
            assert!(
                matches!(error, Error::SimulateDNASeqCIGAREndProblem(_)),
                "rounded deletion consumes the last aligned base"
            );
        }
    }

    #[test]
    fn empty_deletion_at_the_endpoint_keeps_the_sequence() {
        let mut rng = StdRng::seed_from_u64(6);
        // Ordered, bounded fractions guarantee start <= end <= sequence length.
        // An empty slice at len is valid; an out-of-bounds deletion slice is
        // not a constructible state for this small input and needs no injection.
        for bounds in [(0.0, 0.0), (1.0, 1.0)] {
            let (sequence, cigar) = PerfectSeqMatchToNot::seq(b"ACGTACGT".to_vec())
                .unwrap()
                .delete(OrdPair::try_from(bounds).unwrap())
                .build(ReadState::PrimaryFwd, &mut rng)
                .unwrap();
            assert_eq!(sequence, b"ACGTACGT");
            assert_eq!(cigar.unwrap().to_string(), "8M");
        }
    }
}
