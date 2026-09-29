//! Sparse sequence-to-reference coordinate mapping from BAM CIGAR operations.

use crate::Error;
use crate::constants::shared::{MAX_CIGAR_OPERATIONS, MAX_CONTIG_LEN, MAX_SEQ_LEN};

/// Reference coordinates of a record's sequence positions, as used by
/// [`crate::nanalogue_mm_ml_parser`].
///
/// Long reads with indels would otherwise need one entry per base. For well-formed
/// CIGARs the map instead stores one entry per query-consuming operation and looks
/// positions up by binary search.
#[derive(Debug)]
pub(crate) struct SeqToRefMap(Vec<SeqToRefSegment>);

/// A run of sequence positions produced by one query-consuming CIGAR operation.
#[derive(Debug, Clone, Copy)]
struct SeqToRefSegment {
    /// First sequence position in the run.
    query_start: u32,
    /// One past the last sequence position in the run.
    query_end: u32,
    /// Reference position of `query_start` for aligned runs; `None` for insertions and clips.
    ref_start: Option<u32>,
}

impl SeqToRefMap {
    /// Builds an empty map for an unmapped record.
    pub(crate) fn unmapped() -> Self {
        Self(Vec::new())
    }

    /// Returns the reference coordinate of a sequence position that is below the sequence length.
    pub(crate) fn get(&self, position: u32) -> Option<u32> {
        let index = self
            .0
            .partition_point(|segment| segment.query_end <= position);
        let segment = self.0.get(index)?;
        let ref_start = segment.ref_start?;
        let offset = position.checked_sub(segment.query_start)?;
        // Validated when the map was built: every aligned position fits in `u32`.
        ref_start.checked_add(offset)
    }

    /// Builds a segment map from a mapped record's CIGAR.
    ///
    /// Each aligned reference coordinate must fit in `u32`, and the CIGAR must consume
    /// exactly `seq_len` query bases. Padding and hard clipping consume neither coordinate.
    /// Zero-length and unknown operations are rejected instead of relying on dependency
    /// behavior that may panic.
    #[expect(
        clippy::too_many_lines,
        reason = "keeping CIGAR validation and segment construction together makes semantics auditable"
    )]
    pub(crate) fn try_from_raw_cigar(
        raw_cigar: &[u32],
        reference_start: u32,
        seq_len: u32,
    ) -> Result<Self, Error> {
        /// Bits used by the operation code in a BAM CIGAR word.
        const OP_MASK: u32 = 0xf;
        /// Bits the operation length is shifted by in a BAM CIGAR word.
        const LEN_SHIFT: u32 = 4;
        if seq_len > MAX_SEQ_LEN {
            return Err(Error::InvalidSeqLength(format!(
                "sequence length exceeds {MAX_SEQ_LEN}"
            )));
        }
        if reference_start >= MAX_CONTIG_LEN {
            return Err(Error::InvalidModCoords(format!(
                "reference coordinate exceeds maximum contig length {MAX_CONTIG_LEN}"
            )));
        }
        if raw_cigar.is_empty() {
            return Err(Error::InvalidState(
                "mapped record has an empty CIGAR".to_owned(),
            ));
        }
        if raw_cigar.len()
            > usize::try_from(MAX_CIGAR_OPERATIONS).expect("u32 fits in supported usize")
        {
            return Err(Error::InvalidState(format!(
                "max CIGAR operations exceeded {MAX_CIGAR_OPERATIONS}"
            )));
        }
        let mut segments = Vec::with_capacity(
            raw_cigar
                .len()
                .min(usize::try_from(seq_len).expect("u32 fits in supported usize")),
        );
        let mut query = 0u32;
        let mut reference_coordinate = reference_start;
        for &word in raw_cigar {
            let len = word >> LEN_SHIFT;
            if len == 0 {
                return Err(Error::InvalidState(
                    "CIGAR operations must have non-zero length".to_owned(),
                ));
            }
            let query_end = query.saturating_add(len);
            let operation = word & OP_MASK;
            match operation {
                // M, =, X: aligned bases consume query and reference.
                0 | 7 | 8 => {
                    if query_end > seq_len {
                        return Err(Error::InvalidState(
                            "rust_htslib failure! seq coordinates malformed".to_owned(),
                        ));
                    }
                    let last = reference_coordinate.saturating_add(len.saturating_sub(1));
                    if last >= MAX_CONTIG_LEN {
                        return Err(Error::InvalidModCoords(format!(
                            "reference coordinate exceeds maximum contig length {MAX_CONTIG_LEN}"
                        )));
                    }
                    assert!(query_end > query, "non-zero CIGAR operation advances query");
                    segments.push(SeqToRefSegment {
                        query_start: query,
                        query_end,
                        ref_start: Some(reference_coordinate),
                    });
                    query = query_end;
                    reference_coordinate = reference_coordinate.saturating_add(len);
                }
                // I, S: query bases without reference coordinates.
                1 | 4 => {
                    if query_end > seq_len {
                        return Err(Error::InvalidState(
                            "rust_htslib failure! seq coordinates malformed".to_owned(),
                        ));
                    }
                    assert!(query_end > query, "non-zero CIGAR operation advances query");
                    segments.push(SeqToRefSegment {
                        query_start: query,
                        query_end,
                        ref_start: None,
                    });
                    query = query_end;
                }
                // D, N: reference only.
                2 | 3 => {
                    reference_coordinate = reference_coordinate.saturating_add(len);
                }
                // H, P: consume neither query nor reference.
                5 | 6 => {}
                _ => {
                    return Err(Error::InvalidState(format!(
                        "unsupported CIGAR operation code {operation}"
                    )));
                }
            }
        }
        if query == seq_len {
            if reference_coordinate > MAX_CONTIG_LEN {
                return Err(Error::InvalidModCoords(format!(
                    "reference alignment exceeds maximum contig length {MAX_CONTIG_LEN}"
                )));
            }
            assert_eq!(
                segments
                    .iter()
                    .map(|segment| {
                        segment
                            .query_end
                            .checked_sub(segment.query_start)
                            .expect("segment end follows segment start")
                    })
                    .sum::<u32>(),
                seq_len,
                "segment query lengths must sum to the sequence length"
            );
            Ok(Self(segments))
        } else {
            Err(Error::InvalidState(
                "rust_htslib failure! seq coordinates malformed".to_owned(),
            ))
        }
    }
}

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod seq_to_ref_map_tests {
    use super::*;
    use crate::constants::shared::MAX_RECORD_CAPACITY_BYTES;
    use crate::{BaseMods, nanalogue_mm_ml_parser};
    use rust_htslib::bam;
    use rust_htslib::bam::ext::BamRecordExtensions as _;
    use rust_htslib::bam::record::{Aux, Cigar, CigarString};

    /// Builds a mapped record with the given start and CIGAR.
    fn record(pos: i64, operations: Vec<Cigar>) -> bam::Record {
        let cigar = CigarString::from(operations);
        let query_len = cigar
            .iter()
            .map(|op| match *op {
                Cigar::Match(len)
                | Cigar::Ins(len)
                | Cigar::SoftClip(len)
                | Cigar::Equal(len)
                | Cigar::Diff(len) => usize::try_from(len).expect("fits"),
                Cigar::Del(_) | Cigar::RefSkip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => 0,
            })
            .sum::<usize>();
        let mut record = bam::Record::new();
        record.set_tid(0);
        record.unset_unmapped();
        record.set(
            b"read",
            Some(&cigar),
            &vec![b'A'; query_len],
            &vec![30; query_len],
        );
        record.set_pos(pos);
        record
    }

    /// Collects reference coordinates the way the parser did before segment maps.
    fn dense(record: &bam::Record, seq_len: usize) -> Result<Vec<Option<u32>>, String> {
        let temp = record
            .aligned_pairs_full()
            .filter(|x| x[0].is_some())
            .take(seq_len.saturating_add(1))
            .map(|x| match x[1] {
                None => Ok(None),
                Some(v) => u32::try_from(v).map(Some).map_err(|e| {
                    Error::InvalidModCoords(format!(
                        "reference coordinate from aligned_pairs_full is invalid: {e}"
                    ))
                    .to_string()
                }),
            })
            .collect::<Result<Vec<_>, _>>()?;
        if temp.len() == seq_len {
            Ok(temp)
        } else {
            Err(
                Error::InvalidState("rust_htslib failure! seq coordinates malformed".to_owned())
                    .to_string(),
            )
        }
    }

    /// Asserts that the segment map agrees with the dense collection for every position.
    fn assert_matches_dense(record: &bam::Record, seq_len: usize) {
        let expected = dense(record, seq_len);
        let seq_len_u32 = u32::try_from(seq_len).expect("test sequence length fits u32");
        let reference_start = u32::try_from(record.pos()).map_err(|error| {
            Error::InvalidModCoords(format!("reference start coordinate is invalid: {error}"))
                .to_string()
        });
        let actual = reference_start
            .and_then(|validated_reference_start| {
                SeqToRefMap::try_from_raw_cigar(
                    record.raw_cigar(),
                    validated_reference_start,
                    seq_len_u32,
                )
                .map_err(|error| error.to_string())
            })
            .map(|map| {
                assert!(map.0.len() <= usize::try_from(seq_len_u32).expect("u32 fits in usize"));
                assert_eq!(
                    map.0
                        .iter()
                        .map(|segment| {
                            segment
                                .query_end
                                .checked_sub(segment.query_start)
                                .expect("segment end follows segment start")
                        })
                        .sum::<u32>(),
                    seq_len_u32
                );
                (0..seq_len_u32)
                    .map(|index| map.get(index))
                    .collect::<Vec<_>>()
            });
        assert_eq!(actual, expected);
    }

    /// Adds calls at every sequence position and parses without filtering.
    fn parse_all_positions(record: &mut bam::Record) -> Result<BaseMods, Error> {
        let distances = std::iter::repeat_n("0", record.seq_len())
            .collect::<Vec<_>>()
            .join(",");
        record.push_aux(b"MM", Aux::String(&format!("N+n?,{distances};")))?;
        let probabilities = vec![200; record.seq_len()];
        record.push_aux(b"ML", Aux::ArrayU8((&probabilities).into()))?;
        nanalogue_mm_ml_parser(record, |&_| true, |&_| true, |&_, &_, &_| true, 0)
    }

    /// Replaces the first raw CIGAR operation code while retaining its length.
    #[expect(
        clippy::cast_ptr_alignment,
        reason = "BAM storage guarantees that raw CIGAR words are aligned to four bytes"
    )]
    fn replace_first_operation_code(record: &mut bam::Record, operation: u32) {
        const OP_MASK: u32 = 0xf;
        let qname_capacity = usize::from(record.inner().core.l_qname);
        let data = record.inner_mut().data;
        // SAFETY: `Record::set` allocated at least one aligned CIGAR word after the query name.
        let word = unsafe { data.add(qname_capacity) }.cast::<u32>();
        // SAFETY: the pointer addresses the initialized first CIGAR word.
        let current = unsafe { word.read() };
        // SAFETY: the record is exclusively borrowed, and the write retains the word's length.
        unsafe {
            word.write((current & !OP_MASK) | operation);
        }
    }

    /// Encodes a raw BAM CIGAR word.
    fn word(len: u32, operation: u32) -> u32 {
        (len << 4) | operation
    }

    /// Deterministic xorshift generator for reproducible differential tests.
    struct XorShift(u64);

    impl XorShift {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }

        fn below(&mut self, upper_bound: u64) -> u64 {
            self.next()
                .checked_rem(upper_bound)
                .expect("upper bound is non-zero")
        }
    }

    #[test]
    fn segments_match_aligned_pairs_for_every_operation() {
        let cigar = vec![
            Cigar::HardClip(4),
            Cigar::SoftClip(3),
            Cigar::Match(4),
            Cigar::Ins(2),
            Cigar::Del(3),
            Cigar::RefSkip(5),
            Cigar::Equal(2),
            Cigar::Diff(1),
            Cigar::Ins(1),
            Cigar::Match(3),
            Cigar::SoftClip(2),
            Cigar::HardClip(1),
        ];
        let record = record(17, cigar);
        assert_matches_dense(&record, record.seq_len());
    }

    #[test]
    fn segments_match_aligned_pairs_errors() {
        let near_end = i64::from(MAX_CONTIG_LEN) - 3;
        let at_limit = record(near_end, vec![Cigar::SoftClip(2), Cigar::Match(3)]);
        assert_matches_dense(&at_limit, at_limit.seq_len());

        let spliced = record(5, vec![Cigar::Match(3), Cigar::Del(2), Cigar::Match(3)]);
        let seq_len = spliced.seq_len();
        assert_matches_dense(&spliced, seq_len.saturating_sub(1));
        assert_matches_dense(&spliced, seq_len.saturating_add(1));
    }

    #[test]
    fn contiguous_start_above_u32_max_retains_its_error() {
        let mut invalid = record(i64::from(u32::MAX) + 1, vec![Cigar::Match(1)]);
        let error = parse_all_positions(&mut invalid)
            .expect_err("a contiguous alignment start must fit in u32");
        assert!(matches!(error, Error::InvalidModCoords(message)
            if message.starts_with("reference start coordinate is invalid:")));
    }

    #[test]
    fn reference_start_at_or_above_contig_limit_is_rejected() {
        for start in [MAX_CONTIG_LEN, u32::MAX] {
            let mut invalid = record(i64::from(start), vec![Cigar::Match(1)]);
            let error = parse_all_positions(&mut invalid)
                .expect_err("the reference start must be below the contig limit");
            assert!(matches!(error, Error::InvalidModCoords(message)
            if message == format!(
                "reference coordinate exceeds maximum contig length {MAX_CONTIG_LEN}"
            )));
        }
    }

    #[test]
    fn mapped_cigar_without_aligned_query_bases_has_no_reference_positions() {
        for operation in [Cigar::SoftClip(1), Cigar::Ins(1)] {
            let mut record = record(17, vec![operation]);
            let parsed = parse_all_positions(&mut record)
                .expect("the parser retains the previous permissive behavior");
            let annotation = parsed
                .base_mods
                .first()
                .expect("one modification group")
                .ranges
                .annotations
                .first()
                .expect("one annotation");
            assert_eq!(annotation.ref_pos, None);
        }
    }

    #[test]
    fn unaligned_cigar_still_requires_start_below_contig_limit() {
        // Only the start check catches this: with no M/=/X/D/N operation, the final
        // half-open end stays equal to the start, which is not above `MAX_CONTIG_LEN`.
        let mut invalid = record(i64::from(MAX_CONTIG_LEN), vec![Cigar::SoftClip(1)]);
        let error = parse_all_positions(&mut invalid)
            .expect_err("an unaligned record cannot start at the contig limit");
        assert!(matches!(error, Error::InvalidModCoords(message)
        if message == format!(
            "reference coordinate exceeds maximum contig length {MAX_CONTIG_LEN}"
        )));
    }

    #[test]
    fn sequence_length_limit_returns_a_controlled_error() {
        let record = record(17, vec![Cigar::Match(1)]);
        let error = SeqToRefMap::try_from_raw_cigar(record.raw_cigar(), 17, MAX_SEQ_LEN + 1)
            .expect_err("the sequence length limit must be enforced");
        assert!(matches!(error, Error::InvalidSeqLength(message)
            if message == format!("sequence length exceeds {MAX_SEQ_LEN}")));
    }

    #[test]
    fn excessive_query_consumption_returns_a_controlled_error() {
        let mut operations = vec![Cigar::Match(1)];
        operations.extend(std::iter::repeat_n(Cigar::Ins(MAX_SEQ_LEN), 18));
        let mut malformed = bam::Record::new();
        malformed.set(b"read", Some(&CigarString::from(operations)), b"A", &[30]);
        malformed.set_tid(0);
        malformed.unset_unmapped();
        malformed.set_pos(17);

        let error = parse_all_positions(&mut malformed)
            .expect_err("excessive CIGAR query consumption must be rejected");
        assert!(matches!(error, Error::InvalidState(message)
            if message == "rust_htslib failure! seq coordinates malformed"));
    }

    #[test]
    fn cigar_operation_limit_is_a_record_capacity_fail_safe() {
        assert_eq!(MAX_CIGAR_OPERATIONS, MAX_RECORD_CAPACITY_BYTES);
    }

    #[test]
    fn padding_preserves_forward_and_reverse_reference_coordinates() -> Result<(), Error> {
        for flags in [0, 16] {
            let mut padded = record(17, vec![Cigar::Match(1), Cigar::Pad(4), Cigar::Match(2)]);
            padded.set_flags(flags);
            let parsed = parse_all_positions(&mut padded)?;
            let annotations = &parsed
                .base_mods
                .first()
                .expect("one N+n group")
                .ranges
                .annotations;
            assert_eq!(
                annotations
                    .iter()
                    .map(|annotation| annotation.ref_pos)
                    .collect::<Vec<_>>(),
                [Some(17), Some(18), Some(19)]
            );
        }
        Ok(())
    }

    #[test]
    fn maximum_sequence_length_maps_first_and_last_base() {
        let map = SeqToRefMap::try_from_raw_cigar(&[word(MAX_SEQ_LEN, 0)], 0, MAX_SEQ_LEN)
            .expect("the maximum supported sequence length is valid");
        assert_eq!(map.get(0), Some(0));
        assert_eq!(map.get(MAX_SEQ_LEN - 1), Some(MAX_SEQ_LEN - 1));
        assert_eq!(map.get(MAX_SEQ_LEN), None);
    }

    #[test]
    fn lookups_are_exact_at_segment_boundaries() {
        // 2S 3M 2I 3M => query 0..2 clip, 2..5 aligned, 5..7 insertion, 7..10 aligned.
        let raw = [word(2, 4), word(3, 0), word(2, 1), word(3, 0)];
        let map = SeqToRefMap::try_from_raw_cigar(&raw, 100, 10).expect("valid CIGAR");
        let positions: Vec<_> = (0..12).map(|index| map.get(index)).collect();
        assert_eq!(
            positions,
            [
                None,
                None,
                Some(100),
                Some(101),
                Some(102),
                None,
                None,
                Some(103),
                Some(104),
                Some(105),
                None,
                None,
            ]
        );
    }

    #[test]
    fn padding_and_hard_clips_do_not_shift_positions() {
        let raw = [word(2, 0), word(3, 6), word(1, 5), word(2, 0)];
        let map = SeqToRefMap::try_from_raw_cigar(&raw, 7, 4).expect("valid CIGAR");
        let positions: Vec<_> = (0..4).map(|index| map.get(index)).collect();
        assert_eq!(positions, [Some(7), Some(8), Some(9), Some(10)]);
    }

    #[test]
    fn reference_coordinate_saturation_does_not_wrap() {
        let mut valid: Vec<u32> = std::iter::repeat_n(word(MAX_SEQ_LEN, 2), 16).collect();
        valid.push(word(1, 0));
        let map = SeqToRefMap::try_from_raw_cigar(&valid, 0, 1)
            .expect("the aligned base remains below the contig limit");
        assert_eq!(map.get(0), Some(16 * MAX_SEQ_LEN));

        let mut overflowing: Vec<u32> = std::iter::repeat_n(word(MAX_SEQ_LEN, 2), 17).collect();
        overflowing.push(word(1, 0));
        let error = SeqToRefMap::try_from_raw_cigar(&overflowing, 0, 1)
            .expect_err("a saturated reference coordinate must be rejected");
        assert!(matches!(error, Error::InvalidModCoords(_)));
    }

    #[test]
    fn trailing_reference_consumption_obeys_contig_limit() {
        let start = MAX_CONTIG_LEN - 3;
        for (deleted, expected_ok) in [(2, true), (3, false)] {
            let raw = [word(1, 0), word(deleted, 2)];
            let result = SeqToRefMap::try_from_raw_cigar(&raw, start, 1);
            assert_eq!(result.is_ok(), expected_ok, "trailing deletion {deleted}");
        }
    }

    #[test]
    fn reverse_reads_leave_clipped_and_inserted_bases_unmapped() {
        let mut read = record(
            17,
            vec![
                Cigar::SoftClip(2),
                Cigar::Match(3),
                Cigar::Ins(1),
                Cigar::Match(1),
            ],
        );
        read.set_flags(16);
        let parsed = parse_all_positions(&mut read).expect("valid record");
        let mut positions: Vec<(u32, Option<u32>)> = parsed
            .base_mods
            .first()
            .expect("one modification group")
            .ranges
            .annotations
            .iter()
            .map(|annotation| (annotation.pos, annotation.ref_pos))
            .collect();
        positions.sort_unstable();
        assert_eq!(
            positions,
            [
                (0, None),
                (1, None),
                (2, Some(17)),
                (3, Some(18)),
                (4, Some(19)),
                (5, None),
                (6, Some(20)),
            ]
        );
    }

    #[test]
    fn random_cigars_match_htslib_aligned_pairs() {
        let mut rng = XorShift(0x9E37_79B9_7F4A_7C15);
        for _ in 0..5_000 {
            let operation_count = usize::try_from(1 + rng.below(8)).expect("small value");
            let mut operations: Vec<Cigar> = std::iter::repeat_with(|| {
                let len = u32::try_from(1 + rng.below(6)).expect("small value");
                match rng.below(9) {
                    0 => Cigar::Match(len),
                    1 => Cigar::Ins(len),
                    2 => Cigar::Del(len),
                    3 => Cigar::RefSkip(len),
                    4 => Cigar::SoftClip(len),
                    5 => Cigar::HardClip(len),
                    6 => Cigar::Pad(len),
                    7 => Cigar::Equal(len),
                    _ => Cigar::Diff(len),
                }
            })
            .take(operation_count)
            .collect();
            operations.push(Cigar::Match(1));

            let start = if rng.below(2) == 0 {
                i64::try_from(rng.below(1_000)).expect("small value")
            } else {
                i64::from(MAX_CONTIG_LEN) - i64::try_from(rng.below(12)).expect("small value")
            };
            let record_with_padding = record(start, operations.clone());
            // rust-htslib's iterator panics on padding, whose correct effect is no coordinate
            // consumption, so compare against an otherwise identical pad-free record.
            let oracle_record = record(
                start,
                operations
                    .iter()
                    .filter(|operation| !matches!(operation, Cigar::Pad(_)))
                    .copied()
                    .collect(),
            );
            let seq_len_usize = record_with_padding.seq_len();
            let expected = dense(&oracle_record, seq_len_usize);
            let within_limit = expected.as_ref().is_ok_and(|coordinates| {
                coordinates
                    .iter()
                    .flatten()
                    .all(|coordinate| *coordinate < MAX_CONTIG_LEN)
            }) && oracle_record.reference_end() <= i64::from(MAX_CONTIG_LEN);
            let seq_len_u32 = u32::try_from(seq_len_usize).expect("small value");
            let actual = SeqToRefMap::try_from_raw_cigar(
                record_with_padding.raw_cigar(),
                u32::try_from(start).expect("start fits u32"),
                seq_len_u32,
            );

            match (within_limit, actual) {
                (true, Ok(map)) => {
                    let positions: Vec<_> = (0..seq_len_u32).map(|index| map.get(index)).collect();
                    assert_eq!(Ok(positions), expected, "{operations:?} at {start}");
                }
                (false, Err(_)) => {}
                (expected_ok, result) => assert_eq!(
                    expected_ok,
                    result.is_ok(),
                    "{operations:?} at {start}: expected ok={expected_ok}, got {result:?}"
                ),
            }
        }
    }

    #[test]
    fn zero_length_operations_return_controlled_errors() {
        let operations = [
            Cigar::Match(0),
            Cigar::Ins(0),
            Cigar::Del(0),
            Cigar::RefSkip(0),
            Cigar::SoftClip(0),
            Cigar::HardClip(0),
            Cigar::Pad(0),
            Cigar::Equal(0),
            Cigar::Diff(0),
        ];
        for flags in [0, 16] {
            for operation in operations {
                let mut malformed = record(17, vec![Cigar::Match(3), operation]);
                malformed.set_flags(flags);
                let error = parse_all_positions(&mut malformed)
                    .expect_err("zero-length CIGAR operations must be rejected");
                assert!(
                    matches!(&error, Error::InvalidState(message)
                        if message == "CIGAR operations must have non-zero length"),
                    "unexpected error for {operation:?} on flags {flags}: {error}"
                );
            }
        }
    }

    #[test]
    fn mapped_empty_cigar_returns_a_controlled_error() {
        let mut malformed = bam::Record::new();
        malformed.set(
            b"read",
            Some(&CigarString::from(Vec::<Cigar>::new())),
            b"A",
            &[30],
        );
        malformed.set_tid(0);
        malformed.unset_unmapped();
        malformed.set_pos(17);
        let error = parse_all_positions(&mut malformed)
            .expect_err("a mapped record must describe its alignment");
        assert!(matches!(error, Error::InvalidState(message)
            if message == "mapped record has an empty CIGAR"));
    }

    #[test]
    fn every_unknown_raw_operation_returns_a_controlled_error() {
        for operation in 9..=15 {
            let mut malformed = record(17, vec![Cigar::Match(3)]);
            replace_first_operation_code(&mut malformed, operation);
            let error = parse_all_positions(&mut malformed)
                .expect_err("unknown raw CIGAR operations must be rejected");
            assert!(matches!(error, Error::InvalidState(message)
                if message == format!("unsupported CIGAR operation code {operation}")));
        }
    }
}
