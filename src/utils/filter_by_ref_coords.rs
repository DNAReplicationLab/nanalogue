//! `FilterByRefCoords` trait for filtering by coordinates on the reference genome
//! Provides interface for coordinate-based filtering operations

use fibertools_rs::utils::bamranges::Ranges;

/// Implements filter by coordinates on the reference genome.
pub trait FilterByRefCoords {
    /// filters by reference position i.e. all pos such that start <= pos < end
    /// are retained. does not use contig in filtering.
    fn filter_by_ref_pos(&mut self, _: i64, _: i64) {
        unimplemented!()
    }
}

/// Implements filter by reference coordinates for the Ranges
/// struct that contains our modification information.
/// NOTE: Ranges does not contain contig information, so we cannot
/// filter by that here.
impl FilterByRefCoords for Ranges {
    /// filters by reference position i.e. all pos such that start <= pos < end
    /// are retained. does not use contig in filtering.
    fn filter_by_ref_pos(&mut self, start: i64, end: i64) {
        let (start_index, stop_index) = {
            let mut last_invalid_window: Option<usize> = None;
            let mut last_valid_window: Option<usize> = None;
            let mut previous_start: Option<i64> = None;
            let mut previous_end: Option<i64> = None;
            let mut is_possible_invalid: bool;
            for k in self
                .reference_starts
                .iter()
                .zip(self.reference_ends.iter())
                .enumerate()
            {
                // ensure start <= end for each interval and intervals are sorted.
                assert!(
                    ((*k.1.0).is_none() || (*k.1.1).is_none()) || *(k.1.0) <= *(k.1.1),
                    "start {:?} <= end {:?} expected",
                    *(k.1.0),
                    *(k.1.1)
                );
                assert!(
                    (*k.1.0).is_none() || (*k.1.0) > previous_start,
                    "start {:?} > previous_start {:?} expected",
                    *(k.1.0),
                    previous_start
                );
                assert!(
                    (*k.1.1).is_none() || (*k.1.1) > previous_start,
                    "end {:?} > previous_start {:?} expected",
                    *(k.1.1),
                    previous_start
                );
                assert!(
                    (*k.1.0).is_none() || (*k.1.0) >= previous_end,
                    "start {:?} >= previous_end {:?} expected",
                    *(k.1.0),
                    previous_end
                );
                assert!(
                    (*k.1.1).is_none() || (*k.1.1) >= previous_end,
                    "end {:?} >= previous_end {:?} expected",
                    *(k.1.1),
                    previous_end
                );

                if let Some(v) = *(k.1.0) {
                    previous_start = Some(v);
                    if v < end {
                        last_valid_window = Some(k.0);
                    }
                    is_possible_invalid = v < start;
                } else {
                    is_possible_invalid = false;
                }

                if let Some(v) = *(k.1.1) {
                    previous_end = Some(v);
                    if is_possible_invalid && v <= start {
                        last_invalid_window = Some(k.0);
                    }
                }
            }
            (
                last_invalid_window.map_or(0, |x| x.checked_add(1).expect("overflow error")),
                last_valid_window.map_or(0, |x| x.checked_add(1).expect("overflow error")),
            )
        };

        for k in [
            &mut self.starts,
            &mut self.ends,
            &mut self.lengths,
            &mut self.reference_starts,
            &mut self.reference_ends,
            &mut self.reference_lengths,
        ] {
            k.truncate(stop_index);
            k.drain(0..start_index).for_each(drop);
        }
        self.qual.truncate(stop_index);
        self.qual.drain(0..start_index).for_each(drop);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn direct_ranges_filter_by_ref_pos() {
        // Create a Ranges object with multiple ranges
        // All vectors have the same length as required
        let mut ranges = Ranges {
            // Each entry represents a genomic range with different properties
            starts: vec![Some(10), Some(20), Some(30), Some(40)],
            ends: vec![Some(15), Some(25), Some(35), Some(45)],
            lengths: vec![Some(5), Some(5), Some(5), Some(5)],
            qual: vec![100, 120, 140, 160],
            reference_starts: vec![Some(10), Some(20), Some(30), Some(40)],
            reference_ends: vec![Some(15), Some(25), Some(35), Some(45)],
            reference_lengths: vec![Some(5), Some(5), Some(5), Some(5)],
            seq_len: 50,
            reverse: false,
        };

        // Filter ranges to keep only those overlapping with reference region 18-32
        ranges.filter_by_ref_pos(18, 32);

        // Verify that only the ranges that overlap with [18, 32) are kept
        // Should keep ranges [20-25] and [30-35] (indexes 1 and 2)
        assert_eq!(ranges.starts.len(), 2);
        assert_eq!(ranges.ends.len(), 2);
        assert_eq!(ranges.lengths.len(), 2);
        assert_eq!(ranges.qual.len(), 2);
        assert_eq!(ranges.reference_starts.len(), 2);
        assert_eq!(ranges.reference_ends.len(), 2);
        assert_eq!(ranges.reference_lengths.len(), 2);

        // Check the specific values were retained correctly
        assert_eq!(ranges.starts, vec![Some(20), Some(30)]);
        assert_eq!(ranges.ends, vec![Some(25), Some(35)]);
        assert_eq!(ranges.reference_starts, vec![Some(20), Some(30)]);
        assert_eq!(ranges.reference_ends, vec![Some(25), Some(35)]);
        assert_eq!(ranges.qual, vec![120, 140]);

        // Verify seq_len and reverse flag are preserved
        assert_eq!(ranges.seq_len, 50);
        assert!(!ranges.reverse);
    }

    #[test]
    fn ranges_filter_by_ref_pos_no_overlap() {
        // Create a Ranges object with ranges that don't overlap the target region
        let mut ranges = Ranges {
            starts: vec![Some(10), Some(20), Some(60), Some(70)],
            ends: vec![Some(15), Some(25), Some(65), Some(75)],
            lengths: vec![Some(5), Some(5), Some(5), Some(5)],
            qual: vec![100, 120, 140, 160],
            reference_starts: vec![Some(10), Some(20), Some(60), Some(70)],
            reference_ends: vec![Some(15), Some(25), Some(65), Some(75)],
            reference_lengths: vec![Some(5), Some(5), Some(5), Some(5)],
            seq_len: 80,
            reverse: true,
        };

        // Filter ranges for region [30, 50) - none should match
        ranges.filter_by_ref_pos(30, 50);

        // Verify that all ranges were filtered out
        assert_eq!(ranges.starts.len(), 0);
        assert_eq!(ranges.ends.len(), 0);
        assert_eq!(ranges.lengths.len(), 0);
        assert_eq!(ranges.qual.len(), 0);
        assert_eq!(ranges.reference_starts.len(), 0);
        assert_eq!(ranges.reference_ends.len(), 0);
        assert_eq!(ranges.reference_lengths.len(), 0);

        // Verify metadata is preserved
        assert_eq!(ranges.seq_len, 80);
        assert!(ranges.reverse);
    }

    #[test]
    fn ranges_filter_by_ref_pos_with_none_values() {
        // Create a Ranges object with some None values
        let mut ranges = Ranges {
            starts: vec![Some(10), Some(20), None, Some(40)],
            ends: vec![Some(15), Some(25), None, Some(45)],
            lengths: vec![Some(5), Some(5), None, Some(5)],
            qual: vec![100, 120, 140, 160],
            reference_starts: vec![Some(10), Some(20), None, Some(40)],
            reference_ends: vec![Some(15), Some(25), None, Some(45)],
            reference_lengths: vec![Some(5), Some(5), None, Some(5)],
            seq_len: 50,
            reverse: false,
        };

        // Filter ranges to keep only those overlapping with reference region 18-22
        // Only the range at index 1 should be kept
        ranges.filter_by_ref_pos(18, 22);

        // Verify that only the range that overlaps with [18, 22) is kept
        assert_eq!(ranges.starts.len(), 1);
        assert_eq!(ranges.reference_starts, vec![Some(20)]);
        assert_eq!(ranges.reference_ends, vec![Some(25)]);
        assert_eq!(ranges.qual, vec![120]);

        // Verify metadata is preserved
        assert_eq!(ranges.seq_len, 50);
        assert!(!ranges.reverse);
    }

    #[test]
    fn ranges_filter_by_ref_pos_with_none_values_2() {
        // Create a Ranges object with some None values
        let mut ranges = Ranges {
            starts: vec![Some(10), Some(20), None, Some(21), Some(40)],
            ends: vec![Some(15), Some(21), None, Some(22), Some(45)],
            lengths: vec![Some(5), Some(1), None, Some(1), Some(5)],
            qual: vec![100, 120, 140, 150, 160],
            reference_starts: vec![Some(10), Some(20), None, Some(21), Some(40)],
            reference_ends: vec![Some(15), Some(21), None, Some(22), Some(45)],
            reference_lengths: vec![Some(5), Some(1), None, Some(1), Some(5)],
            seq_len: 50,
            reverse: true,
        };

        // Filter ranges to keep only those overlapping with reference region 18-22
        ranges.filter_by_ref_pos(18, 22);

        // Verify that only the ranges that overlaps with [18, 22) are kept
        assert_eq!(ranges.starts.len(), 3);
        assert_eq!(ranges.reference_starts, vec![Some(20), None, Some(21)]);
        assert_eq!(ranges.reference_ends, vec![Some(21), None, Some(22)]);
        assert_eq!(ranges.qual, vec![120, 140, 150]);

        // Verify metadata is preserved
        assert_eq!(ranges.seq_len, 50);
        assert!(ranges.reverse);
    }
}
