//! `FilterByRefCoords` trait for filtering by coordinates on the reference genome
//! Provides interface for coordinate-based filtering operations

use fibertools_rs::utils::bamannotations::Ranges;

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
                .reference_starts()
                .iter()
                .zip(self.reference_ends().iter())
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
                    (*k.1.0).is_none() || (*k.1.0) >= previous_end,
                    "start {:?} >= previous_end {:?} expected",
                    *(k.1.0),
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

        self.annotations.truncate(stop_index);
        self.annotations.drain(0..start_index).for_each(drop);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fibertools_rs::utils::bamannotations::FiberAnnotation;

    #[test]
    fn direct_ranges_filter_by_ref_pos() {
        // Create a Ranges object with multiple ranges
        // All vectors have the same length as required
        let mut ranges = Ranges {
            // Each entry represents a genomic range with different properties
            annotations: vec![
                FiberAnnotation {
                    start: 10,
                    end: 15,
                    length: 5,
                    qual: 100,
                    reference_start: Some(10),
                    reference_end: Some(15),
                    reference_length: Some(5),
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 20,
                    end: 25,
                    length: 5,
                    qual: 120,
                    reference_start: Some(20),
                    reference_end: Some(25),
                    reference_length: Some(5),
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 30,
                    end: 35,
                    length: 5,
                    qual: 140,
                    reference_start: Some(30),
                    reference_end: Some(35),
                    reference_length: Some(5),
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 40,
                    end: 45,
                    length: 5,
                    qual: 160,
                    reference_start: Some(40),
                    reference_end: Some(45),
                    reference_length: Some(5),
                    extra_columns: None,
                },
            ],
            seq_len: 50,
            reverse: false,
        };

        // Filter ranges to keep only those overlapping with reference region 18-32
        ranges.filter_by_ref_pos(18, 32);

        // Verify that only the ranges that overlap with [18, 32) are kept
        // Should keep ranges [20-25] and [30-35] (indexes 1 and 2)
        assert_eq!(ranges.lengths().len(), 2);
        assert_eq!(ranges.reference_lengths().len(), 2);

        // Check the specific values were retained correctly
        assert_eq!(ranges.starts(), vec![20, 30]);
        assert_eq!(ranges.ends(), vec![25, 35]);
        assert_eq!(ranges.reference_starts(), vec![Some(20), Some(30)]);
        assert_eq!(ranges.reference_ends(), vec![Some(25), Some(35)]);
        assert_eq!(ranges.qual(), vec![120, 140]);

        // Verify seq_len and reverse flag are preserved
        assert_eq!(ranges.seq_len, 50);
        assert!(!ranges.reverse);
    }

    #[test]
    fn ranges_filter_by_ref_pos_no_overlap() {
        // Create a Ranges object with ranges that don't overlap the target region
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation {
                    start: 10,
                    end: 15,
                    length: 5,
                    qual: 100,
                    reference_start: Some(10),
                    reference_end: Some(15),
                    reference_length: Some(5),
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 20,
                    end: 25,
                    length: 5,
                    qual: 120,
                    reference_start: Some(20),
                    reference_end: Some(25),
                    reference_length: Some(5),
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 60,
                    end: 65,
                    length: 5,
                    qual: 140,
                    reference_start: Some(60),
                    reference_end: Some(65),
                    reference_length: Some(5),
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 70,
                    end: 75,
                    length: 5,
                    qual: 160,
                    reference_start: Some(70),
                    reference_end: Some(75),
                    reference_length: Some(5),
                    extra_columns: None,
                },
            ],
            seq_len: 80,
            reverse: true,
        };

        // Filter ranges for region [30, 50) - none should match
        ranges.filter_by_ref_pos(30, 50);

        // Verify that all ranges were filtered out
        assert_eq!(ranges.annotations.len(), 0);

        // Verify metadata is preserved
        assert_eq!(ranges.seq_len, 80);
        assert!(ranges.reverse);
    }

    #[test]
    fn ranges_filter_by_ref_pos_with_none_values() {
        // Create a Ranges object with some None values
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation {
                    start: 10,
                    end: 15,
                    length: 5,
                    qual: 100,
                    reference_start: Some(10),
                    reference_end: Some(15),
                    reference_length: Some(5),
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 20,
                    end: 25,
                    length: 5,
                    qual: 120,
                    reference_start: Some(20),
                    reference_end: Some(25),
                    reference_length: Some(5),
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 30,
                    end: 35,
                    length: 5,
                    qual: 140,
                    reference_start: None,
                    reference_end: None,
                    reference_length: None,
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 40,
                    end: 45,
                    length: 5,
                    qual: 160,
                    reference_start: Some(40),
                    reference_end: Some(45),
                    reference_length: Some(5),
                    extra_columns: None,
                },
            ],
            seq_len: 50,
            reverse: false,
        };

        // Filter ranges to keep only those overlapping with reference region 18-22
        // Only the range at index 1 should be kept
        ranges.filter_by_ref_pos(18, 22);

        // Verify that only the range that overlaps with [18, 22) is kept
        assert_eq!(ranges.annotations.len(), 1);
        assert_eq!(ranges.reference_starts(), vec![Some(20)]);
        assert_eq!(ranges.reference_ends(), vec![Some(25)]);
        assert_eq!(ranges.qual(), vec![120]);

        // Verify metadata is preserved
        assert_eq!(ranges.seq_len, 50);
        assert!(!ranges.reverse);
    }

    #[test]
    fn ranges_filter_by_ref_pos_with_none_values_2() {
        // Create a Ranges object with some None values
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation {
                    start: 10,
                    end: 15,
                    length: 5,
                    qual: 100,
                    reference_start: Some(10),
                    reference_end: Some(15),
                    reference_length: Some(5),
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 19,
                    end: 20,
                    length: 1,
                    qual: 120,
                    reference_start: Some(20),
                    reference_end: Some(21),
                    reference_length: Some(1),
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 20,
                    end: 21,
                    length: 1,
                    qual: 140,
                    reference_start: None,
                    reference_end: None,
                    reference_length: None,
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 21,
                    end: 22,
                    length: 1,
                    qual: 150,
                    reference_start: Some(21),
                    reference_end: Some(22),
                    reference_length: Some(1),
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 40,
                    end: 45,
                    length: 5,
                    qual: 160,
                    reference_start: Some(40),
                    reference_end: Some(45),
                    reference_length: Some(5),
                    extra_columns: None,
                },
            ],
            seq_len: 50,
            reverse: true,
        };

        // Filter ranges to keep only those overlapping with reference region 18-22
        ranges.filter_by_ref_pos(18, 22);

        // Verify that only the ranges that overlaps with [18, 22) are kept
        assert_eq!(ranges.annotations.len(), 3);
        assert_eq!(ranges.reference_starts(), vec![Some(20), None, Some(21)]);
        assert_eq!(ranges.reference_ends(), vec![Some(21), None, Some(22)]);
        assert_eq!(ranges.qual(), vec![120, 140, 150]);

        // Verify metadata is preserved
        assert_eq!(ranges.seq_len, 50);
        assert!(ranges.reverse);
    }

    #[test]
    fn ranges_filter_by_ref_pos_with_start_equals_end() {
        // Test the edge case where start == end.
        // Since filtering keeps positions where start <= pos < end,
        // when start == end, no positions satisfy this condition, so no data should be retained.
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation {
                    start: 10,
                    end: 15,
                    length: 5,
                    qual: 100,
                    reference_start: Some(10),
                    reference_end: Some(15),
                    reference_length: Some(5),
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 20,
                    end: 25,
                    length: 5,
                    qual: 120,
                    reference_start: Some(20),
                    reference_end: Some(25),
                    reference_length: Some(5),
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 30,
                    end: 35,
                    length: 5,
                    qual: 140,
                    reference_start: Some(30),
                    reference_end: Some(35),
                    reference_length: Some(5),
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 40,
                    end: 45,
                    length: 5,
                    qual: 160,
                    reference_start: Some(40),
                    reference_end: Some(45),
                    reference_length: Some(5),
                    extra_columns: None,
                },
            ],
            seq_len: 50,
            reverse: false,
        };

        // Filter with start == end (e.g., [25, 25))
        ranges.filter_by_ref_pos(25, 25);

        // Verify that all data was filtered out
        assert_eq!(ranges.annotations.len(), 0);

        // Verify metadata is preserved
        assert_eq!(ranges.seq_len, 50);
        assert!(!ranges.reverse);
    }

    #[test]
    #[should_panic(expected = "start")]
    fn ranges_filter_panics_when_start_greater_than_end() {
        // Create a Ranges object where start > end, violating the first assertion
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation {
                    start: 10,
                    end: 15,
                    length: 5,
                    qual: 100,
                    reference_start: Some(10),
                    reference_end: Some(15),
                    reference_length: Some(5),
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 30,
                    end: 25,
                    length: 5,
                    qual: 120,
                    reference_start: Some(30),
                    reference_end: Some(25), // 30 > 25, this should panic
                    reference_length: Some(5),
                    extra_columns: None,
                },
            ],
            seq_len: 50,
            reverse: false,
        };

        // This should panic because start (30) > end (25)
        ranges.filter_by_ref_pos(0, 50);
    }

    #[test]
    #[should_panic(expected = "start")]
    fn ranges_filter_panics_when_starts_not_increasing() {
        // Create a Ranges object where starts are not strictly increasing
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation {
                    start: 10,
                    end: 15,
                    length: 5,
                    qual: 100,
                    reference_start: Some(10),
                    reference_end: Some(15),
                    reference_length: Some(5),
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 10,
                    end: 20,
                    length: 10,
                    qual: 120,
                    reference_start: Some(10), // Same start value, should panic
                    reference_end: Some(20),
                    reference_length: Some(10),
                    extra_columns: None,
                },
            ],
            seq_len: 50,
            reverse: false,
        };

        // This should panic because start values are not strictly increasing
        ranges.filter_by_ref_pos(0, 50);
    }

    #[test]
    #[should_panic(expected = "start")]
    fn ranges_filter_panics_when_starts_decreasing() {
        // Create a Ranges object where starts are decreasing
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation {
                    start: 20,
                    end: 25,
                    length: 5,
                    qual: 100,
                    reference_start: Some(20),
                    reference_end: Some(25),
                    reference_length: Some(5),
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 10,
                    end: 15,
                    length: 5,
                    qual: 120,
                    reference_start: Some(10), // Decreasing start values, should panic
                    reference_end: Some(15),
                    reference_length: Some(5),
                    extra_columns: None,
                },
            ],
            seq_len: 50,
            reverse: false,
        };

        // This should panic because starts are decreasing
        ranges.filter_by_ref_pos(0, 50);
    }

    #[test]
    #[should_panic(expected = "start")]
    fn ranges_filter_panics_when_start_less_than_previous_end() {
        // Create a Ranges object where a start is less than previous end (overlapping ranges)
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation {
                    start: 10,
                    end: 15,
                    length: 5,
                    qual: 100,
                    reference_start: Some(10),
                    reference_end: Some(15),
                    reference_length: Some(5),
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 14,
                    end: 20,
                    length: 6,
                    qual: 120,
                    reference_start: Some(14), // start (14) < previous_end (15)
                    reference_end: Some(20),
                    reference_length: Some(6),
                    extra_columns: None,
                },
            ],
            seq_len: 50,
            reverse: false,
        };

        // This should panic because start (14) < previous_end (15)
        ranges.filter_by_ref_pos(0, 50);
    }
}
