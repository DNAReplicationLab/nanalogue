//! `FilterModsByRefCoords` trait for filtering by coordinates on the reference genome
//! Provides interface for coordinate-based filtering operations

use crate::{Error, FiberAnnotation, OrdPair, Ranges};

/// Implements filter by coordinates on the reference genome.
pub trait FilterModsByRefCoords {
    /// filters by reference position i.e. all pos such that start <= pos < end
    /// are retained. does not use contig in filtering.
    ///
    /// # Errors
    /// Up to the user to set errors accordingly
    fn filter_mods_by_ref_pos(&mut self, _: u32, _: u32) -> Result<(), Error>;
}

/// Implements filter by reference coordinates for the Ranges
/// struct that contains our modification information.
/// NOTE: Ranges does not contain contig information, so we cannot
/// filter by that here.
impl FilterModsByRefCoords for Ranges {
    /// filters by reference position i.e. all pos such that start <= pos < end
    /// are retained. does not use contig in filtering.
    fn filter_mods_by_ref_pos(&mut self, start: u32, end: u32) -> Result<(), Error> {
        let interval = OrdPair::new(start, end)?;
        let matching_indices = self
            .annotations
            .iter()
            .enumerate()
            .filter_map(|(idx, ann)| {
                ann.ref_pos()
                    .filter(|pos| (interval.low()..interval.high()).contains(pos))
                    .map(|_| idx)
            });
        let (start_index, stop_index) = matching_indices
            .clone()
            .next()
            .zip(matching_indices.clone().next_back())
            .map_or((0, 0), |(first, last)| {
                (
                    first,
                    last.checked_add(1)
                        .expect("an existing vector index can be incremented"),
                )
            });

        self.annotations.truncate(stop_index);
        self.annotations.drain(..start_index).for_each(drop);

        if !self
            .annotations
            .iter()
            .filter_map(FiberAnnotation::ref_pos)
            .is_sorted_by(|previous, next| previous < next)
        {
            return Err(Error::WrongOrder(
                "reference positions are not strictly increasing".to_owned(),
            ));
        }
        Ok(())
    }
}

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use super::*;
    use crate::constants::shared::MAX_CONTIG_LEN;

    #[test]
    fn direct_ranges_filter_mods_by_ref_pos() {
        // Create a Ranges object with multiple per-base annotations.
        // All vectors have the same length as required.
        let mut ranges = Ranges {
            // Each entry represents a single base position with different properties.
            annotations: vec![
                FiberAnnotation::try_new(10, 100, Some(10)).unwrap(),
                FiberAnnotation::try_new(20, 120, Some(20)).unwrap(),
                FiberAnnotation::try_new(30, 140, Some(30)).unwrap(),
                FiberAnnotation::try_new(40, 160, Some(40)).unwrap(),
            ],
            seq_len: 50,
            reverse: false,
        };

        // Filter annotations to keep only those overlapping with reference region 18-32.
        ranges.filter_mods_by_ref_pos(18, 32).unwrap();

        // Verify that only the annotations that overlap with [18, 32) are kept.
        // Should keep the points at ref_pos 20 and 30 (indexes 1 and 2).
        assert_eq!(ranges.annotations.len(), 2);

        // Check the specific values were retained correctly.
        assert_eq!(ranges.pos().collect::<Vec<_>>(), vec![20, 30]);
        assert_eq!(
            ranges.ref_pos().collect::<Vec<_>>(),
            vec![Some(20), Some(30)]
        );
        assert_eq!(ranges.qual().collect::<Vec<_>>(), vec![120, 140]);

        // Verify seq_len and reverse flag are preserved
        assert_eq!(ranges.seq_len, 50);
        assert!(!ranges.reverse);
    }

    #[test]
    fn ranges_filter_mods_by_ref_pos_no_overlap() {
        // Create a Ranges object with annotations that don't overlap the target region.
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation::try_new(10, 100, Some(10)).unwrap(),
                FiberAnnotation::try_new(20, 120, Some(20)).unwrap(),
                FiberAnnotation::try_new(60, 140, Some(60)).unwrap(),
                FiberAnnotation::try_new(70, 160, Some(70)).unwrap(),
            ],
            seq_len: 80,
            reverse: true,
        };

        // Filter annotations for region [30, 50) - none should match.
        ranges.filter_mods_by_ref_pos(30, 50).unwrap();

        // Verify that all annotations were filtered out.
        assert_eq!(ranges.annotations.len(), 0);

        // Verify metadata is preserved
        assert_eq!(ranges.seq_len, 80);
        assert!(ranges.reverse);
    }

    #[test]
    fn point_annotation_is_included_at_its_own_position() {
        // A point at ref_pos 20 occupies [20, 21). Filtering [20, 21) must keep it;
        // filtering [21, 22) must not.
        let base = Ranges {
            annotations: vec![FiberAnnotation::try_new(5, 100, Some(20)).unwrap()],
            seq_len: 50,
            reverse: false,
        };

        let mut included = base.clone();
        included.filter_mods_by_ref_pos(20, 21).unwrap();
        assert_eq!(included.annotations.len(), 1);

        let mut excluded = base.clone();
        excluded.filter_mods_by_ref_pos(21, 22).unwrap();
        assert!(excluded.annotations.is_empty());
    }

    #[test]
    fn maximum_reference_position_does_not_overflow() {
        // The largest storable reference position is `MAX_CONTIG_LEN - 1`.
        let max_ref_pos = MAX_CONTIG_LEN - 1;
        let base = Ranges {
            annotations: vec![FiberAnnotation::try_new(5, 100, Some(max_ref_pos)).unwrap()],
            seq_len: 10,
            reverse: false,
        };

        let mut excluded = base.clone();
        excluded
            .filter_mods_by_ref_pos(MAX_CONTIG_LEN, u32::MAX)
            .unwrap();
        assert!(excluded.annotations.is_empty());

        let mut included = base;
        included
            .filter_mods_by_ref_pos(max_ref_pos, u32::MAX)
            .unwrap();
        assert_eq!(included.annotations.len(), 1);
    }

    #[test]
    fn ranges_filter_mods_by_ref_pos_with_none_values() {
        // Create a Ranges object with some None reference positions.
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation::try_new(10, 100, Some(10)).unwrap(),
                FiberAnnotation::try_new(20, 120, Some(20)).unwrap(),
                FiberAnnotation::try_new(30, 140, None).unwrap(),
                FiberAnnotation::try_new(40, 160, Some(40)).unwrap(),
            ],
            seq_len: 50,
            reverse: false,
        };

        // Filter annotations to keep only those overlapping with reference region 18-22.
        // Only the annotation at index 1 should be kept.
        ranges.filter_mods_by_ref_pos(18, 22).unwrap();

        // Verify that only the annotation that overlaps with [18, 22) is kept.
        assert_eq!(ranges.annotations.len(), 1);
        assert_eq!(ranges.ref_pos().collect::<Vec<_>>(), vec![Some(20)]);
        assert_eq!(ranges.qual().collect::<Vec<_>>(), vec![120]);

        // Verify metadata is preserved
        assert_eq!(ranges.seq_len, 50);
        assert!(!ranges.reverse);
    }

    #[test]
    fn ranges_filter_mods_by_ref_pos_with_none_values_2() {
        // Create a Ranges object with some None reference positions.
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation::try_new(10, 100, Some(10)).unwrap(),
                FiberAnnotation::try_new(19, 120, Some(20)).unwrap(),
                FiberAnnotation::try_new(20, 140, None).unwrap(),
                FiberAnnotation::try_new(21, 150, Some(21)).unwrap(),
                FiberAnnotation::try_new(40, 160, Some(40)).unwrap(),
            ],
            seq_len: 50,
            reverse: true,
        };

        // Filter annotations to keep only those overlapping with reference region 18-22.
        ranges.filter_mods_by_ref_pos(18, 22).unwrap();

        // Verify that only the annotations that overlap with [18, 22) are kept.
        assert_eq!(ranges.annotations.len(), 3);
        assert_eq!(
            ranges.ref_pos().collect::<Vec<_>>(),
            vec![Some(20), None, Some(21)]
        );
        assert_eq!(ranges.qual().collect::<Vec<_>>(), vec![120, 140, 150]);

        // Verify metadata is preserved
        assert_eq!(ranges.seq_len, 50);
        assert!(ranges.reverse);
    }

    #[test]
    fn ranges_filter_mods_by_ref_pos_with_start_equals_end() {
        // Test the edge case where start == end.
        // As this is a 0-bp interval, there must be no data remaining after filtering.
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation::try_new(10, 100, Some(10)).unwrap(),
                FiberAnnotation::try_new(20, 120, Some(20)).unwrap(),
                FiberAnnotation::try_new(30, 140, Some(30)).unwrap(),
                FiberAnnotation::try_new(40, 160, Some(40)).unwrap(),
            ],
            seq_len: 50,
            reverse: false,
        };

        // Filter with start == end (e.g., [20, 20))
        ranges.filter_mods_by_ref_pos(20, 20).unwrap();

        // Verify that no data remains
        assert!(ranges.annotations.is_empty());
    }

    #[test]
    #[should_panic(expected = "WrongOrder")]
    fn ranges_filter_panics_when_start_greater_than_end() {
        // Create a Ranges object with decreasing reference positions (wrong order).
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation::try_new(10, 100, Some(10)).unwrap(),
                // Decreasing reference positions should panic.
                FiberAnnotation::try_new(30, 120, Some(9)).unwrap(),
            ],
            seq_len: 50,
            reverse: false,
        };

        // This should panic because reference positions are not strictly increasing.
        ranges.filter_mods_by_ref_pos(0, 50).unwrap();
    }

    #[test]
    #[should_panic(expected = "WrongOrder")]
    fn ranges_filter_panics_when_starts_not_increasing() {
        // Create a Ranges object where reference positions are not strictly increasing.
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation::try_new(10, 100, Some(10)).unwrap(),
                // Same reference position should panic.
                FiberAnnotation::try_new(11, 120, Some(10)).unwrap(),
            ],
            seq_len: 50,
            reverse: false,
        };

        // This should panic because reference positions are not strictly increasing.
        ranges.filter_mods_by_ref_pos(0, 50).unwrap();
    }

    #[test]
    #[should_panic(expected = "WrongOrder")]
    fn ranges_filter_panics_when_starts_decreasing() {
        // Create a Ranges object where reference positions are decreasing.
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation::try_new(20, 100, Some(20)).unwrap(),
                // Decreasing reference positions should panic.
                FiberAnnotation::try_new(10, 120, Some(10)).unwrap(),
            ],
            seq_len: 50,
            reverse: false,
        };

        // This should panic because reference positions are decreasing.
        ranges.filter_mods_by_ref_pos(0, 50).unwrap();
    }

    #[test]
    #[should_panic(expected = "WrongOrder")]
    fn ranges_filter_panics_when_start_less_than_previous_end() {
        // Create a Ranges object where reference windows overlap (duplicate ref_pos).
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation::try_new(10, 100, Some(10)).unwrap(),
                // Overlapping reference windows (duplicate position) should panic.
                FiberAnnotation::try_new(11, 120, Some(10)).unwrap(),
            ],
            seq_len: 50,
            reverse: false,
        };

        // This should panic because reference windows are not strictly ordered (they overlap).
        ranges.filter_mods_by_ref_pos(0, 50).unwrap();
    }

    #[test]
    fn ignores_wrong_order_outside_retained_span() {
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation::try_new(10, 100, Some(40)).unwrap(),
                FiberAnnotation::try_new(20, 120, Some(20)).unwrap(),
            ],
            seq_len: 50,
            reverse: false,
        };

        ranges.filter_mods_by_ref_pos(15, 25).unwrap();

        assert_eq!(ranges.ref_pos().collect::<Vec<_>>(), vec![Some(20)]);
    }

    #[test]
    fn wrong_order_error_leaves_retained_span() {
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation::try_new(10, 100, Some(10)).unwrap(),
                FiberAnnotation::try_new(20, 120, Some(40)).unwrap(),
                FiberAnnotation::try_new(30, 140, Some(20)).unwrap(),
                FiberAnnotation::try_new(40, 160, Some(50)).unwrap(),
            ],
            seq_len: 60,
            reverse: false,
        };

        let error = ranges.filter_mods_by_ref_pos(15, 45).unwrap_err();

        assert!(matches!(error, Error::WrongOrder(_)));
        assert_eq!(
            ranges.ref_pos().collect::<Vec<_>>(),
            vec![Some(40), Some(20)]
        );
    }

    #[test]
    fn ranges_filter_mods_by_ref_pos_first_few_entries_none() {
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation::try_new(1, 100, None).unwrap(),
                FiberAnnotation::try_new(2, 120, None).unwrap(),
                FiberAnnotation::try_new(41, 140, Some(50)).unwrap(),
            ],
            seq_len: 96,
            reverse: false,
        };

        ranges.filter_mods_by_ref_pos(40, 61).unwrap();

        // Verify that only one point comes through
        assert_eq!(ranges.annotations.len(), 1);
        assert_eq!(ranges.qual().collect::<Vec<_>>(), vec![140u8]);

        ranges.filter_mods_by_ref_pos(70, 91).unwrap();

        // Verify that no point comes through
        assert!(ranges.annotations.is_empty());
        assert!(ranges.qual().next().is_none());
    }
}
