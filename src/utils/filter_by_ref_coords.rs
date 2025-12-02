//! `FilterByRefCoords` trait for filtering by coordinates on the reference genome
//! Provides interface for coordinate-based filtering operations

use crate::{Error, Intersects as _, OrdPair};
use fibertools_rs::utils::bamannotations::Ranges;
use serde::{Deserialize, Serialize};
use std::{cmp::Ordering, cmp::max, fmt};

/// Categorizes the types of windows into two possibilities
#[derive(Clone, Copy, Default, Debug, Serialize, Deserialize)]
pub struct WindowState(Option<OrdPair<u64>>);

impl fmt::Display for WindowState {
    /// converts to string for display i.e. "low, high" or "undefined"
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            WindowState(None) => "undefined".to_owned(),
            WindowState(Some(v)) => v.to_string(),
        }
        .fmt(f)
    }
}

impl WindowState {
    /// Construct a new `WindowState` given `Options` for `start` and `end`
    ///
    /// # Errors
    /// If open windows are encountered i.e. `start` or `end` are `None` but not
    /// both (both are `None` is fine), or if `start` or `end` are negative or
    /// ordered incorrectly, which will lead to an `OrdPair` error i.e. `start > end`
    pub fn new(start: Option<i64>, end: Option<i64>) -> Result<Self, Error> {
        Ok(match (start, end) {
            (None, None) => WindowState(None),
            (Some(v), None) | (None, Some(v)) => {
                return Err(Error::NotImplemented(format!(
                    "window was (None, {v}) or ({v}, None) - we cannot deal with these"
                )));
            }
            (Some(v), Some(w)) => WindowState(Some(OrdPair::<u64>::new(
                u64::try_from(v)?,
                u64::try_from(w)?,
            )?)),
        })
    }

    /// Intersects with a genomic region
    #[must_use]
    #[expect(
        clippy::arithmetic_side_effects,
        reason = "genomic coords << 2^64, so overflow on `lo + 1` unlikely"
    )]
    pub fn intersects(&self, interval: OrdPair<u64>) -> bool {
        match *self {
            WindowState(None) => false,
            WindowState(Some(v)) => {
                // For various reasons, we want to permit 0-bp intevals on the `WindowState` here
                // (but not on the genomic interval) and treat them like 1-bp intervals.
                // This is because `fibertools-rs` uses 0 bp intervals in the `FiberAnnotation`
                // struct (as that crate is evolving (still in version 0.x at time of writing),
                // this may change in the future).

                let lo = v.get_low();
                let hi = v.get_high();

                (interval.get_low()..interval.get_high()).intersects(&(lo..max(lo + 1, hi)))
            }
        }
    }
}

impl PartialEq for WindowState {
    fn eq(&self, other: &Self) -> bool {
        match (*self, *other) {
            (WindowState(None), WindowState(None)) => true,
            (WindowState(Some(v)), WindowState(Some(w))) if v == w => true,
            _ => false,
        }
    }
}

impl PartialOrd for WindowState {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (*self, *other) {
            (WindowState(None), WindowState(None)) => Some(Ordering::Equal),
            (WindowState(None), WindowState(Some(_))) => Some(Ordering::Less),
            (WindowState(Some(_)), WindowState(None)) => Some(Ordering::Greater),
            (WindowState(Some(v)), WindowState(Some(w))) if v == w => Some(Ordering::Equal),
            (WindowState(Some(v)), WindowState(Some(w))) if v.get_high() <= w.get_low() => {
                Some(Ordering::Less)
            }
            (WindowState(Some(v)), WindowState(Some(w))) if w.get_high() <= v.get_low() => {
                Some(Ordering::Greater)
            }
            _ => None,
        }
    }
}

/// Implements filter by coordinates on the reference genome.
pub trait FilterByRefCoords {
    /// filters by reference position i.e. all pos such that start <= pos < end
    /// are retained. does not use contig in filtering.
    ///
    /// # Errors
    /// Up to the user to set errors accordingly
    fn filter_by_ref_pos(&mut self, _: i64, _: i64) -> Result<(), Error> {
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
    fn filter_by_ref_pos(&mut self, start: i64, end: i64) -> Result<(), Error> {
        let genomic_interval = format!("{start}-{end}");
        let interval = OrdPair::<u64>::from_interval(&genomic_interval)?;

        let (start_index, stop_index_plus_one) = {
            let mut coord_limits: Option<OrdPair<usize>> = None;
            let mut previous_window = WindowState(None);
            for k in self
                .reference_starts()
                .iter()
                .zip(self.reference_ends().iter())
                .enumerate()
            {
                let window_state = WindowState::new(*(k.1.0), *(k.1.1))?;
                match window_state {
                    WindowState(None) => {}
                    w @ WindowState(Some(_)) if previous_window < w => previous_window = w,
                    v => {
                        return Err(Error::WrongOrder(format!(
                            "windows are not ordered, previous window {previous_window} is > or overlaps with {v}!"
                        )));
                    }
                }
                if window_state.intersects(interval) {
                    if coord_limits.is_none() {
                        coord_limits = Some(OrdPair::<usize>::new(k.0, k.0)?);
                    } else {
                        let Some(ref mut v) = coord_limits else {
                            unreachable!("we've checked for the `Some` variant already")
                        };
                        v.update_high(k.0)?;
                    }
                }
            }
            if let Some(item) = coord_limits {
                (
                    item.get_low(),
                    item.get_high().checked_add(1).ok_or(Error::Arithmetic(
                        "overflow error in coordinates while filtering by ref".to_owned(),
                    ))?,
                )
            } else {
                (0usize, 0usize)
            }
            /*
            if coord_limits.is_none() {
                (0usize, 0usize)
            } else {
                (
                    coord_limits
                        .expect("no error as we've checked if `Some`")
                        .get_low(),
                    coord_limits
                        .expect("no error as we've checked if `Some`")
                        .get_high()
                        .checked_add(1)
                        .ok_or(Error::Arithmetic("overflow error in coordinates while filtering by ref".to_owned()))?
                )
            }*/
        };
        self.annotations.truncate(stop_index_plus_one);
        self.annotations.drain(0..start_index).for_each(drop);
        Ok(())
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
        ranges.filter_by_ref_pos(18, 32).unwrap();

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
        ranges.filter_by_ref_pos(30, 50).unwrap();

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
        ranges.filter_by_ref_pos(18, 22).unwrap();

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
        ranges.filter_by_ref_pos(18, 22).unwrap();

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
    #[should_panic(expected = "OrdPairConversion")]
    fn ranges_filter_by_ref_pos_with_start_equals_end() {
        // Test the edge case where start == end.
        // As this is a 0-bp interval, we should just get an error.
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
        ranges.filter_by_ref_pos(25, 25).unwrap();
    }

    #[test]
    #[should_panic(expected = "WrongOrder")]
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
        ranges.filter_by_ref_pos(0, 50).unwrap();
    }

    #[test]
    #[should_panic(expected = "WrongOrder")]
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
        ranges.filter_by_ref_pos(0, 50).unwrap();
    }

    #[test]
    #[should_panic(expected = "WrongOrder")]
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
        ranges.filter_by_ref_pos(0, 50).unwrap();
    }

    #[test]
    #[should_panic(expected = "WrongOrder")]
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
        ranges.filter_by_ref_pos(0, 50).unwrap();
    }

    #[test]
    fn ranges_filter_by_ref_pos_first_few_entries_none() {
        let mut ranges = Ranges {
            annotations: vec![
                FiberAnnotation {
                    start: 1,
                    end: 2,
                    length: 1,
                    qual: 100,
                    reference_start: None,
                    reference_end: None,
                    reference_length: None,
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 2,
                    end: 3,
                    length: 1,
                    qual: 120,
                    reference_start: None,
                    reference_end: None,
                    reference_length: None,
                    extra_columns: None,
                },
                FiberAnnotation {
                    start: 41,
                    end: 42,
                    length: 1,
                    qual: 140,
                    reference_start: Some(50),
                    reference_end: Some(51),
                    reference_length: Some(1),
                    extra_columns: None,
                },
            ],
            seq_len: 96,
            reverse: false,
        };

        ranges.filter_by_ref_pos(40, 61).unwrap();

        // Verify that only one point comes through
        assert_eq!(ranges.annotations.len(), 1);
        assert_eq!(ranges.qual(), vec![140u8]);
    }
}
