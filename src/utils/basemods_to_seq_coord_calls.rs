//! `BaseModsToSeqCoordCalls` trait, converts mod data in the fibertools
//! `BaseMods` struct into essentially an array of length = sequence length where every
//! entry is a `Vec` of `u8` containing modification probability per mod.

use crate::{Error, ModChar};
use fibertools_rs::utils::basemods::BaseMods;

/// Holds mod calls for every position along a sequence
/// (irrespective of whether every position has an associated mod call).
///
/// * `mod_types` is a `Vec` of tuples of [`ModChar`], `bool`,
///   which mean modification type, and whether the calls are on the basecalled
///   strand or the opposite.
/// * `mod calls` is a `Vec` holding a `Vec` of mod calls per position,
///   and whose length equals length of the sequence.
#[derive(Debug)]
#[non_exhaustive]
pub struct SeqCoordCalls {
    /// Types of modification + strand combination (max 8)
    mod_types: Vec<(ModChar, bool)>,

    /// Mod calls per position
    mod_calls: Vec<Vec<u8>>,
}

impl SeqCoordCalls {
    /// Returns a boolean vector indicating whether any modification is present at each position.
    ///
    /// The returned vector has the same length as the DNA sequence, with `true` at positions
    /// where at least one modification type has a non-zero quality score, and `false` otherwise.
    ///
    /// # Example
    ///
    /// ```
    /// # use nanalogue_core::utils::basemods_to_seq_coord_calls::SeqCoordCalls;
    /// # use fibertools_rs::utils::basemods::{BaseMods, BaseMod};
    /// # use fibertools_rs::utils::bamannotations::{FiberAnnotation, Ranges};
    /// // Given BaseMods with T+ at position 0 (qual=100) and C+m at position 2 (qual=200)
    /// // for a sequence of length 5:
    /// let base_mods = BaseMods {
    ///     base_mods: vec![
    ///         BaseMod {
    ///             modified_base: b'T',
    ///             strand: '+',
    ///             modification_type: 'T',
    ///             ranges: Ranges {
    ///                 annotations: vec![FiberAnnotation {
    ///                     start: 0, end: 1, length: 1, qual: 100,
    ///                     reference_start: Some(0), reference_end: Some(0),
    ///                     reference_length: Some(0), extra_columns: None,
    ///                 }],
    ///                 seq_len: 5, reverse: false,
    ///             },
    ///             record_is_reverse: false,
    ///         },
    ///         BaseMod {
    ///             modified_base: b'C',
    ///             strand: '+',
    ///             modification_type: 'm',
    ///             ranges: Ranges {
    ///                 annotations: vec![FiberAnnotation {
    ///                     start: 2, end: 3, length: 1, qual: 200,
    ///                     reference_start: None, reference_end: None,
    ///                     reference_length: None, extra_columns: None,
    ///                 }],
    ///                 seq_len: 5, reverse: false,
    ///             },
    ///             record_is_reverse: false,
    ///         },
    ///     ],
    /// };
    ///
    /// let seq_coord_calls = SeqCoordCalls::try_from(&base_mods).unwrap();
    /// let collapsed = seq_coord_calls.collapse_mod_calls();
    ///
    /// // Returns [true, false, true, false, false]
    /// // Position 0: true (T+ modified)
    /// // Position 1: false (no modifications)
    /// // Position 2: true (C+m modified)
    /// // Positions 3-4: false (no modifications)
    /// assert_eq!(collapsed, vec![true, false, true, false, false]);
    /// ```
    #[must_use]
    pub fn collapse_mod_calls(&self) -> Vec<bool> {
        self.mod_calls
            .iter()
            .map(|k| k.iter().any(|x| *x > 0))
            .collect()
    }

    /// Returns the modification types and their strand information.
    ///
    /// Each entry in the returned vector is a tuple `(ModChar, bool)` where:
    /// - The first element is the modification type (e.g., 'T', 'm', 'h')
    /// - The second element indicates strand: `true` for '+' strand, `false` for '-' strand
    ///
    /// The order of mod types in this vector corresponds to the order of quality scores
    /// in the arrays returned by [`mod_calls()`](Self::mod_calls).
    ///
    /// # Example
    ///
    /// ```
    /// # use nanalogue_core::utils::basemods_to_seq_coord_calls::SeqCoordCalls;
    /// # use nanalogue_core::ModChar;
    /// # use fibertools_rs::utils::basemods::{BaseMods, BaseMod};
    /// # use fibertools_rs::utils::bamannotations::{FiberAnnotation, Ranges};
    /// // Given BaseMods with T+ and C+m modifications:
    /// let base_mods = BaseMods {
    ///     base_mods: vec![
    ///         BaseMod {
    ///             modified_base: b'T', strand: '+', modification_type: 'T',
    ///             ranges: Ranges {
    ///                 annotations: vec![FiberAnnotation {
    ///                     start: 0, end: 1, length: 1, qual: 100,
    ///                     reference_start: Some(0), reference_end: Some(0),
    ///                     reference_length: Some(0), extra_columns: None,
    ///                 }],
    ///                 seq_len: 5, reverse: false,
    ///             },
    ///             record_is_reverse: false,
    ///         },
    ///         BaseMod {
    ///             modified_base: b'C', strand: '+', modification_type: 'm',
    ///             ranges: Ranges {
    ///                 annotations: vec![FiberAnnotation {
    ///                     start: 2, end: 3, length: 1, qual: 200,
    ///                     reference_start: None, reference_end: None,
    ///                     reference_length: None, extra_columns: None,
    ///                 }],
    ///                 seq_len: 5, reverse: false,
    ///             },
    ///             record_is_reverse: false,
    ///         },
    ///     ],
    /// };
    ///
    /// let seq_coord_calls = SeqCoordCalls::try_from(&base_mods).unwrap();
    /// let mod_types = seq_coord_calls.mod_types();
    ///
    /// // Returns 2 mod types: T+ and C+m
    /// assert_eq!(mod_types.len(), 2);
    /// assert_eq!(mod_types[0].0, ModChar::new('T'));
    /// assert_eq!(mod_types[0].1, true);  // '+' strand
    /// assert_eq!(mod_types[1].0, ModChar::new('m'));
    /// assert_eq!(mod_types[1].1, true);  // '+' strand
    /// ```
    #[must_use]
    pub fn mod_types(&self) -> &Vec<(ModChar, bool)> {
        &self.mod_types
    }

    /// Returns modification quality scores at a given sequence position.
    ///
    /// Returns an array of `u8` quality scores, one for each modification type.
    /// The order corresponds to the modification types returned by [`mod_types()`](Self::mod_types).
    /// A value of `0` indicates no modification detected at that position for that type.
    ///
    /// Returns `None` if the position is out of bounds.
    ///
    /// # Arguments
    ///
    /// * `pos` - The sequence position (0-based index)
    ///
    /// # Example
    ///
    /// ```
    /// # use nanalogue_core::utils::basemods_to_seq_coord_calls::SeqCoordCalls;
    /// # use fibertools_rs::utils::basemods::{BaseMods, BaseMod};
    /// # use fibertools_rs::utils::bamannotations::{FiberAnnotation, Ranges};
    /// // Given BaseMods with T+ at position 0 (qual=100) and C+m at position 2 (qual=200):
    /// let base_mods = BaseMods {
    ///     base_mods: vec![
    ///         BaseMod {
    ///             modified_base: b'T', strand: '+', modification_type: 'T',
    ///             ranges: Ranges {
    ///                 annotations: vec![FiberAnnotation {
    ///                     start: 0, end: 1, length: 1, qual: 100,
    ///                     reference_start: Some(0), reference_end: Some(0),
    ///                     reference_length: Some(0), extra_columns: None,
    ///                 }],
    ///                 seq_len: 5, reverse: false,
    ///             },
    ///             record_is_reverse: false,
    ///         },
    ///         BaseMod {
    ///             modified_base: b'C', strand: '+', modification_type: 'm',
    ///             ranges: Ranges {
    ///                 annotations: vec![FiberAnnotation {
    ///                     start: 2, end: 3, length: 1, qual: 200,
    ///                     reference_start: None, reference_end: None,
    ///                     reference_length: None, extra_columns: None,
    ///                 }],
    ///                 seq_len: 5, reverse: false,
    ///             },
    ///             record_is_reverse: false,
    ///         },
    ///     ],
    /// };
    ///
    /// let seq_coord_calls = SeqCoordCalls::try_from(&base_mods).unwrap();
    ///
    /// // Position 0: T+ modified (100), C+m not modified (0)
    /// assert_eq!(seq_coord_calls.mod_calls(0), Some(&[100u8, 0u8][..]));
    ///
    /// // Position 1: Neither modification present
    /// assert_eq!(seq_coord_calls.mod_calls(1), Some(&[0u8, 0u8][..]));
    ///
    /// // Position 2: T+ not modified (0), C+m modified (200)
    /// assert_eq!(seq_coord_calls.mod_calls(2), Some(&[0u8, 200u8][..]));
    ///
    /// // Out of bounds returns None
    /// assert_eq!(seq_coord_calls.mod_calls(5), None);
    /// ```
    #[must_use]
    pub fn mod_calls(&self, pos: usize) -> Option<&[u8]> {
        self.mod_calls.get(pos).map(|v| &**v)
    }
}

impl TryFrom<&BaseMods> for SeqCoordCalls {
    type Error = Error;

    /// Converts `BaseMods` data into sequence-coordinate-indexed modification calls.
    ///
    /// Transforms modification data from the fibertools `BaseMods` format (which stores
    /// per-modification-type annotations with positions and quality scores) into a
    /// sequence-indexed format where each position has quality scores for all modification types.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The input `BaseMods` is empty (no sequence length can be determined)
    /// - Different `BaseMod` entries have inconsistent sequence lengths
    /// - Any strand character is not '+' or '-'
    /// - Any annotation is not per-base (e.g., multi-base ranges are not supported)
    /// - Position coordinates are out of bounds
    ///
    /// # Example
    ///
    /// ```
    /// # use nanalogue_core::utils::basemods_to_seq_coord_calls::SeqCoordCalls;
    /// # use fibertools_rs::utils::basemods::{BaseMods, BaseMod};
    /// # use fibertools_rs::utils::bamannotations::{FiberAnnotation, Ranges};
    /// // Create BaseMods with multiple modification types (T+ and C+m):
    /// let base_mods = BaseMods {
    ///     base_mods: vec![
    ///         BaseMod {
    ///             modified_base: b'T',
    ///             strand: '+',
    ///             modification_type: 'T',
    ///             ranges: Ranges {
    ///                 annotations: vec![FiberAnnotation {
    ///                     start: 0, end: 1, length: 1, qual: 100,
    ///                     reference_start: Some(0), reference_end: Some(0),
    ///                     reference_length: Some(0), extra_columns: None,
    ///                 }],
    ///                 seq_len: 5, reverse: false,
    ///             },
    ///             record_is_reverse: false,
    ///         },
    ///         BaseMod {
    ///             modified_base: b'C',
    ///             strand: '+',
    ///             modification_type: 'm',
    ///             ranges: Ranges {
    ///                 annotations: vec![FiberAnnotation {
    ///                     start: 2, end: 3, length: 1, qual: 200,
    ///                     reference_start: None, reference_end: None,
    ///                     reference_length: None, extra_columns: None,
    ///                 }],
    ///                 seq_len: 5, reverse: false,
    ///             },
    ///             record_is_reverse: false,
    ///         },
    ///     ],
    /// };
    ///
    /// // Convert to sequence-indexed format
    /// let seq_coord_calls = SeqCoordCalls::try_from(&base_mods).unwrap();
    ///
    /// // Now we can access modification data by sequence position:
    /// // - Position 0: T+ has qual=100, C+m has qual=0
    /// // - Position 2: T+ has qual=0, C+m has qual=200
    /// // - Other positions: both have qual=0
    /// assert_eq!(seq_coord_calls.mod_calls(0), Some(&[100u8, 0u8][..]));
    /// assert_eq!(seq_coord_calls.mod_calls(2), Some(&[0u8, 200u8][..]));
    /// ```
    fn try_from(value: &BaseMods) -> Result<Self, Error> {
        let (seq_lengths, mod_type_collection, starts_qual) = value
            .base_mods
            .iter()
            .map(|x| {
                let seq_len = x.ranges.seq_len;
                Ok((
                    usize::try_from(seq_len)?,
                    (ModChar::new(x.modification_type),
                    match x.strand{
                        '+' => true,
                        '-' => false,
                        v => return Err(Error::InvalidState(format!("{v} is not a valid strand!"))),
                    }),
                    x.ranges
                        .annotations
                        .iter()
                        .map(|y| {
                            if (0..seq_len).contains(&y.start) && y.end.checked_sub(y.start) == Some(1) {
                                Ok((usize::try_from(y.start)?, y.qual))
                            } else {
                                Err(Error::NotImplemented(
                                    "cannot use a non-per-base annotation and/or malformed coordinates in `BaseMods`"
                                        .to_owned(),
                                ))
                            }
                        })
                        .collect::<Result<Vec<(usize, u8)>, Error>>()?,
                ))
            })
            .collect::<Result<(Vec<usize>, Vec<(ModChar, bool)>, Vec<Vec<(usize, u8)>>), Error>>(
            )?;

        let Some(seq_len) = seq_lengths.first() else {
            return Err(Error::InvalidSeqLength(
                "In an empty `BaseMods`, we don't know sequence length!".to_owned(),
            ));
        };

        if seq_lengths.iter().any(|x| *x != *seq_len) {
            return Err(Error::NotImplemented(
                "malformed `seq_len` in `BaseMods`".to_owned(),
            ));
        }

        let mut mod_list = vec![vec![0u8; seq_lengths.len()]; *seq_len];

        for (k, item) in starts_qual.into_iter().enumerate() {
            for l in item {
                *mod_list
                    .get_mut(l.0)
                    .expect("we've checked all pos < seq_len and all seq_len equal")
                    .get_mut(k)
                    .expect("we've set vector length to `seq_lengths.len()`") = l.1;
            }
        }

        Ok(SeqCoordCalls {
            mod_types: mod_type_collection,
            mod_calls: mod_list,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fibertools_rs::utils::bamannotations::{FiberAnnotation, Ranges};
    use fibertools_rs::utils::basemods::BaseMod;

    #[test]
    fn seq_coord_calls_from_basemods_example_1() -> Result<(), Error> {
        // Example from src/lib.rs nanalogue_mm_ml_parser docstring (count 0)
        let base_mods = BaseMods {
            base_mods: vec![BaseMod {
                modified_base: b'T',
                strand: '+',
                modification_type: 'T',
                ranges: Ranges {
                    annotations: vec![
                        FiberAnnotation {
                            start: 0,
                            end: 1,
                            length: 1,
                            qual: 4,
                            reference_start: Some(9),
                            reference_end: Some(9),
                            reference_length: Some(0),
                            extra_columns: None,
                        },
                        FiberAnnotation {
                            start: 3,
                            end: 4,
                            length: 1,
                            qual: 7,
                            reference_start: Some(12),
                            reference_end: Some(12),
                            reference_length: Some(0),
                            extra_columns: None,
                        },
                        FiberAnnotation {
                            start: 4,
                            end: 5,
                            length: 1,
                            qual: 9,
                            reference_start: Some(13),
                            reference_end: Some(13),
                            reference_length: Some(0),
                            extra_columns: None,
                        },
                        FiberAnnotation {
                            start: 7,
                            end: 8,
                            length: 1,
                            qual: 6,
                            reference_start: Some(16),
                            reference_end: Some(16),
                            reference_length: Some(0),
                            extra_columns: None,
                        },
                    ],
                    seq_len: 8,
                    reverse: false,
                },
                record_is_reverse: false,
            }],
        };

        let seq_coord_calls = SeqCoordCalls::try_from(&base_mods)?;

        // Test mod_types
        let mod_types = seq_coord_calls.mod_types();
        assert_eq!(mod_types.len(), 1);
        assert_eq!(
            mod_types.first().expect("has 1 element").0,
            ModChar::new('T')
        );
        assert!(mod_types.first().expect("has 1 element").1); // + strand

        // Test mod_calls for each position
        assert_eq!(seq_coord_calls.mod_calls(0), Some(&[4u8][..]));
        assert_eq!(seq_coord_calls.mod_calls(1), Some(&[0u8][..]));
        assert_eq!(seq_coord_calls.mod_calls(2), Some(&[0u8][..]));
        assert_eq!(seq_coord_calls.mod_calls(3), Some(&[7u8][..]));
        assert_eq!(seq_coord_calls.mod_calls(4), Some(&[9u8][..]));
        assert_eq!(seq_coord_calls.mod_calls(5), Some(&[0u8][..]));
        assert_eq!(seq_coord_calls.mod_calls(6), Some(&[0u8][..]));
        assert_eq!(seq_coord_calls.mod_calls(7), Some(&[6u8][..]));
        assert_eq!(seq_coord_calls.mod_calls(8), None); // Out of bounds

        // Test collapse_mod_calls
        let collapsed = seq_coord_calls.collapse_mod_calls();
        assert_eq!(
            collapsed,
            vec![true, false, false, true, true, false, false, true]
        );

        Ok(())
    }

    #[test]
    fn seq_coord_calls_from_basemods_example_2() -> Result<(), Error> {
        // Example from src/lib.rs nanalogue_mm_ml_parser docstring (count 2)
        let base_mods = BaseMods {
            base_mods: vec![BaseMod {
                modified_base: b'T',
                strand: '+',
                modification_type: 'T',
                ranges: Ranges {
                    annotations: vec![
                        FiberAnnotation {
                            start: 12,
                            end: 13,
                            length: 1,
                            qual: 3,
                            reference_start: Some(15),
                            reference_end: Some(15),
                            reference_length: Some(0),
                            extra_columns: None,
                        },
                        FiberAnnotation {
                            start: 13,
                            end: 14,
                            length: 1,
                            qual: 3,
                            reference_start: Some(16),
                            reference_end: Some(16),
                            reference_length: Some(0),
                            extra_columns: None,
                        },
                        FiberAnnotation {
                            start: 16,
                            end: 17,
                            length: 1,
                            qual: 4,
                            reference_start: Some(19),
                            reference_end: Some(19),
                            reference_length: Some(0),
                            extra_columns: None,
                        },
                        FiberAnnotation {
                            start: 19,
                            end: 20,
                            length: 1,
                            qual: 3,
                            reference_start: Some(22),
                            reference_end: Some(22),
                            reference_length: Some(0),
                            extra_columns: None,
                        },
                        FiberAnnotation {
                            start: 20,
                            end: 21,
                            length: 1,
                            qual: 182,
                            reference_start: Some(23),
                            reference_end: Some(23),
                            reference_length: Some(0),
                            extra_columns: None,
                        },
                    ],
                    seq_len: 33,
                    reverse: true,
                },
                record_is_reverse: true,
            }],
        };

        let seq_coord_calls = SeqCoordCalls::try_from(&base_mods)?;

        // Test mod_types
        let mod_types = seq_coord_calls.mod_types();
        assert_eq!(mod_types.len(), 1);
        assert_eq!(
            mod_types.first().expect("has 1 element").0,
            ModChar::new('T')
        );
        assert!(mod_types.first().expect("has 1 element").1); // + strand

        // Test mod_calls for positions with modifications
        assert_eq!(seq_coord_calls.mod_calls(12), Some(&[3u8][..]));
        assert_eq!(seq_coord_calls.mod_calls(13), Some(&[3u8][..]));
        assert_eq!(seq_coord_calls.mod_calls(16), Some(&[4u8][..]));
        assert_eq!(seq_coord_calls.mod_calls(19), Some(&[3u8][..]));
        assert_eq!(seq_coord_calls.mod_calls(20), Some(&[182u8][..]));

        // Test mod_calls for positions without modifications
        for k in 0..33 {
            if [12, 13, 16, 19, 20].iter().all(|x| *x != k) {
                assert_eq!(seq_coord_calls.mod_calls(k), Some(&[0u8][..]));
            }
        }
        assert_eq!(seq_coord_calls.mod_calls(33), None); // Out of bounds

        // Test collapse_mod_calls
        let collapsed = seq_coord_calls.collapse_mod_calls();
        assert_eq!(collapsed.len(), 33);
        assert!(collapsed.get(12).expect("in range"));
        assert!(collapsed.get(13).expect("in range"));
        assert!(collapsed.get(16).expect("in range"));
        assert!(collapsed.get(19).expect("in range"));
        assert!(collapsed.get(20).expect("in range"));

        // Test collapsed_mod_calls for positions without modifications
        for k in 0..33 {
            if [12, 13, 16, 19, 20].iter().all(|x| *x != k) {
                assert!(!collapsed.get(k).expect("in range"));
            }
        }
        assert!(collapsed.get(33).is_none()); // Out of bounds

        Ok(())
    }

    #[test]
    fn seq_coord_calls_multiple_mod_types() -> Result<(), Error> {
        // Test with multiple modification types (both T+ and C+m)
        let base_mods = BaseMods {
            base_mods: vec![
                BaseMod {
                    modified_base: b'T',
                    strand: '+',
                    modification_type: 'T',
                    ranges: Ranges {
                        annotations: vec![FiberAnnotation {
                            start: 0,
                            end: 1,
                            length: 1,
                            qual: 100,
                            reference_start: Some(0),
                            reference_end: Some(0),
                            reference_length: Some(0),
                            extra_columns: None,
                        }],
                        seq_len: 5,
                        reverse: false,
                    },
                    record_is_reverse: false,
                },
                BaseMod {
                    modified_base: b'C',
                    strand: '+',
                    modification_type: 'm',
                    ranges: Ranges {
                        annotations: vec![FiberAnnotation {
                            start: 2,
                            end: 3,
                            length: 1,
                            qual: 200,
                            reference_start: None,
                            reference_end: None,
                            reference_length: None,
                            extra_columns: None,
                        }],
                        seq_len: 5,
                        reverse: false,
                    },
                    record_is_reverse: false,
                },
            ],
        };

        let seq_coord_calls = SeqCoordCalls::try_from(&base_mods)?;

        // Test mod_types - should have 2 types
        let mod_types = seq_coord_calls.mod_types();
        assert_eq!(mod_types.len(), 2);
        assert_eq!(
            mod_types.first().expect("has 2 elements").0,
            ModChar::new('T')
        );
        assert!(mod_types.first().expect("has 2 elements").1);
        assert_eq!(
            mod_types.get(1).expect("has 2 elements").0,
            ModChar::new('m')
        );
        assert!(mod_types.get(1).expect("has 2 elements").1);

        // Test mod_calls - each position should have 2 values (one for each mod type)
        assert_eq!(seq_coord_calls.mod_calls(0), Some(&[100u8, 0u8][..])); // T modified, C not
        assert_eq!(seq_coord_calls.mod_calls(1), Some(&[0u8, 0u8][..])); // Neither modified
        assert_eq!(seq_coord_calls.mod_calls(2), Some(&[0u8, 200u8][..])); // C modified, T not
        assert_eq!(seq_coord_calls.mod_calls(3), Some(&[0u8, 0u8][..])); // Neither modified
        assert_eq!(seq_coord_calls.mod_calls(4), Some(&[0u8, 0u8][..])); // Neither modified

        // Test collapse_mod_calls
        let collapsed = seq_coord_calls.collapse_mod_calls();
        assert_eq!(collapsed, vec![true, false, true, false, false]);

        Ok(())
    }

    #[test]
    fn seq_coord_calls_negative_strand() -> Result<(), Error> {
        // Test with negative strand
        let base_mods = BaseMods {
            base_mods: vec![BaseMod {
                modified_base: b'G',
                strand: '-',
                modification_type: 'h',
                ranges: Ranges {
                    annotations: vec![FiberAnnotation {
                        start: 1,
                        end: 2,
                        length: 1,
                        qual: 50,
                        reference_start: Some(1),
                        reference_end: Some(1),
                        reference_length: Some(0),
                        extra_columns: None,
                    }],
                    seq_len: 3,
                    reverse: false,
                },
                record_is_reverse: false,
            }],
        };

        let seq_coord_calls = SeqCoordCalls::try_from(&base_mods)?;

        // Test mod_types
        let mod_types = seq_coord_calls.mod_types();
        assert_eq!(mod_types.len(), 1);
        assert_eq!(
            mod_types.first().expect("has 1 element").0,
            ModChar::new('h')
        );
        assert!(!mod_types.first().expect("has 1 element").1); // - strand

        Ok(())
    }

    #[test]
    #[should_panic(expected = "InvalidSeqLength")]
    fn seq_coord_calls_empty_basemods() {
        // Test with empty BaseMods - should fail
        let base_mods = BaseMods { base_mods: vec![] };
        let _: SeqCoordCalls = SeqCoordCalls::try_from(&base_mods).unwrap();
    }

    #[test]
    #[should_panic(expected = "NotImplemented")]
    fn seq_coord_calls_mismatched_seq_len() {
        // Test with mismatched seq_len in different BaseMod entries
        let base_mods = BaseMods {
            base_mods: vec![
                BaseMod {
                    modified_base: b'T',
                    strand: '+',
                    modification_type: 'T',
                    ranges: Ranges {
                        annotations: vec![],
                        seq_len: 5,
                        reverse: false,
                    },
                    record_is_reverse: false,
                },
                BaseMod {
                    modified_base: b'C',
                    strand: '+',
                    modification_type: 'm',
                    ranges: Ranges {
                        annotations: vec![],
                        seq_len: 10, // Different seq_len!
                        reverse: false,
                    },
                    record_is_reverse: false,
                },
            ],
        };
        let _: SeqCoordCalls = SeqCoordCalls::try_from(&base_mods).unwrap();
    }

    #[test]
    #[should_panic(expected = "InvalidState")]
    fn seq_coord_calls_invalid_strand() {
        // Test with invalid strand character
        let base_mods = BaseMods {
            base_mods: vec![BaseMod {
                modified_base: b'T',
                strand: 'x', // Invalid strand
                modification_type: 'T',
                ranges: Ranges {
                    annotations: vec![],
                    seq_len: 5,
                    reverse: false,
                },
                record_is_reverse: false,
            }],
        };
        let _: SeqCoordCalls = SeqCoordCalls::try_from(&base_mods).unwrap();
    }

    #[test]
    #[should_panic(expected = "NotImplemented")]
    fn seq_coord_calls_non_per_base_annotation() {
        // Test with multi-base annotation (not per-base)
        let base_mods = BaseMods {
            base_mods: vec![BaseMod {
                modified_base: b'T',
                strand: '+',
                modification_type: 'T',
                ranges: Ranges {
                    annotations: vec![FiberAnnotation {
                        start: 0,
                        end: 3, // Multi-base annotation (length 3, not 1)
                        length: 3,
                        qual: 100,
                        reference_start: Some(0),
                        reference_end: Some(2),
                        reference_length: Some(2),
                        extra_columns: None,
                    }],
                    seq_len: 5,
                    reverse: false,
                },
                record_is_reverse: false,
            }],
        };
        let _: SeqCoordCalls = SeqCoordCalls::try_from(&base_mods).unwrap();
    }

    #[test]
    #[should_panic(expected = "NotImplemented")]
    fn seq_coord_calls_start_beyond_seq_len() {
        // Test with start coordinate beyond seq_len
        let base_mods = BaseMods {
            base_mods: vec![BaseMod {
                modified_base: b'T',
                strand: '+',
                modification_type: 'T',
                ranges: Ranges {
                    annotations: vec![FiberAnnotation {
                        start: 10, // Start beyond seq_len (5)
                        end: 11,
                        length: 1,
                        qual: 100,
                        reference_start: Some(10),
                        reference_end: Some(10),
                        reference_length: Some(0),
                        extra_columns: None,
                    }],
                    seq_len: 5,
                    reverse: false,
                },
                record_is_reverse: false,
            }],
        };
        let _: SeqCoordCalls = SeqCoordCalls::try_from(&base_mods).unwrap();
    }
}
