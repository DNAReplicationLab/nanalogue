//! Integration tests for `FilterByRefCoords` functionality
//! Tests focus on filtering reads by reference coordinates

use indoc::indoc;
use nanalogue_core::{CurrRead, FilterByRefCoords, read_utils::AlignAndModData};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ranges_filter_by_ref_coords_via_curr_read() {
        // Input JSON with modification data spanning reference positions 10-30
        let input_json = indoc! {r#"
        {
          "alignment_type": "primary_forward",
          "alignment": {
            "start": 5,
            "end": 35
          },
          "mod_table": [
            {
              "base": "T",
              "mod_code": "T",
              "data": [
                [0, 10, 200],
                [1, 15, 180],
                [2, 20, 220],
                [3, 25, 190],
                [4, 30, 210]
              ]
            },
            {
              "base": "C",
              "is_strand_plus": true,
              "mod_code": "m",
              "implicit": false,
              "data": [
                [5, 12, 150],
                [6, 18, 160],
                [7, 28, 170]
              ]
            }
          ],
          "seq_len": 8
        }"#};

        // Expected JSON after filtering for range 15-25
        let expected_json = indoc! {r#"
        {
          "alignment_type": "primary_forward",
          "alignment": {
            "start": 5,
            "end": 35
          },
          "mod_table": [
            {
              "base": "T",
              "mod_code": "T",
              "data": [
                [1, 15, 180],
                [2, 20, 220]
              ]
            },
            {
              "base": "C",
              "is_strand_plus": true,
              "mod_code": "m",
              "data": [
                [6, 18, 160]
              ]
            }
          ],
          "seq_len": 8
        }"#};

        // Deserialize input, apply filter, and compare with expected
        let mut curr_read: CurrRead<AlignAndModData> = serde_json::from_str(input_json).unwrap();
        curr_read.filter_by_ref_pos(15, 25);

        let expected_curr_read: CurrRead<AlignAndModData> =
            serde_json::from_str(expected_json).unwrap();
        assert_eq!(curr_read, expected_curr_read);
    }

    #[test]
    fn ranges_filter_no_overlap_1() {
        // Input JSON with modification data
        let input_json = indoc! {r#"
        {
          "alignment_type": "primary_forward",
          "alignment": {
            "start": 10,
            "end": 30
          },
          "mod_table": [
            {
              "base": "T",
              "mod_code": "T",
              "data": [
                [0, 15, 200],
                [1, 20, 180],
                [2, 25, 220]
              ]
            }
          ],
          "seq_len": 3
        }"#};

        // Expected JSON after filtering with no overlap (26-60)
        let expected_json = indoc! {r#"
        {
          "alignment_type": "primary_forward",
          "alignment": {
            "start": 10,
            "end": 30
          },
          "mod_table": [
            {
              "base": "T",
              "mod_code": "T",
              "data": [
              ]
            }
          ],
          "seq_len": 3
        }"#};

        // Deserialize input, apply filter, and compare with expected
        let mut curr_read: CurrRead<AlignAndModData> = serde_json::from_str(input_json).unwrap();
        curr_read.filter_by_ref_pos(26, 60);

        let expected_curr_read: CurrRead<AlignAndModData> =
            serde_json::from_str(expected_json).unwrap();
        assert_eq!(curr_read, expected_curr_read);
    }

    #[test]
    fn ranges_filter_no_overlap_2() {
        // Input JSON with modification data (same as test_1)
        let input_json = indoc! {r#"
        {
          "alignment_type": "primary_forward",
          "alignment": {
            "start": 10,
            "end": 30
          },
          "mod_table": [
            {
              "base": "T",
              "mod_code": "T",
              "data": [
                [0, 15, 200],
                [1, 20, 180],
                [2, 25, 220]
              ]
            }
          ],
          "seq_len": 3
        }"#};

        // Expected JSON after filtering with no overlap (0-15)
        let expected_json = indoc! {r#"
        {
          "alignment_type": "primary_forward",
          "alignment": {
            "start": 10,
            "end": 30
          },
          "mod_table": [
            {
              "base": "T",
              "mod_code": "T",
              "data": []
            }
          ],
          "seq_len": 3
        }"#};

        // Deserialize input, apply filter, and compare with expected
        let mut curr_read: CurrRead<AlignAndModData> = serde_json::from_str(input_json).unwrap();
        curr_read.filter_by_ref_pos(0, 15);

        let expected_curr_read: CurrRead<AlignAndModData> =
            serde_json::from_str(expected_json).unwrap();
        assert_eq!(curr_read, expected_curr_read);
    }

    #[test]
    fn ranges_filter_partial_overlap() {
        // Input JSON with modification data
        let input_json = indoc! {r#"
        {
          "alignment_type": "primary_forward",
          "alignment": {
            "start": 10,
            "end": 40
          },
          "mod_table": [
            {
              "base": "T",
              "mod_code": "T",
              "data": [
                [0, 10, 200],
                [1, 15, 180],
                [2, 20, 220],
                [3, 25, 190],
                [4, 30, 210],
                [5, 35, 170]
              ]
            }
          ],
          "seq_len": 6
        }"#};

        // Expected JSON after filtering for range 20-30
        let expected_json = indoc! {r#"
        {
          "alignment_type": "primary_forward",
          "alignment": {
            "start": 10,
            "end": 40
          },
          "mod_table": [
            {
              "base": "T",
              "mod_code": "T",
              "data": [
                [2, 20, 220],
                [3, 25, 190]
              ]
            }
          ],
          "seq_len": 6
        }"#};

        // Deserialize input, apply filter, and compare with expected
        let mut curr_read: CurrRead<AlignAndModData> = serde_json::from_str(input_json).unwrap();
        curr_read.filter_by_ref_pos(20, 30);

        let expected_curr_read: CurrRead<AlignAndModData> =
            serde_json::from_str(expected_json).unwrap();
        assert_eq!(curr_read, expected_curr_read);
    }

    #[test]
    fn ranges_filter_reverse_strand() {
        // Input JSON with reverse strand modification data
        let input_json = indoc! {r#"
        {
          "alignment_type": "primary_reverse",
          "alignment": {
            "start": 10,
            "end": 40
          },
          "mod_table": [
            {
              "base": "T",
              "mod_code": "T",
              "data": [
                [0, 15, 170],
                [1, 20, 210],
                [2, 25, 190],
                [3, 30, 220],
                [4, 35, 180]
              ]
            }
          ],
          "seq_len": 5
        }"#};

        // Expected JSON after filtering for range 20-30
        let expected_json = indoc! {r#"
        {
          "alignment_type": "primary_reverse",
          "alignment": {
            "start": 10,
            "end": 40
          },
          "mod_table": [
            {
              "base": "T",
              "mod_code": "T",
              "data": [
                [1, 20, 210],
                [2, 25, 190]
              ]
            }
          ],
          "seq_len": 5
        }"#};

        // Deserialize input, apply filter, and compare with expected
        let mut curr_read: CurrRead<AlignAndModData> = serde_json::from_str(input_json).unwrap();
        curr_read.filter_by_ref_pos(20, 30);

        let expected_curr_read: CurrRead<AlignAndModData> =
            serde_json::from_str(expected_json).unwrap();
        assert_eq!(curr_read, expected_curr_read);
    }
}
