//! MM-tag parsing helpers.

use crate::{
    Error, ModChar,
    constants::shared::{MAX_MM_GAP, MAX_MM_TAG_LENGTH, MAX_MOD_TYPES},
};
use std::{collections::HashSet, str::FromStr as _};

/// Base, strand, longest numeric modification code, and optional `?` or `.` suffix.
const MAX_MM_HEADER_LENGTH: usize = ModChar::MAX_NUMERIC_CODE_LENGTH.saturating_add(3);

/// Parsed representation of a single MM-tag group.
#[derive(Debug)]
#[non_exhaustive]
pub struct ParsedMmGroup {
    /// Modified base code.
    pub mod_base: u8,
    /// Modification strand.
    pub mod_strand: char,
    /// Parsed modification type.
    pub modification_type: ModChar,
    /// Whether omitted positions are implicitly unmodified.
    pub is_implicit: bool,
    /// Distances between successive modified bases.
    pub mod_dists: Vec<u32>,
}

/// Parses the comma-separated distances that follow an MM group's header.
///
/// Each field must be a non-empty sequence of ASCII decimal digits whose value does not
/// exceed [`MAX_MM_GAP`].
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the digit match bounds subtraction, and MAX_MM_GAP bounds decimal accumulation"
)]
fn parse_mm_distances(distances: &str) -> Result<Vec<u32>, Error> {
    if distances.len()
        >= usize::try_from(MAX_MM_TAG_LENGTH).expect("no error on 32-bit platforms and above")
    {
        return Err(Error::InvalidModCoords(
            "MM distance list is too long to process".to_owned(),
        ));
    }
    // Every field holds at least one byte, and fields are separated by one comma.
    let mut mod_dists = Vec::with_capacity(distances.len().div_ceil(2).max(1));
    let mut field_start = 0usize;
    let mut value = 0u32;
    // A synthetic trailing comma finalizes the last field through the same path as all others.
    for (index, byte) in distances.bytes().chain(std::iter::once(b',')).enumerate() {
        match byte {
            b',' => {
                if field_start == index {
                    return Err(Error::InvalidModCoords(
                        "invalid MM distance ``: expected ASCII decimal digits".to_owned(),
                    ));
                }
                mod_dists.push(value);
                field_start = index.saturating_add(1);
                value = 0;
            }
            b'0'..=b'9' => {
                let digit = u32::from(byte - b'0');
                value = value * 10 + digit;
                if value > MAX_MM_GAP {
                    let field = distances
                        .get(field_start..)
                        .unwrap_or_default()
                        .split(',')
                        .next()
                        .unwrap_or_default();
                    return Err(Error::InvalidModCoords(format!(
                        "invalid MM distance `{field}`: value exceeds maximum gap {MAX_MM_GAP}"
                    )));
                }
            }
            _ => {
                let field = distances
                    .get(field_start..)
                    .unwrap_or_default()
                    .split(',')
                    .next()
                    .unwrap_or_default();
                return Err(Error::InvalidModCoords(format!(
                    "invalid MM distance `{field}`: expected ASCII decimal digits"
                )));
            }
        }
    }
    Ok(mod_dists)
}

/// Parse semicolon-delimited MM-tag text into groups.
///
/// This parser intentionally requires every MM group, including the final
/// one, to be terminated with `;` as specified by the SAM tag format.
/// Unterminated trailing groups are rejected rather than accepted leniently.
///
/// # Errors
/// Returns an error when any MM-tag group is malformed.
#[expect(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::string_slice,
    clippy::missing_asserts_for_indexing,
    reason = "bounds are checked before slicing and index arithmetic is tightly controlled"
)]
#[expect(
    clippy::too_many_lines,
    reason = "keeping validation in the single-pass MM group parser makes its state explicit"
)]
#[expect(
    clippy::missing_panics_doc,
    reason = "u32 -> usize conversion will not fail"
)]
pub fn mm_groups(group: &str) -> Result<Vec<ParsedMmGroup>, Error> {
    let mut groups = Vec::<ParsedMmGroup>::new();
    let mut seen_combinations = HashSet::new();
    let mut group_start = 0usize;
    let group_bytes = group.as_bytes();

    // Checks `mm_text` is not too long.
    // Other checks such as record size checks may stop this from ever triggering.
    if group_bytes.len()
        > usize::try_from(MAX_MM_TAG_LENGTH).expect("no error on 32-bit platforms and above")
    {
        return Err(Error::InvalidState(
            "MM tag is too long to process".to_owned(),
        ));
    }

    for (index, byte) in group_bytes.iter().copied().enumerate() {
        if byte != b';' {
            continue;
        }

        // `group` is valid UTF-8 and `;` is ASCII, so these are character boundaries.
        let raw_group = &group[group_start..index];
        if raw_group.len()
            >= usize::try_from(MAX_MM_TAG_LENGTH).expect("no error on 32-bit platforms and above")
        {
            return Err(Error::InvalidModCoords(
                "MM group is too long to process".to_owned(),
            ));
        }
        let (header, distances) = match raw_group.split_once(',') {
            Some((header, distances)) => (header, Some(distances)),
            None => (raw_group, None),
        };
        if header.len() > MAX_MM_HEADER_LENGTH {
            return Err(Error::InvalidModType(format!(
                "MM group header exceeds {MAX_MM_HEADER_LENGTH} bytes"
            )));
        }
        let header_bytes = header.as_bytes();
        // smallest valid length is 3 e.g. something like "C+m"
        if header_bytes.len() < 3 {
            return Err(Error::InvalidModType(
                "malformed MM group encountered while parsing MM tag".to_owned(),
            ));
        }
        let mod_base = match header_bytes[0] {
            v @ (b'A' | b'C' | b'G' | b'T' | b'U' | b'N') => v,
            v => {
                return Err(Error::InvalidBase(format!(
                    "invalid MM base `{}`",
                    char::from(v)
                )));
            }
        };
        let mod_strand = match &header_bytes[1] {
            &b'+' => '+',
            &b'-' => '-',
            v => return Err(Error::InvalidModType(format!("invalid MM strand `{v}`"))),
        };

        let (modification_type, is_implicit) = {
            let mod_type_and_implicit_flag = &header[2..];
            let n = mod_type_and_implicit_flag.len();
            if mod_type_and_implicit_flag.ends_with('?') {
                (
                    ModChar::from_str(&mod_type_and_implicit_flag[0..n - 1])?,
                    false,
                )
            } else if mod_type_and_implicit_flag.ends_with('.') {
                (
                    ModChar::from_str(&mod_type_and_implicit_flag[0..n - 1])?,
                    true,
                )
            } else {
                (ModChar::from_str(&mod_type_and_implicit_flag[0..n])?, true)
            }
        };

        // Check for duplicate strand, modification_type combinations
        // NOTE: If MM data has mods like C+m, A+m, this will fire and
        // reject the data. Whereas, if we do this check after rejecting
        // a base as specified by the user, this data would be o.k.
        // We let this be as C+m, A+m is very unlikely i.e. an experimental
        // scenario where both cytosines and adenines are replaced by a 5-methyl cytosine.
        // Two bases replaced by the _same_ kind of modification is very unlikely.
        if !seen_combinations.insert((mod_strand, modification_type)) {
            return Err(Error::InvalidDuplicates(format!(
                "Duplicate strand '{mod_strand}' and modification_type '{modification_type}' combination found",
            )));
        }

        // Check we don't process too many mod types
        if seen_combinations.len() > usize::from(MAX_MOD_TYPES) {
            return Err(Error::InvalidState(format!(
                "max types of mods exceeded {MAX_MOD_TYPES}"
            )));
        }

        let mod_dists = match distances {
            Some(distance_text) => parse_mm_distances(distance_text)?,
            None => Vec::new(),
        };

        assert!(
            mod_dists.len()
                < usize::try_from(MAX_MM_TAG_LENGTH)
                    .expect("no error on 32-bit platforms and above"),
            "no error as we have checked mm text does not exceed this length so no way the array exceeds the length"
        );

        assert!(
            matches!(mod_base, b'A' | b'C' | b'G' | b'T' | b'U' | b'N'),
            "MM base was validated before constructing the parsed group"
        );
        assert!(
            matches!(mod_strand, '+' | '-'),
            "MM strand was validated before constructing the parsed group"
        );
        groups.push(ParsedMmGroup {
            mod_base,
            mod_strand,
            modification_type,
            is_implicit,
            mod_dists,
        });
        group_start = index + 1;
    }

    if group_start != group.len() {
        return Err(Error::InvalidModType(
            "MM tag must end with `;` terminators for each group".to_owned(),
        ));
    }

    Ok(groups)
}

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use super::*;

    #[test]
    fn mm_group_headers_respect_the_derived_length_limit() {
        assert_eq!(MAX_MM_HEADER_LENGTH, 10);
        let maximum = mm_groups("C+1114111?,0;").expect("maximum header should parse");
        assert_eq!(maximum.len(), 1);
        let overlong = mm_groups("C+12345678?,0;").expect_err("overlong header should fail");
        assert!(matches!(overlong, Error::InvalidModType(_)));
    }

    #[test]
    fn mm_distance_list_must_be_shorter_than_the_tag_limit() {
        let limit = usize::try_from(MAX_MM_TAG_LENGTH).expect("supported platform");
        let distances = "0".repeat(limit);
        let error = parse_mm_distances(&distances).expect_err("distance list at limit");
        assert!(matches!(error, Error::InvalidModCoords(_)));
    }

    #[test]
    fn mm_distances_accept_only_ascii_decimal_values_within_the_gap_limit() {
        let maximum = MAX_MM_GAP.to_string();
        let excessive_value = u64::from(MAX_MM_GAP)
            .checked_add(1)
            .expect("MM gap limit fits below u64::MAX");
        let excessive = excessive_value.to_string();
        for (entry, expected) in [
            ("0".to_owned(), 0),
            ("7".to_owned(), 7),
            ("42".to_owned(), 42),
            (maximum.clone(), MAX_MM_GAP),
            ("0000000000012".to_owned(), 12),
        ] {
            assert_eq!(
                parse_mm_distances(&entry).expect("decimal u32 value should parse"),
                vec![expected],
                "entry `{entry}`"
            );
        }

        for entry in [
            excessive,
            "4294967295".to_owned(),
            "4294967296".to_owned(),
            "+5".to_owned(),
            "+".to_owned(),
            String::new(),
            "-1".to_owned(),
            "1a".to_owned(),
            " 1".to_owned(),
            "\u{661}".to_owned(),
        ] {
            let message = parse_mm_distances(&entry)
                .expect_err("invalid or excessive decimal field should fail")
                .to_string();
            assert!(
                message.contains(&format!("invalid MM distance `{entry}`")),
                "entry `{entry}` gave `{message}`"
            );
        }

        assert_eq!(
            parse_mm_distances(&format!("{maximum},1,{maximum}"))
                .expect("the gap cap applies to each field, not their sum"),
            [MAX_MM_GAP, 1, MAX_MM_GAP]
        );
        for group in [
            format!("{excessive_value},1"),
            format!("1,{excessive_value},2"),
        ] {
            let _error = parse_mm_distances(&group).expect_err("every gap must respect the cap");
        }
    }

    #[test]
    fn mm_distances_preserve_order_and_reject_invalid_fields() {
        assert_eq!(
            parse_mm_distances("3,0,12,4").expect("distances should parse"),
            vec![3, 0, 12, 4]
        );
        let _signed = parse_mm_distances("1,+4,2").expect_err("signed field");
        let _empty_middle = parse_mm_distances("1,,2").expect_err("empty field");
        let _empty_list = parse_mm_distances("").expect_err("empty distance list");
        let _trailing_comma = parse_mm_distances("3,0,12,4,").expect_err("trailing comma");
    }

    #[test]
    fn mm_groups_allows_empty_input() {
        let groups = mm_groups("").expect("empty MM tag input should parse successfully");
        assert!(
            groups.is_empty(),
            "empty MM tag input should produce no groups"
        );
    }

    #[test]
    fn mm_groups_rejects_just_semicolon() {
        let result = mm_groups(";");
        assert!(result.is_err(), "just a semicolon should fail");
    }

    #[test]
    fn mm_groups_rejects_missing_semicolon_terminator() {
        let result = mm_groups("C+m,0,1");
        assert!(result.is_err(), "unterminated MM group should fail");
    }

    #[test]
    fn mm_groups_parses_question_flag_without_distances() {
        let groups = mm_groups("C+m?;").expect("MM group with `?` and no distances should parse");
        assert_eq!(groups.len(), 1);
        let group = groups.first().expect("one MM group should be present");
        assert_eq!(group.mod_base, b'C');
        assert_eq!(group.mod_strand, '+');
        assert_eq!(group.modification_type.val(), 'm');
        assert!(!group.is_implicit);
        assert_eq!(group.mod_dists, Vec::<u32>::new());
    }

    #[test]
    fn mm_groups_parses_period_flag_without_distances() {
        let groups = mm_groups("C+m.;").expect("MM group with `.` and no distances should parse");
        assert_eq!(groups.len(), 1);
        let group = groups.first().expect("one MM group should be present");
        assert_eq!(group.mod_base, b'C');
        assert_eq!(group.mod_strand, '+');
        assert_eq!(group.modification_type.val(), 'm');
        assert!(group.is_implicit);
        assert_eq!(group.mod_dists, Vec::<u32>::new());
    }

    #[test]
    fn mm_groups_parses_no_flag_without_distances() {
        let groups = mm_groups("C+m;").expect("MM group with no distances should parse");
        assert_eq!(groups.len(), 1);
        let group = groups.first().expect("one MM group should be present");
        assert_eq!(group.mod_base, b'C');
        assert_eq!(group.mod_strand, '+');
        assert_eq!(group.modification_type.val(), 'm');
        assert!(group.is_implicit);
        assert_eq!(group.mod_dists, Vec::<u32>::new());
    }

    #[test]
    fn mm_groups_rejects_just_semicolon_after_one_valid_mod() {
        let result = mm_groups("C+m?;;");
        assert!(
            result.is_err(),
            "just a semicolon present after a valid mod should fail"
        );
    }

    #[test]
    fn mm_groups_fails_question_flag_without_distances_but_with_comma() {
        let result = mm_groups("C+m?,;");
        assert!(
            result.is_err(),
            "question flag without distances but with comma should fail"
        );
    }

    #[test]
    fn mm_groups_parses_multiple_groups() {
        let groups =
            mm_groups("C+m,0,1;A-a?,2;").expect("multiple MM groups should parse successfully");
        assert_eq!(groups.len(), 2);

        let first_group = groups.first().expect("first MM group should be present");
        assert_eq!(first_group.mod_base, b'C');
        assert_eq!(first_group.mod_strand, '+');
        assert_eq!(first_group.modification_type.val(), 'm');
        assert!(first_group.is_implicit);
        assert_eq!(first_group.mod_dists, vec![0, 1]);

        let second_group = groups.get(1).expect("second MM group should be present");
        assert_eq!(second_group.mod_base, b'A');
        assert_eq!(second_group.mod_strand, '-');
        assert_eq!(second_group.modification_type.val(), 'a');
        assert!(!second_group.is_implicit);
        assert_eq!(second_group.mod_dists, vec![2]);
    }

    #[test]
    fn mm_groups_rejects_too_many_mods() {
        let mut long_mm_tag = String::with_capacity(800);
        for k in 0..101 {
            let tag = format!("C+{k},0;");
            long_mm_tag.push_str(&tag);
        }
        let result = mm_groups(&long_mm_tag);
        assert!(result.is_err(), "too many mods should fail");
    }

    #[test]
    fn mm_groups_rejects_different_base_same_mod_type_strand() {
        let result = mm_groups("C+m,0;A+m,0;");
        assert!(result.is_err(), "different base same mod type should fail");
    }

    #[test]
    fn mm_groups_rejects_same_base_same_mod_type_strand() {
        let result = mm_groups("C+m,0;C+m,0;");
        assert!(result.is_err(), "same base same mod type should fail");
    }

    #[test]
    fn mm_groups_rejects_really_long_mod_tag() {
        let result = mm_groups("C+123456789,0;");
        assert!(result.is_err(), "really long mod tag should fail");
    }

    #[test]
    fn mm_groups_accepts_char_limit() {
        let result = mm_groups("C+1114111?,0;");
        assert!(
            result.is_ok(),
            "really long mod tag within char limit should pass"
        );
    }

    #[test]
    fn mm_groups_rejects_above_char_limit() {
        let result = mm_groups("C+1114112?,0;");
        assert!(
            result.is_err(),
            "really long mod tag above char limit should fail"
        );
    }
}
