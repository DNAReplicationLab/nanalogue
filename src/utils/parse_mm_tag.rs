//! MM-tag parsing helpers.

use crate::{Error, ModChar, constants::shared::MAX_MOD_TYPES};
use std::{collections::HashSet, str, str::FromStr as _};

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
    clippy::missing_panics_doc,
    reason = "u32 -> usize conversion will not fail"
)]
pub fn mm_groups(group: &str) -> Result<Vec<ParsedMmGroup>, Error> {
    let mut groups = Vec::<ParsedMmGroup>::new();
    let mut seen_combinations = HashSet::new();
    let mut group_start = 0usize;
    let group_bytes = group.as_bytes();

    for (index, byte) in group_bytes.iter().copied().enumerate() {
        if byte != b';' {
            continue;
        }

        let raw_group = &group_bytes[group_start..index];
        // smallest valid length is 3 e.g. something like "C+m"
        if raw_group.len() < 3 {
            return Err(Error::InvalidModType(
                "malformed MM group encountered while parsing MM tag".to_owned(),
            ));
        }
        let mod_base = match raw_group[0] {
            v @ (b'A' | b'C' | b'G' | b'T' | b'U' | b'N') => v,
            v => {
                return Err(Error::InvalidBase(format!(
                    "invalid MM base `{}`",
                    char::from(v)
                )));
            }
        };
        let mod_strand = match &raw_group[1] {
            &b'+' => '+',
            &b'-' => '-',
            v => return Err(Error::InvalidModType(format!("invalid MM strand `{v}`"))),
        };

        let (mod_type, is_implicit) = {
            let mod_type_and_implicit_flag_str =
                str::from_utf8(&raw_group[2..raw_group.len().min(11)])?;
            // 2..11 because max char is 1114111 (7) + an optional ?/. (1) + comma (1) = 9
            let mod_type_and_implicit_flag = mod_type_and_implicit_flag_str
                .split(',')
                .take(1)
                .collect::<Vec<_>>()[0];
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
        if !seen_combinations.insert((mod_strand, mod_type)) {
            return Err(Error::InvalidDuplicates(format!(
                "Duplicate strand '{mod_strand}' and modification_type '{mod_type}' combination found",
            )));
        }

        // Check we don't process too many mod types
        if seen_combinations.len() > usize::from(MAX_MOD_TYPES) {
            return Err(Error::InvalidState(format!(
                "max types of mods exceeded {MAX_MOD_TYPES}"
            )));
        }

        let mod_dists = str::from_utf8(raw_group)?
            .split(',')
            .skip(1)
            .take(usize::try_from(u32::MAX).expect("no error on 32-bit platforms and higher"))
            .map(|entry| {
                entry.parse::<u32>().map_err(|err| {
                    Error::InvalidModCoords(format!("invalid MM distance `{entry}`: {err}"))
                })
            })
            .collect::<Result<Vec<_>, Error>>()?;

        if mod_dists.len()
            == usize::try_from(u32::MAX).expect("no error on 32-bit platforms and higher")
        {
            return Err(Error::InvalidState("mod dists too long!".to_owned()));
        }

        groups.push(ParsedMmGroup {
            mod_base,
            mod_strand,
            modification_type: mod_type,
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
mod tests {
    use super::*;

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
        assert!(group.mod_dists.is_empty());
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
        assert!(group.mod_dists.is_empty());
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
        assert!(group.mod_dists.is_empty());
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
