//! MM-tag parsing helpers.

use crate::{Error, ModChar};
use std::str;

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
    pub mod_dists: Vec<usize>,
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
    clippy::string_from_utf8_as_bytes,
    clippy::string_slice,
    reason = "bounds are checked before slicing and index arithmetic is tightly controlled"
)]
pub fn mm_groups(group: &str) -> Result<Vec<ParsedMmGroup>, Error> {
    let mut groups = Vec::<ParsedMmGroup>::new();
    let mut group_start = 0usize;

    for (index, byte) in group.as_bytes().iter().copied().enumerate() {
        if byte != b';' {
            continue;
        }

        let raw_group = &group[group_start..index];
        if raw_group.is_empty() {
            return Err(Error::InvalidModType(
                "empty MM group encountered while parsing MM tag".to_owned(),
            ));
        }

        let mut parsed_mod_base: Option<u8> = None;
        let mut parsed_mod_strand: Option<char> = None;
        let mut parsed_mod_code_end: Option<usize> = None;

        for (local_index, local_byte) in raw_group.as_bytes().iter().copied().enumerate() {
            match local_index {
                0 => {
                    if !matches!(local_byte, b'A' | b'C' | b'G' | b'T' | b'U' | b'N') {
                        return Err(Error::InvalidBase(format!(
                            "invalid MM base `{}`",
                            char::from(local_byte)
                        )));
                    }
                    parsed_mod_base = Some(local_byte);
                }
                1 => {
                    let strand = char::from(local_byte);
                    if !matches!(strand, '+' | '-') {
                        return Err(Error::InvalidModType(format!(
                            "invalid MM strand `{strand}` in group `{raw_group}`"
                        )));
                    }
                    parsed_mod_strand = Some(strand);
                }
                _ => {
                    if !local_byte.is_ascii_alphanumeric() {
                        parsed_mod_code_end = Some(local_index);
                        break;
                    }
                }
            }
        }

        let Some(mod_base) = parsed_mod_base else {
            return Err(Error::InvalidModType("empty MM group".to_owned()));
        };
        let Some(mod_strand) = parsed_mod_strand else {
            return Err(Error::InvalidModType(format!(
                "missing MM strand for group `{raw_group}`"
            )));
        };

        let mod_code_end = parsed_mod_code_end.unwrap_or(raw_group.len());
        if mod_code_end == 2 {
            return Err(Error::EmptyModType(format!(
                "missing MM modification type in group `{raw_group}`"
            )));
        }
        let modification_type: ModChar =
            str::from_utf8(&raw_group.as_bytes()[2..mod_code_end])?.parse()?;

        let punctuation = raw_group.as_bytes().get(mod_code_end).copied();
        let is_implicit = match punctuation {
            Some(b'?') => false,
            // A comma here is not a punctuation flag; it means there was no
            // explicit `.`/`?` marker and the distance list starts immediately.
            Some(b'.' | b',') | None => true,
            Some(other) => {
                return Err(Error::InvalidModType(format!(
                    "invalid MM punctuation `{}` in group `{raw_group}`",
                    char::from(other)
                )));
            }
        };
        let distance_start = mod_code_end + usize::from(matches!(punctuation, Some(b'.' | b'?')));

        let mod_dists = str::from_utf8(&raw_group.as_bytes()[distance_start..])?
            .split(',')
            .map(str::trim)
            .filter(|entry| !entry.is_empty())
            .map(|entry| {
                entry.parse::<usize>().map_err(|err| {
                    Error::InvalidModCoords(format!("invalid MM distance `{entry}`: {err}"))
                })
            })
            .collect::<Result<Vec<_>, Error>>()?;

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
    fn mm_groups_rejects_missing_semicolon_terminator() {
        let result = mm_groups("C+m,0,1");
        assert!(result.is_err(), "unterminated MM group should fail");
    }

    #[test]
    fn mm_groups_rejects_consecutive_semicolons() {
        let result = mm_groups("C+m,0;;");
        assert!(result.is_err(), "consecutive semicolons should fail");
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
}
