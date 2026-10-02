//! Test-only compatibility helpers for the existing ANSI viewport goldens.

use super::ansi_screen::{AnsiColor, AnsiScreen, AnsiStyle};
use std::fmt::Write as _;

/// Serializes a parsed screen using the RGB palette already stored in the goldens.
pub(super) fn screen_as_ansi(screen: &AnsiScreen) -> String {
    let mut output = String::new();
    for row in 0..screen.rows() {
        write!(output, "\x1b[{};1H", row.saturating_add(1)).expect("writing to String cannot fail");
        let mut current_style = None;
        for column in 0..screen.cols() {
            let cell = screen
                .cell(row, column)
                .expect("screen contains every in-bounds cell");
            if current_style != Some(cell.style) {
                write_style(cell.style, &mut output);
                current_style = Some(cell.style);
            }
            output.push(cell.character);
        }
        output.push_str("\x1b[0m");
    }
    output
}

/// Parses an existing RGB golden through the production palette-index parser.
pub(super) fn parse_ansi(frame: &str, cols: u16, rows: u16) -> Result<AnsiScreen, String> {
    let normalized = normalize_rgb_sgr(frame)?;
    AnsiScreen::parse(&normalized, cols, rows)
}

/// Rewrites known legacy RGB groups only inside complete SGR sequences.
fn normalize_rgb_sgr(frame: &str) -> Result<String, String> {
    let mut output = String::new();
    let mut offset = 0usize;
    while let Some(relative_start) = frame
        .get(offset..)
        .expect("offset remains on a character boundary")
        .find("\x1b[")
    {
        let start = offset.saturating_add(relative_start);
        output.push_str(
            frame
                .get(offset..start)
                .expect("escape start is on a character boundary"),
        );
        let parameters_start = start.saturating_add(2);
        let Some((relative_command, command)) = frame
            .get(parameters_start..)
            .expect("CSI parameters start on a character boundary")
            .char_indices()
            .find(|(_index, character)| matches!(character, '@'..='~'))
        else {
            output.push_str(
                frame
                    .get(start..)
                    .expect("escape start is on a character boundary"),
            );
            offset = frame.len();
            break;
        };
        let command_start = parameters_start.saturating_add(relative_command);
        let sequence_end = command_start.saturating_add(command.len_utf8());
        if command == 'm' {
            output.push_str("\x1b[");
            output.push_str(&normalize_sgr_parameters(
                frame
                    .get(parameters_start..command_start)
                    .expect("SGR parameters are on character boundaries"),
            )?);
            output.push('m');
        } else {
            output.push_str(
                frame
                    .get(start..sequence_end)
                    .expect("CSI sequence ends on a character boundary"),
            );
        }
        offset = sequence_end;
    }
    output.push_str(
        frame
            .get(offset..)
            .expect("offset remains on a character boundary"),
    );
    Ok(output)
}

/// Rewrites complete known `38;2;r;g;b` and `48;2;r;g;b` parameter groups.
fn normalize_sgr_parameters(parameters: &str) -> Result<String, String> {
    let parts = parameters.split(';').collect::<Vec<_>>();
    let mut normalized = Vec::new();
    let mut position = 0usize;
    while let Some(&part) = parts.get(position) {
        if matches!(part, "38" | "48") {
            let rgb = parts
                .get(position.saturating_add(1)..position.saturating_add(5))
                .ok_or_else(|| String::from("ANSI golden contains incomplete RGB colour"))?;
            let ["2", red, green, blue] = rgb else {
                return Err(String::from(
                    "ANSI golden contains unsupported extended colour",
                ));
            };
            normalized.push(legacy_palette_code(part, red, green, blue).ok_or_else(|| {
                format!("ANSI golden contains unknown RGB colour {red};{green};{blue}")
            })?);
            position = position.saturating_add(5);
        } else {
            normalized.push(part);
            position = position.saturating_add(1);
        }
    }
    Ok(normalized.join(";"))
}

/// Returns the palette SGR corresponding to one legacy golden RGB group.
fn legacy_palette_code(code: &str, red: &str, green: &str, blue: &str) -> Option<&'static str> {
    match (code, red, green, blue) {
        ("38", "181", "189", "104") => Some("32"),
        ("38", "240", "198", "116") => Some("33"),
        ("38", "138", "190", "183") => Some("36"),
        ("38", "102", "102", "102") => Some("90"),
        ("38", "234", "234", "234") => Some("97"),
        ("48", "129", "162", "190") => Some("44"),
        _ => None,
    }
}

/// Writes one canonical style transition.
fn write_style(style: AnsiStyle, output: &mut String) {
    output.push_str("\x1b[0m");
    let mut codes = Vec::new();
    if style.bold {
        codes.push(String::from("1"));
    }
    if style.faint {
        codes.push(String::from("2"));
    }
    if style.underline {
        codes.push(String::from("4"));
    }
    if style.inverse {
        codes.push(String::from("7"));
    }
    if let AnsiColor::Palette(index) = style.foreground {
        codes.push(String::from(canonical_color(index, false)));
    }
    if let AnsiColor::Palette(index) = style.background {
        codes.push(String::from(canonical_color(index, true)));
    }
    if !codes.is_empty() {
        write!(output, "\x1b[{}m", codes.join(";")).expect("writing to String cannot fail");
    }
}

/// Preserves the RGB palette previously used by the ANSI goldens.
fn canonical_color(index: u8, background: bool) -> &'static str {
    match (background, index) {
        (false, 2) => "38;2;181;189;104",
        (false, 3) => "38;2;240;198;116",
        (false, 6) => "38;2;138;190;183",
        (false, 8) => "38;2;102;102;102",
        (false, 15) => "38;2;234;234;234",
        (true, 4) => "48;2;129;162;190",
        _ => unreachable!("ANSI parser only produces viewer palette colours"),
    }
}

mod tests {
    use super::*;

    #[test]
    fn normalization_changes_only_known_complete_sgr_colours() {
        let combined = "\x1b[1;38;2;181;189;104;48;2;129;162;190mX";
        assert_eq!(
            normalize_rgb_sgr(combined),
            Ok(String::from("\x1b[1;32;44mX"))
        );

        let literal = "\x1b[1;1H38;2;181;189;104";
        assert_eq!(normalize_rgb_sgr(literal), Ok(String::from(literal)));

        let malformed_cursor = "\x1b[1;38;2;181;189;104HX";
        assert_eq!(
            normalize_rgb_sgr(malformed_cursor),
            Ok(String::from(malformed_cursor))
        );
        assert_eq!(
            parse_ansi(malformed_cursor, 200, 2).unwrap_err(),
            "viewer cursor position has extra parameters: \"1;38;2;181;189;104\""
        );

        let non_sgr_final = "\x1b[1@;38;2;181;189;104m";
        assert_eq!(
            normalize_rgb_sgr(non_sgr_final),
            Ok(String::from(non_sgr_final))
        );
    }

    #[test]
    fn normalization_rejects_unknown_or_incomplete_rgb() {
        assert_eq!(
            normalize_rgb_sgr("\x1b[38;2;1;2;3mX").unwrap_err(),
            "ANSI golden contains unknown RGB colour 1;2;3"
        );
        assert_eq!(
            normalize_rgb_sgr("\x1b[38;2;181mX").unwrap_err(),
            "ANSI golden contains incomplete RGB colour"
        );
    }
}
