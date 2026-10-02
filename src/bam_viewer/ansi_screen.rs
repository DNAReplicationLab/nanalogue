//! Strict interpretation of the ANSI subset emitted by the BAM viewer.

/// Indexed ANSI colour used by the viewer's generated frames.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum AnsiColor {
    /// Host terminal's default colour.
    Default,
    /// Host terminal palette index from zero through fifteen.
    Palette(u8),
}

/// Supported rendition state for one parsed viewer cell.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[expect(
    clippy::struct_excessive_bools,
    reason = "these are independent SGR attributes in the viewer's restricted ANSI subset"
)]
pub(super) struct AnsiStyle {
    /// Increased intensity.
    pub bold: bool,
    /// Decreased intensity.
    pub faint: bool,
    /// Single underline, used to mark modified bases.
    pub underline: bool,
    /// Reversed foreground and background.
    pub inverse: bool,
    /// Indexed foreground colour.
    pub foreground: AnsiColor,
    /// Indexed background colour.
    pub background: AnsiColor,
}

impl Default for AnsiStyle {
    fn default() -> Self {
        Self {
            bold: false,
            faint: false,
            underline: false,
            inverse: false,
            foreground: AnsiColor::Default,
            background: AnsiColor::Default,
        }
    }
}

/// One interpreted terminal cell.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct AnsiCell {
    /// Visible single-cell character.
    pub character: char,
    /// Rendition active when the character was written.
    pub style: AnsiStyle,
}

impl Default for AnsiCell {
    fn default() -> Self {
        Self {
            character: ' ',
            style: AnsiStyle::default(),
        }
    }
}

/// Strict interpretation of one complete viewer frame.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct AnsiScreen {
    /// Number of terminal columns represented by the cell grid.
    cols: u16,
    /// Number of terminal rows represented by the cell grid.
    rows: u16,
    /// Row-major interpreted terminal cells.
    cells: Vec<AnsiCell>,
}

impl AnsiScreen {
    /// Parses only the absolute cursor moves, SGR codes, and glyphs emitted by this viewer.
    pub(super) fn parse(frame: &str, cols: u16, rows: u16) -> Result<Self, String> {
        let effective_cols = cols.max(1);
        let effective_rows = rows.max(1);
        let cell_count = usize::from(effective_cols)
            .checked_mul(usize::from(effective_rows))
            .ok_or_else(|| String::from("terminal dimensions exceed addressable memory"))?;
        let mut screen = Self {
            cols: effective_cols,
            rows: effective_rows,
            cells: vec![AnsiCell::default(); cell_count],
        };
        let mut row = 0u16;
        let mut column = 0u16;
        let mut style = AnsiStyle::default();
        let mut characters = frame.chars();

        while let Some(character) = characters.next() {
            if character == '\x1b' {
                if characters.next() != Some('[') {
                    return Err(String::from("viewer ANSI contains a non-CSI escape"));
                }
                let mut parameters = String::new();
                let command = loop {
                    let next = characters
                        .next()
                        .ok_or_else(|| String::from("viewer ANSI ends inside a CSI sequence"))?;
                    if next.is_ascii_alphabetic() {
                        break next;
                    }
                    if !next.is_ascii_digit() && next != ';' {
                        return Err(format!(
                            "viewer ANSI contains unsupported CSI byte {next:?}"
                        ));
                    }
                    parameters.push(next);
                };
                match command {
                    'H' => parse_cursor_position(
                        &parameters,
                        &mut row,
                        &mut column,
                        effective_rows,
                        effective_cols,
                    )?,
                    'm' => parse_sgr(&parameters, &mut style)?,
                    _ => return Err(format!("viewer ANSI contains unsupported CSI {command}")),
                }
                continue;
            }
            if !is_supported_character(character) {
                return Err(format!(
                    "viewer ANSI contains unsupported display character {character:?}"
                ));
            }
            if column >= effective_cols {
                return Err(format!(
                    "viewer writes outside {effective_cols}x{effective_rows} terminal"
                ));
            }
            let index = usize::from(row)
                .checked_mul(usize::from(effective_cols))
                .and_then(|offset| offset.checked_add(usize::from(column)))
                .ok_or_else(|| String::from("viewer cursor position exceeds addressable memory"))?;
            let cell = screen.cells.get_mut(index).ok_or_else(|| {
                format!("viewer writes outside {effective_cols}x{effective_rows} terminal")
            })?;
            *cell = AnsiCell { character, style };
            column = column.saturating_add(1);
        }
        Ok(screen)
    }

    /// Returns the number of represented terminal columns.
    pub(super) const fn cols(&self) -> u16 {
        self.cols
    }

    /// Returns the number of represented terminal rows.
    pub(super) const fn rows(&self) -> u16 {
        self.rows
    }

    /// Returns a cell by zero-based terminal coordinates.
    pub(super) fn cell(&self, row: u16, column: u16) -> Option<AnsiCell> {
        if row >= self.rows || column >= self.cols {
            return None;
        }
        usize::from(row)
            .checked_mul(usize::from(self.cols))
            .and_then(|offset| offset.checked_add(usize::from(column)))
            .and_then(|index| self.cells.get(index))
            .copied()
    }
}

/// Parses a one-based absolute cursor position.
fn parse_cursor_position(
    parameters: &str,
    row: &mut u16,
    column: &mut u16,
    rows: u16,
    cols: u16,
) -> Result<(), String> {
    let (parsed_row, parsed_column) = if parameters.is_empty() {
        (1, 1)
    } else {
        let (row_text, column_text) = parameters
            .split_once(';')
            .ok_or_else(|| format!("viewer cursor position is incomplete: {parameters:?}"))?;
        if column_text.contains(';') {
            return Err(format!(
                "viewer cursor position has extra parameters: {parameters:?}"
            ));
        }
        (
            row_text
                .parse::<u16>()
                .map_err(|_parse_error| format!("viewer cursor row is invalid: {row_text:?}"))?,
            column_text.parse::<u16>().map_err(|_parse_error| {
                format!("viewer cursor column is invalid: {column_text:?}")
            })?,
        )
    };
    if parsed_row == 0 || parsed_column == 0 || parsed_row > rows || parsed_column > cols {
        return Err(format!(
            "viewer cursor position {parsed_row};{parsed_column} is outside {cols}x{rows} terminal"
        ));
    }
    *row = parsed_row.saturating_sub(1);
    *column = parsed_column.saturating_sub(1);
    Ok(())
}

/// Applies the restricted SGR vocabulary emitted by the frame builders.
fn parse_sgr(parameters: &str, style: &mut AnsiStyle) -> Result<(), String> {
    let effective_parameters = if parameters.is_empty() {
        "0"
    } else {
        parameters
    };
    for parameter in effective_parameters.split(';') {
        let code = parameter
            .parse::<u8>()
            .map_err(|_parse_error| format!("viewer SGR parameter is invalid: {parameter:?}"))?;
        match code {
            0 => *style = AnsiStyle::default(),
            1 => style.bold = true,
            2 => style.faint = true,
            4 => style.underline = true,
            7 => style.inverse = true,
            22 => {
                style.bold = false;
                style.faint = false;
            }
            24 => style.underline = false,
            32 => style.foreground = AnsiColor::Palette(2),
            33 => style.foreground = AnsiColor::Palette(3),
            36 => style.foreground = AnsiColor::Palette(6),
            39 => style.foreground = AnsiColor::Default,
            44 => style.background = AnsiColor::Palette(4),
            90 => style.foreground = AnsiColor::Palette(8),
            97 => style.foreground = AnsiColor::Palette(15),
            _ => return Err(format!("viewer ANSI contains unsupported SGR code {code}")),
        }
    }
    Ok(())
}

/// Returns whether the frame builders can deliberately emit this one-cell character.
fn is_supported_character(character: char) -> bool {
    character == ' '
        || character.is_ascii_graphic()
        || matches!(
            character,
            '·' | '•'
                | '●'
                | '━'
                | '┃'
                | '┓'
                | '┗'
                | '┛'
                | '┏'
                | '└'
                | '▒'
                | '┬'
                | '─'
                | '┤'
                | '│'
        )
}
