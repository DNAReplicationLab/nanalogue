//! Whole-read plot rasterization and individual-view frame assembly.

use super::{
    cli::ViewMode,
    render::{FrameFooter, fixed_line, position_error_footer, position_prompt_footer},
    state::Viewer,
};
use nanalogue_core::region_sequences::ReadModProfile;
use std::fmt::Write as _;

/// One raster cell in an individual-read plot.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum PlotCell {
    /// Empty background.
    #[default]
    Empty,
    /// One or more raw calls mapped to this cell.
    Dots(u16),
    /// Windowed profile line, which takes precedence over raw calls.
    Line(char),
}

/// Maps a raw ML value onto a plot row using integer round-to-nearest arithmetic.
#[expect(
    clippy::arithmetic_side_effects,
    clippy::integer_division,
    clippy::integer_division_remainder_used,
    reason = "bounded ML values are intentionally mapped with round-to-nearest integer arithmetic"
)]
fn probability_row(probability: u8, plot_rows: usize) -> usize {
    if plot_rows <= 1 {
        return 0;
    }
    let height = u32::try_from(plot_rows.saturating_sub(1)).unwrap_or(u32::MAX);
    usize::try_from((u32::from(255u8.saturating_sub(probability)) * height + 127) / 255)
        .unwrap_or(usize::MAX)
        .min(plot_rows.saturating_sub(1))
}

/// Maps a thresholded window value onto a plot row.
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    reason = "bounded probabilities and terminal dimensions are intentionally rasterised to cells"
)]
fn window_value_row(value: nanalogue_core::F32Bw0and1, plot_rows: usize) -> usize {
    if plot_rows <= 1 {
        return 0;
    }
    let height = plot_rows.saturating_sub(1) as f32;
    ((1.0 - value.val()) * height).round() as usize
}

/// Maps a reference position onto a whole-read plot column.
fn reference_column(profile: &ReadModProfile, position: u32, plot_cols: usize) -> usize {
    reference_column_for_bounds(
        profile.align_start(),
        profile.align_end(),
        position,
        plot_cols,
    )
}

/// Maps a reference position onto a whole-read plot column for explicit alignment bounds.
#[expect(
    clippy::arithmetic_side_effects,
    clippy::integer_division,
    clippy::integer_division_remainder_used,
    reason = "the specified whole-read mapping requires truncating u64 integer arithmetic"
)]
fn reference_column_for_bounds(
    align_start: u32,
    align_end: u32,
    position: u32,
    plot_cols: usize,
) -> usize {
    if plot_cols <= 1 || align_end <= align_start.saturating_add(1) {
        return 0;
    }
    if position >= align_end.saturating_sub(1) {
        return plot_cols.saturating_sub(1);
    }
    let offset = u64::from(position.saturating_sub(align_start));
    let width = u64::try_from(plot_cols).unwrap_or(u64::MAX);
    let span = u64::from(align_end.saturating_sub(align_start));
    usize::try_from(offset.saturating_mul(width) / span)
        .unwrap_or(usize::MAX)
        .min(plot_cols.saturating_sub(1))
}

/// Adds a raw call without replacing a profile line.
fn add_dot(cell: &mut PlotCell) {
    match cell {
        PlotCell::Empty => *cell = PlotCell::Dots(1),
        PlotCell::Dots(count) => *count = count.saturating_add(1),
        PlotCell::Line(_) => {}
    }
}

/// Returns the density glyph for a raw-call count.
fn dot_glyph(count: u16) -> char {
    match count {
        0 | 1 => '·',
        2..=4 => '•',
        _ => '●',
    }
}

/// Rasterises a whole-read profile into plot cells.
fn rasterise_profile(
    profile: &ReadModProfile,
    plot_cols: usize,
    plot_rows: usize,
) -> Vec<Vec<PlotCell>> {
    rasterise_profile_data(
        profile.align_start(),
        profile.align_end(),
        profile.calls(),
        profile.windows(),
        profile.window_series_starts(),
        plot_cols,
        plot_rows,
    )
}

/// Rasterises profile primitives, keeping coordinate mapping independently testable.
#[expect(
    clippy::indexing_slicing,
    clippy::needless_range_loop,
    reason = "all raster indices are clamped and the arguments are the profile's primitive fields"
)]
fn rasterise_profile_data(
    align_start: u32,
    align_end: u32,
    calls: &[(u32, u8)],
    windows: &[(u32, u32, nanalogue_core::F32Bw0and1)],
    window_series_starts: &[usize],
    plot_cols: usize,
    plot_rows: usize,
) -> Vec<Vec<PlotCell>> {
    let mut grid = vec![vec![PlotCell::Empty; plot_cols]; plot_rows];
    if plot_cols == 0 || plot_rows == 0 {
        return grid;
    }
    for &(ref_pos, probability) in calls {
        let row = probability_row(probability, plot_rows);
        let col = reference_column_for_bounds(align_start, align_end, ref_pos, plot_cols);
        add_dot(&mut grid[row][col]);
    }

    let mut previous: Option<(u32, u32, usize, usize)> = None;
    for (window_index, &(ref_start, ref_end, value)) in windows.iter().enumerate() {
        if window_index > 0 && window_series_starts.binary_search(&window_index).is_ok() {
            previous = None;
        }
        if ref_start >= ref_end {
            continue;
        }
        let row = window_value_row(value, plot_rows);
        let start_col = reference_column_for_bounds(align_start, align_end, ref_start, plot_cols);
        let end_col = reference_column_for_bounds(
            align_start,
            align_end,
            ref_end.saturating_sub(1),
            plot_cols,
        )
        .max(start_col);
        let boundary_cell_before = grid[row][start_col];
        for col in start_col..=end_col {
            grid[row][col] = PlotCell::Line('━');
        }

        if let Some((previous_start, previous_end, previous_row, previous_end_col)) = previous
            && ref_start >= previous_start
            && ref_start >= previous_end
        {
            let connector_col = start_col;
            for col in previous_end_col.saturating_add(1)..=connector_col {
                grid[previous_row][col] = PlotCell::Line('━');
            }
            if previous_row == row {
                grid[row][connector_col] = match boundary_cell_before {
                    PlotCell::Line(glyph) if connector_col == previous_end_col && glyph != '━' => {
                        PlotCell::Line(glyph)
                    }
                    PlotCell::Empty | PlotCell::Dots(_) | PlotCell::Line(_) => PlotCell::Line('━'),
                };
                previous = Some((ref_start, ref_end, row, end_col));
                continue;
            }
            let (top, bottom) = if previous_row <= row {
                (previous_row, row)
            } else {
                (row, previous_row)
            };
            for connector_row in top..=bottom {
                grid[connector_row][connector_col] = PlotCell::Line('┃');
            }
            if connector_col == previous_end_col {
                // Multiple windows compressed into one column retain a continuous vertical spine.
            } else if previous_row < row {
                grid[previous_row][connector_col] = PlotCell::Line('┓');
                grid[row][connector_col] = PlotCell::Line('┗');
            } else {
                grid[previous_row][connector_col] = PlotCell::Line('┛');
                grid[row][connector_col] = PlotCell::Line('┏');
            }
        }
        previous = Some((ref_start, ref_end, row, end_col));
    }
    grid
}

/// Truncates and pads Unicode text by terminal cells used by this viewer's single-width glyphs.
fn unicode_line(text: &str, width: u16) -> String {
    let cell_width = usize::from(width);
    let mut line = text.chars().take(cell_width).collect::<String>();
    line.extend(std::iter::repeat_n(
        ' ',
        cell_width.saturating_sub(line.chars().count()),
    ));
    line
}

/// Sanitises and truncates an external label without padding it.
fn compact_label(text: &str, width: usize) -> String {
    text.bytes()
        .take(width)
        .map(|byte| {
            if byte.is_ascii_graphic() || byte == b' ' {
                char::from(byte)
            } else {
                '?'
            }
        })
        .collect()
}

/// Sanitises a label and marks truncation instead of presenting it as complete.
fn abbreviated_label(text: &str, width: usize) -> String {
    let mut label = compact_label(text, width);
    if text.len() > width && width > 0 {
        let _last_character = label.pop();
        label.push('~');
    }
    label
}

/// Returns the first reference base actually plotted in a column.
fn first_base_for_column(profile: &ReadModProfile, column: usize, plot_cols: usize) -> Option<u32> {
    first_base_for_bounds(
        profile.align_start(),
        profile.align_end(),
        column,
        plot_cols,
    )
}

/// Returns the first reference base actually plotted in a column for explicit bounds.
#[expect(
    clippy::arithmetic_side_effects,
    clippy::integer_division,
    clippy::integer_division_remainder_used,
    reason = "ceiling division maps each ruler column to its first reference base"
)]
fn first_base_for_bounds(
    align_start: u32,
    align_end: u32,
    column: usize,
    plot_cols: usize,
) -> Option<u32> {
    if plot_cols == 0 || align_start >= align_end || column >= plot_cols {
        return None;
    }
    let span = u64::from(align_end.saturating_sub(align_start));
    let numerator = u64::try_from(column)
        .unwrap_or(u64::MAX)
        .saturating_mul(span);
    let denominator = u64::try_from(plot_cols).unwrap_or(u64::MAX).max(1);
    let offset = numerator.saturating_add(denominator.saturating_sub(1)) / denominator;
    let candidate = align_start.saturating_add(u32::try_from(offset).unwrap_or(u32::MAX));
    if candidate < align_end
        && reference_column_for_bounds(align_start, align_end, candidate, plot_cols) == column
    {
        Some(candidate)
    } else if column == plot_cols.saturating_sub(1) {
        let last_base = align_end.saturating_sub(1);
        (reference_column_for_bounds(align_start, align_end, last_base, plot_cols) == column)
            .then_some(last_base)
    } else {
        None
    }
}

/// Writes a complete coordinate label when it fits; partial coordinates are never emitted.
fn place_coordinate_label(labels: &mut [char], label_start: usize, coordinate: u32) -> bool {
    let label = coordinate.to_string();
    let label_end = label_start.saturating_add(label.len());
    let Some(destination) = labels.get_mut(label_start..label_end) else {
        return false;
    };
    for (cell, character) in destination.iter_mut().zip(label.chars()) {
        *cell = character;
    }
    true
}

/// Builds the x-axis ruler and its labels for one whole-read profile.
fn plot_ruler(viewer: &Viewer, profile: &ReadModProfile, plot_cols: usize) -> (String, String) {
    let viewport_end = viewer
        .viewport
        .start
        .saturating_add(viewer.current_window_len());
    let shade_start = reference_column(profile, viewer.viewport.start, plot_cols);
    let mut shade_end = reference_column(profile, viewport_end.saturating_sub(1), plot_cols);
    shade_end = shade_end.max(shade_start);
    let mut ruler = String::from("    └");
    for column in 0..plot_cols {
        ruler.push(if (shade_start..=shade_end).contains(&column) {
            '▒'
        } else if column.is_multiple_of(10) {
            '┬'
        } else {
            '─'
        });
    }

    let mut labels = vec![' '; plot_cols.saturating_add(5)];
    let mut previous_coordinate = None;
    for column in (0..plot_cols).step_by(10) {
        let Some(coordinate) = first_base_for_column(profile, column, plot_cols)
            .map(|position| position.saturating_add(1))
        else {
            continue;
        };
        if previous_coordinate == Some(coordinate) {
            continue;
        }
        previous_coordinate = Some(coordinate);
        let label_start = column.saturating_add(5);
        let _label_was_placed = place_coordinate_label(&mut labels, label_start, coordinate);
    }
    (ruler, labels.into_iter().collect())
}

/// Builds a status line without ever partially clipping numeric fields.
pub(super) fn individual_status(
    viewer: &Viewer,
    profiles: &[ReadModProfile],
    width: u16,
) -> String {
    let path_text = viewer
        .path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("?");
    let target_text = viewer.target_name();
    let end = viewer
        .viewport
        .start
        .saturating_add(viewer.current_window_len());
    let mod_label = viewer
        .mod_type
        .map_or_else(|| String::from("?"), |value| value.to_string());
    let win = match viewer.mode {
        ViewMode::Individual { win } => win.get(),
        ViewMode::Table => 0,
    };
    let selected = profiles.get(viewer.viewport.read_offset);
    let read_number = selected.map_or(0, |_profile| viewer.viewport.read_offset.saturating_add(1));
    let middle = format!(
        ":{}-{end} reads {} read {read_number}/{} ",
        viewer.viewport.start.saturating_add(1),
        profiles.len(),
        profiles.len()
    );
    let suffix = selected.map_or_else(
        || format!("mods {mod_label} win {win}"),
        |profile| {
            format!(
                " {} mods {mod_label} win {win}",
                if profile.is_reverse() { '-' } else { '+' }
            )
        },
    );
    let variable_count = if selected.is_some() { 3 } else { 2 };
    let fixed_width = middle.len().saturating_add(suffix.len()).saturating_add(2);
    let available = usize::from(width).saturating_sub(fixed_width);
    if available < variable_count {
        let fallback = format!("win {win}");
        return if fallback.len() <= usize::from(width) {
            fallback
        } else {
            String::new()
        };
    }

    let target_width = target_text
        .len()
        .min(available.saturating_sub(variable_count.saturating_sub(1)));
    let remaining = available.saturating_sub(target_width);
    let path_width = path_text
        .len()
        .min(8)
        .min(remaining.saturating_sub(usize::from(selected.is_some())));
    let id_width = remaining.saturating_sub(path_width);
    let path = abbreviated_label(path_text, path_width);
    let target = abbreviated_label(target_text, target_width);
    selected.map_or_else(
        || format!(" {path} {target}{middle}{suffix}"),
        |profile| {
            format!(
                " {path} {target}{middle}{}{suffix}",
                abbreviated_label(profile.read_id(), id_width)
            )
        },
    )
}

/// Builds an ANSI frame for a whole-read modification probability plot.
#[expect(
    clippy::integer_division,
    clippy::integer_division_remainder_used,
    reason = "frame assembly is cohesive and midpoint arithmetic intentionally uses terminal cells"
)]
pub(super) fn build_individual_frame(
    viewer: &Viewer,
    profiles: &[ReadModProfile],
    cols: u16,
    rows: u16,
    footer_state: FrameFooter<'_>,
) -> String {
    let effective_cols = cols.max(1);
    let mut frame = String::from("\x1b[H");
    frame.push_str("\x1b[1;97;44m");
    frame.push_str(&fixed_line(" nanalogue BAM viewer", effective_cols));
    frame.push_str("\x1b[0m");
    let selected = profiles.get(viewer.viewport.read_offset);

    if rows > 1 {
        let status = individual_status(viewer, profiles, effective_cols);
        frame.push_str("\x1b[2;1H\x1b[36m");
        frame.push_str(&fixed_line(&status, effective_cols));
        frame.push_str("\x1b[0m");
    }

    if rows < 8 {
        // Below eight rows, only title, status, and footer remain legible.
    } else if effective_cols < 15 {
        if rows > 2 {
            frame.push_str("\x1b[3;1H");
            frame.push_str(&fixed_line("terminal narrow", effective_cols));
        }
    } else if let Some(profile) = selected {
        let plot_rows = usize::from(rows.saturating_sub(5));
        let plot_cols = usize::from(effective_cols.saturating_sub(5));
        let grid = rasterise_profile(profile, plot_cols, plot_rows);
        for (row, cells) in grid.iter().enumerate() {
            let terminal_row = row.saturating_add(3);
            write!(frame, "\x1b[{terminal_row};1H").expect("writing to String cannot fail");
            let label = match row {
                0 => "1.0 ┤",
                value if value == plot_rows / 2 && plot_rows % 2 == 1 => "0.5 ┤",
                value if value.saturating_add(1) == plot_rows => "0.0 ┤",
                _ => "    │",
            };
            frame.push_str(label);
            let mut grey = false;
            for cell in cells {
                match cell {
                    PlotCell::Empty => frame.push(' '),
                    PlotCell::Dots(count) => {
                        if !grey {
                            frame.push_str("\x1b[90m");
                            grey = true;
                        }
                        frame.push(dot_glyph(*count));
                    }
                    PlotCell::Line(glyph) => {
                        if grey {
                            frame.push_str("\x1b[39m");
                            grey = false;
                        }
                        frame.push_str("\x1b[1m");
                        frame.push(*glyph);
                        frame.push_str("\x1b[22m");
                    }
                }
            }
            if grey {
                frame.push_str("\x1b[39m");
            }
        }
        if profile.calls().is_empty() && plot_rows > 0 {
            let message = format!(
                "no {} calls in read",
                viewer
                    .mod_type
                    .map_or_else(|| String::from("?"), |value| value.to_string())
            );
            frame.push_str("\x1b[3;7H");
            frame.push_str(&fixed_line(&message, effective_cols.saturating_sub(6)));
        }
        if rows > 4 {
            let (ruler, labels) = plot_ruler(viewer, profile, plot_cols);
            write!(frame, "\x1b[{};1H", rows.saturating_sub(2))
                .expect("writing to String cannot fail");
            frame.push_str(&unicode_line(&ruler, effective_cols));
            write!(frame, "\x1b[{};1H", rows.saturating_sub(1))
                .expect("writing to String cannot fail");
            frame.push_str(&fixed_line(&labels, effective_cols));
        }
    } else if rows > 2 {
        frame.push_str("\x1b[3;1H");
        frame.push_str(&fixed_line(" no reads span this window", effective_cols));
    } else {
        // The title and optional status are the complete degraded view at this height.
    }

    if rows > 3 {
        write!(frame, "\x1b[{rows};1H\x1b[7m").expect("writing to String cannot fail");
        let footer = match footer_state {
            FrameFooter::Controls => format!(
                "j/k read  pgup/dn  home/end  h/l {} bp  g goto  q quit",
                viewer.window_len
            ),
            FrameFooter::PositionPrompt { input, error } => error.map_or_else(
                || position_prompt_footer(input, cols),
                |message| position_error_footer(message, cols),
            ),
        };
        frame.push_str(&fixed_line(&footer, effective_cols));
        frame.push_str("\x1b[0m");
    }
    frame
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plot_probability_and_reference_mapping_cover_boundaries() {
        assert_eq!(probability_row(255, 5), 0);
        assert_eq!(probability_row(0, 5), 4);
        assert_eq!(probability_row(128, 5), 2);
        assert_eq!(
            window_value_row(nanalogue_core::F32Bw0and1::new(0.5).expect("bounded"), 5),
            2
        );

        assert_eq!(reference_column_for_bounds(100, 200, 100, 10), 0);
        assert_eq!(reference_column_for_bounds(100, 200, 199, 10), 9);
        assert_eq!(reference_column_for_bounds(100, 103, 100, 10), 0);
        assert_eq!(reference_column_for_bounds(100, 103, 101, 10), 3);
        assert_eq!(reference_column_for_bounds(100, 103, 102, 10), 9);

        let tick_coordinates = (0..30)
            .step_by(10)
            .map(|column| first_base_for_bounds(100, 103, column, 30))
            .collect::<Vec<_>>();
        assert_eq!(tick_coordinates, [Some(100), Some(101), None]);
        for column in 0..30 {
            if let Some(position) = first_base_for_bounds(100, 103, column, 30) {
                assert!((100..103).contains(&position));
                assert_eq!(reference_column_for_bounds(100, 103, position, 30), column);
            }
        }
        assert_eq!(first_base_for_bounds(100, 101, 0, 30), Some(100));
        assert_eq!(first_base_for_bounds(100, 101, 29, 30), None);

        let mut labels = vec![' '; 80];
        assert!(!place_coordinate_label(&mut labels, 75, 100_935));
        assert!(labels.iter().all(|character| *character == ' '));
        assert!(place_coordinate_label(&mut labels, 74, 100_935));
        assert_eq!(labels.get(74..), Some(&['1', '0', '0', '9', '3', '5'][..]));
    }

    #[test]
    fn plot_dot_density_and_line_precedence_are_unambiguous() {
        assert_eq!(dot_glyph(1), '·');
        assert_eq!(dot_glyph(2), '•');
        assert_eq!(dot_glyph(4), '•');
        assert_eq!(dot_glyph(5), '●');

        let mut cell = PlotCell::Empty;
        add_dot(&mut cell);
        add_dot(&mut cell);
        assert_eq!(cell, PlotCell::Dots(2));
        cell = PlotCell::Line('━');
        add_dot(&mut cell);
        assert_eq!(cell, PlotCell::Line('━'));
    }

    #[test]
    fn reverse_mapped_calls_keep_reference_orientation() {
        let calls = [(25, 0), (26, 255), (27, 0), (28, 255), (29, 0), (30, 255)];
        let grid = rasterise_profile_data(1, 31, &calls, &[], &[], 30, 5);

        for (row, column) in [(4, 24), (0, 25), (4, 26), (0, 27), (4, 28), (0, 29)] {
            assert_eq!(
                grid.get(row).and_then(|cells| cells.get(column)),
                Some(&PlotCell::Dots(1))
            );
        }
    }

    #[test]
    fn plot_windows_form_continuous_steps_but_independent_series_do_not_join() {
        let high = nanalogue_core::F32Bw0and1::new(1.0).expect("bounded");
        let low = nanalogue_core::F32Bw0and1::new(0.0).expect("bounded");
        let windows = [(0, 3, high), (10, 13, low)];
        let connected = rasterise_profile_data(0, 20, &[], &windows, &[0], 20, 3);
        let top = connected.first().expect("top row");
        assert!(
            top.get(3..10)
                .expect("gap columns")
                .iter()
                .all(|cell| *cell == PlotCell::Line('━'))
        );
        assert_eq!(top.get(10), Some(&PlotCell::Line('┓')));
        assert_eq!(
            connected.get(1).and_then(|row| row.get(10)),
            Some(&PlotCell::Line('┃'))
        );
        assert_eq!(
            connected.get(2).and_then(|row| row.get(10)),
            Some(&PlotCell::Line('┗'))
        );

        let separate = rasterise_profile_data(0, 20, &[], &windows, &[0, 1], 20, 3);
        assert!(
            separate
                .first()
                .and_then(|row| row.get(3..10))
                .expect("gap columns")
                .iter()
                .all(|cell| *cell == PlotCell::Empty)
        );

        let three_windows = [(0, 1, high), (10, 11, low), (20, 21, high)];
        let three_step = rasterise_profile_data(0, 30, &[], &three_windows, &[0], 30, 3);
        assert_eq!(
            three_step.first().and_then(|row| row.get(10)),
            Some(&PlotCell::Line('┓')),
            "the third plateau must not erase the middle plateau's incoming corner"
        );

        let collapsed = rasterise_profile_data(
            0,
            100,
            &[],
            &[
                (0, 1, high),
                (1, 2, low),
                (2, 3, nanalogue_core::F32Bw0and1::new(0.5).expect("bounded")),
            ],
            &[0],
            1,
            3,
        );
        assert!(
            collapsed
                .iter()
                .all(|row| row.first() == Some(&PlotCell::Line('┃')))
        );
    }

    #[test]
    fn plot_flat_continuations_preserve_existing_connections() {
        let high = nanalogue_core::F32Bw0and1::new(1.0).expect("bounded");
        let low = nanalogue_core::F32Bw0and1::new(0.0).expect("bounded");
        let flat =
            rasterise_profile_data(0, 100, &[], &[(0, 25, high), (25, 50, high)], &[0], 10, 3);
        assert!(
            flat.first()
                .and_then(|row| row.get(0..5))
                .expect("flat plateau")
                .iter()
                .all(|cell| *cell == PlotCell::Line('━'))
        );

        let collapsed_flat_end = rasterise_profile_data(
            0,
            100,
            &[],
            &[(0, 1, high), (1, 2, low), (2, 3, low)],
            &[0],
            1,
            3,
        );
        assert!(
            collapsed_flat_end
                .iter()
                .all(|row| row.first() == Some(&PlotCell::Line('┃')))
        );

        let descending_then_flat = rasterise_profile_data(
            0,
            100,
            &[],
            &[(0, 10, high), (10, 15, low), (15, 30, low)],
            &[0],
            10,
            3,
        );
        assert_eq!(
            descending_then_flat.get(2).and_then(|row| row.get(1)),
            Some(&PlotCell::Line('┗'))
        );
        let ascending_then_flat = rasterise_profile_data(
            0,
            100,
            &[],
            &[(0, 10, low), (10, 15, high), (15, 30, high)],
            &[0],
            10,
            3,
        );
        assert_eq!(
            ascending_then_flat.first().and_then(|row| row.get(1)),
            Some(&PlotCell::Line('┏'))
        );
    }
}
