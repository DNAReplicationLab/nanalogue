//! ANSI frame assembly for table and shared viewer chrome.

use super::state::Viewer;
use nanalogue_core::region_sequences::RegionSequence;
use std::fmt::Write as _;

/// Converts arbitrary BAM text bytes into fixed-width printable ASCII.
fn printable_label(bytes: &[u8], width: usize) -> String {
    let mut label = bytes
        .iter()
        .take(width)
        .map(|byte| {
            if byte.is_ascii_graphic() || *byte == b' ' {
                char::from(*byte)
            } else {
                '?'
            }
        })
        .collect::<String>();
    while label.len() < width {
        label.push(' ');
    }
    label
}

/// Pads and styles a nanalogue region sequence to the visible genomic screen width.
pub(super) fn sequence_columns(
    record: &RegionSequence,
    width: u16,
    show_insertions: bool,
) -> String {
    let (sequence, modifications) = if show_insertions {
        (
            record.sequence_with_insertions(),
            record.modifications_with_insertions(),
        )
    } else {
        (record.sequence(), record.modifications())
    };
    assert_eq!(
        sequence.len(),
        modifications.len(),
        "region sequence and modification calls must have equal lengths"
    );
    let mut output = String::new();
    let mut bold = false;
    for offset in 0..usize::from(width) {
        let modified = modifications.get(offset).copied().unwrap_or(false);
        if modified != bold {
            output.push_str(if modified { "\x1b[1;4m" } else { "\x1b[22;24m" });
            bold = modified;
        }
        output.push(char::from(
            sequence.as_bytes().get(offset).copied().unwrap_or(b' '),
        ));
    }
    if bold {
        output.push_str("\x1b[22;24m");
    }
    output
}

/// Truncates and pads a line to the terminal width.
pub(super) fn fixed_line(text: &str, width: u16) -> String {
    printable_label(text.as_bytes(), usize::from(width))
}

/// Truncates or pads a label and appends its column separator.
pub(super) fn label_column(label: &str, width: u16) -> String {
    let mut column = printable_label(label.as_bytes(), usize::from(width.saturating_sub(1)));
    column.push(' ');
    column
}

/// Builds the compact footer while preserving its existing action-oriented toggle labels.
fn default_footer(viewer: &Viewer) -> String {
    format!(
        "h/l {} bp  j/k row  pgup/dn  home/end  g goto  r {} IDs  i {} ins  q quit",
        viewer.window_len,
        if viewer.full_read_ids {
            "short"
        } else {
            "full"
        },
        if viewer.show_insertions {
            "hide"
        } else {
            "show"
        }
    )
}

/// Selects the footer state rendered by the viewer.
#[derive(Clone, Copy)]
pub(super) enum FrameFooter<'a> {
    /// Show the normal navigation controls.
    Controls,
    /// Show the active genomic-position prompt or its current error.
    PositionPrompt {
        /// The position text entered so far.
        input: &'a str,
        /// The current validation error, if input submission failed.
        error: Option<&'a str>,
    },
}

/// Builds the ANSI frame that Ghostty parses into a terminal screen.
pub(super) fn build_frame(
    viewer: &Viewer,
    records: &[RegionSequence],
    cols: u16,
    rows: u16,
    footer_state: FrameFooter<'_>,
) -> String {
    let effective_cols = cols.max(1);
    let genome_cols = effective_cols
        .saturating_sub(viewer.read_label_width)
        .max(1)
        .min(u16::try_from(viewer.current_window_len()).expect("maximum region length fits u16"));
    let visible_reads = usize::from(rows.saturating_sub(4));
    let end = viewer
        .viewport
        .start
        .saturating_add(viewer.current_window_len())
        .min(viewer.target_len());
    let mut frame = String::from("\x1b[H");

    let title = fixed_line(" nanalogue BAM viewer", effective_cols);
    frame.push_str("\x1b[1;97;44m");
    frame.push_str(&title);
    frame.push_str("\x1b[0m");

    if rows > 1 {
        let path = viewer
            .path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("?");
        let mod_status = viewer.modification_status();
        let status = format!(
            " {path}  {}:{}-{}  reads {}  row {}{mod_status}",
            viewer.target_name(),
            viewer.viewport.start.saturating_add(1),
            end,
            records.len(),
            viewer.viewport.read_offset.saturating_add(1)
        );
        frame.push_str("\x1b[2;1H\x1b[36m");
        frame.push_str(&fixed_line(&status, effective_cols));
        frame.push_str("\x1b[0m");
    }

    if rows > 2 {
        let mut ruler = label_column("read id", viewer.read_label_width);
        for column in 0..usize::from(genome_cols) {
            let coordinate = viewer
                .viewport
                .start
                .saturating_add(u32::try_from(column).unwrap_or(u32::MAX))
                .saturating_add(1);
            ruler.push(if coordinate.is_multiple_of(10) {
                '|'
            } else {
                '.'
            });
        }
        frame.push_str("\x1b[3;1H\x1b[2m");
        frame.push_str(&fixed_line(&ruler, effective_cols));
        frame.push_str("\x1b[0m");
    }

    for (screen_row, record) in records
        .iter()
        .skip(viewer.viewport.read_offset)
        .take(visible_reads)
        .enumerate()
    {
        let terminal_row = screen_row.saturating_add(4);
        let visible_label_width = viewer.read_label_width.min(effective_cols);
        let visible_sequence_width =
            genome_cols.min(effective_cols.saturating_sub(visible_label_width));
        write!(&mut frame, "\x1b[{terminal_row};1H").expect("writing to String cannot fail");
        frame.push_str(if record.is_reverse() {
            "\x1b[33m"
        } else {
            "\x1b[32m"
        });
        frame.push_str(&fixed_line(
            &label_column(record.read_id(), viewer.read_label_width),
            visible_label_width,
        ));
        frame.push_str(&sequence_columns(
            record,
            visible_sequence_width,
            viewer.show_insertions,
        ));
        frame.push_str("\x1b[0m");
    }

    if rows > 3 {
        write!(&mut frame, "\x1b[{rows};1H\x1b[7m").expect("writing to String cannot fail");
        let footer = match footer_state {
            FrameFooter::Controls => default_footer(viewer),
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

/// Formats a goto prompt whose editable end remains visible at the terminal edge.
pub(super) fn position_prompt_footer(input: &str, cols: u16) -> String {
    const LONG_PREFIX: &str = "Go to CONTIG:START: ";
    let width = usize::from(cols.max(1));
    let prefix = if width > LONG_PREFIX.len() {
        LONG_PREFIX
    } else if width > 3 {
        "g: "
    } else {
        ""
    };
    let available = width.saturating_sub(prefix.len());
    let mut suffix_start = input.len().saturating_sub(available);
    while !input.is_char_boundary(suffix_start) {
        suffix_start = suffix_start.saturating_add(1);
    }
    let suffix = input
        .get(suffix_start..)
        .expect("suffix starts at a checked UTF-8 boundary");
    format!("{prefix}{suffix}")
}

/// Formats a goto error while reserving space for explicit corrective actions.
pub(super) fn position_error_footer(error: &str, cols: u16) -> String {
    const PREFIX: &str = "Error: ";
    const SUFFIX: &str = ". Backspace to correct; Esc to cancel.";
    let width = usize::from(cols.max(1));
    let fixed_width = PREFIX.len().saturating_add(SUFFIX.len());
    if width <= fixed_width {
        return if width >= 28 {
            String::from("Backspace edits; Esc cancels.")
        } else if width >= 19 {
            String::from("Bksp edit; Esc back")
        } else if width >= 8 {
            String::from("Bksp/Esc")
        } else {
            String::from("Esc")
        };
    }
    let detail_width = width.saturating_sub(fixed_width);
    let mut detail = error
        .bytes()
        .take(detail_width)
        .map(|byte| {
            if byte.is_ascii_graphic() || byte == b' ' {
                char::from(byte)
            } else {
                '?'
            }
        })
        .collect::<String>();
    if error.len() > detail_width && detail_width >= 3 {
        detail.truncate(detail_width.saturating_sub(3));
        detail.push_str("...");
    }
    format!("{PREFIX}{detail}{SUFFIX}")
}
