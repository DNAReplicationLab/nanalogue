//! Command-line parsing for the BAM viewer.

use super::{MAX_REGION_LENGTH, READ_LABEL_WIDTH};
use nanalogue_core::ModChar;
use std::{
    ffi::OsString,
    num::{IntErrorKind, NonZeroU32},
    path::PathBuf,
    str::FromStr as _,
};

/// Usage, display conventions, and controls shown for help and argument errors.
pub(super) const USAGE: &str = concat!(
    "Usage: nanalogue_bam_viewer <BAM> <CONTIG:START> [MOD_TYPE [WINDOW_SIZE individual]]\n",
    "START is zero-based; displayed coordinates are one-based.\n",
    "The end coordinate is selected from the terminal width, up to 200 bp.\n",
    "MOD_TYPE is a letter or numeric ChEBI code; calls with probability >= 0.5 are ",
    "bold and underlined.\n",
    "\n",
    "Display:\n",
    "  Green reads are forward; yellow reads are reverse.\n",
    "  Spaces are outside an alignment; dots are deletions or reference skips.\n",
    "  Lowercase bases are insertions.\n",
    "  An asterisk means the BAM alignment has no stored read sequence.\n",
    "  Individual view plots grey raw ML calls over the whole read and a bold ",
    "default-foreground step line of the per-window fraction of calls with ",
    "probability >= 0.5.\n",
    "  WINDOW_SIZE is a positive number of modified bases; windows are non-overlapping.\n",
    "\n",
    "Controls:\n",
    "  Left/Right or h/l move one genomic window.\n",
    "  Horizontal movement truncates read IDs and hides insertions.\n",
    "  Up/Down or k/j move one read; Page Up/Page Down move one read page.\n",
    "  Home/End jump to the first/last read; g prompts for CONTIG:START.\n",
    "  A successful goto also truncates read IDs and hides insertions.\n",
    "  In goto: type CONTIG:START; Backspace edits; Enter submits; Escape cancels.\n",
    "  In table view, r toggles full read IDs, i toggles insertions, and s saves ",
    "the rendered screen as text, plus a modification mask when MOD_TYPE is set.\n",
    "  In individual view, j/k selects one read; r/i/s have no effect.\n",
    "  Outside goto, q or Escape quits; Ctrl-C or Ctrl-D always quits.",
);

/// Viewer presentation selected by command-line arguments.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) enum ViewMode {
    /// Region sequence table.
    #[default]
    Table,
    /// Whole-read modification probability plot with non-overlapping windows.
    Individual {
        /// Number of modification calls per window.
        win: NonZeroU32,
    },
}

/// Returns whether a value contains only decimal digits after an optional leading plus.
fn is_unsigned_decimal(value: &str) -> bool {
    let digits = value.strip_prefix('+').unwrap_or(value);
    !digits.is_empty() && digits.bytes().all(|byte| byte.is_ascii_digit())
}

/// Initial reference and zero-based coordinate supplied on the command line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct InitialPosition {
    /// Reference name.
    pub contig: String,
    /// Zero-based first displayed coordinate.
    pub start: u32,
}

impl InitialPosition {
    /// Parses `CONTIG:START`, splitting at the final colon for colon-containing contig names.
    pub(super) fn parse(value: &str) -> Result<Self, String> {
        let (contig, start_text) = value
            .rsplit_once(':')
            .ok_or_else(|| String::from("position must have the form CONTIG:START"))?;
        if contig.is_empty() || start_text.is_empty() {
            return Err(String::from("position must have the form CONTIG:START"));
        }
        let start = start_text.parse::<u32>().map_err(|error| {
            if error.kind() == &IntErrorKind::PosOverflow && is_unsigned_decimal(start_text) {
                format!("START exceeds the maximum supported value ({})", u32::MAX)
            } else {
                String::from("START must be a non-negative integer")
            }
        })?;
        Ok(Self {
            contig: String::from(contig),
            start,
        })
    }
}

/// Returns the number of sequence columns available at a terminal width.
pub(super) fn window_len_for_columns(cols: u16) -> u32 {
    u32::from(cols.saturating_sub(READ_LABEL_WIDTH).max(1)).min(MAX_REGION_LENGTH)
}

/// Command-line options for the BAM viewer.
#[derive(Debug)]
pub(super) struct Args {
    /// Indexed BAM file to view.
    pub bam: PathBuf,

    /// Initial reference position, such as `chr1:1000`.
    pub position: InitialPosition,

    /// Optional modification type to display in bold.
    pub mod_type: Option<ModChar>,

    /// Requested presentation mode.
    pub mode: ViewMode,
}

impl Args {
    /// Parses the positional arguments, returning `None` for a help request.
    pub(super) fn parse_from<I>(arguments: I) -> Result<Option<Self>, String>
    where
        I: IntoIterator<Item = OsString>,
    {
        let mut argument_iter = arguments.into_iter();
        let Some(bam) = argument_iter.next() else {
            return Err(String::from("missing BAM file"));
        };
        if bam == "-h" || bam == "--help" {
            return Ok(None);
        }
        let position_argument = argument_iter
            .next()
            .ok_or_else(|| String::from("missing CONTIG:START position"))?;
        let position_string = position_argument
            .into_string()
            .map_err(|_position| String::from("position must be valid UTF-8"))?;
        let position = InitialPosition::parse(&position_string)?;
        let mod_type = argument_iter
            .next()
            .map(|mod_type_argument| {
                let mod_type_string = mod_type_argument
                    .into_string()
                    .map_err(|_mod_type| String::from("MOD_TYPE must be valid UTF-8"))?;
                ModChar::from_str(&mod_type_string).map_err(|error| error.to_string())
            })
            .transpose()?;
        let window_option = argument_iter.next();
        let mode_argument = argument_iter.next();
        let mode = match (mod_type, window_option, mode_argument) {
            (None | Some(_), None, None) => ViewMode::Table,
            (Some(_), Some(window_argument), Some(keyword)) => {
                let window_text = window_argument
                    .into_string()
                    .map_err(|_window| String::from("WINDOW_SIZE must be valid UTF-8"))?;
                let window = window_text.parse::<u32>().map_err(|error| {
                    if error.kind() == &IntErrorKind::PosOverflow
                        && is_unsigned_decimal(&window_text)
                    {
                        format!(
                            "WINDOW_SIZE exceeds the maximum supported value ({})",
                            u32::MAX
                        )
                    } else {
                        String::from("WINDOW_SIZE must be a positive integer")
                    }
                })?;
                let win = NonZeroU32::new(window)
                    .ok_or_else(|| String::from("WINDOW_SIZE must be a positive integer"))?;
                if keyword != "individual" {
                    return Err(String::from(
                        "expected literal 'individual' after WINDOW_SIZE",
                    ));
                }
                ViewMode::Individual { win }
            }
            (None, _, Some(_)) | (None, Some(_), None) => {
                return Err(String::from("individual view requires MOD_TYPE"));
            }
            (Some(_), None, Some(_)) | (Some(_), Some(_), None) => {
                return Err(String::from(
                    "WINDOW_SIZE and literal 'individual' must be supplied together",
                ));
            }
        };
        if argument_iter.next().is_some() {
            return Err(String::from("unexpected trailing arguments"));
        }
        Ok(Some(Self {
            bam: PathBuf::from(bam),
            position,
            mod_type,
            mode,
        }))
    }
}
