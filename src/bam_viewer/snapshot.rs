//! Plain-text projection and paired-file persistence for table snapshots.

use super::{render::TableSequenceGeometry, state::Viewer};
use libghostty_vt::{
    render::{CellIterator, RowIterator, Snapshot},
    style::Underline,
};
use std::{
    fs::{self, File, OpenOptions},
    io,
    path::{Path, PathBuf},
};

/// Text and modification-mask projections of one rendered terminal screen.
#[derive(Debug, PartialEq, Eq)]
pub(super) struct TextSnapshot {
    /// Visible character grid with one newline after every terminal row.
    pub text: String,
    /// Matching grid with binary cells in displayed sequence rectangles.
    pub modifications: String,
}

/// Projects Ghostty's interpreted table screen into text and modification grids.
pub(super) fn project_text_snapshot<'alloc>(
    snapshot: &Snapshot<'alloc, '_>,
    row_iterator: &mut RowIterator<'alloc>,
    cell_iterator: &mut CellIterator<'alloc>,
    geometry: TableSequenceGeometry,
) -> Result<TextSnapshot, libghostty_vt::Error> {
    let mut text = String::new();
    let mut modifications = String::new();
    let mut rows = row_iterator.update(snapshot)?;
    let mut row_index = 0u16;
    while let Some(row) = rows.next() {
        let mut cells = cell_iterator.update(row)?;
        let mut column_index = 0u16;
        while let Some(cell) = cells.next() {
            let graphemes = cell.graphemes()?;
            if graphemes.is_empty() {
                text.push(' ');
            } else {
                text.extend(&graphemes);
            }
            if geometry.contains(row_index, column_index) {
                modifications.push(if cell.style()?.underline == Underline::None {
                    '0'
                } else {
                    '1'
                });
            } else if graphemes.is_empty() {
                modifications.push(' ');
            } else {
                modifications.extend(&graphemes);
            }
            column_index = column_index.saturating_add(1);
        }
        text.push('\n');
        modifications.push('\n');
        row_index = row_index.saturating_add(1);
    }
    Ok(TextSnapshot {
        text,
        modifications,
    })
}

/// Converts a filename component into portable ASCII without path separators.
fn safe_component(value: &str) -> String {
    const MAX_COMPONENT_LENGTH: usize = 80;
    let mut safe = String::new();
    let mut replacing = false;
    for character in value.chars() {
        if safe.len() == MAX_COMPONENT_LENGTH {
            break;
        }
        if character.is_ascii_alphanumeric() || matches!(character, '-' | '_' | '.') {
            safe.push(character);
            replacing = false;
        } else {
            if !replacing {
                safe.push('_');
            }
            replacing = true;
        }
    }
    if safe.is_empty() {
        String::from("unknown")
    } else {
        safe
    }
}

/// Generates the shared snapshot prefix from the displayed table position.
pub(super) fn snapshot_prefix(viewer: &Viewer) -> String {
    let bam_name = viewer.path.file_name().map_or_else(
        || String::from("unknown.bam"),
        |name| name.to_string_lossy().into_owned(),
    );
    let displayed_start = viewer.viewport.start.saturating_add(1);
    let displayed_end = viewer
        .viewport
        .start
        .saturating_add(viewer.current_window_len())
        .min(viewer.target_len());
    format!(
        "nanalogue-{}-{}-{displayed_start}-{displayed_end}-row-{}",
        safe_component(&bam_name),
        safe_component(viewer.target_name()),
        viewer.viewport.read_offset.saturating_add(1)
    )
}

/// Paths successfully written for one snapshot pair.
#[derive(Debug, PartialEq, Eq)]
pub(super) struct SnapshotPaths {
    /// Plain screen text path.
    pub text: PathBuf,
    /// Modification mask path.
    pub modifications: PathBuf,
}

/// Opens a destination without replacing any existing filesystem entry.
fn create_new(path: &Path) -> io::Result<File> {
    OpenOptions::new().write(true).create_new(true).open(path)
}

/// Formats an I/O failure with the reason first so narrow footers remain useful.
#[expect(
    clippy::wildcard_enum_match_arm,
    reason = "only two error kinds need concise user-facing wording"
)]
fn io_failure(error: &io::Error) -> String {
    match error.kind() {
        io::ErrorKind::AlreadyExists => String::from("files already exist; nothing replaced"),
        io::ErrorKind::PermissionDenied => String::from("permission denied"),
        _ => error.to_string(),
    }
}

/// Removes files created by a failed paired write and includes cleanup errors in the message.
fn failed_write_error(error: &io::Error, created_paths: &[&Path]) -> String {
    let cleanup_errors = created_paths
        .iter()
        .filter_map(|path| fs::remove_file(path).err().map(|cleanup| (path, cleanup)))
        .map(|(path, cleanup)| format!("{}: {cleanup}", path.display()))
        .collect::<Vec<_>>();
    if cleanup_errors.is_empty() {
        io_failure(error)
    } else {
        format!(
            "{}; cleanup failed for {}",
            io_failure(error),
            cleanup_errors.join(", ")
        )
    }
}

/// Implements paired persistence with an injectable exclusive-create operation.
fn write_snapshot_pair_with<W, F>(
    directory: &Path,
    prefix: &str,
    snapshot: &TextSnapshot,
    mut open: F,
) -> Result<SnapshotPaths, String>
where
    W: io::Write,
    F: FnMut(&Path) -> io::Result<W>,
{
    let text_path = directory.join(format!("{prefix}.txt"));
    let modifications_path = directory.join(format!("{prefix}.mods.txt"));
    let mut text_file = open(&text_path).map_err(|error| {
        format!(
            "{} creating text snapshot ({})",
            io_failure(&error),
            text_path
                .file_name()
                .unwrap_or(text_path.as_os_str())
                .to_string_lossy()
        )
    })?;
    let mut modifications_file = match open(&modifications_path) {
        Ok(file) => file,
        Err(error) => {
            drop(text_file);
            let message = failed_write_error(&error, &[&text_path]);
            return Err(format!(
                "{message} creating modification mask ({})",
                modifications_path
                    .file_name()
                    .unwrap_or(modifications_path.as_os_str())
                    .to_string_lossy()
            ));
        }
    };

    let write_result = (|| -> io::Result<()> {
        text_file.write_all(snapshot.text.as_bytes())?;
        modifications_file.write_all(snapshot.modifications.as_bytes())?;
        text_file.flush()?;
        modifications_file.flush()
    })();
    if let Err(error) = write_result {
        drop(text_file);
        drop(modifications_file);
        return Err(failed_write_error(
            &error,
            &[&text_path, &modifications_path],
        ));
    }

    Ok(SnapshotPaths {
        text: text_path,
        modifications: modifications_path,
    })
}

/// Writes an all-or-neither snapshot pair without overwriting existing entries.
pub(super) fn write_snapshot_pair(
    directory: &Path,
    prefix: &str,
    snapshot: &TextSnapshot,
) -> Result<SnapshotPaths, String> {
    write_snapshot_pair_with(directory, prefix, snapshot, create_new)
}

/// Exercises paired rollback with controlled writers in unit tests.
#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
pub(super) fn write_snapshot_pair_with_test_opener<W, F>(
    directory: &Path,
    prefix: &str,
    snapshot: &TextSnapshot,
    open: F,
) -> Result<SnapshotPaths, String>
where
    W: io::Write,
    F: FnMut(&Path) -> io::Result<W>,
{
    write_snapshot_pair_with(directory, prefix, snapshot, open)
}
