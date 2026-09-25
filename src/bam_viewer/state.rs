//! Viewer state, navigation, and BAM data selection.

use super::{
    MAX_REGION_LENGTH, READ_LABEL_WIDTH,
    cli::{InitialPosition, ViewMode},
};
use crossterm::event::KeyCode;
use nanalogue_core::{
    ModChar,
    region_sequences::{ReadModProfile, RegionSequence, RegionSequenceReader},
};
use std::{error::Error, num::NonZeroU32, path::PathBuf};

/// Position and scroll state of the viewer.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct Viewport {
    /// Numeric BAM target identifier.
    pub tid: u32,
    /// Zero-based first visible reference coordinate.
    pub start: u32,
    /// First visible read row.
    pub read_offset: usize,
}

impl Viewport {
    /// Applies one navigation key while respecting reference and read bounds.
    #[expect(
        clippy::wildcard_enum_match_arm,
        reason = "all other terminal keys intentionally leave the viewport unchanged"
    )]
    pub(super) fn navigate(
        &mut self,
        key: KeyCode,
        contig_len: u32,
        window_len: u32,
        read_count: usize,
        visible_reads: usize,
    ) {
        let max_start = contig_len.saturating_sub(window_len);
        let navigation_page_size = visible_reads.max(1);
        let max_offset = read_count.saturating_sub(navigation_page_size);
        match key {
            KeyCode::Left | KeyCode::Char('h') => {
                self.start = self.start.saturating_sub(window_len);
                self.read_offset = 0;
            }
            KeyCode::Right | KeyCode::Char('l') => {
                self.start = self
                    .start
                    .saturating_add(window_len)
                    .min(max_start)
                    .max(self.start);
                self.read_offset = 0;
            }
            KeyCode::Up | KeyCode::Char('k') => {
                self.read_offset = self.read_offset.saturating_sub(1);
            }
            KeyCode::Down | KeyCode::Char('j') => {
                self.read_offset = self.read_offset.saturating_add(1).min(max_offset);
            }
            KeyCode::PageUp => {
                self.read_offset = self.read_offset.saturating_sub(navigation_page_size);
            }
            KeyCode::PageDown => {
                self.read_offset = self
                    .read_offset
                    .saturating_add(navigation_page_size)
                    .min(max_offset);
            }
            KeyCode::Home => {
                self.read_offset = 0;
            }
            KeyCode::End => {
                self.read_offset = max_offset;
            }
            _ => {}
        }
    }
}

/// Returns enough columns for the longest cached read ID and a separator.
pub(super) fn full_read_label_width(records: &[RegionSequence]) -> u16 {
    records
        .iter()
        .map(|record| {
            u16::try_from(record.read_id().len())
                .unwrap_or(u16::MAX)
                .saturating_add(1)
        })
        .max()
        .unwrap_or(READ_LABEL_WIDTH)
}

/// BAM data and mutable viewport state.
#[derive(Debug)]
pub(super) struct Viewer {
    /// Source BAM path.
    pub path: PathBuf,
    /// Nanalogue reader reused for region sequence tables.
    pub reader: RegionSequenceReader,
    /// Current position and read scroll.
    pub viewport: Viewport,
    /// Number of genomic bases in each fetched window.
    pub window_len: u32,
    /// Width of the read-ID column, including its trailing separator.
    pub read_label_width: u16,
    /// Whether complete read IDs are displayed.
    pub full_read_ids: bool,
    /// Whether lowercase insertion bases are displayed.
    pub show_insertions: bool,
    /// Modification type displayed in bold, if requested.
    pub mod_type: Option<ModChar>,
    /// Active presentation mode.
    pub mode: ViewMode,
}

impl Viewer {
    /// Opens an indexed BAM and chooses the initial viewport.
    #[cfg(test)]
    #[cfg_attr(coverage_nightly, coverage(off))]
    pub(super) fn open(
        path: PathBuf,
        position: &InitialPosition,
        mod_type: Option<ModChar>,
        window_len: u32,
    ) -> Result<Self, Box<dyn Error>> {
        Self::open_mode(path, position, mod_type, ViewMode::Table, window_len)
    }

    /// Opens an indexed BAM with an explicit presentation mode.
    pub(super) fn open_mode(
        path: PathBuf,
        position: &InitialPosition,
        mod_type: Option<ModChar>,
        mode: ViewMode,
        window_len: u32,
    ) -> Result<Self, Box<dyn Error>> {
        let reader = RegionSequenceReader::from_path(&path)?;

        let tid = reader
            .target_id(&position.contig)
            .ok_or_else(|| format!("reference '{}' is not in the BAM header", position.contig))?;
        let contig_len = reader
            .target_len(tid)
            .ok_or("BAM target identifier is invalid")?;
        if position.start >= contig_len {
            return Err(format!(
                "initial position {} is outside reference '{}' (length {contig_len})",
                position.start, position.contig
            )
            .into());
        }
        let viewport = Viewport {
            tid,
            start: position.start,
            read_offset: 0,
        };

        Ok(Self {
            path,
            reader,
            viewport,
            window_len: window_len.clamp(1, MAX_REGION_LENGTH),
            read_label_width: READ_LABEL_WIDTH,
            full_read_ids: false,
            show_insertions: false,
            mod_type,
            mode,
        })
    }

    /// Returns the current target name.
    pub(super) fn target_name(&self) -> &str {
        self.reader.target_name(self.viewport.tid).unwrap_or("?")
    }

    /// Returns the current target length.
    pub(super) fn target_len(&self) -> u32 {
        self.reader.target_len(self.viewport.tid).unwrap_or(0)
    }

    /// Returns the current window length, shortened only at the end of a reference.
    pub(super) fn current_window_len(&self) -> u32 {
        self.window_len
            .min(self.target_len().saturating_sub(self.viewport.start))
    }

    /// Returns the status suffix for the active modification type.
    pub(super) fn modification_status(&self) -> String {
        self.mod_type
            .map(|mod_type| format!("  mods {mod_type}>=0.5"))
            .unwrap_or_default()
    }

    /// Toggles between the default and longest cached read-ID widths.
    pub(super) fn toggle_read_id_width(&mut self, records: &[RegionSequence]) {
        self.full_read_ids = !self.full_read_ids;
        self.read_label_width = if self.full_read_ids {
            full_read_label_width(records)
        } else {
            READ_LABEL_WIDTH
        };
    }

    /// Restores the default read-ID width.
    pub(super) fn reset_read_id_width(&mut self) {
        self.full_read_ids = false;
        self.read_label_width = READ_LABEL_WIDTH;
    }

    /// Restores display options that do not carry across genomic windows.
    fn reset_horizontal_options(&mut self) {
        self.reset_read_id_width();
        self.show_insertions = false;
    }

    /// Moves to a validated genomic position and reports whether new records are needed.
    pub(super) fn go_to(&mut self, position: &InitialPosition) -> Result<bool, String> {
        let tid = self
            .reader
            .target_id(&position.contig)
            .ok_or_else(|| format!("unknown reference '{}'", position.contig))?;
        let contig_len = self
            .reader
            .target_len(tid)
            .ok_or_else(|| String::from("BAM target identifier is invalid"))?;
        if position.start >= contig_len {
            return Err(format!(
                "position {} is outside reference '{}' (length {contig_len})",
                position.start, position.contig
            ));
        }
        let changed = self.viewport.tid != tid || self.viewport.start != position.start;
        self.viewport = Viewport {
            tid,
            start: position.start,
            read_offset: 0,
        };
        self.reset_horizontal_options();
        Ok(changed)
    }

    /// Applies a navigation or display key and reports whether a new region must be fetched.
    pub(super) fn handle_key(
        &mut self,
        key: KeyCode,
        records: &[RegionSequence],
        visible_reads: usize,
    ) -> bool {
        if key == KeyCode::Char('r') && self.mode == ViewMode::Table {
            self.toggle_read_id_width(records);
            return false;
        }
        if key == KeyCode::Char('i') && self.mode == ViewMode::Table {
            self.show_insertions = !self.show_insertions;
            return false;
        }
        if matches!(
            key,
            KeyCode::Left | KeyCode::Right | KeyCode::Char('h' | 'l')
        ) {
            self.reset_horizontal_options();
        }
        let previous_start = self.viewport.start;
        self.viewport.navigate(
            key,
            self.target_len(),
            self.window_len,
            records.len(),
            visible_reads,
        );
        self.viewport.start != previous_start
    }

    /// Fetches records overlapping the visible genomic range.
    pub(super) fn visible_records(&mut self) -> Result<Vec<RegionSequence>, Box<dyn Error>> {
        let end = self
            .viewport
            .start
            .saturating_add(self.current_window_len())
            .min(self.target_len());
        Ok(self
            .reader
            .sequences(self.viewport.tid, self.viewport.start, end, self.mod_type)?)
    }

    /// Fetches modification profiles spanning the visible genomic range.
    pub(super) fn visible_profiles(
        &mut self,
        win: NonZeroU32,
    ) -> Result<Vec<ReadModProfile>, Box<dyn Error>> {
        let end = self
            .viewport
            .start
            .saturating_add(self.current_window_len())
            .min(self.target_len());
        let mod_type = self
            .mod_type
            .ok_or("individual view requires a modification type")?;
        Ok(self
            .reader
            .profiles(self.viewport.tid, self.viewport.start, end, mod_type, win)?)
    }

    /// Applies navigation in individual mode, where exactly one read is visible.
    pub(super) fn handle_individual_key(&mut self, key: KeyCode, read_count: usize) -> bool {
        let previous_start = self.viewport.start;
        let previous_read_offset = self.viewport.read_offset;
        self.viewport
            .navigate(key, self.target_len(), self.window_len, read_count, 1);
        let moved_horizontally = self.viewport.start != previous_start;
        if !moved_horizontally
            && matches!(
                key,
                KeyCode::Left | KeyCode::Right | KeyCode::Char('h' | 'l')
            )
        {
            self.viewport.read_offset = previous_read_offset;
        }
        moved_horizontally
    }
}

/// Cached records for the active view.
#[derive(Debug)]
pub(super) enum ViewerRecords {
    /// Projected sequences for the table view.
    Table(Vec<RegionSequence>),
    /// Whole-read profiles for the individual view.
    Individual(Vec<ReadModProfile>),
}

/// Selected individual alignment retained across a horizontal refetch.
#[derive(Clone, Debug)]
pub(super) struct IndividualSelection(ReadModProfile);

impl ViewerRecords {
    /// Returns the number of cached reads.
    pub(super) fn len(&self) -> usize {
        match self {
            Self::Table(records) => records.len(),
            Self::Individual(records) => records.len(),
        }
    }

    /// Returns a read ID at an index, if present.
    pub(super) fn read_id(&self, index: usize) -> Option<&str> {
        match self {
            Self::Table(records) => records.get(index).map(RegionSequence::read_id),
            Self::Individual(records) => records.get(index).map(ReadModProfile::read_id),
        }
    }

    /// Captures the selected alignment when the viewer is in individual mode.
    pub(super) fn individual_selection(&self, index: usize) -> Option<IndividualSelection> {
        let Self::Individual(profiles) = self else {
            return None;
        };
        profiles.get(index).cloned().map(IndividualSelection)
    }
}

/// Fetches records appropriate for the viewer's active mode.
pub(super) fn fetch_viewer_records(viewer: &mut Viewer) -> Result<ViewerRecords, Box<dyn Error>> {
    match viewer.mode {
        ViewMode::Table => Ok(ViewerRecords::Table(viewer.visible_records()?)),
        ViewMode::Individual { win } => {
            Ok(ViewerRecords::Individual(viewer.visible_profiles(win)?))
        }
    }
}

/// Re-selects a cached table read by ID after a refetch, or clamps the previous index.
pub(super) fn reselect_read(
    records: &ViewerRecords,
    read_id: Option<&str>,
    previous_index: usize,
) -> usize {
    read_id
        .and_then(|selected| {
            (0..records.len()).find(|&index| records.read_id(index) == Some(selected))
        })
        .unwrap_or_else(|| previous_index.min(records.len().saturating_sub(1)))
}

/// Re-selects an individual alignment after refetch, falling back to the first entry.
pub(super) fn reselect_individual_alignment(
    records: &ViewerRecords,
    selection: Option<&IndividualSelection>,
) -> usize {
    let ViewerRecords::Individual(profiles) = records else {
        return 0;
    };
    selection
        .and_then(|IndividualSelection(selected)| {
            profiles
                .iter()
                .position(|profile| profile.is_same_alignment(selected))
        })
        .unwrap_or(0)
}
