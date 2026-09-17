//! Interactive terminal viewer for indexed BAM files.
//!
//! `libghostty-vt` owns the virtual screen and interprets each ANSI frame.
//! Crossterm only handles the host terminal's raw mode, events, and drawing.

#![expect(
    clippy::print_stderr,
    reason = "command-line errors are intentionally reported to stderr"
)]
#![expect(
    clippy::non_ascii_literal,
    reason = "the terminal plot deliberately uses Unicode single-cell plotting glyphs"
)]
#![expect(
    clippy::pattern_type_mismatch,
    reason = "matching references to view enums is clearer without explicit reference patterns"
)]

use crossterm::{
    SynchronizedUpdate as _,
    cursor::{Hide, MoveTo, Show},
    event::{self, Event, KeyCode, KeyEvent, KeyEventKind, KeyModifiers},
    execute, queue,
    style::{Attribute, Color, Print, SetAttribute, SetBackgroundColor, SetForegroundColor},
    terminal::{
        EndSynchronizedUpdate, EnterAlternateScreen, LeaveAlternateScreen, disable_raw_mode,
        enable_raw_mode,
    },
};
use libghostty_vt::{
    RenderState, Terminal, TerminalOptions,
    render::{CellIterator, RowIterator},
    style::Underline,
};
#[cfg(test)]
use nanalogue_core::{
    ModChar,
    region_sequences::{RegionSequence, RegionSequenceReader},
};
use std::{
    env,
    error::Error,
    io::{self, Stdout, Write as _},
    sync::Arc,
};
#[cfg(test)]
use std::{ffi::OsString, fmt::Write as _, num::NonZeroU32, path::PathBuf, str::FromStr as _};

mod cli;
mod plot;
mod render;
mod state;

use cli::{Args, InitialPosition, USAGE, ViewMode, window_len_for_columns};
use plot::build_individual_frame;
use render::{FrameFooter, build_frame};
#[cfg(test)]
use render::{
    fixed_line, label_column, position_error_footer, position_prompt_footer, sequence_columns,
};
#[cfg(test)]
use state::Viewport;
use state::{Viewer, ViewerRecords, fetch_viewer_records, full_read_label_width, reselect_read};

/// Width reserved for read names and the separating space.
const READ_LABEL_WIDTH: u16 = 19;

/// Maximum number of genomic bases displayed regardless of terminal width.
const MAX_REGION_LENGTH: u32 = 200;

/// Type of the process-wide panic callback retained during a TUI session.
type PanicHook = dyn Fn(&std::panic::PanicHookInfo<'_>) + Send + Sync + 'static;

/// Restores host-terminal state when the interactive session ends.
struct TerminalGuard {
    /// Panic callback that was active before this viewer session.
    previous_panic_hook: Arc<PanicHook>,
}

impl std::fmt::Debug for TerminalGuard {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("TerminalGuard").finish_non_exhaustive()
    }
}

/// Best-effort restoration shared by normal cleanup and the panic hook.
fn restore_terminal() {
    let _screen_result = execute!(
        io::stdout(),
        EndSynchronizedUpdate,
        SetAttribute(Attribute::Reset),
        SetForegroundColor(Color::Reset),
        SetBackgroundColor(Color::Reset),
        Show,
        LeaveAlternateScreen
    );
    let _raw_mode_result = disable_raw_mode();
}

impl TerminalGuard {
    /// Enters raw mode and the alternate screen.
    fn enter(stdout: &mut Stdout) -> io::Result<Self> {
        enable_raw_mode()?;
        let previous_panic_hook = Arc::<PanicHook>::from(std::panic::take_hook());
        let panic_hook = Arc::clone(&previous_panic_hook);
        std::panic::set_hook(Box::new(move |info| {
            restore_terminal();
            panic_hook(info);
        }));
        let guard = Self {
            previous_panic_hook,
        };
        if let Err(error) = execute!(stdout, EnterAlternateScreen, Hide) {
            drop(guard);
            return Err(error);
        }
        Ok(guard)
    }
}

impl Drop for TerminalGuard {
    fn drop(&mut self) {
        restore_terminal();
        if !std::thread::panicking() {
            let previous_panic_hook = Arc::clone(&self.previous_panic_hook);
            std::panic::set_hook(Box::new(move |info| previous_panic_hook(info)));
        }
    }
}

/// Bridges ANSI frames through Ghostty's VT model onto the host terminal.
#[derive(Debug)]
struct GhosttyRenderer {
    /// Ghostty virtual terminal.
    terminal: Terminal<'static, 'static>,
    /// Snapshot state reused between frames.
    render_state: RenderState<'static>,
    /// Reusable row iterator.
    rows: RowIterator<'static>,
    /// Reusable cell iterator.
    cells: CellIterator<'static>,
}

impl GhosttyRenderer {
    /// Creates the VT model at the requested size.
    fn new(cols: u16, rows: u16) -> Result<Self, Box<dyn Error>> {
        Ok(Self {
            terminal: Terminal::new(TerminalOptions {
                cols: cols.max(1),
                rows: rows.max(1),
                max_scrollback: 0,
            })?,
            render_state: RenderState::new()?,
            rows: RowIterator::new()?,
            cells: CellIterator::new()?,
        })
    }

    /// Parses a complete ANSI frame with Ghostty and paints its resulting cells.
    fn draw(
        &mut self,
        stdout: &mut Stdout,
        frame: &str,
        cols: u16,
        rows: u16,
    ) -> Result<(), Box<dyn Error>> {
        let effective_cols = cols.max(1);
        let effective_rows = rows.max(1);
        self.terminal.resize(effective_cols, effective_rows, 1, 1)?;
        self.terminal.vt_write(b"\x1bc\x1b[2J\x1b[H\x1b[?25l");
        self.terminal.vt_write(frame.as_bytes());

        let snapshot = self.render_state.update(&self.terminal)?;
        let mut row_iter = self.rows.update(&snapshot)?;
        stdout.sync_update(|output| -> Result<(), Box<dyn Error>> {
            let mut y = 0u16;
            queue!(output, MoveTo(0, 0))?;
            while let Some(row) = row_iter.next() {
                let mut cell_iter = self.cells.update(row)?;
                let mut x = 0u16;
                while let Some(cell) = cell_iter.next() {
                    let style = cell.style()?;
                    let fg = cell.fg_color()?.map_or(Color::Reset, |color| Color::Rgb {
                        r: color.r,
                        g: color.g,
                        b: color.b,
                    });
                    let bg = cell.bg_color()?.map_or(Color::Reset, |color| Color::Rgb {
                        r: color.r,
                        g: color.g,
                        b: color.b,
                    });
                    queue!(
                        output,
                        MoveTo(x, y),
                        SetAttribute(Attribute::Reset),
                        SetForegroundColor(fg),
                        SetBackgroundColor(bg)
                    )?;
                    if style.bold {
                        queue!(output, SetAttribute(Attribute::Bold))?;
                    }
                    if style.italic {
                        queue!(output, SetAttribute(Attribute::Italic))?;
                    }
                    if style.underline != Underline::None {
                        queue!(output, SetAttribute(Attribute::Underlined))?;
                    }
                    if style.inverse {
                        queue!(output, SetAttribute(Attribute::Reverse))?;
                    }
                    let graphemes = cell.graphemes()?;
                    if graphemes.is_empty() {
                        queue!(output, Print(' '))?;
                    } else {
                        for grapheme in graphemes {
                            queue!(output, Print(grapheme))?;
                        }
                    }
                    x = x.saturating_add(1);
                }
                y = y.saturating_add(1);
            }
            queue!(
                output,
                SetAttribute(Attribute::Reset),
                SetForegroundColor(Color::Reset),
                SetBackgroundColor(Color::Reset)
            )?;
            Ok(())
        })??;
        Ok(())
    }
}

/// Returns whether a key is one of the viewer's standard exit sequences.
fn should_quit(key: KeyEvent) -> bool {
    matches!(key.code, KeyCode::Char('q') | KeyCode::Esc)
        || (key.modifiers.contains(KeyModifiers::CONTROL)
            && matches!(key.code, KeyCode::Char('c' | 'd')))
}

/// Result of editing a genomic position in the footer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PositionPromptOutcome {
    /// The user cancelled the prompt.
    Cancelled,
    /// The user requested that the viewer close.
    Quit,
    /// A valid position was accepted; the value reports whether records changed.
    Navigated(bool),
}

/// Builds the frame appropriate for the active cached-record type.
fn build_viewer_frame(
    viewer: &Viewer,
    cached_records: &ViewerRecords,
    cols: u16,
    rows: u16,
    footer_state: FrameFooter<'_>,
) -> String {
    match cached_records {
        ViewerRecords::Table(records) => build_frame(viewer, records, cols, rows, footer_state),
        ViewerRecords::Individual(profiles) => {
            build_individual_frame(viewer, profiles, cols, rows, footer_state)
        }
    }
}

/// Applies one key press to the goto prompt.
#[expect(
    clippy::wildcard_enum_match_arm,
    reason = "unrelated terminal keys intentionally leave the position prompt unchanged"
)]
fn handle_position_prompt_key(
    input: &mut String,
    input_error: &mut Option<String>,
    key: KeyEvent,
    viewer: &mut Viewer,
) -> Option<PositionPromptOutcome> {
    if key.modifiers.contains(KeyModifiers::CONTROL) && matches!(key.code, KeyCode::Char('c' | 'd'))
    {
        return Some(PositionPromptOutcome::Quit);
    }
    match key.code {
        KeyCode::Esc => Some(PositionPromptOutcome::Cancelled),
        KeyCode::Enter => {
            let parsed_position = InitialPosition::parse(input);
            match parsed_position.and_then(|position| viewer.go_to(&position)) {
                Ok(changed) => Some(PositionPromptOutcome::Navigated(changed)),
                Err(error) => {
                    *input_error = Some(error);
                    None
                }
            }
        }
        KeyCode::Backspace => {
            let _removed_character = input.pop();
            *input_error = None;
            None
        }
        KeyCode::Char(character)
            if !key
                .modifiers
                .intersects(KeyModifiers::CONTROL | KeyModifiers::ALT) =>
        {
            input.push(character);
            *input_error = None;
            None
        }
        _ => None,
    }
}

/// Reads a genomic position while using only the existing footer row as a prompt.
fn prompt_for_position(
    renderer: &mut GhosttyRenderer,
    stdout: &mut Stdout,
    viewer: &mut Viewer,
    records: &ViewerRecords,
) -> Result<PositionPromptOutcome, Box<dyn Error>> {
    let mut input = String::new();
    let mut input_error: Option<String> = None;
    loop {
        let (cols, rows) = crossterm::terminal::size()?;
        let frame = build_viewer_frame(
            viewer,
            records,
            cols,
            rows,
            FrameFooter::PositionPrompt {
                input: &input,
                error: input_error.as_deref(),
            },
        );
        renderer.draw(stdout, &frame, cols, rows)?;

        let Event::Key(key) = event::read()? else {
            continue;
        };
        if key.kind != KeyEventKind::Press {
            continue;
        }
        if let Some(outcome) = handle_position_prompt_key(&mut input, &mut input_error, key, viewer)
        {
            return Ok(outcome);
        }
    }
}

/// Refetches width-dependent table records after a terminal resize.
///
/// Individual profiles contain the complete selected read and are redrawn from the existing cache.
fn handle_terminal_resize(
    viewer: &mut Viewer,
    records: &mut ViewerRecords,
    resized_window_len: u32,
) -> Result<(), Box<dyn Error>> {
    if viewer.mode != ViewMode::Table || resized_window_len == viewer.window_len {
        return Ok(());
    }

    let selected_read_id = records
        .read_id(viewer.viewport.read_offset)
        .map(String::from);
    viewer.window_len = resized_window_len;
    *records = fetch_viewer_records(viewer)?;
    viewer.viewport.read_offset = reselect_read(
        records,
        selected_read_id.as_deref(),
        viewer.viewport.read_offset,
    );
    if viewer.full_read_ids
        && let ViewerRecords::Table(table_records) = records
    {
        viewer.read_label_width = full_read_label_width(table_records);
    }
    Ok(())
}

/// Runs the interactive event loop.
fn run(args: Args) -> Result<(), Box<dyn Error>> {
    let (initial_cols, initial_rows) = crossterm::terminal::size()?;
    let mut viewer = Viewer::open_mode(
        args.bam,
        &args.position,
        args.mod_type,
        args.mode,
        window_len_for_columns(initial_cols),
    )?;
    let mut renderer = GhosttyRenderer::new(initial_cols, initial_rows)?;
    let mut stdout = io::stdout();
    let _guard = TerminalGuard::enter(&mut stdout)?;
    let mut records = fetch_viewer_records(&mut viewer)?;

    loop {
        let (cols, rows) = crossterm::terminal::size()?;
        let resized_window_len = window_len_for_columns(cols);
        handle_terminal_resize(&mut viewer, &mut records, resized_window_len)?;
        let visible_reads = if viewer.mode == ViewMode::Table {
            usize::from(rows.saturating_sub(4))
        } else {
            1
        };
        viewer.viewport.read_offset = viewer
            .viewport
            .read_offset
            .min(records.len().saturating_sub(visible_reads));
        let frame = build_viewer_frame(&viewer, &records, cols, rows, FrameFooter::Controls);
        renderer.draw(&mut stdout, &frame, cols, rows)?;

        let Event::Key(key) = event::read()? else {
            continue;
        };
        if key.kind != KeyEventKind::Press {
            continue;
        }
        if should_quit(key) {
            break;
        }
        if key.code == KeyCode::Char('g') {
            match prompt_for_position(&mut renderer, &mut stdout, &mut viewer, &records)? {
                PositionPromptOutcome::Cancelled => {}
                PositionPromptOutcome::Quit => break,
                PositionPromptOutcome::Navigated(changed) => {
                    if changed {
                        records = fetch_viewer_records(&mut viewer)?;
                    }
                }
            }
            continue;
        }
        let refetch = match &records {
            ViewerRecords::Table(table_records) => {
                viewer.handle_key(key.code, table_records, visible_reads)
            }
            ViewerRecords::Individual(profiles) => {
                viewer.handle_individual_key(key.code, profiles.len())
            }
        };
        if refetch {
            records = fetch_viewer_records(&mut viewer)?;
        }
    }
    Ok(())
}

/// Parses arguments and reports errors after terminal cleanup.
fn main() {
    let args = match Args::parse_from(env::args_os().skip(1)) {
        Ok(Some(args)) => args,
        Ok(None) => {
            let _write_result = writeln!(io::stdout().lock(), "{USAGE}");
            return;
        }
        Err(error) => {
            eprintln!("Error: {error}\n{USAGE}");
            std::process::exit(2);
        }
    };
    if let Err(error) = run(args) {
        eprintln!("Error: {error}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plot::{abbreviated_individual_read_label, individual_status};
    use nanalogue_core::{
        simulate_mod_bam::{AlignmentFormat, SimulationConfig, TempBamSimulation},
        uuid, write_bam_denovo,
    };
    use rust_htslib::bam::{
        self,
        record::{Aux, Cigar, CigarString},
    };
    use std::path::Path;

    const DEMO_GOLDEN_COLS: u16 = 90;
    const DEMO_GOLDEN_ROWS: u16 = 20;
    const DEMO_VISIBLE_READS: usize = 16;

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct CapturedStyle {
        style: libghostty_vt::style::Style,
        foreground: Option<libghostty_vt::style::RgbColor>,
        background: Option<libghostty_vt::style::RgbColor>,
    }

    impl CapturedStyle {
        fn write_ansi(self, output: &mut String) {
            output.push_str("\x1b[0m");
            let mut codes = Vec::new();
            codes.extend(
                [
                    (self.style.bold, "1"),
                    (self.style.faint, "2"),
                    (self.style.italic, "3"),
                    (self.style.underline != Underline::None, "4"),
                    (self.style.blink, "5"),
                    (self.style.inverse, "7"),
                    (self.style.invisible, "8"),
                    (self.style.strikethrough, "9"),
                    (self.style.overline, "53"),
                ]
                .into_iter()
                .filter(|&(enabled, _code)| enabled)
                .map(|(_enabled, code)| String::from(code)),
            );
            if let Some(color) = self.foreground {
                codes.push(format!("38;2;{};{};{}", color.r, color.g, color.b));
            }
            if let Some(color) = self.background {
                codes.push(format!("48;2;{};{};{}", color.r, color.g, color.b));
            }
            if !codes.is_empty() {
                write!(output, "\x1b[{}m", codes.join(";")).expect("writing to String cannot fail");
            }
        }
    }

    fn viewport_as_ansi<'alloc>(
        snapshot: &libghostty_vt::render::Snapshot<'alloc, '_>,
        row_iterator: &mut RowIterator<'alloc>,
        cell_iterator: &mut CellIterator<'alloc>,
    ) -> Result<String, Box<dyn Error>> {
        let mut output = String::new();
        let mut rows = row_iterator.update(snapshot)?;
        let mut row_number = 1usize;
        let mut loop_iterations = 0usize;
        while let Some(row) = rows.next() {
            loop_iterations = loop_iterations.saturating_add(1);
            assert!(
                loop_iterations <= 5_000_000,
                "ANSI viewport capture exceeds five million loop iterations"
            );
            write!(output, "\x1b[{row_number};1H").expect("writing to String cannot fail");
            let mut cells = cell_iterator.update(row)?;
            let mut current_style = None;
            while let Some(cell) = cells.next() {
                loop_iterations = loop_iterations.saturating_add(1);
                assert!(
                    loop_iterations <= 5_000_000,
                    "ANSI viewport capture exceeds five million loop iterations"
                );
                let captured_style = CapturedStyle {
                    style: cell.style()?,
                    foreground: cell.fg_color()?,
                    background: cell.bg_color()?,
                };
                if current_style != Some(captured_style) {
                    captured_style.write_ansi(&mut output);
                    current_style = Some(captured_style);
                }
                let graphemes = cell.graphemes()?;
                if graphemes.is_empty() {
                    output.push(' ');
                } else {
                    output.extend(graphemes);
                }
            }
            output.push_str("\x1b[0m");
            row_number = row_number.saturating_add(1);
        }
        Ok(output)
    }

    fn render_ansi_viewport(
        viewer: &Viewer,
        records: &[RegionSequence],
        cols: u16,
        rows: u16,
        footer_state: FrameFooter<'_>,
    ) -> Result<String, Box<dyn Error>> {
        let frame = build_frame(viewer, records, cols, rows, footer_state);
        render_frame_as_ansi(&frame, cols, rows)
    }

    fn render_frame_as_ansi(frame: &str, cols: u16, rows: u16) -> Result<String, Box<dyn Error>> {
        let mut terminal = Terminal::new(TerminalOptions {
            cols,
            rows,
            max_scrollback: 0,
        })?;
        terminal.vt_write(frame.as_bytes());
        let mut render_state = RenderState::new()?;
        let snapshot = render_state.update(&terminal)?;
        let mut row_iterator = RowIterator::new()?;
        let mut cell_iterator = CellIterator::new()?;
        let actual = viewport_as_ansi(&snapshot, &mut row_iterator, &mut cell_iterator)?;

        let mut replayed_terminal = Terminal::new(TerminalOptions {
            cols,
            rows,
            max_scrollback: 0,
        })?;
        replayed_terminal.vt_write(actual.as_bytes());
        let mut replayed_render_state = RenderState::new()?;
        let replayed_snapshot = replayed_render_state.update(&replayed_terminal)?;
        let mut replayed_row_iterator = RowIterator::new()?;
        let mut replayed_cell_iterator = CellIterator::new()?;
        let replayed = viewport_as_ansi(
            &replayed_snapshot,
            &mut replayed_row_iterator,
            &mut replayed_cell_iterator,
        )?;
        assert_eq!(actual, replayed, "ANSI golden must reproduce its viewport");
        Ok(actual)
    }

    fn assert_ansi_golden(name: &str, actual: &str) -> Result<(), Box<dyn Error>> {
        let golden_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/goldens")
            .join(name);
        let update_golden = env::var("NANALOGUE_UPDATE_GOLDENS").as_deref() == Ok("1");
        assert!(
            !update_golden || env::var_os("CI").is_none(),
            "golden files must not be updated in CI"
        );
        if update_golden {
            std::fs::write(&golden_path, actual.as_bytes())?;
        }
        let expected = std::fs::read(&golden_path)?;
        assert_eq!(actual.as_bytes(), expected);
        Ok(())
    }

    fn render_demo_viewport(
        viewer: &Viewer,
        records: &[RegionSequence],
        footer_state: FrameFooter<'_>,
    ) -> Result<String, Box<dyn Error>> {
        render_ansi_viewport(
            viewer,
            records,
            DEMO_GOLDEN_COLS,
            DEMO_GOLDEN_ROWS,
            footer_state,
        )
    }

    fn handle_demo_key(viewer: &mut Viewer, key: KeyCode, records: &[RegionSequence]) -> bool {
        viewer.handle_key(key, records, DEMO_VISIBLE_READS)
    }

    fn write_zero_sequence_viewer_bam(
        simulation: &TempBamSimulation,
    ) -> Result<PathBuf, Box<dyn Error>> {
        let mut record = bam::Record::new();
        record.set_tid(0);
        record.set_pos(0);
        record.set_mapq(60);
        record.unset_unmapped();
        record.set(
            b"sequence-not-stored",
            Some(&CigarString::from(vec![Cigar::Match(5)])),
            b"",
            &[],
        );
        record.push_aux(b"MM", Aux::String("A+a?,0;"))?;
        record.push_aux(b"ML", Aux::ArrayU8((&[200][..]).into()))?;
        let path = Path::new(simulation.bam_path())
            .parent()
            .expect("simulation BAM has an owning directory")
            .join("zero-sequence.bam");
        write_bam_denovo(
            [record],
            [(String::from("chr1"), 20)],
            [String::from("rg1")],
            Vec::<String>::new(),
            &path,
        )?;
        Ok(path)
    }

    fn profile_record(
        read_id: &[u8],
        tid: i32,
        start: i64,
        length: u32,
        reverse: bool,
        probabilities: Option<&[u8]>,
    ) -> Result<bam::Record, Box<dyn Error>> {
        let mut record = bam::Record::new();
        record.set_tid(tid);
        record.set_pos(start);
        record.set_mapq(60);
        record.unset_unmapped();
        if reverse {
            record.set_reverse();
        }
        let base = if reverse { b'T' } else { b'A' };
        record.set(
            read_id,
            Some(&CigarString::from(vec![Cigar::Match(length)])),
            &vec![base; usize::try_from(length)?],
            &vec![30; usize::try_from(length)?],
        );
        if let Some(values) = probabilities {
            let deltas = std::iter::repeat_n("0", values.len())
                .collect::<Vec<_>>()
                .join(",");
            record.push_aux(b"MM", Aux::String(&format!("A+a?,{deltas};")))?;
            record.push_aux(b"ML", Aux::ArrayU8(values.into()))?;
        }
        Ok(record)
    }

    fn write_individual_semantics_bam() -> Result<PathBuf, Box<dyn Error>> {
        let mut reverse_duplicate = profile_record(
            b"duplicate",
            0,
            1,
            30,
            true,
            Some(&[255, 0, 255, 0, 255, 0]),
        )?;
        reverse_duplicate.set_secondary();
        let records = [
            profile_record(
                b"duplicate",
                0,
                0,
                30,
                false,
                Some(&[255, 0, 255, 0, 255, 0]),
            )?,
            reverse_duplicate,
            profile_record(b"no-calls", 0, 2, 30, false, None)?,
            profile_record(b"goto-target", 1, 4, 20, false, Some(&[0, 255, 128]))?,
        ];
        let path = env::temp_dir().join(format!("{}.bam", uuid::v4_random()));
        write_bam_denovo(
            records,
            [(String::from("first"), 50), (String::from("second"), 50)],
            [String::from("rg1")],
            Vec::<String>::new(),
            &path,
        )?;
        Ok(path)
    }

    fn remove_viewer_test_bam(path: &Path) -> Result<(), Box<dyn Error>> {
        std::fs::remove_file(path)?;
        std::fs::remove_file(format!("{}.bai", path.display()))?;
        Ok(())
    }

    fn individual_semantics_viewer(path: &Path) -> Result<Viewer, Box<dyn Error>> {
        Viewer::open_mode(
            path.to_path_buf(),
            &InitialPosition {
                contig: String::from("first"),
                start: 5,
            },
            Some(ModChar::new('a')),
            ViewMode::Individual {
                win: NonZeroU32::new(3).expect("non-zero"),
            },
            5,
        )
    }

    fn goto_test_viewer() -> Viewer {
        let position = InitialPosition {
            contig: String::from("dummyIII"),
            start: 23,
        };
        let mut viewer = Viewer::open(PathBuf::from("examples/example_1.bam"), &position, None, 7)
            .expect("position should open");
        viewer.viewport.read_offset = 1;
        viewer.full_read_ids = true;
        viewer.read_label_width = 40;
        viewer.show_insertions = true;
        viewer
    }

    fn prompt_key(
        input: &mut String,
        input_error: &mut Option<String>,
        code: KeyCode,
        viewer: &mut Viewer,
    ) -> Option<PositionPromptOutcome> {
        handle_position_prompt_key(
            input,
            input_error,
            KeyEvent::new(code, KeyModifiers::NONE),
            viewer,
        )
    }

    #[test]
    fn navigation_clamps_each_axis() {
        let mut viewport = Viewport {
            tid: 0,
            start: 0,
            read_offset: 4,
        };
        viewport.navigate(KeyCode::Down, 100, 20, 7, 3);
        assert_eq!(viewport.start, 0);
        assert_eq!(viewport.read_offset, 4);

        viewport.navigate(KeyCode::Right, 100, 20, 7, 3);
        assert_eq!(viewport.start, 20);
        assert_eq!(viewport.read_offset, 0);

        viewport.navigate(KeyCode::Left, 100, 20, 7, 3);
        viewport.navigate(KeyCode::Up, 100, 20, 7, 3);
        assert_eq!(viewport.start, 0);
        assert_eq!(viewport.read_offset, 0);

        viewport.navigate(KeyCode::Right, 1_000, 20, 7, 3);
        assert_eq!(viewport.start, 20);

        viewport.read_offset = 7;
        viewport.navigate(KeyCode::PageUp, 1_000, 20, 20, 4);
        assert_eq!(viewport.read_offset, 3);
        viewport.navigate(KeyCode::PageUp, 1_000, 20, 20, 4);
        assert_eq!(viewport.read_offset, 0);
        viewport.navigate(KeyCode::PageDown, 1_000, 20, 20, 4);
        assert_eq!(viewport.read_offset, 4);
        viewport.navigate(KeyCode::PageDown, 1_000, 20, 20, 4);
        assert_eq!(viewport.read_offset, 8);
        viewport.read_offset = 15;
        viewport.navigate(KeyCode::PageDown, 1_000, 20, 20, 4);
        assert_eq!(viewport.read_offset, 16);
        viewport.navigate(KeyCode::Home, 1_000, 20, 20, 4);
        assert_eq!(viewport.read_offset, 0);
        viewport.navigate(KeyCode::End, 1_000, 20, 20, 4);
        assert_eq!(viewport.read_offset, 16);
    }

    #[test]
    fn viewer_uses_nanalogue_region_sequences() {
        let mut reader = RegionSequenceReader::from_path("examples/example_1.bam")
            .expect("open indexed example");
        assert!(
            reader
                .sequences(2, 20, 30, None)
                .expect("retrieve a partially covered region")
                .is_empty(),
            "a read that does not span the full region must be excluded"
        );
        let rows = reader
            .sequences(2, 23, 30, None)
            .expect("retrieve nanalogue region sequences");
        let row = rows.first().expect("one overlapping read");
        assert_eq!(row.read_id(), "a4f36092-b4d5-47a9-813e-c22c3b477a0c");
        assert_eq!(sequence_columns(row, 10, false), "ACATCAA   ");
    }

    #[test]
    fn table_reselection_clamps_when_qname_disappears() {
        let position = InitialPosition {
            contig: String::from("dummyIII"),
            start: 23,
        };
        let mut viewer = Viewer::open(PathBuf::from("examples/example_1.bam"), &position, None, 7)
            .expect("position should open");
        let record = viewer
            .visible_records()
            .expect("records should load")
            .into_iter()
            .next()
            .expect("one record");
        let records = ViewerRecords::Table(vec![record.clone(), record]);

        assert_eq!(
            reselect_read(&records, Some("read-that-no-longer-spans"), usize::MAX),
            1
        );
    }

    #[test]
    fn viewer_bolds_only_requested_high_probability_modifications() {
        let mut reader = RegionSequenceReader::from_path("examples/example_1.bam")
            .expect("open indexed example");
        let unstyled_rows = reader
            .sequences(2, 23, 30, None)
            .expect("retrieve sequences without modification parsing");
        assert!(
            unstyled_rows
                .first()
                .expect("one overlapping read")
                .modifications()
                .iter()
                .all(|modified| !modified)
        );

        let styled_rows = reader
            .sequences(2, 23, 30, Some(ModChar::new('T')))
            .expect("retrieve thresholded T modifications");
        let row = styled_rows.first().expect("one overlapping read");
        assert_eq!(
            row.modifications(),
            [false, false, false, true, false, false, false]
        );
        assert_eq!(
            sequence_columns(row, 10, false),
            "ACA\x1b[1;4mT\x1b[22;24mCAA   "
        );
    }

    #[test]
    fn individual_resize_reuses_cached_duplicate_qname_profiles() -> Result<(), Box<dyn Error>> {
        let path = write_individual_semantics_bam()?;
        let mut viewer = individual_semantics_viewer(&path)?;
        let mut records = fetch_viewer_records(&mut viewer)?;
        let ViewerRecords::Individual(profiles) = &records else {
            return Err("individual viewer must fetch profiles".into());
        };
        assert_eq!(profiles.len(), 3);
        assert_eq!(
            abbreviated_individual_read_label(profiles, 0, usize::MAX).as_deref(),
            Some("duplicate#1")
        );
        assert_eq!(
            abbreviated_individual_read_label(profiles, 1, usize::MAX).as_deref(),
            Some("duplicate#2")
        );
        assert_eq!(
            abbreviated_individual_read_label(profiles, 1, 8).as_deref(),
            Some("dupli~#2")
        );
        viewer.viewport.read_offset = 1;
        assert!(individual_status(&viewer, profiles, 100).contains("duplicate#2"));
        let ten_duplicates =
            std::iter::repeat_n(profiles.first().expect("first duplicate").clone(), 10)
                .collect::<Vec<_>>();
        viewer.viewport.read_offset = 9;
        for width in 1..=100 {
            let status = individual_status(&viewer, &ten_duplicates, width);
            assert!(
                !status.contains("mods") || status.contains("#10"),
                "a detailed width-{width} status must show the complete ordinal: {status:?}"
            );
        }
        viewer.viewport.read_offset = 1;
        // An individual fetch now fails, so a successful resize proves that it reused the cache.
        viewer.mod_type = None;
        handle_terminal_resize(&mut viewer, &mut records, 10)?;
        assert_eq!(viewer.mod_type, None);
        assert_eq!(viewer.window_len, 5);
        assert_eq!(viewer.viewport.read_offset, 1);
        let ViewerRecords::Individual(resized_profiles) = &records else {
            return Err("individual resize must preserve cached profiles".into());
        };
        assert_eq!(resized_profiles.len(), 3);
        let selected_calls = resized_profiles
            .get(viewer.viewport.read_offset)
            .expect("selected profile")
            .calls();
        assert_eq!(
            selected_calls,
            [(25, 0), (26, 255), (27, 0), (28, 255), (29, 0), (30, 255)],
            "the reverse duplicate record must remain selected"
        );

        remove_viewer_test_bam(&path)
    }

    #[test]
    fn ghostty_preserves_modified_base_columns_and_styles() -> Result<(), Box<dyn Error>> {
        let position = InitialPosition {
            contig: String::from("dummyIII"),
            start: 23,
        };
        let mut viewer = Viewer::open(
            PathBuf::from("examples/example_1.bam"),
            &position,
            Some(ModChar::new('T')),
            7,
        )?;
        let records = viewer.visible_records()?;
        let frame = build_frame(&viewer, &records, 26, 6, FrameFooter::Controls);
        let mut terminal = Terminal::new(TerminalOptions {
            cols: 26,
            rows: 6,
            max_scrollback: 0,
        })?;
        terminal.vt_write(frame.as_bytes());
        let mut render_state = RenderState::new()?;
        let snapshot = render_state.update(&terminal)?;
        let mut row_iterator = RowIterator::new()?;
        let mut rows = row_iterator.update(&snapshot)?;
        for _row_index in 0..4 {
            let _row = rows.next().expect("read row should exist");
        }
        let mut cell_iterator = CellIterator::new()?;
        let mut cells = cell_iterator.update(&rows)?;

        cells.select(21)?;
        assert_eq!(cells.graphemes()?, ['A']);
        assert!(!cells.style()?.bold);
        assert_eq!(cells.style()?.underline, Underline::None);
        cells.select(22)?;
        assert_eq!(cells.graphemes()?, ['T']);
        assert!(cells.style()?.bold);
        assert_eq!(cells.style()?.underline, Underline::Single);
        cells.select(23)?;
        assert_eq!(cells.graphemes()?, ['C']);
        assert!(!cells.style()?.bold);
        assert_eq!(cells.style()?.underline, Underline::None);
        Ok(())
    }

    #[test]
    fn visible_viewport_matches_ansi_golden() -> Result<(), Box<dyn Error>> {
        let position = InitialPosition {
            contig: String::from("dummyIII"),
            start: 23,
        };
        let mut viewer = Viewer::open(
            PathBuf::from("examples/example_1.bam"),
            &position,
            Some(ModChar::new('T')),
            7,
        )?;
        let records = viewer.visible_records()?;
        let actual = render_ansi_viewport(&viewer, &records, 26, 6, FrameFooter::Controls)?;
        assert_ansi_golden("bam_viewer_visible.ansi", &actual)?;
        Ok(())
    }

    fn assert_end_and_goto_goldens(simulation: &TempBamSimulation) -> Result<(), Box<dyn Error>> {
        let initial_position = InitialPosition {
            contig: String::from("contig_00000"),
            start: 25,
        };

        for (name_suffix, mod_type) in [("mods", Some(ModChar::new('m'))), ("no_mods", None)] {
            let mut viewer = Viewer::open(
                PathBuf::from(simulation.bam_path()),
                &initial_position,
                mod_type,
                window_len_for_columns(DEMO_GOLDEN_COLS),
            )?;
            viewer.path = PathBuf::from("nanalogue-viewer-demo.bam");
            let initial_records = viewer.visible_records()?;
            assert_eq!(initial_records.len(), 187);
            assert_eq!(
                initial_records.iter().any(|record| {
                    record
                        .modifications_with_insertions()
                        .iter()
                        .any(|modified| *modified)
                }),
                mod_type.is_some()
            );

            assert!(!handle_demo_key(
                &mut viewer,
                KeyCode::Char('r'),
                &initial_records
            ));
            assert!(!handle_demo_key(
                &mut viewer,
                KeyCode::Char('i'),
                &initial_records
            ));
            assert!(!handle_demo_key(
                &mut viewer,
                KeyCode::End,
                &initial_records
            ));
            assert_eq!(viewer.viewport.read_offset, 171);
            let end_viewport =
                render_demo_viewport(&viewer, &initial_records, FrameFooter::Controls)?;
            assert_ansi_golden(
                &format!("bam_viewer_end_key_{name_suffix}.ansi"),
                &end_viewport,
            )?;

            if mod_type.is_some() {
                assert!(!handle_demo_key(
                    &mut viewer,
                    KeyCode::Home,
                    &initial_records
                ));
                assert_eq!(viewer.viewport.read_offset, 0);
                assert!(viewer.full_read_ids);
                assert!(viewer.show_insertions);
                let home_viewport =
                    render_demo_viewport(&viewer, &initial_records, FrameFooter::Controls)?;
                assert_ansi_golden("bam_viewer_key_home.ansi", &home_viewport)?;
            }

            assert!(viewer.go_to(&InitialPosition {
                contig: String::from("contig_00000"),
                start: 45,
            })?);
            assert!(!viewer.full_read_ids);
            assert!(!viewer.show_insertions);
            let goto_records = viewer.visible_records()?;
            assert_eq!(goto_records.len(), 187);
            assert_eq!(
                goto_records
                    .iter()
                    .any(|record| { record.modifications().iter().any(|modified| *modified) }),
                mod_type.is_some()
            );
            let goto_viewport =
                render_demo_viewport(&viewer, &goto_records, FrameFooter::Controls)?;
            assert_ansi_golden(
                &format!("bam_viewer_goto_{name_suffix}.ansi"),
                &goto_viewport,
            )?;
        }
        Ok(())
    }

    fn assert_display_key_goldens(simulation: &TempBamSimulation) -> Result<(), Box<dyn Error>> {
        let transition_position = InitialPosition {
            contig: String::from("contig_00000"),
            start: 365,
        };
        let mut viewer = Viewer::open(
            PathBuf::from(simulation.bam_path()),
            &transition_position,
            Some(ModChar::new('m')),
            window_len_for_columns(DEMO_GOLDEN_COLS),
        )?;
        viewer.path = PathBuf::from("nanalogue-viewer-demo.bam");
        let records = viewer.visible_records()?;
        assert_eq!(records.len(), 187);
        assert!(records.iter().any(|record| {
            record
                .sequence_with_insertions()
                .bytes()
                .any(|base| base.is_ascii_lowercase())
        }));

        let default_viewport = render_demo_viewport(&viewer, &records, FrameFooter::Controls)?;
        assert_ansi_golden("bam_viewer_key_default.ansi", &default_viewport)?;

        assert!(!handle_demo_key(&mut viewer, KeyCode::Char('i'), &records));
        assert!(!handle_demo_key(&mut viewer, KeyCode::Char('i'), &records));
        assert!(!viewer.show_insertions);
        assert_eq!(
            render_demo_viewport(&viewer, &records, FrameFooter::Controls)?,
            default_viewport
        );

        assert!(!handle_demo_key(&mut viewer, KeyCode::Char('r'), &records));
        assert!(!handle_demo_key(&mut viewer, KeyCode::Char('r'), &records));
        assert!(!viewer.full_read_ids);
        assert_eq!(
            render_demo_viewport(&viewer, &records, FrameFooter::Controls)?,
            default_viewport
        );

        assert!(!handle_demo_key(&mut viewer, KeyCode::Down, &records));
        assert_eq!(viewer.viewport.read_offset, 1);
        assert!(!handle_demo_key(&mut viewer, KeyCode::Up, &records));
        assert_eq!(viewer.viewport.read_offset, 0);
        assert_eq!(
            render_demo_viewport(&viewer, &records, FrameFooter::Controls)?,
            default_viewport
        );

        assert!(!handle_demo_key(&mut viewer, KeyCode::Char('i'), &records));
        let insertion_viewport = render_demo_viewport(&viewer, &records, FrameFooter::Controls)?;
        assert_ne!(insertion_viewport, default_viewport);
        assert_ansi_golden("bam_viewer_key_insertions.ansi", &insertion_viewport)?;

        assert!(!handle_demo_key(&mut viewer, KeyCode::Char('r'), &records));
        let full_id_viewport = render_demo_viewport(&viewer, &records, FrameFooter::Controls)?;
        assert_ne!(full_id_viewport, insertion_viewport);
        assert_ansi_golden("bam_viewer_key_full_ids.ansi", &full_id_viewport)?;

        assert!(!handle_demo_key(&mut viewer, KeyCode::PageDown, &records));
        assert_eq!(viewer.viewport.read_offset, 16);
        assert!(viewer.full_read_ids);
        assert!(viewer.show_insertions);
        let page_down_viewport = render_demo_viewport(&viewer, &records, FrameFooter::Controls)?;
        assert_ansi_golden("bam_viewer_key_page_down.ansi", &page_down_viewport)?;

        assert!(!handle_demo_key(&mut viewer, KeyCode::PageUp, &records));
        assert_eq!(viewer.viewport.read_offset, 0);
        assert!(viewer.full_read_ids);
        assert!(viewer.show_insertions);
        let page_up_viewport = render_demo_viewport(&viewer, &records, FrameFooter::Controls)?;
        assert_eq!(page_up_viewport, full_id_viewport);
        assert_ansi_golden("bam_viewer_key_full_ids.ansi", &page_up_viewport)?;

        assert!(!handle_demo_key(&mut viewer, KeyCode::PageDown, &records));
        assert_eq!(viewer.viewport.read_offset, 16);
        assert!(handle_demo_key(&mut viewer, KeyCode::Right, &records));
        assert_eq!(viewer.viewport.start, 436);
        assert_eq!(viewer.viewport.read_offset, 0);
        assert!(!viewer.full_read_ids);
        assert!(!viewer.show_insertions);
        let right_records = viewer.visible_records()?;
        let right_viewport = render_demo_viewport(&viewer, &right_records, FrameFooter::Controls)?;
        assert_ansi_golden("bam_viewer_key_right.ansi", &right_viewport)?;

        assert!(handle_demo_key(&mut viewer, KeyCode::Left, &right_records));
        assert_eq!(viewer.viewport.start, 365);
        assert_eq!(viewer.viewport.read_offset, 0);
        assert!(!viewer.full_read_ids);
        assert!(!viewer.show_insertions);
        let left_records = viewer.visible_records()?;
        let left_viewport = render_demo_viewport(&viewer, &left_records, FrameFooter::Controls)?;
        assert_eq!(left_viewport, default_viewport);
        assert_ansi_golden("bam_viewer_key_default.ansi", &left_viewport)?;
        Ok(())
    }

    fn assert_successful_prompt_goto(
        viewer: &mut Viewer,
        records: &[RegionSequence],
    ) -> Result<(), Box<dyn Error>> {
        assert!(!handle_demo_key(viewer, KeyCode::Char('r'), records));
        assert!(!handle_demo_key(viewer, KeyCode::Char('i'), records));
        let mut input = String::new();
        let mut input_error = None;
        for character in "contig_00001:120".chars() {
            assert_eq!(
                prompt_key(
                    &mut input,
                    &mut input_error,
                    KeyCode::Char(character),
                    viewer,
                ),
                None
            );
        }
        assert_eq!(
            prompt_key(&mut input, &mut input_error, KeyCode::Enter, viewer),
            Some(PositionPromptOutcome::Navigated(true))
        );
        assert_eq!(viewer.target_name(), "contig_00001");
        assert_eq!(viewer.viewport.start, 120);
        assert_eq!(viewer.viewport.read_offset, 0);
        assert!(!viewer.full_read_ids);
        assert!(!viewer.show_insertions);
        let goto_records = viewer.visible_records()?;
        let goto_viewport = render_demo_viewport(viewer, &goto_records, FrameFooter::Controls)?;
        assert_ansi_golden("bam_viewer_goto_success.ansi", &goto_viewport)?;
        Ok(())
    }

    fn assert_goto_prompt_goldens(simulation: &TempBamSimulation) -> Result<(), Box<dyn Error>> {
        let prompt_position = InitialPosition {
            contig: String::from("contig_00000"),
            start: 25,
        };
        let mut viewer = Viewer::open(
            PathBuf::from(simulation.bam_path()),
            &prompt_position,
            Some(ModChar::new('m')),
            window_len_for_columns(DEMO_GOLDEN_COLS),
        )?;
        viewer.path = PathBuf::from("nanalogue-viewer-demo.bam");
        let records = viewer.visible_records()?;
        let before_prompt = render_demo_viewport(&viewer, &records, FrameFooter::Controls)?;
        let mut input = String::new();
        let mut input_error = None;
        let prompt_viewport = render_demo_viewport(
            &viewer,
            &records,
            FrameFooter::PositionPrompt {
                input: &input,
                error: input_error.as_deref(),
            },
        )?;
        assert_ansi_golden("bam_viewer_goto_prompt.ansi", &prompt_viewport)?;

        for character in "missing:0".chars() {
            assert_eq!(
                prompt_key(
                    &mut input,
                    &mut input_error,
                    KeyCode::Char(character),
                    &mut viewer,
                ),
                None
            );
        }
        assert_eq!(
            prompt_key(&mut input, &mut input_error, KeyCode::Enter, &mut viewer),
            None
        );
        assert!(input_error.is_some());
        let error_viewport = render_demo_viewport(
            &viewer,
            &records,
            FrameFooter::PositionPrompt {
                input: &input,
                error: input_error.as_deref(),
            },
        )?;
        assert_ansi_golden("bam_viewer_goto_error.ansi", &error_viewport)?;
        let narrow_error_viewport = render_ansi_viewport(
            &viewer,
            &records,
            27,
            10,
            FrameFooter::PositionPrompt {
                input: &input,
                error: input_error.as_deref(),
            },
        )?;
        assert_ansi_golden("bam_viewer_narrow_error.ansi", &narrow_error_viewport)?;

        assert_eq!(
            prompt_key(
                &mut input,
                &mut input_error,
                KeyCode::Backspace,
                &mut viewer,
            ),
            None
        );
        assert_eq!(input, "missing:");
        assert!(input_error.is_none());
        let correcting_viewport = render_demo_viewport(
            &viewer,
            &records,
            FrameFooter::PositionPrompt {
                input: &input,
                error: input_error.as_deref(),
            },
        )?;
        assert_ansi_golden("bam_viewer_goto_correcting.ansi", &correcting_viewport)?;

        assert_eq!(
            prompt_key(&mut input, &mut input_error, KeyCode::Esc, &mut viewer),
            Some(PositionPromptOutcome::Cancelled)
        );
        let after_cancel = render_demo_viewport(&viewer, &records, FrameFooter::Controls)?;
        assert_eq!(after_cancel, before_prompt);
        assert_successful_prompt_goto(&mut viewer, &records)
    }

    fn assert_empty_and_contig_end_goldens() -> Result<(), Box<dyn Error>> {
        let mut empty_viewer = Viewer::open(
            PathBuf::from("examples/example_1.bam"),
            &InitialPosition {
                contig: String::from("dummyIII"),
                start: 20,
            },
            None,
            10,
        )?;
        let empty_records = empty_viewer.visible_records()?;
        assert!(empty_records.is_empty());
        let empty_viewport =
            render_demo_viewport(&empty_viewer, &empty_records, FrameFooter::Controls)?;
        assert_ansi_golden("bam_viewer_no_reads.ansi", &empty_viewport)?;

        let mut end_viewer = Viewer::open(
            PathBuf::from("examples/example_1.bam"),
            &InitialPosition {
                contig: String::from("dummyI"),
                start: 20,
            },
            None,
            10,
        )?;
        assert_eq!(end_viewer.current_window_len(), 2);
        let end_records = end_viewer.visible_records()?;
        let end_viewport = render_demo_viewport(&end_viewer, &end_records, FrameFooter::Controls)?;
        assert_ansi_golden("bam_viewer_contig_end.ansi", &end_viewport)?;
        Ok(())
    }

    fn assert_zero_sequence_golden(simulation: &TempBamSimulation) -> Result<(), Box<dyn Error>> {
        let path = write_zero_sequence_viewer_bam(simulation)?;
        let mut viewer = Viewer::open(
            path.clone(),
            &InitialPosition {
                contig: String::from("chr1"),
                start: 0,
            },
            Some(ModChar::new('a')),
            5,
        )?;
        viewer.path = PathBuf::from("zero-sequence.bam");
        let records = viewer.visible_records()?;
        assert_eq!(records.len(), 1);
        assert_eq!(records.first().expect("one record").sequence(), "*");
        let viewport = render_ansi_viewport(&viewer, &records, 50, 8, FrameFooter::Controls)?;
        assert_ansi_golden("bam_viewer_zero_sequence.ansi", &viewport)?;
        Ok(())
    }

    fn assert_terminal_size_goldens(simulation: &TempBamSimulation) -> Result<(), Box<dyn Error>> {
        let mut viewer = Viewer::open(
            PathBuf::from(simulation.bam_path()),
            &InitialPosition {
                contig: String::from("contig_00000"),
                start: 365,
            },
            Some(ModChar::new('m')),
            window_len_for_columns(DEMO_GOLDEN_COLS),
        )?;
        viewer.path = PathBuf::from("nanalogue-viewer-demo.bam");
        let records = viewer.visible_records()?;
        let narrow_prompt = render_ansi_viewport(
            &viewer,
            &records,
            27,
            10,
            FrameFooter::PositionPrompt {
                input: "contig_name_that_is_far_too_long:123456789",
                error: None,
            },
        )?;
        assert_ansi_golden("bam_viewer_narrow_prompt.ansi", &narrow_prompt)?;

        let short_terminal = render_ansi_viewport(
            &viewer,
            &records,
            DEMO_GOLDEN_COLS,
            4,
            FrameFooter::Controls,
        )?;
        assert_ansi_golden("bam_viewer_short_terminal.ansi", &short_terminal)?;
        Ok(())
    }

    #[test]
    fn demo_navigation_viewports_match_ansi_goldens() -> Result<(), Box<dyn Error>> {
        let config: SimulationConfig = serde_json::from_str(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/examples/bam_viewer_demo.json"
        )))?;
        let simulation = TempBamSimulation::new(config, AlignmentFormat::Bam)?;
        assert_end_and_goto_goldens(&simulation)?;
        assert_display_key_goldens(&simulation)?;
        assert_goto_prompt_goldens(&simulation)?;
        assert_terminal_size_goldens(&simulation)?;
        assert_empty_and_contig_end_goldens()?;
        assert_zero_sequence_golden(&simulation)?;
        Ok(())
    }

    #[test]
    fn individual_reverse_duplicate_alignments_are_independently_selectable()
    -> Result<(), Box<dyn Error>> {
        let path = write_individual_semantics_bam()?;
        let mut viewer = individual_semantics_viewer(&path)?;
        let profiles = viewer.visible_profiles(NonZeroU32::new(3).expect("non-zero"))?;
        assert_eq!(profiles.len(), 3);

        let forward = profiles.first().expect("forward duplicate profile");
        assert_eq!(forward.read_id(), "duplicate");
        assert!(!forward.is_reverse());
        assert_eq!((forward.align_start(), forward.align_end()), (0, 30));

        assert!(!viewer.handle_individual_key(KeyCode::Char('j'), profiles.len()));
        assert_eq!(viewer.viewport.read_offset, 1);
        let reverse = profiles
            .get(viewer.viewport.read_offset)
            .expect("selected reverse duplicate profile");
        assert_eq!(reverse.read_id(), "duplicate");
        assert!(reverse.is_reverse());
        assert_eq!((reverse.align_start(), reverse.align_end()), (1, 31));
        assert_eq!(
            reverse.calls(),
            [(25, 0), (26, 255), (27, 0), (28, 255), (29, 0), (30, 255)]
        );
        let first_window = reverse.windows().first().expect("first reverse window");
        let second_window = reverse.windows().get(1).expect("second reverse window");
        assert_eq!(first_window.0..first_window.1, 25..28);
        assert!((first_window.2.val() - 1.0 / 3.0).abs() < f32::EPSILON);
        assert_eq!(second_window.0..second_window.1, 28..31);
        assert!((second_window.2.val() - 2.0 / 3.0).abs() < f32::EPSILON);
        assert!(individual_status(&viewer, &profiles, 80).contains("duplicate#2 - mods a win 3"));

        assert!(!viewer.handle_individual_key(KeyCode::Char('k'), profiles.len()));
        assert_eq!(viewer.viewport.read_offset, 0);
        assert!(
            !profiles
                .first()
                .expect("forward duplicate profile")
                .is_reverse()
        );
        remove_viewer_test_bam(&path)
    }

    #[test]
    fn individual_empty_and_no_call_states_match_goldens() -> Result<(), Box<dyn Error>> {
        let path = write_individual_semantics_bam()?;
        let mut viewer = individual_semantics_viewer(&path)?;
        viewer.path = PathBuf::from("individual-states.bam");
        let profiles = viewer.visible_profiles(NonZeroU32::new(3).expect("non-zero"))?;
        let no_calls_index = profiles
            .iter()
            .position(|profile| profile.read_id() == "no-calls")
            .expect("no-call profile");
        let no_calls = profiles.get(no_calls_index).expect("no-call profile");
        assert!(!no_calls.is_reverse());
        assert_eq!((no_calls.align_start(), no_calls.align_end()), (2, 32));
        assert!(no_calls.calls().is_empty());
        assert!(no_calls.windows().is_empty());
        viewer.viewport.read_offset = no_calls_index;
        let no_calls_frame =
            build_individual_frame(&viewer, &profiles, 60, 16, FrameFooter::Controls);
        assert!(no_calls_frame.contains("no a calls in read"));
        assert_ansi_golden(
            "bam_viewer_individual_no_calls.ansi",
            &render_frame_as_ansi(&no_calls_frame, 60, 16)?,
        )?;

        assert!(viewer.go_to(&InitialPosition {
            contig: String::from("first"),
            start: 40,
        })?);
        let empty_profiles = viewer.visible_profiles(NonZeroU32::new(3).expect("non-zero"))?;
        assert!(empty_profiles.is_empty());
        assert_eq!(viewer.viewport.read_offset, 0);
        assert!(individual_status(&viewer, &empty_profiles, 60).contains("read 0/0"));
        let empty_frame =
            build_individual_frame(&viewer, &empty_profiles, 60, 16, FrameFooter::Controls);
        assert!(empty_frame.contains("no reads span this window"));
        assert_ansi_golden(
            "bam_viewer_individual_no_reads.ansi",
            &render_frame_as_ansi(&empty_frame, 60, 16)?,
        )?;
        remove_viewer_test_bam(&path)
    }

    #[test]
    fn individual_goto_resets_selection_and_fetches_target() -> Result<(), Box<dyn Error>> {
        let path = write_individual_semantics_bam()?;
        let mut viewer = individual_semantics_viewer(&path)?;
        viewer.viewport.read_offset = 2;

        assert!(viewer.go_to(&InitialPosition {
            contig: String::from("second"),
            start: 5,
        })?);
        assert_eq!(viewer.target_name(), "second");
        assert_eq!(viewer.viewport.start, 5);
        assert_eq!(viewer.viewport.read_offset, 0);
        let profiles = viewer.visible_profiles(NonZeroU32::new(3).expect("non-zero"))?;
        assert_eq!(profiles.len(), 1);
        let profile = profiles.first().expect("goto target profile");
        assert_eq!(profile.read_id(), "goto-target");
        assert!(!profile.is_reverse());
        assert_eq!((profile.align_start(), profile.align_end()), (4, 24));
        assert_eq!(profile.calls(), [(4, 0), (5, 255), (6, 128)]);
        let window = profile.windows().first().expect("one goto profile window");
        assert_eq!(window.0..window.1, 4..7);
        assert!((window.2.val() - 2.0 / 3.0).abs() < f32::EPSILON);
        assert!(individual_status(&viewer, &profiles, 80).contains("goto-target + mods a win 3"));
        remove_viewer_test_bam(&path)
    }

    #[test]
    fn individual_navigation_viewports_match_ansi_goldens() -> Result<(), Box<dyn Error>> {
        let config: SimulationConfig = serde_json::from_str(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/examples/bam_viewer_individual_demo.json"
        )))?;
        let simulation = TempBamSimulation::new(config, AlignmentFormat::Bam)?;
        let mode = ViewMode::Individual {
            win: NonZeroU32::new(300).expect("non-zero"),
        };
        let mut viewer = Viewer::open_mode(
            PathBuf::from(simulation.bam_path()),
            &InitialPosition {
                contig: String::from("contig_00000"),
                start: 30_000,
            },
            Some(ModChar::new('T')),
            mode,
            window_len_for_columns(80),
        )?;
        viewer.path = PathBuf::from("individual-demo.bam");
        let mut profiles = viewer.visible_profiles(NonZeroU32::new(300).expect("non-zero"))?;
        assert!((15..=30).contains(&profiles.len()));

        let default_frame =
            build_individual_frame(&viewer, &profiles, 80, 24, FrameFooter::Controls);
        assert!(default_frame.contains("contig_00000:30001"));
        assert!(default_frame.contains("mods T win 300"));
        viewer.mod_type = Some(ModChar::from_str("472232").expect("numeric modification code"));
        let numeric_status = individual_status(&viewer, &profiles, 80);
        assert!(numeric_status.len() <= 80);
        assert!(numeric_status.ends_with("mods 472232 win 300"));
        viewer.mod_type = Some(ModChar::new('T'));
        assert_eq!(individual_status(&viewer, &profiles, 6), "");
        assert_eq!(individual_status(&viewer, &profiles, 7), "win 300");
        viewer.mode = ViewMode::Individual {
            win: NonZeroU32::new(u32::MAX).expect("non-zero"),
        };
        assert_eq!(individual_status(&viewer, &profiles, 13), "");
        assert_eq!(individual_status(&viewer, &profiles, 14), "win 4294967295");
        viewer.mode = mode;
        assert_ansi_golden(
            "bam_viewer_individual_default.ansi",
            &render_frame_as_ansi(&default_frame, 80, 24)?,
        )?;

        assert!(!viewer.handle_individual_key(KeyCode::Char('j'), profiles.len()));
        let after_j_frame =
            build_individual_frame(&viewer, &profiles, 80, 24, FrameFooter::Controls);
        assert_ansi_golden(
            "bam_viewer_individual_after_j.ansi",
            &render_frame_as_ansi(&after_j_frame, 80, 24)?,
        )?;

        assert!(viewer.handle_individual_key(KeyCode::Char('l'), profiles.len()));
        profiles = viewer.visible_profiles(NonZeroU32::new(300).expect("non-zero"))?;
        let horizontal_frame =
            build_individual_frame(&viewer, &profiles, 80, 24, FrameFooter::Controls);
        assert_ansi_golden(
            "bam_viewer_individual_after_l.ansi",
            &render_frame_as_ansi(&horizontal_frame, 80, 24)?,
        )?;

        let narrow = build_individual_frame(&viewer, &profiles, 14, 24, FrameFooter::Controls);
        assert_ansi_golden(
            "bam_viewer_individual_narrow.ansi",
            &render_frame_as_ansi(&narrow, 14, 24)?,
        )?;
        let short = build_individual_frame(&viewer, &profiles, 80, 6, FrameFooter::Controls);
        assert!(!short.contains("1.0 ┤"));
        assert!(!short.contains('└'));
        assert_ansi_golden(
            "bam_viewer_individual_short.ansi",
            &render_frame_as_ansi(&short, 80, 6)?,
        )?;

        assert!(viewer.go_to(&InitialPosition {
            contig: String::from("contig_00001"),
            start: 30_000,
        })?);
        let second_contig_profiles =
            viewer.visible_profiles(NonZeroU32::new(300).expect("non-zero"))?;
        let second_contig_frame = build_individual_frame(
            &viewer,
            &second_contig_profiles,
            100,
            24,
            FrameFooter::Controls,
        );
        assert!(second_contig_frame.contains("contig_00001:30001"));
        assert!(!second_contig_frame.contains("contig_000:30001"));
        Ok(())
    }

    #[test]
    fn read_id_width_toggles_to_the_longest_cached_id() {
        let position = InitialPosition {
            contig: String::from("dummyIII"),
            start: 23,
        };
        let mut viewer = Viewer::open(PathBuf::from("examples/example_1.bam"), &position, None, 7)
            .expect("position should open");
        let records = viewer.visible_records().expect("records should load");
        let longest_id_width = records
            .iter()
            .map(|record| record.read_id().len())
            .max()
            .expect("one record")
            .saturating_add(1);

        viewer.toggle_read_id_width(&records);
        assert!(viewer.full_read_ids);
        assert_eq!(usize::from(viewer.read_label_width), longest_id_width);
        let first_read_id = records.first().expect("one record").read_id();
        assert!(
            build_frame(&viewer, &records, 80, 10, FrameFooter::Controls).contains(first_read_id)
        );

        viewer.viewport.navigate(KeyCode::Down, 76, 7, 10, 1);
        assert_eq!(usize::from(viewer.read_label_width), longest_id_width);

        viewer.toggle_read_id_width(&records);
        assert!(!viewer.full_read_ids);
        assert_eq!(viewer.read_label_width, READ_LABEL_WIDTH);

        viewer.toggle_read_id_width(&records);
        viewer.reset_read_id_width();
        assert!(!viewer.full_read_ids);
        assert_eq!(viewer.read_label_width, READ_LABEL_WIDTH);
    }

    #[test]
    fn display_keys_do_not_fetch_and_horizontal_keys_reset_options() {
        let position = InitialPosition {
            contig: String::from("dummyIII"),
            start: 23,
        };
        let mut viewer = Viewer::open(PathBuf::from("examples/example_1.bam"), &position, None, 7)
            .expect("position should open");
        let records = viewer.visible_records().expect("records should load");

        assert!(!viewer.handle_key(KeyCode::Char('i'), &records, 1));
        assert!(viewer.show_insertions);
        for key in [
            KeyCode::Up,
            KeyCode::Down,
            KeyCode::Char('k'),
            KeyCode::Char('j'),
            KeyCode::PageUp,
            KeyCode::PageDown,
            KeyCode::Home,
            KeyCode::End,
        ] {
            assert!(!viewer.handle_key(key, &records, 1));
            assert!(viewer.show_insertions);
        }
        assert!(!viewer.handle_key(KeyCode::Char('i'), &records, 1));
        assert!(!viewer.show_insertions);

        assert!(!viewer.handle_key(KeyCode::Char('r'), &records, 1));
        assert!(viewer.full_read_ids);
        assert!(!viewer.handle_key(KeyCode::Char('r'), &records, 1));
        assert!(!viewer.full_read_ids);

        for key in [
            KeyCode::Left,
            KeyCode::Right,
            KeyCode::Char('h'),
            KeyCode::Char('l'),
        ] {
            viewer.show_insertions = true;
            viewer.full_read_ids = true;
            viewer.read_label_width = full_read_label_width(&records);
            assert!(viewer.handle_key(key, &records, 1));
            assert!(!viewer.show_insertions);
            assert!(!viewer.full_read_ids);
            assert_eq!(viewer.read_label_width, READ_LABEL_WIDTH);
        }

        viewer.viewport.start = 0;
        viewer.show_insertions = true;
        viewer.full_read_ids = true;
        assert!(!viewer.handle_key(KeyCode::Left, &records, 1));
        assert!(!viewer.show_insertions);
        assert!(!viewer.full_read_ids);
    }

    #[test]
    fn goto_validates_position_and_resets_horizontal_options() {
        let initial_position = InitialPosition {
            contig: String::from("dummyIII"),
            start: 23,
        };
        let mut viewer = Viewer::open(
            PathBuf::from("examples/example_1.bam"),
            &initial_position,
            None,
            7,
        )
        .expect("position should open");
        viewer.viewport.read_offset = 1;
        viewer.full_read_ids = true;
        viewer.read_label_width = 40;
        viewer.show_insertions = true;

        let changed = viewer
            .go_to(&InitialPosition {
                contig: String::from("dummyI"),
                start: 10,
            })
            .expect("valid position");
        assert!(changed);
        assert_eq!(viewer.target_name(), "dummyI");
        assert_eq!(viewer.viewport.start, 10);
        assert_eq!(viewer.viewport.read_offset, 0);
        assert!(!viewer.full_read_ids);
        assert_eq!(viewer.read_label_width, READ_LABEL_WIDTH);
        assert!(!viewer.show_insertions);

        assert!(
            !viewer
                .go_to(&InitialPosition {
                    contig: String::from("dummyI"),
                    start: 10,
                })
                .expect("unchanged position remains valid")
        );
        let unchanged_viewport = viewer.viewport;
        let _unknown_reference_error = viewer
            .go_to(&InitialPosition {
                contig: String::from("missing"),
                start: 0,
            })
            .expect_err("unknown reference should fail");
        assert_eq!(viewer.viewport, unchanged_viewport);
        let _out_of_range_error = viewer
            .go_to(&InitialPosition {
                contig: String::from("dummyI"),
                start: 22,
            })
            .expect_err("out-of-range position should fail");
        assert_eq!(viewer.viewport, unchanged_viewport);
    }

    #[test]
    fn invalid_goto_can_be_corrected_and_submitted() {
        let mut viewer = goto_test_viewer();
        let mut input = String::from("missing:0");
        let mut input_error = None;

        assert_eq!(
            prompt_key(&mut input, &mut input_error, KeyCode::Enter, &mut viewer),
            None
        );
        assert!(
            input_error
                .as_deref()
                .is_some_and(|error| error.contains("unknown reference"))
        );
        assert_eq!(viewer.target_name(), "dummyIII");
        assert_eq!(viewer.viewport.read_offset, 1);
        assert!(viewer.full_read_ids);
        assert!(viewer.show_insertions);

        for _ in 0..input.len() {
            assert_eq!(
                prompt_key(
                    &mut input,
                    &mut input_error,
                    KeyCode::Backspace,
                    &mut viewer,
                ),
                None
            );
        }
        for character in "dummyI:10".chars() {
            assert_eq!(
                prompt_key(
                    &mut input,
                    &mut input_error,
                    KeyCode::Char(character),
                    &mut viewer,
                ),
                None
            );
        }
        assert_eq!(
            prompt_key(&mut input, &mut input_error, KeyCode::Enter, &mut viewer),
            Some(PositionPromptOutcome::Navigated(true))
        );
        assert_eq!(viewer.target_name(), "dummyI");
        assert_eq!(viewer.viewport.start, 10);
        assert_eq!(viewer.viewport.read_offset, 0);
        assert!(!viewer.full_read_ids);
        assert!(!viewer.show_insertions);
    }

    #[test]
    fn corrected_goto_can_remain_invalid_without_changing_viewer() {
        let mut viewer = goto_test_viewer();
        let original_viewport = viewer.viewport;
        let mut input = String::from("missing:0");
        let mut input_error = None;

        assert_eq!(
            prompt_key(&mut input, &mut input_error, KeyCode::Enter, &mut viewer),
            None
        );
        assert_eq!(
            prompt_key(
                &mut input,
                &mut input_error,
                KeyCode::Backspace,
                &mut viewer,
            ),
            None
        );
        for character in "nonsense".chars() {
            let _outcome = prompt_key(
                &mut input,
                &mut input_error,
                KeyCode::Char(character),
                &mut viewer,
            );
        }
        assert_eq!(
            prompt_key(&mut input, &mut input_error, KeyCode::Enter, &mut viewer),
            None
        );
        assert!(
            input_error
                .as_deref()
                .is_some_and(|error| error.contains("non-negative integer"))
        );
        assert_eq!(viewer.target_name(), "dummyIII");
        assert_eq!(viewer.viewport, original_viewport);
        assert!(viewer.full_read_ids);
        assert!(viewer.show_insertions);
    }

    #[test]
    fn escape_cancels_invalid_goto_and_cached_navigation_still_works() {
        let mut viewer = goto_test_viewer();
        let records = viewer.visible_records().expect("records should load");
        let original_viewport = viewer.viewport;
        let mut input = String::from("missing:0");
        let mut input_error = None;

        assert_eq!(
            prompt_key(&mut input, &mut input_error, KeyCode::Enter, &mut viewer),
            None
        );
        assert_eq!(
            prompt_key(&mut input, &mut input_error, KeyCode::Esc, &mut viewer),
            Some(PositionPromptOutcome::Cancelled)
        );
        assert_eq!(viewer.viewport, original_viewport);
        assert!(viewer.full_read_ids);
        assert!(viewer.show_insertions);

        assert!(!viewer.handle_key(KeyCode::Home, &records, 1));
        assert_eq!(viewer.viewport.read_offset, 0);
        assert!(viewer.full_read_ids);
        assert!(viewer.show_insertions);
    }

    #[test]
    fn ruler_and_read_labels_share_the_same_dynamic_column() {
        assert_eq!(label_column("read id", 5), "read ");
        assert_eq!(label_column("read", 12), "read        ");
    }

    #[test]
    fn plain_argument_parser_accepts_bam_and_position() {
        let args = Args::parse_from([OsString::from("reads.bam"), OsString::from("chr1:10")])
            .expect("valid arguments")
            .expect("not a help request");
        assert_eq!(args.bam, PathBuf::from("reads.bam"));
        assert_eq!(
            args.position,
            InitialPosition {
                contig: String::from("chr1"),
                start: 10,
            }
        );
        assert_eq!(args.mod_type, None);
    }

    #[test]
    fn usage_describes_display_and_controls() {
        assert!(USAGE.contains("bold and underlined"));
        assert!(
            USAGE.contains("Lowercase bases are insertions; dots are deletions or reference skips")
        );
        assert!(USAGE.contains("asterisk means the BAM alignment has no stored read sequence"));
        assert!(USAGE.contains("Horizontal movement truncates read IDs and hides insertions"));
        assert!(USAGE.contains("Page Up/Page Down"));
        assert!(USAGE.contains("Home/End jump to the first/last read"));
        assert!(USAGE.contains("g prompts for CONTIG:START"));
        assert!(USAGE.contains("A successful goto also truncates read IDs and hides insertions"));
        assert!(USAGE.contains("Backspace edits; Enter submits; Escape cancels"));
        assert!(USAGE.contains("r toggles full read IDs and i toggles insertions"));
        assert!(USAGE.contains("grey raw ML calls"));
        assert!(USAGE.contains("WINDOW_SIZE is a positive number of modified bases"));
        assert!(USAGE.contains("In individual view, j/k selects one read; r/i have no effect"));
        assert!(USAGE.contains("Ctrl-C or Ctrl-D always quits"));
    }

    #[test]
    fn goto_prompt_keeps_the_editable_suffix_visible() {
        let mut input = format!("{}:123", "long-contig".repeat(8));
        let footer = position_prompt_footer(&input, 80);
        assert_eq!(footer.len(), 80);
        assert!(footer.ends_with("long-contig:123"));

        let _removed_character = input.pop();
        let edited_footer = position_prompt_footer(&input, 80);
        assert_eq!(edited_footer.len(), 80);
        assert!(edited_footer.ends_with("long-contig:12"));
    }

    #[test]
    fn goto_error_explains_how_to_correct_or_cancel() {
        assert_eq!(
            position_error_footer("unknown reference 'missing'", 80),
            "Error: unknown reference 'missing'. Backspace to correct; Esc to cancel."
        );
        let long_error = position_error_footer(&"x".repeat(100), 80);
        assert_eq!(long_error.len(), 80);
        assert!(long_error.contains("..."));
        assert!(long_error.ends_with(". Backspace to correct; Esc to cancel."));
        assert_eq!(
            position_error_footer("unknown reference 'missing'", 45),
            "Backspace edits; Esc cancels."
        );
        assert_eq!(
            position_error_footer("unknown reference 'missing'", 27),
            "Bksp edit; Esc back"
        );
        assert_eq!(
            position_error_footer("unknown reference 'missing'", 18),
            "Bksp/Esc"
        );
        assert_eq!(
            fixed_line(&position_error_footer("unknown reference 'missing'", 7), 7),
            "Esc    "
        );
    }

    #[test]
    fn plain_argument_parser_accepts_a_modification_type() {
        let args = Args::parse_from([
            OsString::from("reads.bam"),
            OsString::from("chr1:10"),
            OsString::from("472232"),
        ])
        .expect("valid arguments")
        .expect("not a help request");
        assert_eq!(
            args.mod_type.expect("modification type").to_string(),
            "472232"
        );
        assert_eq!(args.mode, ViewMode::Table);
    }

    #[test]
    fn argument_parser_accepts_individual_mode() {
        let args = Args::parse_from([
            OsString::from("reads.bam"),
            OsString::from("chr1:10"),
            OsString::from("m"),
            OsString::from("300"),
            OsString::from("individual"),
        ])
        .expect("valid individual arguments")
        .expect("not help");
        assert_eq!(args.mod_type, Some(ModChar::new('m')));
        assert_eq!(
            args.mode,
            ViewMode::Individual {
                win: NonZeroU32::new(300).expect("non-zero")
            }
        );
    }

    #[test]
    fn argument_parser_rejects_incomplete_or_invalid_individual_mode() {
        for arguments in [
            vec!["reads.bam", "chr1:10", "m", "individual"],
            vec!["reads.bam", "chr1:10", "m", "300"],
            vec!["reads.bam", "chr1:10", "300", "individual"],
            vec!["reads.bam", "chr1:10", "m", "0", "individual"],
            vec!["reads.bam", "chr1:10", "m", "abc", "individual"],
            vec!["reads.bam", "chr1:10", "m", "300", "individual", "extra"],
        ] {
            let _error = Args::parse_from(arguments.into_iter().map(OsString::from))
                .expect_err("invalid individual argument combination must fail");
        }
    }

    #[test]
    fn argument_parser_requires_a_position() {
        let error = Args::parse_from([OsString::from("reads.bam")])
            .expect_err("position should be compulsory");
        assert!(error.contains("missing CONTIG:START position"));
    }

    #[test]
    fn argument_parser_rejects_an_end_coordinate() {
        let error = Args::parse_from([OsString::from("reads.bam"), OsString::from("chr1:10-20")])
            .expect_err("an end coordinate should not be accepted");
        assert!(error.contains("START must be a non-negative integer"));
    }

    #[test]
    fn position_parser_allows_colons_in_contig_names() {
        assert_eq!(
            InitialPosition::parse("chr1:alternate:10").expect("valid position"),
            InitialPosition {
                contig: String::from("chr1:alternate"),
                start: 10,
            }
        );
    }

    #[test]
    fn terminal_width_defines_the_genomic_window() {
        assert_eq!(window_len_for_columns(10), 1);
        assert_eq!(window_len_for_columns(100), 81);
        assert_eq!(window_len_for_columns(500), 200);
    }

    #[test]
    fn viewer_starts_at_the_requested_position() {
        let position = InitialPosition {
            contig: String::from("dummyIII"),
            start: 10,
        };
        let mut viewer = Viewer::open(PathBuf::from("examples/example_1.bam"), &position, None, 40)
            .expect("position should open");
        assert_eq!(viewer.viewport.start, 10);
        assert_eq!(viewer.current_window_len(), 40);
        let frame = build_frame(&viewer, &[], 80, 10, FrameFooter::Controls);
        assert!(frame.contains("h/l 40 bp"));
        assert!(frame.contains("pgup/dn"));
        assert!(frame.contains("home/end"));
        assert!(frame.contains("g goto"));
        assert!(frame.contains("r full IDs"));
        assert!(frame.contains("i show ins"));
        assert!(frame.contains("q quit"));

        viewer.window_len = 200;
        viewer.full_read_ids = true;
        viewer.show_insertions = true;
        let toggled_frame = build_frame(&viewer, &[], 80, 10, FrameFooter::Controls);
        assert!(toggled_frame.contains("r short IDs"));
        assert!(toggled_frame.contains("i hide ins"));
        assert!(toggled_frame.contains("q quit"));

        let prompt_frame = build_frame(
            &viewer,
            &[],
            80,
            10,
            FrameFooter::PositionPrompt {
                input: "dummyI:10",
                error: None,
            },
        );
        assert!(prompt_frame.contains("Go to CONTIG:START: dummyI:10"));
        assert!(!prompt_frame.contains("q quit"));
    }

    #[test]
    fn window_is_shortened_at_the_end_of_a_contig() {
        let position = InitialPosition {
            contig: String::from("dummyI"),
            start: 10,
        };
        let viewer = Viewer::open(
            PathBuf::from("examples/example_1.bam"),
            &position,
            None,
            200,
        )
        .expect("position should open");
        assert_eq!(viewer.viewport.start, 10);
        assert_eq!(viewer.current_window_len(), 12);
    }

    #[test]
    fn last_contig_base_is_a_valid_one_base_window() {
        let position = InitialPosition {
            contig: String::from("dummyI"),
            start: 21,
        };
        let viewer = Viewer::open(
            PathBuf::from("examples/example_1.bam"),
            &position,
            None,
            200,
        )
        .expect("last base should open");
        assert_eq!(viewer.viewport.start, 21);
        assert_eq!(viewer.current_window_len(), 1);
    }

    #[test]
    fn position_at_contig_end_is_rejected() {
        let position = InitialPosition {
            contig: String::from("dummyI"),
            start: 22,
        };
        let error = Viewer::open(
            PathBuf::from("examples/example_1.bam"),
            &position,
            None,
            200,
        )
        .expect_err("position at contig end should fail");
        assert!(error.to_string().contains("initial position 22"));
    }

    #[test]
    fn right_navigation_never_moves_a_shortened_window_backward() {
        let mut viewport = Viewport {
            tid: 0,
            start: 90,
            read_offset: 4,
        };
        viewport.navigate(KeyCode::Right, 100, 20, 7, 3);
        assert_eq!(viewport.start, 90);
        assert_eq!(viewport.read_offset, 0);
    }

    #[test]
    fn standard_exit_keys_are_recognized() {
        assert!(should_quit(KeyEvent::new(
            KeyCode::Char('q'),
            KeyModifiers::NONE
        )));
        assert!(should_quit(KeyEvent::new(KeyCode::Esc, KeyModifiers::NONE)));
        assert!(should_quit(KeyEvent::new(
            KeyCode::Char('c'),
            KeyModifiers::CONTROL
        )));
        assert!(should_quit(KeyEvent::new(
            KeyCode::Char('d'),
            KeyModifiers::CONTROL
        )));
        assert!(!should_quit(KeyEvent::new(
            KeyCode::Char('c'),
            KeyModifiers::NONE
        )));
    }

    #[test]
    fn ghostty_parses_a_complete_frame() -> Result<(), Box<dyn Error>> {
        let mut terminal = Terminal::new(TerminalOptions {
            cols: 20,
            rows: 3,
            max_scrollback: 0,
        })?;
        terminal.vt_write(b"\x1b[H\x1b[1;32mBAM\x1b[0m\x1b[3;1Hq quit");
        let mut state = RenderState::new()?;
        let snapshot = state.update(&terminal)?;
        assert_eq!(snapshot.cols()?, 20);
        assert_eq!(snapshot.rows()?, 3);
        Ok(())
    }
}
