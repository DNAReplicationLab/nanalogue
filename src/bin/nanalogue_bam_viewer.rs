//! Interactive terminal viewer for indexed BAM files.
//!
//! `libghostty-vt` owns the virtual screen and interprets each ANSI frame.
//! Crossterm only handles the host terminal's raw mode, events, and drawing.

#![expect(
    clippy::print_stderr,
    reason = "command-line errors are intentionally reported to stderr"
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
use nanalogue_core::reads_table::{RegionSequence, RegionSequenceReader};
use std::{
    env,
    error::Error,
    ffi::OsString,
    fmt::Write as _,
    io::{self, Stdout, Write as _},
    path::PathBuf,
    sync::Arc,
};

/// Width reserved for read names and the separating space.
const READ_LABEL_WIDTH: u16 = 19;

/// Maximum number of genomic bases displayed regardless of terminal width.
const MAX_REGION_LENGTH: u32 = 200;

/// Short usage text for this deliberately minimal two-argument program.
const USAGE: &str = "Usage: nanalogue_bam_viewer <BAM> [CONTIG:START]\n\
START is zero-based; displayed coordinates are one-based.\n\
The end coordinate is selected from the terminal width, up to 200 bp.";

/// Initial reference and zero-based coordinate supplied on the command line.
#[derive(Debug, Clone, PartialEq, Eq)]
struct InitialPosition {
    /// Reference name.
    contig: String,
    /// Zero-based first displayed coordinate.
    start: u32,
}

impl InitialPosition {
    /// Parses `CONTIG:START`, splitting at the final colon for colon-containing contig names.
    fn parse(value: &str) -> Result<Self, String> {
        let (contig, start_text) = value
            .rsplit_once(':')
            .ok_or_else(|| String::from("position must have the form CONTIG:START"))?;
        if contig.is_empty() || start_text.is_empty() {
            return Err(String::from("position must have the form CONTIG:START"));
        }
        let start = start_text
            .parse::<u32>()
            .map_err(|_error| String::from("START must be a non-negative integer"))?;
        Ok(Self {
            contig: String::from(contig),
            start,
        })
    }
}

/// Returns the number of sequence columns available at a terminal width.
fn window_len_for_columns(cols: u16) -> u32 {
    u32::from(cols.saturating_sub(READ_LABEL_WIDTH).max(1)).min(MAX_REGION_LENGTH)
}

/// Command-line options for the BAM viewer.
#[derive(Debug)]
struct Args {
    /// Indexed BAM file to view.
    bam: PathBuf,

    /// Initial reference position, such as `chr1:1000`.
    position: Option<InitialPosition>,
}

impl Args {
    /// Parses the two positional arguments, returning `None` for a help request.
    fn parse_from<I>(arguments: I) -> Result<Option<Self>, String>
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
        let position = argument_iter
            .next()
            .map(|position_argument| {
                let position_string = position_argument
                    .into_string()
                    .map_err(|_position| String::from("position must be valid UTF-8"))?;
                InitialPosition::parse(&position_string)
            })
            .transpose()?;
        if argument_iter.next().is_some() {
            return Err(String::from(
                "expected only BAM and optional CONTIG:START arguments",
            ));
        }
        Ok(Some(Self {
            bam: PathBuf::from(bam),
            position,
        }))
    }
}

/// Position and scroll state of the viewer.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Viewport {
    /// Numeric BAM target identifier.
    tid: u32,
    /// Zero-based first visible reference coordinate.
    start: u32,
    /// First visible read row.
    read_offset: usize,
}

impl Viewport {
    /// Applies one navigation key while respecting reference and read bounds.
    #[expect(
        clippy::wildcard_enum_match_arm,
        reason = "all other terminal keys intentionally leave the viewport unchanged"
    )]
    fn navigate(
        &mut self,
        key: KeyCode,
        contig_len: u32,
        window_len: u32,
        read_count: usize,
        visible_reads: usize,
    ) {
        let max_start = contig_len.saturating_sub(window_len);
        let max_offset = read_count.saturating_sub(visible_reads);
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
            _ => {}
        }
    }
}

/// BAM data and mutable viewport state.
#[derive(Debug)]
struct Viewer {
    /// Source BAM path.
    path: PathBuf,
    /// Nanalogue reader reused for region sequence tables.
    reader: RegionSequenceReader,
    /// Current position and read scroll.
    viewport: Viewport,
    /// Number of genomic bases in each fetched window.
    window_len: u32,
    /// Width of the read-ID column, including its trailing separator.
    read_label_width: u16,
    /// Whether complete read IDs are displayed.
    full_read_ids: bool,
}

impl Viewer {
    /// Opens an indexed BAM and chooses the initial viewport.
    fn open(
        path: PathBuf,
        position: Option<InitialPosition>,
        window_len: u32,
    ) -> Result<Self, Box<dyn Error>> {
        let reader = RegionSequenceReader::from_path(&path)?;

        let viewport = if let Some(initial) = position {
            let tid = reader.target_id(&initial.contig).ok_or_else(|| {
                format!("reference '{}' is not in the BAM header", initial.contig)
            })?;
            let contig_len = reader
                .target_len(tid)
                .ok_or("BAM target identifier is invalid")?;
            if initial.start >= contig_len {
                return Err(format!(
                    "initial position {} is outside reference '{}' (length {contig_len})",
                    initial.start, initial.contig
                )
                .into());
            }
            Viewport {
                tid,
                start: initial.start,
                read_offset: 0,
            }
        } else {
            let viewport = Viewport::default();
            if reader.target_len(viewport.tid).unwrap_or(0) == 0 {
                return Err("first BAM reference has zero length".into());
            }
            viewport
        };

        Ok(Self {
            path,
            reader,
            viewport,
            window_len: window_len.clamp(1, MAX_REGION_LENGTH),
            read_label_width: READ_LABEL_WIDTH,
            full_read_ids: false,
        })
    }

    /// Returns the current target name.
    fn target_name(&self) -> &str {
        self.reader.target_name(self.viewport.tid).unwrap_or("?")
    }

    /// Returns the current target length.
    fn target_len(&self) -> u32 {
        self.reader.target_len(self.viewport.tid).unwrap_or(0)
    }

    /// Returns the current window length, shortened only at the end of a reference.
    fn current_window_len(&self) -> u32 {
        self.window_len
            .min(self.target_len().saturating_sub(self.viewport.start))
    }

    /// Toggles between the default and longest cached read-ID widths.
    fn toggle_read_id_width(&mut self, records: &[RegionSequence]) {
        self.full_read_ids = !self.full_read_ids;
        self.read_label_width = if self.full_read_ids {
            full_read_label_width(records)
        } else {
            READ_LABEL_WIDTH
        };
    }

    /// Restores the default read-ID width.
    fn reset_read_id_width(&mut self) {
        self.full_read_ids = false;
        self.read_label_width = READ_LABEL_WIDTH;
    }

    /// Fetches records spanning the visible genomic range.
    fn visible_records(&mut self) -> Result<Vec<RegionSequence>, Box<dyn Error>> {
        let end = self
            .viewport
            .start
            .saturating_add(self.current_window_len())
            .min(self.target_len());
        Ok(self
            .reader
            .sequences(self.viewport.tid, self.viewport.start, end)?)
    }
}

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

/// Pads a nanalogue region sequence to the visible genomic screen width.
fn sequence_columns(record: &RegionSequence, width: u16) -> String {
    let mut output = vec![b' '; usize::from(width)];
    for (offset, base) in record.sequence().bytes().enumerate() {
        if let Some(cell) = output.get_mut(offset) {
            *cell = base;
        }
    }
    String::from_utf8(output).expect("nanalogue region sequences contain ASCII")
}

/// Truncates and pads a line to the terminal width.
fn fixed_line(text: &str, width: u16) -> String {
    printable_label(text.as_bytes(), usize::from(width))
}

/// Truncates or pads a label and appends its column separator.
fn label_column(label: &str, width: u16) -> String {
    let mut column = printable_label(label.as_bytes(), usize::from(width.saturating_sub(1)));
    column.push(' ');
    column
}

/// Returns enough columns for the longest cached read ID and a separator.
fn full_read_label_width(records: &[RegionSequence]) -> u16 {
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

/// Builds the ANSI frame that Ghostty parses into a terminal screen.
fn build_frame(viewer: &Viewer, records: &[RegionSequence], cols: u16, rows: u16) -> String {
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
        let status = format!(
            " {path}  {}:{}-{}  reads {}  row {}",
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
        let mut ruler = label_column("reference", viewer.read_label_width);
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
        let mut line = label_column(record.read_id(), viewer.read_label_width);
        line.push_str(&sequence_columns(record, genome_cols));
        let terminal_row = screen_row.saturating_add(4);
        write!(&mut frame, "\x1b[{terminal_row};1H").expect("writing to String cannot fail");
        frame.push_str(if record.is_reverse() {
            "\x1b[33m"
        } else {
            "\x1b[32m"
        });
        frame.push_str(&fixed_line(&line, effective_cols));
        frame.push_str("\x1b[0m");
    }

    if rows > 3 {
        write!(&mut frame, "\x1b[{rows};1H\x1b[7m").expect("writing to String cannot fail");
        let footer = format!(
            " left/h right/l move {} bp  up/k down/j scroll reads  i {} IDs  q quit",
            viewer.window_len,
            if viewer.full_read_ids {
                "short"
            } else {
                "full"
            }
        );
        frame.push_str(&fixed_line(&footer, effective_cols));
        frame.push_str("\x1b[0m");
    }
    frame
}

/// Returns whether a key is one of the viewer's standard exit sequences.
fn should_quit(key: KeyEvent) -> bool {
    matches!(key.code, KeyCode::Char('q') | KeyCode::Esc)
        || (key.modifiers.contains(KeyModifiers::CONTROL)
            && matches!(key.code, KeyCode::Char('c' | 'd')))
}

/// Runs the interactive event loop.
fn run(args: Args) -> Result<(), Box<dyn Error>> {
    let (initial_cols, initial_rows) = crossterm::terminal::size()?;
    let mut viewer = Viewer::open(
        args.bam,
        args.position,
        window_len_for_columns(initial_cols),
    )?;
    let mut renderer = GhosttyRenderer::new(initial_cols, initial_rows)?;
    let mut stdout = io::stdout();
    let _guard = TerminalGuard::enter(&mut stdout)?;
    let mut records = viewer.visible_records()?;

    loop {
        let (cols, rows) = crossterm::terminal::size()?;
        let resized_window_len = window_len_for_columns(cols);
        if resized_window_len != viewer.window_len {
            viewer.window_len = resized_window_len;
            records = viewer.visible_records()?;
            if viewer.full_read_ids {
                viewer.read_label_width = full_read_label_width(&records);
            }
        }
        let visible_reads = usize::from(rows.saturating_sub(4));
        viewer.viewport.read_offset = viewer
            .viewport
            .read_offset
            .min(records.len().saturating_sub(visible_reads));
        let frame = build_frame(&viewer, &records, cols, rows);
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
        if key.code == KeyCode::Char('i') {
            viewer.toggle_read_id_width(&records);
            continue;
        }
        if matches!(
            key.code,
            KeyCode::Left | KeyCode::Right | KeyCode::Char('h' | 'l')
        ) {
            viewer.reset_read_id_width();
        }
        let previous_start = viewer.viewport.start;
        viewer.viewport.navigate(
            key.code,
            viewer.target_len(),
            viewer.window_len,
            records.len(),
            visible_reads,
        );
        if viewer.viewport.start != previous_start {
            records = viewer.visible_records()?;
        }
    }
    Ok(())
}

/// Parses arguments and reports errors after terminal cleanup.
fn main() {
    // SAFETY: this is called before HTSlib work or additional threads begin.
    unsafe {
        nanalogue_core::init_ssl_certificates();
    }
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
    }

    #[test]
    fn viewer_uses_nanalogue_region_sequences() {
        let mut reader = RegionSequenceReader::from_path("examples/example_1.bam")
            .expect("open indexed example");
        assert!(
            reader
                .sequences(2, 20, 30)
                .expect("retrieve a partially covered region")
                .is_empty(),
            "a read that does not span the full region must be excluded"
        );
        let rows = reader
            .sequences(2, 23, 30)
            .expect("retrieve nanalogue region sequences");
        let row = rows.first().expect("one overlapping read");
        assert_eq!(row.read_id(), "a4f36092-b4d5-47a9-813e-c22c3b477a0c");
        assert_eq!(sequence_columns(row, 10), "ACATCAA   ");
    }

    #[test]
    fn read_id_width_toggles_to_the_longest_cached_id() {
        let position = InitialPosition {
            contig: String::from("dummyIII"),
            start: 23,
        };
        let mut viewer = Viewer::open(PathBuf::from("examples/example_1.bam"), Some(position), 7)
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
        assert!(build_frame(&viewer, &records, 80, 10).contains(first_read_id));

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
    fn ruler_and_read_labels_share_the_same_dynamic_column() {
        assert_eq!(label_column("reference", 5), "refe ");
        assert_eq!(label_column("read", 12), "read        ");
    }

    #[test]
    fn plain_argument_parser_accepts_bam_and_position() {
        let args = Args::parse_from([OsString::from("reads.bam"), OsString::from("chr1:10")])
            .expect("valid arguments")
            .expect("not a help request");
        assert_eq!(args.bam, PathBuf::from("reads.bam"));
        assert_eq!(
            args.position.expect("position"),
            InitialPosition {
                contig: String::from("chr1"),
                start: 10,
            }
        );
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
        let viewer = Viewer::open(PathBuf::from("examples/example_1.bam"), Some(position), 40)
            .expect("position should open");
        assert_eq!(viewer.viewport.start, 10);
        assert_eq!(viewer.current_window_len(), 40);
        let frame = build_frame(&viewer, &[], 80, 10);
        assert!(frame.contains("move 40 bp"));
        assert!(frame.contains("i full IDs"));
    }

    #[test]
    fn window_is_shortened_at_the_end_of_a_contig() {
        let position = InitialPosition {
            contig: String::from("dummyI"),
            start: 10,
        };
        let viewer = Viewer::open(PathBuf::from("examples/example_1.bam"), Some(position), 200)
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
        let viewer = Viewer::open(PathBuf::from("examples/example_1.bam"), Some(position), 200)
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
        let error = Viewer::open(PathBuf::from("examples/example_1.bam"), Some(position), 200)
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
    fn omitted_position_starts_on_the_first_contig() {
        let viewer = Viewer::open(PathBuf::from("examples/example_1.bam"), None, 200)
            .expect("first contig should open");
        assert_eq!(viewer.viewport.start, 0);
        assert_eq!(viewer.current_window_len(), 22);
    }

    #[test]
    fn standard_exit_keys_are_recognized() {
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
