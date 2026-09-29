//! Verifies the complete on-screen lifecycle of simulated reads while scrolling right.
//!
//! A read can have four observable relationships to a fixed-width table viewport:
//!
//! ```text
//! ABSENT                       "          "
//! RIGHT_BOUNDARY_HAS_SEQUENCE  "     ACGTA"  (the read is entering)
//! FULL                         "ACGTACGTAC"
//! LEFT_BOUNDARY_HAS_SEQUENCE   "ACGTA     "  (the read is leaving)
//! ```
//!
//! `ABSENT` has two distinct meanings in a lifecycle: the read has not appeared yet, or
//! the read has already disappeared. Keeping those phases separate prevents a pairwise
//! transition check from accidentally allowing a read to reappear. For monotonically
//! increasing, contiguous 10 bp viewports and this fixture's 50 bp alignments, the complete
//! state machine is:
//!
//! ```text
//! +---------------+   +-----------------------------+   +------+   +----------------------------+   +--------------+
//! | ABSENT_BEFORE |-->| RIGHT_BOUNDARY_HAS_SEQUENCE |-->| FULL |-->| LEFT_BOUNDARY_HAS_SEQUENCE |-->| ABSENT_AFTER |
//! +-------+-------+   +-----------------------------+   +--+---+   +----------------------------+   +-------+------+
//!         |                                                   |                                               |
//!         +-- remains absent                                 +-- additional full windows                     +-- remains absent
//!         |                                                                                                  ^
//!         +---------------- boundary-aligned start ----------> FULL -------- boundary-aligned end ------------+
//! ```
//!
//! The boundary-aligned shortcuts skip the partial frame: `ABSENT_BEFORE -> FULL`
//! when a read starts at a viewport boundary, and `FULL -> ABSENT_AFTER` when it ends
//! at one. The complete transition table is:
//!
//! | Current state | Allowed next states |
//! | --- | --- |
//! | `ABSENT_BEFORE` | `ABSENT_BEFORE`, `RIGHT_BOUNDARY_HAS_SEQUENCE`, `FULL` |
//! | `RIGHT_BOUNDARY_HAS_SEQUENCE` | `FULL` |
//! | `FULL` | `FULL`, `LEFT_BOUNDARY_HAS_SEQUENCE`, `ABSENT_AFTER` |
//! | `LEFT_BOUNDARY_HAS_SEQUENCE` | `ABSENT_AFTER` |
//! | `ABSENT_AFTER` | `ABSENT_AFTER` |
//!
//! This test generates randomly positioned 50 bp alignments with a centered 20 bp deletion
//! along a fixed-length contig. The alignments are longer than one viewport, and deletion
//! dots count as alignment content rather than whitespace. The test starts at the first
//! reference base and uses the same `l` key handling and record refetch path as the
//! interactive viewer until it reaches the final viewport. Every read is classified from
//! its rendered boundary cells in every frame, then checked against the table above. It
//! also checks that every read contributes exactly 20 deletion columns across the frames.
//! Finally, the test checks that the simulated population contains the expected number of
//! complete non-boundary lifecycles, including leading whitespace on entry, several full
//! frames, trailing whitespace on exit, and absence afterwards.

use super::{
    handle_viewer_key,
    render::sequence_columns,
    state::{Viewer, ViewerRecords, fetch_viewer_records},
};
use crate::cli::InitialPosition;
use crossterm::event::KeyCode;
use nanalogue_core::{
    region_sequences::RegionSequence,
    simulate_mod_bam::{AlignmentFormat, SimulationConfig, TempBamSimulation},
};
use rust_htslib::bam::{self, Read as _};
use std::{
    collections::{BTreeMap, BTreeSet},
    error::Error,
    path::{Path, PathBuf},
};

const WINDOW_LEN: u32 = 10;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum VisibleState {
    RightBoundaryHasSequence,
    Full,
    LeftBoundaryHasSequence,
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
enum LifecycleState {
    AbsentBefore,
    RightBoundaryHasSequence,
    Full,
    LeftBoundaryHasSequence,
    AbsentAfter,
}

#[derive(Debug)]
struct Frame {
    start: u32,
    reads: BTreeMap<String, VisibleObservation>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct VisibleObservation {
    state: VisibleState,
    deletion_columns: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct AlignmentSpan {
    start: u32,
    end: u32,
}

fn classify_visible(record: &RegionSequence) -> Result<VisibleObservation, Box<dyn Error>> {
    let rendered = sequence_columns(
        record,
        u16::try_from(WINDOW_LEN).expect("test window length fits u16"),
        false,
    );
    assert_eq!(
        rendered.len(),
        usize::try_from(WINDOW_LEN).expect("test window length fits usize")
    );
    let has_left_boundary = rendered.as_bytes().first().copied() != Some(b' ');
    let has_right_boundary = rendered.as_bytes().last().copied() != Some(b' ');

    let state = match (has_left_boundary, has_right_boundary) {
        (false, true) => VisibleState::RightBoundaryHasSequence,
        (true, true) => VisibleState::Full,
        (true, false) => VisibleState::LeftBoundaryHasSequence,
        (false, false) => {
            return Err(format!(
                "read {} is contained within one viewport; the fixture promises longer reads",
                record.read_id()
            )
            .into());
        }
    };
    Ok(VisibleObservation {
        state,
        deletion_columns: rendered.bytes().filter(|&base| base == b'.').count(),
    })
}

fn capture_frame(viewer: &Viewer, cached_records: &ViewerRecords) -> Result<Frame, Box<dyn Error>> {
    let ViewerRecords::Table(table_records) = cached_records else {
        return Err("lifecycle test must remain in table mode".into());
    };
    let mut reads = BTreeMap::new();
    for record in table_records {
        let previous = reads.insert(String::from(record.read_id()), classify_visible(record)?);
        assert!(previous.is_none(), "visible read IDs must be unique");
    }
    Ok(Frame {
        start: viewer.viewport.start,
        reads,
    })
}

fn expected_visible_state(
    span: AlignmentSpan,
    viewport_start: u32,
) -> Result<Option<VisibleState>, Box<dyn Error>> {
    let viewport_end = viewport_start.saturating_add(WINDOW_LEN);
    if span.end <= viewport_start || span.start >= viewport_end {
        return Ok(None);
    }
    Ok(
        match (span.start <= viewport_start, span.end >= viewport_end) {
            (false, true) => Some(VisibleState::RightBoundaryHasSequence),
            (true, true) => Some(VisibleState::Full),
            (true, false) => Some(VisibleState::LeftBoundaryHasSequence),
            (false, false) => {
                return Err(format!(
                    "alignment {span:?} is contained within viewport {viewport_start}..{viewport_end}; the fixture promises longer reads"
                )
                .into());
            }
        },
    )
}

const fn lifecycle_state(current: LifecycleState, visible: Option<VisibleState>) -> LifecycleState {
    match visible {
        Some(VisibleState::RightBoundaryHasSequence) => LifecycleState::RightBoundaryHasSequence,
        Some(VisibleState::Full) => LifecycleState::Full,
        Some(VisibleState::LeftBoundaryHasSequence) => LifecycleState::LeftBoundaryHasSequence,
        None if matches!(current, LifecycleState::AbsentBefore) => LifecycleState::AbsentBefore,
        None => LifecycleState::AbsentAfter,
    }
}

const fn transition_is_allowed(from: LifecycleState, to: LifecycleState) -> bool {
    matches!(
        (from, to),
        (
            LifecycleState::AbsentBefore,
            LifecycleState::AbsentBefore
                | LifecycleState::RightBoundaryHasSequence
                | LifecycleState::Full
        ) | (
            LifecycleState::RightBoundaryHasSequence,
            LifecycleState::Full
        ) | (
            LifecycleState::Full,
            LifecycleState::Full
                | LifecycleState::LeftBoundaryHasSequence
                | LifecycleState::AbsentAfter
        ) | (
            LifecycleState::LeftBoundaryHasSequence | LifecycleState::AbsentAfter,
            LifecycleState::AbsentAfter,
        )
    )
}

fn mapped_alignment_spans(path: &Path) -> Result<BTreeMap<String, AlignmentSpan>, Box<dyn Error>> {
    let mut bam_reader = bam::Reader::from_path(path)?;
    let mut alignment_spans = BTreeMap::new();
    for record_result in bam_reader.records() {
        let record = record_result?;
        if record.is_unmapped() {
            continue;
        }
        let read_id = String::from_utf8(record.qname().to_vec())?;
        let span = AlignmentSpan {
            start: u32::try_from(record.pos())?,
            end: u32::try_from(record.cigar().end_pos())?,
        };
        assert_eq!(span.end.saturating_sub(span.start), 50);
        assert!(
            alignment_spans.insert(read_id, span).is_none(),
            "simulated read IDs must be unique"
        );
    }
    Ok(alignment_spans)
}

fn capture_scroll_frames(viewer: &mut Viewer) -> Result<Vec<Frame>, Box<dyn Error>> {
    let mut records = fetch_viewer_records(viewer)?;
    let mut frames = Vec::new();

    loop {
        frames.push(capture_frame(viewer, &records)?);
        if viewer.viewport.start.saturating_add(WINDOW_LEN) >= viewer.target_len() {
            break;
        }
        assert!(
            handle_viewer_key(viewer, &mut records, KeyCode::Char('l'), usize::MAX)?,
            "right navigation must advance before the final viewport"
        );
    }
    Ok(frames)
}

fn lifecycle_is_complete(
    read_id: &str,
    span: AlignmentSpan,
    frames: &[Frame],
) -> Result<bool, Box<dyn Error>> {
    let mut current = LifecycleState::AbsentBefore;
    let mut observed = BTreeSet::new();
    let mut deletion_columns = 0usize;
    let mut full_frames = 0usize;

    for frame in frames {
        let visible = frame.reads.get(read_id).copied();
        let expected = expected_visible_state(span, frame.start)?;
        assert_eq!(
            visible.map(|observation| observation.state),
            expected,
            "read {read_id} with alignment {span:?} has incorrect visibility in viewport {}..{}",
            frame.start,
            frame.start.saturating_add(WINDOW_LEN)
        );
        let next = lifecycle_state(current, visible.map(|observation| observation.state));
        assert!(
            transition_is_allowed(current, next),
            "read {read_id} made forbidden transition {current:?} -> {next:?} at viewport {}..{}",
            frame.start,
            frame.start.saturating_add(WINDOW_LEN)
        );
        let _new_state = observed.insert(next);
        deletion_columns = deletion_columns
            .saturating_add(visible.map_or(0, |observation| observation.deletion_columns));
        full_frames = full_frames.saturating_add(usize::from(matches!(next, LifecycleState::Full)));
        current = next;
    }

    assert_eq!(
        deletion_columns, 20,
        "read {read_id} must expose its complete centered 20 bp deletion"
    );
    let complete = [
        LifecycleState::AbsentBefore,
        LifecycleState::RightBoundaryHasSequence,
        LifecycleState::Full,
        LifecycleState::LeftBoundaryHasSequence,
        LifecycleState::AbsentAfter,
    ]
    .into_iter()
    .all(|state| observed.contains(&state));
    if complete {
        assert!(
            full_frames > 1,
            "complete lifecycle for read {read_id} must include several full frames"
        );
    }
    Ok(complete)
}

#[test]
fn simulated_reads_follow_the_full_scroll_lifecycle() -> Result<(), Box<dyn Error>> {
    let config: SimulationConfig = serde_json::from_str(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/examples/bam_viewer_full_lifecycle_demo.json"
    )))?;
    let simulation = TempBamSimulation::new(config, AlignmentFormat::Bam)?;
    let alignment_spans = mapped_alignment_spans(Path::new(simulation.bam_path()))?;
    let mut viewer = Viewer::open(
        PathBuf::from(simulation.bam_path()),
        &InitialPosition {
            contig: String::from("contig_00000"),
            start: 0,
        },
        None,
        WINDOW_LEN,
    )?;
    let frames = capture_scroll_frames(&mut viewer)?;

    assert_eq!(
        frames.iter().map(|frame| frame.start).collect::<Vec<_>>(),
        (0..viewer.target_len())
            .step_by(WINDOW_LEN as usize)
            .collect::<Vec<_>>()
    );
    let read_ids = frames
        .iter()
        .flat_map(|frame| frame.reads.keys().cloned())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        read_ids,
        alignment_spans.keys().cloned().collect(),
        "the viewer must show every mapped simulated read and no others"
    );

    let mut complete_non_boundary_lifecycles = 0usize;
    for read_id in read_ids {
        let span = *alignment_spans
            .get(&read_id)
            .expect("viewer read ID came from the mapped BAM records");
        if lifecycle_is_complete(&read_id, span, &frames)? {
            complete_non_boundary_lifecycles = complete_non_boundary_lifecycles.saturating_add(1);
        }
    }

    // Six of the simulator's seven equally likely read states are mapped. A 50 bp read has
    // 51 possible leftmost starts on a 100 bp contig: 0..=50. Of those, 24 cannot show all
    // five lifecycle stages: 0..=10 lack an observed absent-before or partial-entry stage;
    // 20 and 30 are boundary-aligned; and 40..=50 lack a partial-exit or observed
    // absent-after stage. Thus 10_000 * (6 / 7) * (27 / 51) is about 4_538 reads.
    assert!(
        (4_488..=4_588).contains(&complete_non_boundary_lifecycles),
        "expected 4,488 to 4,588 complete non-boundary lifecycles, observed {complete_non_boundary_lifecycles}"
    );
    Ok(())
}
