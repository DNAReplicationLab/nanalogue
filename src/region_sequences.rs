//! Sequence retrieval for the interactive BAM viewer.

use crate::{
    BamPreFilt as _, CurrRead, Error, F32Bw0and1, ModChar, SeqCoordCalls, ThresholdState,
    analysis::threshold_and_mean,
    constants::shared::{MAX_RECORD_CAPACITY_BYTES, MAX_RECORDS},
    ensure_bounded_counter, ensure_record_data_capacity, nanalogue_indexed_bam_reader,
    read_utils::{CurrReadState, CurrReadStateWithAlign},
};
use rust_htslib::bam;
use std::{collections::HashMap, num::NonZeroU32};

/// A read ID and sequence projected onto a requested reference region.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct RegionSequence {
    /// Read identifier.
    read_id: String,
    /// Number of requested reference positions before this alignment begins.
    region_offset: u32,
    /// Region sequence with insertions omitted and deletions represented by `.`.
    sequence: String,
    /// Region sequence with insertions represented by lowercase bases.
    sequence_with_insertions: String,
    /// Modification calls corresponding to `sequence`.
    modifications: Vec<bool>,
    /// Modification calls corresponding to `sequence_with_insertions`.
    modifications_with_insertions: Vec<bool>,
    /// Whether the alignment is on the reverse strand.
    reverse: bool,
}

impl RegionSequence {
    /// Returns the read identifier.
    #[must_use]
    pub fn read_id(&self) -> &str {
        &self.read_id
    }

    /// Returns the number of requested reference positions before this alignment begins.
    #[must_use]
    pub fn region_offset(&self) -> u32 {
        self.region_offset
    }

    /// Returns the sequence intersecting the requested region.
    #[must_use]
    pub fn sequence(&self) -> &str {
        &self.sequence
    }

    /// Returns the sequence intersecting the requested region with lowercase insertions.
    #[must_use]
    pub fn sequence_with_insertions(&self) -> &str {
        &self.sequence_with_insertions
    }

    /// Returns modification calls corresponding to [`Self::sequence`].
    #[must_use]
    pub fn modifications(&self) -> &[bool] {
        &self.modifications
    }

    /// Returns modification calls corresponding to [`Self::sequence_with_insertions`].
    #[must_use]
    pub fn modifications_with_insertions(&self) -> &[bool] {
        &self.modifications_with_insertions
    }

    /// Returns whether the alignment is on the reverse strand.
    #[must_use]
    pub fn is_reverse(&self) -> bool {
        self.reverse
    }
}

/// Raw and windowed modification calls for one aligned read.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct ReadModProfile {
    /// BAM alignment identity used when a viewer window is refetched.
    identity: AlignmentIdentity,
    /// Reference position and raw ML probability for mapped calls.
    calls: Vec<(u32, u8)>,
    /// Reference bounds and thresholded mean for each complete call window.
    windows: Vec<(u32, u32, F32Bw0and1)>,
    /// Indices in `windows` at which independent modification series begin.
    window_series_starts: Vec<usize>,
}

/// BAM fields used to identify one alignment across viewer window refetches.
#[derive(Debug, Clone, PartialEq, Eq)]
struct AlignmentIdentity {
    /// Read identifier.
    read_id: String,
    /// Supported alignment and strand state.
    read_state: crate::ReadState,
    /// Target identifier.
    tid: i32,
    /// Zero-based alignment start.
    pos: u32,
    /// Mapping quality.
    mapq: u8,
    /// Zero-based exclusive alignment end on the reference.
    reference_end: u32,
    /// Occurrence among records with the same preceding identity fields.
    duplicate_index: u32,
}

/// BAM fields shared by alignments that require a duplicate index.
type AlignmentIdentityBase = (String, crate::ReadState, i32, u32, u8, u32);

/// Returns one read's alignment identity with an index among equal base identities.
fn alignment_identity<S>(
    read: &CurrRead<S>,
    occurrences: &mut HashMap<AlignmentIdentityBase, u32>,
) -> Result<AlignmentIdentity, Error>
where
    S: CurrReadState + CurrReadStateWithAlign,
{
    let read_id = String::from(read.read_id());
    let read_state = read.read_state();
    let (tid, pos) = read.contig_id_and_start()?;
    let mapq = read.mapq();
    let reference_end = pos.checked_add(read.align_len()?).ok_or_else(|| {
        Error::InvalidState(String::from("alignment reference end exceeds u32::MAX"))
    })?;
    let next_duplicate_index = occurrences
        .entry((read_id.clone(), read_state, tid, pos, mapq, reference_end))
        .or_insert(0u32);
    let duplicate_index = *next_duplicate_index;
    *next_duplicate_index = next_duplicate_index.checked_add(1).ok_or_else(|| {
        Error::InvalidState(String::from(
            "too many identical alignments in modification profiles",
        ))
    })?;
    Ok(AlignmentIdentity {
        read_id,
        read_state,
        tid,
        pos,
        mapq,
        reference_end,
        duplicate_index,
    })
}

impl ReadModProfile {
    /// Returns the read identifier.
    #[must_use]
    pub fn read_id(&self) -> &str {
        &self.identity.read_id
    }

    /// Returns whether two profiles represent the same BAM alignment.
    ///
    /// Unlike comparing read identifiers, this distinguishes alignment state,
    /// coordinates, mapping quality, and otherwise equal duplicate-name records.
    #[must_use]
    pub fn is_same_alignment(&self, other: &Self) -> bool {
        self.identity == other.identity
    }

    /// Returns whether the alignment is on the reverse strand.
    #[must_use]
    pub fn is_reverse(&self) -> bool {
        self.identity.read_state.strand() == '-'
    }

    /// Returns the zero-based inclusive alignment start.
    #[must_use]
    pub fn align_start(&self) -> u32 {
        self.identity.pos
    }

    /// Returns the zero-based exclusive alignment end.
    #[must_use]
    pub fn align_end(&self) -> u32 {
        self.identity.reference_end
    }

    /// Returns mapped raw calls as `(reference position, ML probability)`.
    #[must_use]
    pub fn calls(&self) -> &[(u32, u8)] {
        &self.calls
    }

    /// Returns complete call windows as `(reference start, reference end, value)`.
    #[must_use]
    pub fn windows(&self) -> &[(u32, u32, F32Bw0and1)] {
        &self.windows
    }

    /// Returns the window indices at which independent mod-strand series begin.
    #[must_use]
    pub fn window_series_starts(&self) -> &[usize] {
        &self.window_series_starts
    }
}

/// Indexed reader that exposes read-table region sequences without leaking BAM types.
#[derive(Debug)]
pub struct RegionSequenceReader {
    /// Indexed alignment reader.
    reader: bam::IndexedReader,
    /// Reference names copied from the alignment header.
    target_names: Vec<String>,
    /// Reference lengths copied from the alignment header.
    target_lengths: Vec<u32>,
}

/// Formats nanalogue's reference-coordinate calls with and without insertion entries.
fn format_region_sequences<I>(
    coordinates: I,
    read_sequence: &[u8],
    read_modifications: &[bool],
) -> Result<(String, String, Vec<bool>, Vec<bool>), Error>
where
    I: IntoIterator<Item = Option<(bool, u32)>>,
{
    let mut sequence = Vec::new();
    let mut sequence_with_insertions = Vec::new();
    let mut modifications = Vec::new();
    let mut modifications_with_insertions = Vec::new();
    for base in coordinates {
        match base {
            Some((true, sequence_position)) => {
                let index = usize::try_from(sequence_position)?;
                let uppercase_nucleotide = read_sequence
                    .get(index)
                    .ok_or_else(|| {
                        Error::UnavailableData(String::from(
                            "sequence coordinate is outside the read sequence",
                        ))
                    })?
                    .to_ascii_uppercase();
                let modified = *read_modifications.get(index).ok_or_else(|| {
                    Error::UnavailableData(String::from(
                        "sequence coordinate is outside the modification calls",
                    ))
                })?;
                sequence.push(uppercase_nucleotide);
                sequence_with_insertions.push(uppercase_nucleotide);
                modifications.push(modified);
                modifications_with_insertions.push(modified);
            }
            None => {
                sequence.push(b'.');
                sequence_with_insertions.push(b'.');
                modifications.push(false);
                modifications_with_insertions.push(false);
            }
            Some((false, sequence_position)) => {
                let index = usize::try_from(sequence_position)?;
                let insertion = read_sequence
                    .get(index)
                    .ok_or_else(|| {
                        Error::UnavailableData(String::from(
                            "insertion coordinate is outside the read sequence",
                        ))
                    })?
                    .to_ascii_lowercase();
                let modified = *read_modifications.get(index).ok_or_else(|| {
                    Error::UnavailableData(String::from(
                        "insertion coordinate is outside the modification calls",
                    ))
                })?;
                sequence_with_insertions.push(insertion);
                modifications_with_insertions.push(modified);
            }
        }
    }
    Ok((
        String::from_utf8(sequence)?,
        String::from_utf8(sequence_with_insertions)?,
        modifications,
        modifications_with_insertions,
    ))
}

impl RegionSequenceReader {
    /// Opens an indexed BAM file.
    ///
    /// # Errors
    /// Returns an error if the BAM, its index, or its reference header is invalid.
    pub fn from_path<P>(path: &P) -> Result<Self, Error>
    where
        P: AsRef<std::path::Path> + ?Sized,
    {
        let reader = nanalogue_indexed_bam_reader(path, bam::FetchDefinition::All)?;
        let header = bam::Read::header(&reader);
        if header.target_count() == 0 {
            return Err(Error::UnavailableData(String::from(
                "BAM header has no reference sequences",
            )));
        }
        let target_names = header
            .target_names()
            .into_iter()
            .map(|name| String::from_utf8(name.to_vec()).map_err(Error::from))
            .collect::<Result<Vec<_>, _>>()?;
        let target_lengths = (0..header.target_count())
            .map(|tid| {
                header
                    .target_len(tid)
                    .ok_or_else(|| {
                        Error::UnavailableData(String::from("BAM header target has no length"))
                    })?
                    .try_into()
                    .map_err(Error::from)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            reader,
            target_names,
            target_lengths,
        })
    }

    /// Returns the numeric identifier for a reference name.
    #[must_use]
    pub fn target_id(&self, name: &str) -> Option<u32> {
        let tid = self.target_names.iter().position(|target| target == name)?;
        u32::try_from(tid).ok()
    }

    /// Returns a reference name by numeric identifier.
    #[must_use]
    pub fn target_name(&self, tid: u32) -> Option<&str> {
        self.target_names
            .get(usize::try_from(tid).ok()?)
            .map(String::as_str)
    }

    /// Returns a reference length by numeric identifier.
    #[must_use]
    pub fn target_len(&self, tid: u32) -> Option<u32> {
        self.target_lengths.get(usize::try_from(tid).ok()?).copied()
    }

    /// Validates a half-open reference interval against the BAM header.
    fn validate_region(&self, tid: u32, start: u32, end: u32) -> Result<(), Error> {
        if start >= end {
            return Err(Error::InvalidAlignCoords(format!("{tid}:{start}-{end}")));
        }
        let target_len = self.target_len(tid).ok_or_else(|| {
            Error::InvalidAlignCoords(format!("target {tid} is not in the BAM header"))
        })?;
        if start >= target_len || end > target_len {
            let target_name = self.target_name(tid).unwrap_or("unknown");
            return Err(Error::InvalidAlignCoords(format!(
                "{tid}:{start}-{end} is outside reference '{target_name}' (length {target_len})"
            )));
        }
        Ok(())
    }

    /// Retrieves reads overlapping a reference interval and their projected sequences.
    ///
    /// Insertions are omitted and deletions or reference skips are represented by `.`.
    /// When `mod_type` is present, matching calls with probability at least 128 are marked.
    /// Modification tags are not parsed when `mod_type` is absent.
    /// Rows are sorted by read ID. [`RegionSequence::region_offset`] reports any requested
    /// reference positions before an alignment starts. A mapped record whose sequence is omitted
    /// is returned as `*`.
    ///
    /// # Errors
    /// Returns an error if the interval or any fetched alignment record is invalid.
    pub fn sequences(
        &mut self,
        tid: u32,
        start: u32,
        end: u32,
        mod_type: Option<ModChar>,
    ) -> Result<Vec<RegionSequence>, Error> {
        self.validate_region(tid, start, end)?;
        self.reader.fetch((tid, i64::from(start), i64::from(end)))?;
        let region = crate::GenomicBed3::new(i32::try_from(tid)?, start, end);
        let mut rows = Vec::new();
        let mut record_count = 0u32;
        for record_result in bam::Read::records(&mut self.reader) {
            let record = record_result?;
            if !record.filt_by_region(&region, false) {
                continue;
            }
            ensure_bounded_counter(&mut record_count, MAX_RECORDS, "region sequence table")?;
            ensure_record_data_capacity(
                record.inner().m_data,
                MAX_RECORD_CAPACITY_BYTES,
                "region sequence table",
            )?;
            let (curr_read, has_sequence) =
                match CurrRead::default().try_from_only_alignment(&record) {
                    Ok(curr_read) => (curr_read, true),
                    Err(Error::ZeroSeqLen(_)) => (
                        CurrRead::default().try_from_only_alignment_zero_seq_len(&record)?,
                        false,
                    ),
                    Err(error) => return Err(error),
                };
            let read_id = String::from(curr_read.read_id());
            let reverse = curr_read.strand() == '-';
            let (_alignment_tid, alignment_start) = curr_read.contig_id_and_start()?;
            let region_offset = alignment_start.saturating_sub(start);
            let (sequence, sequence_with_insertions, modifications, modifications_with_insertions) =
                if has_sequence {
                    let read_sequence = record.seq().as_bytes();
                    let coordinates = curr_read.seq_coords_from_ref_coords(&record, &region)?;
                    let read_modifications = if let Some(requested_mod_type) = mod_type {
                        let curr_read_with_mods = curr_read.set_mod_data_restricted(
                            &record,
                            ThresholdState::GtEq(128),
                            |_| true,
                            |_, _, observed_mod_type| *observed_mod_type == requested_mod_type,
                            0,
                        )?;
                        match SeqCoordCalls::try_from(&curr_read_with_mods.mod_data().0) {
                            Ok(calls) => calls.collapse_mod_calls(),
                            Err(Error::UnavailableData(_)) => vec![false; read_sequence.len()],
                            Err(error) => return Err(error),
                        }
                    } else {
                        vec![false; read_sequence.len()]
                    };
                    format_region_sequences(coordinates, &read_sequence, &read_modifications)?
                } else {
                    (
                        String::from("*"),
                        String::from("*"),
                        vec![false],
                        vec![false],
                    )
                };
            rows.push(RegionSequence {
                read_id,
                region_offset,
                sequence,
                sequence_with_insertions,
                modifications,
                modifications_with_insertions,
                reverse,
            });
        }
        rows.sort_by(|left, right| left.read_id.cmp(&right.read_id));
        Ok(rows)
    }

    /// Retrieves raw and non-overlapping windowed modification profiles for reads spanning a
    /// complete reference interval.
    ///
    /// Raw insertion calls participate in windows but are omitted from [`ReadModProfile::calls`].
    /// Windows containing no mapped calls are omitted. Rows are sorted by read ID.
    ///
    /// # Errors
    /// Returns an error if the interval, an alignment, or its modification tags are invalid.
    pub fn profiles(
        &mut self,
        tid: u32,
        start: u32,
        end: u32,
        mod_type: ModChar,
        win: NonZeroU32,
    ) -> Result<Vec<ReadModProfile>, Error> {
        self.validate_region(tid, start, end)?;
        self.reader.fetch((tid, i64::from(start), i64::from(end)))?;
        let region = crate::GenomicBed3::new(i32::try_from(tid)?, start, end);
        let mut profiles = Vec::new();
        let mut alignment_occurrences = HashMap::new();
        let mut record_count = 0u32;
        let win_size = usize::try_from(win.get())?;

        for record_result in bam::Read::records(&mut self.reader) {
            let record = record_result?;
            if !record.filt_by_region(&region, true) {
                continue;
            }
            ensure_bounded_counter(
                &mut record_count,
                MAX_RECORDS,
                "region modification profiles",
            )?;
            ensure_record_data_capacity(
                record.inner().m_data,
                MAX_RECORD_CAPACITY_BYTES,
                "region modification profiles",
            )?;
            let (curr_read, has_sequence) =
                match CurrRead::default().try_from_only_alignment(&record) {
                    Ok(curr_read) => (curr_read, true),
                    Err(Error::ZeroSeqLen(_)) => (
                        CurrRead::default().try_from_only_alignment_zero_seq_len(&record)?,
                        false,
                    ),
                    Err(error) => return Err(error),
                };
            let identity = alignment_identity(&curr_read, &mut alignment_occurrences)?;
            let mut calls = Vec::new();
            let mut windows = Vec::new();
            let mut window_series_starts = Vec::new();

            if has_sequence {
                let read_with_mods = curr_read.set_mod_data_restricted(
                    &record,
                    ThresholdState::GtEq(0),
                    |_| true,
                    |_, _, observed_mod_type| *observed_mod_type == mod_type,
                    0,
                )?;
                for base_mod in &read_with_mods.mod_data().0.base_mods {
                    calls.extend(base_mod.ranges.annotations.iter().filter_map(|annotation| {
                        annotation.ref_pos.map(|ref_pos| (ref_pos, annotation.qual))
                    }));
                    let series_start = windows.len();
                    for chunk in base_mod.ranges.annotations.chunks_exact(win_size) {
                        let mut reference_positions = chunk.iter().filter_map(|item| item.ref_pos);
                        let Some(first_reference_position) = reference_positions.next() else {
                            continue;
                        };
                        let (ref_win_start, ref_win_max) = reference_positions.fold(
                            (first_reference_position, first_reference_position),
                            |(minimum, maximum), position| {
                                (minimum.min(position), maximum.max(position))
                            },
                        );
                        let ref_win_end = ref_win_max.checked_add(1).ok_or_else(|| {
                            Error::InvalidState(String::from(
                                "reference modification window ends at u32::MAX",
                            ))
                        })?;
                        let probabilities = chunk.iter().map(|item| item.qual).collect::<Vec<_>>();
                        windows.push((
                            ref_win_start,
                            ref_win_end,
                            threshold_and_mean(&probabilities)?,
                        ));
                    }
                    if windows.len() > series_start {
                        window_series_starts.push(series_start);
                    }
                }
                calls.sort_unstable_by_key(|&(ref_pos, _probability)| ref_pos);
            }

            profiles.push(ReadModProfile {
                identity,
                calls,
                windows,
                window_series_starts,
            });
        }
        profiles.sort_by(|left, right| left.read_id().cmp(right.read_id()));
        Ok(profiles)
    }
}

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use super::*;
    use crate::{file_utils::write_bam_denovo, uuid};
    use rust_htslib::bam::record::{Aux, Cigar, CigarString};

    fn write_test_bam(mm: &str, ml: &[u8]) -> Result<std::path::PathBuf, Error> {
        let mut record = bam::Record::new();
        record.set_tid(0);
        record.set_pos(0);
        record.set_mapq(60);
        record.unset_unmapped();
        record.set(
            b"read",
            Some(&CigarString::from(vec![Cigar::Match(5)])),
            b"ACGTA",
            &[30; 5],
        );
        record.push_aux(b"MM", Aux::String(mm))?;
        record.push_aux(b"ML", Aux::ArrayU8(ml.into()))?;
        let path = std::env::temp_dir().join(format!("{}.bam", uuid::v4_random()));
        write_bam_denovo(
            [record],
            [(String::from("chr1"), 5)],
            [String::from("rg1")],
            Vec::<String>::new(),
            &path,
        )?;
        Ok(path)
    }

    fn write_zero_sequence_test_bam(
        cigar: Vec<Cigar>,
        mod_tags: Option<(&str, &[u8])>,
    ) -> Result<std::path::PathBuf, Error> {
        let mut record = bam::Record::new();
        record.set_tid(0);
        record.set_pos(0);
        record.set_mapq(60);
        record.unset_unmapped();
        record.set(
            b"sequence-not-stored",
            Some(&CigarString::from(cigar)),
            b"",
            &[],
        );
        if let Some((mm, ml)) = mod_tags {
            record.push_aux(b"MM", Aux::String(mm))?;
            record.push_aux(b"ML", Aux::ArrayU8(ml.into()))?;
        }
        let path = std::env::temp_dir().join(format!("{}.bam", uuid::v4_random()));
        write_bam_denovo(
            [record],
            [(String::from("chr1"), 20)],
            [String::from("rg1")],
            Vec::<String>::new(),
            &path,
        )?;
        Ok(path)
    }

    fn assert_zero_sequence_placeholder(row: &RegionSequence) {
        assert_eq!(row.read_id(), "sequence-not-stored");
        assert_eq!(row.region_offset(), 0);
        assert_eq!(row.sequence(), "*");
        assert_eq!(row.sequence_with_insertions(), "*");
        assert_eq!(row.modifications(), [false]);
        assert_eq!(row.modifications_with_insertions(), [false]);
    }

    fn remove_test_bam(path: &std::path::Path) {
        std::fs::remove_file(path).expect("remove test BAM");
        std::fs::remove_file(format!("{}.bai", path.display())).expect("remove test BAM index");
    }

    #[test]
    fn alignment_identity_uses_read_state_end_and_duplicate_index() -> Result<(), Error> {
        let read = |end, read_state| {
            crate::CurrReadBuilder::default()
                .alignment_type(read_state)
                .alignment(
                    crate::AlignmentInfoBuilder::default()
                        .start(5)
                        .end(end)
                        .contig(String::from("chr1"))
                        .contig_id(0)
                        .build()?,
                )
                .read_id(String::from("shared-name"))
                .mapq(60)
                .seq_len(20)
                .build()
        };
        let mut occurrences = HashMap::new();
        let first = read(25, crate::ReadState::PrimaryFwd)?;
        let first_identity = alignment_identity(&first, &mut occurrences)?;
        let same_bounds = read(25, crate::ReadState::PrimaryFwd)?;
        let same_bounds_identity = alignment_identity(&same_bounds, &mut occurrences)?;
        assert_eq!(first_identity.duplicate_index, 0);
        assert_eq!(same_bounds_identity.duplicate_index, 1);

        let longer = read(26, crate::ReadState::PrimaryFwd)?;
        let longer_identity = alignment_identity(&longer, &mut occurrences)?;
        assert_eq!(longer_identity.reference_end, 26);
        assert_eq!(longer_identity.duplicate_index, 0);

        let reverse = read(25, crate::ReadState::PrimaryRev)?;
        let reverse_identity = alignment_identity(&reverse, &mut occurrences)?;
        assert_eq!(reverse_identity.duplicate_index, 0);
        Ok(())
    }

    #[test]
    fn profiles_keep_mapped_calls_and_window_in_read_order() -> Result<(), Error> {
        let mut record = bam::Record::new();
        record.set_tid(0);
        record.set_pos(10);
        record.set_mapq(60);
        record.unset_unmapped();
        record.set(
            b"profile",
            Some(&CigarString::from(vec![
                Cigar::Match(2),
                Cigar::Ins(1),
                Cigar::Match(2),
                Cigar::Del(1),
                Cigar::Match(3),
            ])),
            b"AAAAAAAA",
            &[30; 8],
        );
        record.push_aux(b"MM", Aux::String("A+a?,0,0,0,0,0,0,0,0;"))?;
        record.push_aux(
            b"ML",
            Aux::ArrayU8((&[0, 127, 255, 128, 255, 0, 255, 128][..]).into()),
        )?;
        let path = std::env::temp_dir().join(format!("{}.bam", uuid::v4_random()));
        write_bam_denovo(
            [record],
            [(String::from("chr1"), 30)],
            [String::from("rg1")],
            Vec::<String>::new(),
            &path,
        )?;
        let mut reader = RegionSequenceReader::from_path(&path)?;

        let profiles = reader.profiles(
            0,
            10,
            18,
            ModChar::new('a'),
            NonZeroU32::new(3).expect("non-zero"),
        )?;
        let profile = profiles.first().expect("one spanning read");
        assert_eq!(profile.align_start(), 10);
        assert_eq!(profile.align_end(), 18);
        assert_eq!(profile.calls().len(), 7, "the insertion call is not a dot");
        assert_eq!(
            profile
                .calls()
                .iter()
                .map(|call| call.0)
                .collect::<Vec<_>>(),
            [10, 11, 12, 13, 15, 16, 17]
        );
        assert_eq!(
            profile.windows().len(),
            2,
            "the trailing two calls are unused"
        );
        let first_window = profile.windows().first().expect("first window");
        let second_window = profile.windows().get(1).expect("second window");
        assert_eq!(first_window.0..first_window.1, 10..12);
        assert!((first_window.2.val() - 1.0 / 3.0).abs() < f32::EPSILON);
        assert_eq!(second_window.0..second_window.1, 12..16);
        assert!((second_window.2.val() - 2.0 / 3.0).abs() < f32::EPSILON);

        remove_test_bam(&path);
        Ok(())
    }

    #[test]
    fn profiles_return_empty_data_for_an_omitted_sequence() -> Result<(), Error> {
        let path = write_zero_sequence_test_bam(vec![Cigar::Match(5)], None)?;
        let mut reader = RegionSequenceReader::from_path(&path)?;

        let profiles = reader.profiles(
            0,
            0,
            5,
            ModChar::new('a'),
            NonZeroU32::new(3).expect("non-zero"),
        )?;
        let profile = profiles.first().expect("one spanning read");
        assert!(profile.calls().is_empty());
        assert!(profile.windows().is_empty());

        remove_test_bam(&path);
        Ok(())
    }

    #[test]
    fn profiles_reject_empty_and_reversed_ranges() -> Result<(), Error> {
        let path = write_test_bam("A+a?,0;", &[255])?;
        let mut reader = RegionSequenceReader::from_path(&path)?;
        let win = NonZeroU32::new(3).expect("non-zero");

        for (start, end) in [(2, 2), (4, 1)] {
            let error = reader
                .profiles(0, start, end, ModChar::new('a'), win)
                .expect_err("invalid ranges must be rejected");
            assert!(matches!(error, Error::InvalidAlignCoords(_)));
        }

        remove_test_bam(&path);
        Ok(())
    }

    #[test]
    fn profiles_preserve_independent_mod_strand_window_series() -> Result<(), Error> {
        let path = write_test_bam("A+a?,0,0;T-a?,0;", &[255, 255, 0])?;
        let mut reader = RegionSequenceReader::from_path(&path)?;

        let profiles = reader.profiles(
            0,
            0,
            5,
            ModChar::new('a'),
            NonZeroU32::new(1).expect("non-zero"),
        )?;
        let profile = profiles.first().expect("one spanning read");
        assert_eq!(profile.windows().len(), 3);
        assert_eq!(profile.window_series_starts(), [0, 2]);

        remove_test_bam(&path);
        Ok(())
    }

    #[test]
    fn region_sequences_optionally_show_modified_lowercase_insertions() {
        let (sequence, sequence_with_insertions, modifications, modifications_with_insertions) =
            format_region_sequences(
                [Some((true, 0)), None, Some((false, 1)), Some((true, 2))],
                b"aGt",
                &[true, true, false],
            )
            .expect("valid sequence calls");
        assert_eq!(sequence, "A.T");
        assert_eq!(sequence_with_insertions, "A.gT");
        assert_eq!(modifications, [true, false, false]);
        assert_eq!(modifications_with_insertions, [true, false, true, false]);
    }

    #[test]
    fn region_sequences_filter_mod_type_at_inclusive_threshold() -> Result<(), Error> {
        let path = write_test_bam("A+a?,0,0;C+m?,0;", &[127, 128, 200])?;
        let mut reader = RegionSequenceReader::from_path(&path)?;

        let no_mods = reader.sequences(0, 0, 5, None)?;
        assert_eq!(
            no_mods.first().expect("one read").modifications(),
            [false; 5]
        );
        let a_mods = reader.sequences(0, 0, 5, Some(ModChar::new('a')))?;
        assert_eq!(
            a_mods.first().expect("one read").modifications(),
            [false, false, false, false, true]
        );
        let m_mods = reader.sequences(0, 0, 5, Some(ModChar::new('m')))?;
        assert_eq!(
            m_mods.first().expect("one read").modifications(),
            [false, true, false, false, false]
        );
        let absent_mods = reader.sequences(0, 0, 5, Some(ModChar::new('h')))?;
        assert_eq!(
            absent_mods.first().expect("one read").modifications(),
            [false; 5]
        );

        remove_test_bam(&path);
        Ok(())
    }

    #[test]
    fn omitted_mod_type_does_not_parse_malformed_tags() -> Result<(), Error> {
        let path = write_test_bam("A+a?,0,0;", &[200])?;
        let mut reader = RegionSequenceReader::from_path(&path)?;

        let no_mods = reader.sequences(0, 0, 5, None)?;
        assert_eq!(
            no_mods.first().expect("one read").modifications(),
            [false; 5]
        );
        let _error = reader
            .sequences(0, 0, 5, Some(ModChar::new('a')))
            .expect_err("malformed modification tags should fail when parsed");

        remove_test_bam(&path);
        Ok(())
    }

    #[test]
    fn zero_length_sequence_is_returned_as_an_asterisk() -> Result<(), Error> {
        let path = write_zero_sequence_test_bam(vec![Cigar::Match(5)], None)?;
        let mut reader = RegionSequenceReader::from_path(&path)?;

        let rows = reader.sequences(0, 0, 5, None)?;
        assert_zero_sequence_placeholder(rows.first().expect("one read"));
        let clipped_rows = reader.sequences(0, 2, 5, None)?;
        assert_zero_sequence_placeholder(clipped_rows.first().expect("one clipped read"));

        remove_test_bam(&path);
        Ok(())
    }

    #[test]
    fn zero_length_sequence_ignores_complex_cigar_operations() -> Result<(), Error> {
        let path = write_zero_sequence_test_bam(
            vec![
                Cigar::Match(2),
                Cigar::Ins(2),
                Cigar::Del(1),
                Cigar::RefSkip(1),
                Cigar::Match(3),
            ],
            None,
        )?;
        let mut reader = RegionSequenceReader::from_path(&path)?;

        let rows = reader.sequences(0, 0, 7, None)?;
        assert_zero_sequence_placeholder(rows.first().expect("one read"));

        remove_test_bam(&path);
        Ok(())
    }

    #[test]
    fn zero_length_sequence_skips_malformed_modification_tags() -> Result<(), Error> {
        let path =
            write_zero_sequence_test_bam(vec![Cigar::Match(5)], Some(("A+a?,0,0;", &[200])))?;
        let mut reader = RegionSequenceReader::from_path(&path)?;

        let rows = reader.sequences(0, 0, 5, Some(ModChar::new('a')))?;
        assert_zero_sequence_placeholder(rows.first().expect("one read"));

        remove_test_bam(&path);
        Ok(())
    }

    #[test]
    fn zero_length_sequence_skips_unpaired_modification_probabilities() -> Result<(), Error> {
        let path = write_zero_sequence_test_bam(vec![Cigar::Match(5)], Some(("", &[200])))?;
        let mut reader = RegionSequenceReader::from_path(&path)?;

        let rows = reader.sequences(0, 0, 5, Some(ModChar::new('a')))?;
        assert_zero_sequence_placeholder(rows.first().expect("one read"));

        remove_test_bam(&path);
        Ok(())
    }

    #[test]
    fn zero_length_sequence_skips_paired_modification_tags() -> Result<(), Error> {
        let path = write_zero_sequence_test_bam(vec![Cigar::Match(5)], Some(("A+a?,0;", &[200])))?;
        let mut reader = RegionSequenceReader::from_path(&path)?;

        let rows = reader.sequences(0, 0, 5, Some(ModChar::new('a')))?;
        assert_zero_sequence_placeholder(rows.first().expect("one read"));

        remove_test_bam(&path);
        Ok(())
    }
}
