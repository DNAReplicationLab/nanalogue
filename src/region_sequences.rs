//! Sequence retrieval for the interactive BAM viewer.

use crate::{
    BamPreFilt as _, CurrRead, Error, ModChar, SeqCoordCalls, ThresholdState,
    constants::shared::{MAX_RECORD_CAPACITY_BYTES, MAX_RECORDS},
    ensure_bounded_counter, ensure_record_data_capacity, nanalogue_indexed_bam_reader,
};
use rust_htslib::bam;

/// A read ID and sequence projected onto a requested reference region.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub struct RegionSequence {
    /// Read identifier.
    read_id: String,
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

    /// Retrieves reads spanning a complete reference interval and their projected sequences.
    ///
    /// This uses the same full-region filter and coordinate conversion as the read-table command.
    /// Insertions are omitted and deletions or reference skips are represented by `.`.
    /// When `mod_type` is present, matching calls with probability at least 128 are marked.
    /// Modification tags are not parsed when `mod_type` is absent.
    /// Rows are sorted by read ID. A mapped record whose sequence is omitted is returned as `*`.
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
        if start >= end {
            return Err(Error::InvalidAlignCoords(format!("{tid}:{start}-{end}")));
        }
        self.reader.fetch((tid, i64::from(start), i64::from(end)))?;
        let region = crate::GenomicBed3::new(i32::try_from(tid)?, start, end);
        let mut rows = Vec::new();
        let mut record_count = 0u32;
        for record_result in bam::Read::records(&mut self.reader) {
            let record = record_result?;
            if !record.filt_by_region(&region, true) {
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
}

#[cfg(test)]
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

    fn remove_test_bam(path: &std::path::Path) {
        std::fs::remove_file(path).expect("remove test BAM");
        std::fs::remove_file(format!("{}.bai", path.display())).expect("remove test BAM index");
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
}
