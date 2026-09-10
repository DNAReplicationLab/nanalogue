//! Sequence retrieval for the interactive BAM viewer.

use crate::{
    BamPreFilt as _, CurrRead, Error,
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
fn format_region_sequences<I>(coordinates: I) -> Result<(String, String), Error>
where
    I: IntoIterator<Item = Option<(bool, u8, u8)>>,
{
    let mut sequence = Vec::new();
    let mut sequence_with_insertions = Vec::new();
    for base in coordinates {
        match base {
            Some((true, nucleotide, _quality)) => {
                let uppercase_nucleotide = nucleotide.to_ascii_uppercase();
                sequence.push(uppercase_nucleotide);
                sequence_with_insertions.push(uppercase_nucleotide);
            }
            None => {
                sequence.push(b'.');
                sequence_with_insertions.push(b'.');
            }
            Some((false, insertion, _quality)) => {
                sequence_with_insertions.push(insertion.to_ascii_lowercase());
            }
        }
    }
    Ok((
        String::from_utf8(sequence)?,
        String::from_utf8(sequence_with_insertions)?,
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
    /// Rows are sorted by read ID. A mapped record whose sequence is omitted is returned as `*`.
    ///
    /// # Errors
    /// Returns an error if the interval or any fetched alignment record is invalid.
    pub fn sequences(
        &mut self,
        tid: u32,
        start: u32,
        end: u32,
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
            let (sequence, sequence_with_insertions) = if has_sequence {
                format_region_sequences(curr_read.seq_and_qual_on_ref_coords(&record, &region)?)?
            } else {
                (String::from("*"), String::from("*"))
            };
            rows.push(RegionSequence {
                read_id: String::from(curr_read.read_id()),
                sequence,
                sequence_with_insertions,
                reverse: curr_read.strand() == '-',
            });
        }
        rows.sort_by(|left, right| left.read_id.cmp(&right.read_id));
        Ok(rows)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn region_sequences_optionally_show_lowercase_insertions() {
        let (sequence, sequence_with_insertions) = format_region_sequences([
            Some((true, b'a', 30)),
            None,
            Some((false, b'G', 30)),
            Some((true, b't', 30)),
        ])
        .expect("valid sequence calls");
        assert_eq!(sequence, "A.T");
        assert_eq!(sequence_with_insertions, "A.gT");
    }
}
