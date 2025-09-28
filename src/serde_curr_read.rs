//! Serialization and deserialization implementations for CurrRead structs.

use crate::Error;
use crate::read_utils::{AlignAndModData, CurrRead};
use crate::utils::{ModChar, ReadState};
use fibertools_rs::utils::basemods::BaseMods;
use serde::{Deserialize, Serialize};

/// Serialized representation of CurrRead with condensed JSON format
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct SerializedCurrRead {
    /// The type of alignment (primary, secondary, supplementary, unmapped)
    pub alignment_type: ReadState,
    /// Alignment information, None for unmapped reads
    #[serde(skip_serializing_if = "Option::is_none")]
    pub alignment: Option<AlignmentInfo>,
    /// Condensed modification data table
    pub mod_table: Vec<ModTableEntry>,
    /// Read identifier
    pub read_id: String,
    /// Sequence length
    pub seq_len: u64,
}

impl Default for SerializedCurrRead {
    fn default() -> Self {
        Self {
            alignment_type: ReadState::Unmapped, // note that default is unmapped now, not primary
            alignment: None,
            mod_table: Vec::new(),
            read_id: String::new(),
            seq_len: 0,
        }
    }
}

/// Alignment information for mapped reads
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default)]
pub struct AlignmentInfo {
    /// Start position on reference
    pub start: u64,
    /// End position on reference
    pub end: u64,
    /// Contig/chromosome name
    pub contig: String,
    /// Contig/chromosome ID
    pub contig_id: i32,
}

/// Individual modification table entry
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default)]
pub struct ModTableEntry {
    /// Base that is modified (A, C, G, T, etc.)
    pub base: char,
    /// Whether this is on the plus strand
    pub is_strand_plus: bool,
    /// Modification code (character or numeric)
    pub mod_code: ModChar,
    /// Whether this modification data is implicit
    pub implicit: bool,
    /// Modification data as [start, ref_start, qual] tuples
    pub data: Vec<(u64, i64, u8)>,
}

impl SerializedCurrRead {
    /// Validates that coordinates in modification data are within expected ranges
    pub fn check_coordinates(&self) -> Result<(), Error> {
        for entry in &self.mod_table {
            // Check that implicit is false
            if entry.implicit {
                return Err(Error::DeSerializeImplicit);
            }

            for &(start, ref_start, _qual) in &entry.data {
                // Check that sequence coordinates are in range [0, seq_len)
                if start >= self.seq_len {
                    return Err(Error::InvalidSeqLength);
                }

                // Check reference coordinates based on alignment status
                match &self.alignment {
                    Some(alignment) => {
                        // For aligned reads, ref_start should be in [start, end) or -1
                        let align_range = alignment.start as i64..alignment.end as i64;
                        if ref_start != -1 && !align_range.contains(&ref_start) {
                            return Err(Error::InvalidAlignCoords);
                        }
                    }
                    None => {
                        // For unmapped reads, all ref_start values must be -1
                        if ref_start != -1 {
                            return Err(Error::InvalidAlignCoords);
                        }
                    }
                }
            }
        }

        Ok(())
    }
}

impl Serialize for CurrRead<AlignAndModData> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let serialized: SerializedCurrRead =
            self.clone().try_into().map_err(serde::ser::Error::custom)?;
        serialized.serialize(serializer)
    }
}

impl TryFrom<CurrRead<AlignAndModData>> for SerializedCurrRead {
    type Error = Error;

    fn try_from(curr_read: CurrRead<AlignAndModData>) -> Result<Self, Self::Error> {
        let alignment_type = curr_read.read_state();

        let alignment = match curr_read.read_state() {
            ReadState::Unmapped => None,
            _ => {
                let (contig_id, start) = curr_read.contig_id_and_start()?;
                let align_len = curr_read.align_len()?;
                let contig = curr_read.contig_name()?.to_string();
                let end = start + align_len;

                Some(AlignmentInfo {
                    start,
                    end,
                    contig,
                    contig_id,
                })
            }
        };

        let mod_table = condense_base_mods(&curr_read.mod_data().0)?;

        let read_id = curr_read.read_id()?.to_string();
        let seq_len = curr_read.seq_len()?;

        Ok(SerializedCurrRead {
            alignment_type,
            alignment,
            mod_table,
            read_id,
            seq_len,
        })
    }
}

impl<'de> Deserialize<'de> for CurrRead<AlignAndModData> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let serialized = SerializedCurrRead::deserialize(deserializer)?;
        serialized.try_into().map_err(serde::de::Error::custom)
    }
}

impl TryFrom<SerializedCurrRead> for CurrRead<AlignAndModData> {
    type Error = Error;

    fn try_from(serialized: SerializedCurrRead) -> Result<Self, Self::Error> {
        // Validate coordinates before proceeding with deserialization
        serialized.check_coordinates()?;

        // Reconstruct BaseMods from mod_table
        let base_mods = reconstruct_base_mods(
            &serialized.mod_table,
            serialized.alignment_type,
            serialized.seq_len,
        )?;

        // Extract alignment information
        let (align_len, contig_id_and_start, contig_name) = match &serialized.alignment {
            Some(alignment) => {
                let align_len = Some(alignment.end - alignment.start);
                let contig_id_and_start = Some((alignment.contig_id, alignment.start));
                let contig_name = Some(alignment.contig.clone());
                (align_len, contig_id_and_start, contig_name)
            }
            None => (None, None, None),
        };

        // Create CurrRead using the new constructor
        Ok(CurrRead::from_deserialized_data(
            serialized.alignment_type,
            serialized.read_id,
            serialized.seq_len,
            align_len,
            base_mods,
            contig_id_and_start,
            contig_name,
        ))
    }
}

/// Convert BaseMods to condensed mod_table format
fn condense_base_mods(base_mods: &BaseMods) -> Result<Vec<ModTableEntry>, Error> {
    let mut mod_table = Vec::new();

    for base_mod in &base_mods.base_mods {
        let entries: Result<Vec<_>, Error> = base_mod
            .ranges
            .starts
            .iter()
            .zip(base_mod.ranges.reference_starts.iter())
            .zip(base_mod.ranges.qual.iter())
            .map(|((start_opt, ref_start_opt), &qual)| {
                let start = u64::try_from(start_opt.ok_or(Error::UnavailableData)?)?;
                let ref_start = ref_start_opt.unwrap_or(-1);
                Ok((start, ref_start, qual))
            })
            .collect();
        let entries = entries?;

        mod_table.push(ModTableEntry {
            base: base_mod.modified_base as char,
            is_strand_plus: base_mod.strand == '+',
            mod_code: ModChar::new(base_mod.modification_type),
            implicit: false, // Always false to handle modBAM format requirements
            data: entries,
        });
    }

    Ok(mod_table)
}

/// Reconstruct BaseMods from condensed mod_table format
fn reconstruct_base_mods(
    mod_table: &[ModTableEntry],
    alignment_type: ReadState,
    seq_len: u64,
) -> Result<BaseMods, Error> {
    use fibertools_rs::utils::bamranges::Ranges;
    use fibertools_rs::utils::basemods::BaseMod;

    let mut base_mods = Vec::new();

    for entry in mod_table {
        let mut starts = Vec::new();
        let mut reference_starts = Vec::new();
        let mut qual = Vec::new();

        for &(start, ref_start, q) in &entry.data {
            starts.push(Some(i64::try_from(start)?));
            reference_starts.push(if ref_start == -1 {
                None
            } else {
                Some(ref_start)
            });
            qual.push(q);
        }

        // Calculate ends: starts + 1 where available, None otherwise
        let ends: Vec<Option<i64>> = starts
            .iter()
            .map(|&start_opt| start_opt.map(|start| start + 1))
            .collect();

        // Calculate lengths: Some(1) where starts available, None otherwise
        let lengths: Vec<Option<i64>> = starts
            .iter()
            .map(|&start_opt| start_opt.map(|_| 1))
            .collect();

        // Set reference_ends identical to reference_starts
        let reference_ends = reference_starts.clone();

        // Calculate reference_lengths: Some(0) where reference_starts available, None otherwise
        let reference_lengths: Vec<Option<i64>> = reference_starts
            .iter()
            .map(|&ref_start_opt| ref_start_opt.map(|_| 0))
            .collect();

        let ranges = Ranges {
            starts,
            ends,
            lengths,
            reference_starts,
            reference_ends,
            reference_lengths,
            qual,
            seq_len: seq_len.try_into()?,
            reverse: matches!(
                alignment_type,
                ReadState::PrimaryRev | ReadState::SecondaryRev | ReadState::SupplementaryRev
            ),
        };

        let strand = if entry.is_strand_plus { '+' } else { '-' };

        base_mods.push(BaseMod {
            modified_base: entry.base as u8,
            strand,
            modification_type: entry.mod_code.val(),
            ranges,
            record_is_reverse: false, // Default value
        });
    }

    Ok(BaseMods { base_mods })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ThresholdState, nanalogue_bam_reader};
    use indoc::indoc;
    use rust_htslib::bam::Read;

    // Simple test to verify the module compiles and basic structures work
    #[test]
    fn test_serialized_curr_read_creation() {
        let alignment = AlignmentInfo {
            start: 1000,
            end: 2000,
            contig: "chr1".to_string(),
            contig_id: 0,
        };

        let mod_entry = ModTableEntry {
            base: 'A',
            is_strand_plus: true,
            mod_code: ModChar::new('a'),
            implicit: false,
            data: vec![(10, 1100, 50), (20, 1200, 75)],
        };

        let serialized = SerializedCurrRead {
            alignment_type: ReadState::PrimaryFwd,
            alignment: Some(alignment),
            mod_table: vec![mod_entry],
            read_id: "test_read".to_string(),
            seq_len: 1000,
        };

        // Test JSON serialization
        let json = serde_json::to_string(&serialized).expect("Serialization should work");

        // Test JSON deserialization back to SerializedCurrRead
        let _deserialized: SerializedCurrRead =
            serde_json::from_str(&json).expect("Deserialization should work");
    }

    #[test]
    fn test_first_record_serde() -> Result<(), Error> {
        // Read the first record from the example BAM file
        let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
        let first_record = reader.records().next().unwrap()?;

        // Create CurrRead with alignment and modification data
        let curr_read = CurrRead::default()
            .try_from_only_alignment(&first_record)?
            .set_mod_data(&first_record, ThresholdState::GtEq(0), 0)?;

        let actual_json: serde_json::Value = serde_json::to_value(&curr_read)?;

        let expected_json_str = indoc! {r#"
            {
              "alignment_type": "primary_forward",
              "alignment": {
                "start": 9,
                "end": 17,
                "contig": "dummyI",
                "contig_id": 0
              },
              "mod_table": [
                {
                  "base": "T",
                  "is_strand_plus": true,
                  "mod_code": "T",
                  "implicit": false,
                  "data": [
                    [0, 9, 4],
                    [3, 12, 7],
                    [4, 13, 9],
                    [7, 16, 6]
                  ]
                }
              ],
              "read_id": "5d10eb9a-aae1-4db8-8ec6-7ebb34d32575",
              "seq_len": 8
            }"#};

        let expected_json: serde_json::Value = serde_json::from_str(expected_json_str)?;

        // Compare expected and actual outputs
        assert_eq!(actual_json, expected_json);

        // Also test deserialization: deserialize the expected JSON and compare with original CurrRead
        let deserialized_curr_read: CurrRead<AlignAndModData> =
            serde_json::from_str(expected_json_str)?;
        assert_eq!(deserialized_curr_read, curr_read);

        Ok(())
    }

    #[test]
    fn test_first_record_roundtrip() -> Result<(), Error> {
        // Read the first record from the example BAM file (same as serialization test)
        let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
        let first_record = reader.records().next().unwrap()?;

        // Create the original CurrRead with alignment and modification data
        let original_curr_read = CurrRead::default()
            .try_from_only_alignment(&first_record)?
            .set_mod_data(&first_record, ThresholdState::GtEq(0), 0)?;

        // Serialize to JSON
        let json_str = serde_json::to_string_pretty(&original_curr_read)?;

        // Deserialize back to CurrRead
        let deserialized_curr_read: CurrRead<AlignAndModData> = serde_json::from_str(&json_str)?;

        // The deserialized CurrRead should be equal to the original
        assert_eq!(deserialized_curr_read, original_curr_read);

        Ok(())
    }

    #[test]
    fn test_blank_json_record_roundtrip() -> Result<(), Error> {
        let json_str = indoc! {r#"
            {
            }"#};

        // Deserialize JSON to CurrRead
        let curr_read: CurrRead<AlignAndModData> = serde_json::from_str(json_str)?;

        // Serialize back to JSON
        let serialized_json = serde_json::to_string_pretty(&curr_read)?;

        // Deserialize again
        let roundtrip_curr_read: CurrRead<AlignAndModData> =
            serde_json::from_str(&serialized_json)?;

        // Check that the roundtrip preserves equality
        assert_eq!(curr_read, roundtrip_curr_read);

        Ok(())
    }
}
