//! # ReadUtils
//!
//! Implements CurrRead Struct for processing information retrieved from BAM files
//! and the mod information in the BAM file using a parser implemented in
//! another module.

use bedrs::prelude::StrandedBed3;
use bedrs::{Coordinates, Strand};
use fibertools_rs::utils::basemods::{BaseMod, BaseMods};
use rust_htslib::{bam::ext::BamRecordExtensions, bam::record::Record};
use std::collections::HashMap;
use std::fmt;
use std::num::NonZeroU64;
use std::rc::Rc;

// Import from our crate
use crate::Error;
use crate::F32Bw0and1;
use crate::ModChar;
use crate::OrdPair;
use crate::nanalogue_mm_ml_parser;

/// Alignment state of a read; seven possibilities + one unknown state
#[derive(Debug, Clone, Default, Copy, PartialEq)]
pub enum ReadState {
    /// Unknown alignment
    #[default]
    Unknown,
    /// Primary alignment to the reference strand
    PrimaryFwd,
    /// Primary alignment opposite the reference strand
    PrimaryRev,
    /// Secondary alignment to the reference strand
    SecondaryFwd,
    /// Secondary alignment opposite the reference strand
    SecondaryRev,
    /// Supplementary alignment to the reference strand
    SupplementaryFwd,
    /// Supplementary alignment opposite the reference strand
    SupplementaryRev,
    /// Marked as unmapped in the BAM file
    Unmapped,
}

/// Implements printing of read state
impl fmt::Display for ReadState {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let printable = match *self {
            ReadState::Unknown => "unknown",
            ReadState::PrimaryFwd => "primary_forward",
            ReadState::SecondaryFwd => "secondary_forward",
            ReadState::SupplementaryFwd => "supplementary_forward",
            ReadState::PrimaryRev => "primary_reverse",
            ReadState::SecondaryRev => "secondary_reverse",
            ReadState::SupplementaryRev => "supplementary_reverse",
            ReadState::Unmapped => "unmapped",
        };
        write!(f, "{printable}")
    }
}

/// Types of thresholds on modification level that can be applied to modification data
#[derive(Debug, Clone, PartialEq)]
pub enum ThresholdState {
    /// modification probability >= this value
    GtEq(u8),
    /// modification probability <= this value
    LtEq(u8),
    /// modification probability >= low and <= high of this ordered pair
    GtEqLtEq(OrdPair<u8>),
}

/// default threshold is >= 0 i.e. all mods are allowed
impl Default for ThresholdState {
    fn default() -> Self {
        ThresholdState::GtEq(0)
    }
}

/// Our main struct that receives and stores from one BAM record.
/// Also has methods for processing this information.
/// The information within the struct can only be manipulated by
/// the methods defined here. This is to ensure the struct
/// doesn't fall into an invalid state, which could cause mistakes
/// in calculations associated with the struct. For example:
/// if I want to measure mean modification density along windows
/// of the raw modification data, I need a guarantee that the
/// modification data is stored in a sorted way. We can guarantee
/// this when the modification data is parsed, but we cannot
/// guarantee this if we allow free access to the internals
/// of the struct, as a user can mess with this, making results
/// of the windowing function invalid. To prevent these kinds
/// of problems, all internals of this struct are private.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct CurrRead {
    /// Stores alignment type
    state: ReadState,

    /// Read ID of molecule, also called query name in some contexts
    read_id: Option<String>,

    /// Length of the stored sequence. This is usually the basecalled
    /// sequence but is not guaranteed to be so.
    seq_len: Option<u64>,

    /// Length of the segment on the reference genome the molecule maps to
    align_len: Option<u64>,

    /// Stores modification information along with any applied thresholds
    mods: Option<(BaseMods, ThresholdState)>,

    /// Name of the reference genome contig and the starting position on
    /// contig that the molecule maps or aligns to.
    contig_and_start: Option<(i32, u64)>,
}

impl CurrRead {
    /// sets the alignment of the read using BAM record
    pub fn set_read_state(&mut self, record: &Record) -> Result<bool, Error> {
        // set read state
        match &self.state {
            ReadState::Unknown => {
                match (
                    record.is_reverse(),
                    record.is_unmapped(),
                    record.is_secondary(),
                    record.is_supplementary(),
                ) {
                    (false, true, false, false) => {
                        self.state = ReadState::Unmapped;
                        Ok(true)
                    }
                    (true, false, false, false) => {
                        self.state = ReadState::PrimaryRev;
                        Ok(true)
                    }
                    (true, false, true, false) => {
                        self.state = ReadState::SecondaryRev;
                        Ok(true)
                    }
                    (false, false, true, false) => {
                        self.state = ReadState::SecondaryFwd;
                        Ok(true)
                    }
                    (true, false, false, true) => {
                        self.state = ReadState::SupplementaryRev;
                        Ok(true)
                    }
                    (false, false, false, true) => {
                        self.state = ReadState::SupplementaryFwd;
                        Ok(true)
                    }
                    (false, false, false, false) => {
                        self.state = ReadState::PrimaryFwd;
                        Ok(true)
                    }
                    _ => Err(Error::UnknownAlignState),
                }
            }
            _ => Err(Error::InvalidDuplicates(
                "cannot set align state again!".to_string(),
            )),
        }
    }
    /// gets the read state
    #[must_use]
    pub fn get_read_state(&self) -> ReadState {
        self.state
    }
    /// set length of sequence from BAM record
    pub fn set_seq_len(&mut self, record: &Record) -> Result<Option<u64>, Error> {
        match &self.seq_len {
            Some(_) => Err(Error::InvalidDuplicates(
                "cannot set sequence length again!".to_string(),
            )),
            None => {
                self.seq_len = Some(
                    NonZeroU64::new(record.seq_len().try_into()?)
                        .ok_or(Error::InvalidSeqLength)?
                        .get(),
                );
                Ok(self.seq_len)
            }
        }
    }
    /// set alignment length from BAM record if available
    pub fn set_align_len(&mut self, record: &Record) -> Result<Option<u64>, Error> {
        match &self.align_len {
            Some(_) => Err(Error::InvalidDuplicates(
                "cannot set alignment length again!".to_string(),
            )),
            None => match self.state {
                ReadState::Unknown => Err(Error::UnknownAlignState),
                ReadState::Unmapped => Ok(None),
                _ => {
                    let st: i64 = record.pos();
                    let en: i64 = record.reference_end();
                    if en > st && st >= 0 {
                        self.align_len = Some((en - st).try_into()?);
                        Ok(self.align_len)
                    } else {
                        Err(Error::InvalidAlignLength)
                    }
                }
            },
        }
    }
    /// sets contig and start from BAM record if available
    pub fn set_contig_and_start(&mut self, record: &Record) -> Result<bool, Error> {
        match &self.contig_and_start {
            Some(_) => Err(Error::InvalidDuplicates(
                "cannot set contig and start again!".to_string(),
            )),
            None => match self.state {
                ReadState::Unknown => Err(Error::UnknownAlignState),
                ReadState::Unmapped => Ok(false),
                _ => {
                    self.contig_and_start = Some((record.tid(), record.pos().try_into()?));
                    Ok(true)
                }
            },
        }
    }
    /// sets read ID (also called query name) from BAM record
    pub fn set_read_id(&mut self, record: &Record) -> Result<&str, Error> {
        match &self.read_id {
            Some(_) => Err(Error::InvalidDuplicates(
                "cannot set read id again!".to_string(),
            )),
            None => match str::from_utf8(record.qname()) {
                Ok(v) => {
                    self.read_id = Some(v.to_string());
                    self.get_read_id()
                }
                Err(_) => Err(Error::InvalidReadID),
            },
        }
    }
    /// gets read id
    #[must_use]
    pub fn get_read_id(&self) -> Result<&str, Error> {
        match &self.read_id {
            None => Err(Error::InvalidState("read id not available".to_string())),
            Some(v) => Ok(v.as_str()),
        }
    }
    /// sets modification data using the BAM record
    pub fn set_mod_data(&mut self, record: &Record, mod_thres: u8) {
        self.mods = Some((
            nanalogue_mm_ml_parser(record, mod_thres, None),
            ThresholdState::GtEq(mod_thres),
        ));
    }
    /// sets modification data using BAM record but restricted to the specified
    /// type of modification
    pub fn set_mod_data_one_tag(&mut self, record: &Record, mod_thres: u8, mod_tag: ModChar) {
        self.mods = Some((
            nanalogue_mm_ml_parser(record, mod_thres, Some(mod_tag)),
            ThresholdState::GtEq(mod_thres),
        ));
    }
    /// window modification data
    #[must_use]
    pub fn windowed_mod_data(
        &self,
        win_size: usize,
        slide_size: usize,
        tag: ModChar,
    ) -> Result<Option<Vec<F32Bw0and1>>, Error> {
        let mut result = Vec::<F32Bw0and1>::new();
        let mut is_track_seen: bool = false;
        let tag_char = tag.get_val();
        if let Some((BaseMods { base_mods: v }, _)) = &self.mods {
            for k in v {
                match k {
                    BaseMod {
                        modified_base: _,
                        strand: _,
                        record_is_reverse: _,
                        modification_type: x,
                        ranges: track,
                    } if *x == tag_char && !is_track_seen => {
                        is_track_seen = true;
                        let data = &track.qual;
                        if win_size > data.len() {
                            continue;
                        }
                        result = (0..=data.len() - win_size)
                            .step_by(slide_size)
                            .map(|i| {
                                let window_slice = &data[i..i + win_size];
                                let sum: f32 = window_slice.iter().map(|&val| val as f32).sum();
                                F32Bw0and1::new(sum / (256.0 * win_size as f32))
                            })
                            .collect::<Result<Vec<F32Bw0and1>, _>>()?;
                    }
                    BaseMod {
                        modified_base: _,
                        strand: _,
                        record_is_reverse: _,
                        modification_type: x,
                        ranges: _,
                    } if *x == tag_char && is_track_seen => {
                        return Err(Error::NotImplementedError(
                            "cannot window on data on multiple tracks".to_string(),
                        ));
                    }
                    _ => {}
                }
            }
        }
        if !result.is_empty() {
            Ok(Some(result))
        } else {
            Ok(None)
        }
    }
    /// Performs a count of number of modified bases per modified type.
    /// Note that this result depends on the type of filtering done
    /// while the struct was created e.g. by modification threshold.
    #[must_use]
    pub fn mod_count_per_mod(&self) -> Option<HashMap<ModChar, u32>> {
        let mut output = HashMap::<ModChar, u32>::new();
        match &self.mods {
            Some((BaseMods { base_mods: v }, _)) => {
                if v.is_empty() {
                    None
                } else {
                    for k in v {
                        let mod_count = k.ranges.qual.len() as u32;
                        output
                            .entry(ModChar::new(k.modification_type))
                            .and_modify(|e| *e += mod_count)
                            .or_insert(mod_count);
                    }
                    Some(output)
                }
            }
            None => None,
        }
    }
    /// filter modification data so that only data corresponding to the given
    /// range of positions is retained. We have copied and adapted code
    /// from the fibertools_rs repository here.
    fn filter_by_ref_pos(&mut self, start_end: OrdPair<i64>) -> Result<bool, Error> {
        macro_rules! subset {
            ( $to_be_subset: expr, $vec_indices:expr ) => {
                $to_be_subset = $vec_indices.iter().map(|&i| $to_be_subset[i]).collect();
            };
        }
        let start = start_end.get_low();
        let end = start_end.get_high();
        match (self.state, &mut self.mods) {
            (ReadState::Unknown, _) => Err(Error::UnknownAlignState),
            (ReadState::Unmapped, _) | (_, None) => Ok(false),
            (_, Some((BaseMods { base_mods: v }, _))) => {
                for k in v {
                    let w = &mut k.ranges;
                    let to_keep = w
                        .reference_starts
                        .iter()
                        .enumerate()
                        .filter_map(|(i, &s)| {
                            if s >= Some(start) && s <= Some(end) {
                                Some(i)
                            } else {
                                None
                            }
                        })
                        .collect::<Vec<_>>();
                    subset!(w.qual, to_keep);
                    subset!(w.starts, to_keep);
                    subset!(w.ends, to_keep);
                    subset!(w.lengths, to_keep);
                    subset!(w.reference_starts, to_keep);
                    subset!(w.reference_ends, to_keep);
                    subset!(w.reference_lengths, to_keep);
                }
                Ok(true)
            }
        }
    }
    /// Uses only alignment information and no modification information to
    /// create the struct. Use this if you want to perform operations that
    /// do not involve reading or manipulating the modification data.
    pub fn try_from_only_alignment(record: Record) -> Result<Self, Error> {
        let mut curr_read_state = CurrRead::default();
        curr_read_state.set_read_state(&record)?;
        curr_read_state.set_align_len(&record)?;
        curr_read_state.set_contig_and_start(&record)?;
        Ok(curr_read_state)
    }
}

impl fmt::Display for CurrRead {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut output_string = String::from("");

        if let Some(v) = &self.read_id {
            output_string = output_string + "\t\"read_id\": \"" + v + "\",\n";
        }

        if let Some(v) = self.seq_len {
            output_string = output_string + "\t\"sequence_length\": " + &v.to_string() + ",\n";
        }

        if let Some(v) = self.align_len {
            output_string = output_string + "\t\"alignment_length\": " + &v.to_string() + ",\n";
        }

        if let Some((v, w)) = &self.contig_and_start {
            output_string = output_string + "\t\"contig\": \"" + &v.to_string() + "\",\n";
            output_string = output_string + "\t\"reference_start\": " + &w.to_string() + ",\n";
        }

        output_string =
            output_string + "\t\"alignment_type\": \"" + &self.state.to_string() + "\",\n";

        output_string += "\t\"mod_count\": ";
        if let Some((BaseMods { base_mods: v }, _)) = &self.mods {
            if !v.is_empty() {
                for k in v {
                    output_string += format!(
                        "\"{}{}{}:{};",
                        k.modified_base as char,
                        k.strand,
                        ModChar::new(k.modification_type),
                        k.ranges.qual.len()
                    )
                    .as_str();
                }
                output_string.pop();
                output_string += "\"\n";
            } else {
                output_string += " \"0\"\n";
            }
        } else {
            output_string += " \"0\"\n";
        }

        write!(f, "{{\n{output_string}}}")
    }
}

impl TryFrom<CurrRead> for StrandedBed3<i32, u64> {
    type Error = crate::Error;

    fn try_from(value: CurrRead) -> Result<Self, Self::Error> {
        match (value.state, value.align_len, value.contig_and_start) {
            (ReadState::Unknown, _, _) => Err(Error::UnknownAlignState),
            (ReadState::Unmapped, _, _) => Ok(StrandedBed3::empty()),
            (_, None, _) => Err(Error::InvalidAlignLength),
            (_, _, None) => Err(Error::InvalidContigAndStart),
            (
                ReadState::PrimaryFwd | ReadState::SecondaryFwd | ReadState::SupplementaryFwd,
                Some(al),
                Some((cg, st)),
            ) => Ok(StrandedBed3::new(cg, st, st + al, Strand::Forward)),
            (
                ReadState::PrimaryRev | ReadState::SecondaryRev | ReadState::SupplementaryRev,
                Some(al),
                Some((cg, st)),
            ) => Ok(StrandedBed3::new(cg, st, st + al, Strand::Reverse)),
        }
    }
}

impl TryFrom<Record> for CurrRead {
    type Error = crate::Error;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        let mut curr_read_state = CurrRead::default();
        curr_read_state.set_read_state(&record)?;
        curr_read_state.set_read_id(&record)?;
        curr_read_state.set_seq_len(&record)?;
        curr_read_state.set_align_len(&record)?;
        curr_read_state.set_mod_data(&record, 128);
        curr_read_state.set_contig_and_start(&record)?;
        Ok(curr_read_state)
    }
}

impl TryFrom<Rc<Record>> for CurrRead {
    type Error = crate::Error;

    fn try_from(record: Rc<Record>) -> Result<Self, Self::Error> {
        let mut curr_read_state = CurrRead::default();
        curr_read_state.set_read_state(&record)?;
        curr_read_state.set_read_id(&record)?;
        curr_read_state.set_seq_len(&record)?;
        curr_read_state.set_align_len(&record)?;
        curr_read_state.set_mod_data(&record, 128);
        curr_read_state.set_contig_and_start(&record)?;
        Ok(curr_read_state)
    }
}
