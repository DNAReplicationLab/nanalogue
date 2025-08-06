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
use std::ops::RangeInclusive;
use std::rc::Rc;

// Import from our crate
use crate::Error;
use crate::F32Bw0and1;
use crate::ModChar;
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
    /// Marked as unmapped in the BAM file. We are assuming
    /// that unmapped sequences will not be stored as reversed
    /// complements, as what would be the point of that?
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

/// Types of thresholds on modification level that can be applied to modification data.
/// Two possible use cases: (1) to specify that reading mod data should be restricted
/// to bases at least this level of modified, or (2) to specify that only bases
/// in this range should be regarded as modified.
/// Values are 0 to 255 below as that's how they are stored in a modBAM file and
/// this struct is expected to be used in contexts dealing directly with this data.
#[derive(Debug, Clone, PartialEq)]
pub enum ThresholdState {
    /// modification probability >= this value, values are 0 to 255
    GtEq(u8),
    /// modification probability <= this value, values are 0 to 255
    LtEq(u8),
    /// modification probability within this range, values are 0 to 255
    GtEqLtEq(RangeInclusive<u8>),
}

/// default threshold is >= 0 i.e. all mods are allowed
impl Default for ThresholdState {
    fn default() -> Self {
        ThresholdState::GtEq(0)
    }
}

impl fmt::Display for ThresholdState {
    /// display the u8 thresholds as a floating point number between 0 and 1
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let printable = match &self {
            ThresholdState::GtEq(v) => format!("probabilities >= {}", F32Bw0and1::from(*v)),
            ThresholdState::LtEq(v) => format!("probabilities <= {}", F32Bw0and1::from(*v)),
            ThresholdState::GtEqLtEq(v) => {
                format!(
                    "{} <= probabilities <= {}",
                    F32Bw0and1::from(*v.start()),
                    F32Bw0and1::from(*v.end())
                )
            }
        };
        write!(f, "{printable}")
    }
}

impl From<ThresholdState> for RangeInclusive<u8> {
    // convert threshold state into an inclusive range
    fn from(value: ThresholdState) -> Self {
        match value {
            ThresholdState::GtEq(v) => v..=u8::MAX,
            ThresholdState::LtEq(v) => 0..=v,
            ThresholdState::GtEqLtEq(w) => w,
        }
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
/// modification data is sorted by position. We can guarantee
/// this when the modification data is parsed, but we cannot
/// guarantee this if we allow free access to the struct,
/// as a user can mess with this. To prevent these kinds
/// of problems, all internals of this struct are private.
/// NOTE: we could have implemented these as a trait extension
/// to the rust htslib Record struct, but we have chosen not to,
/// as we may want to persist data like modifications and do
/// multiple operations on them. And Record has inconvenient
/// return types like i64 instead of u64 for positions along the
/// genome.
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
    /// NOTE: the contig here is numeric and refers to an index on the BAM
    /// header. To convert this into an alphanumeric string, you have to
    /// process the header and store it in `contig_name` below.
    /// We have left it this way as it is easier to store and process integers.
    contig_and_start: Option<(i32, u64)>,

    /// Contig name.
    contig_name: Option<String>,
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
    pub fn read_state(&self) -> ReadState {
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
    /// gets length of sequence
    pub fn seq_len(&self) -> Result<Option<u64>, Error> {
        match self.read_state() {
            ReadState::Unknown => Err(Error::UnknownAlignState),
            _ => Ok(self.seq_len),
        }
    }
    /// set alignment length from BAM record if available
    pub fn set_align_len(&mut self, record: &Record) -> Result<Option<u64>, Error> {
        match &self.align_len {
            Some(_) => Err(Error::InvalidDuplicates(
                "cannot set alignment length again!".to_string(),
            )),
            None => match self.read_state() {
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
    /// gets alignment length
    pub fn align_len(&self) -> Result<Option<u64>, Error> {
        match self.read_state() {
            ReadState::Unknown => Err(Error::UnknownAlignState),
            ReadState::Unmapped => Ok(None),
            _ => Ok(self.align_len),
        }
    }
    /// sets contig and start from BAM record if available
    pub fn set_contig_and_start(&mut self, record: &Record) -> Result<bool, Error> {
        match &self.contig_and_start {
            Some(_) => Err(Error::InvalidDuplicates(
                "cannot set contig and start again!".to_string(),
            )),
            None => match self.read_state() {
                ReadState::Unknown => Err(Error::UnknownAlignState),
                ReadState::Unmapped => Ok(false),
                _ => {
                    self.contig_and_start = Some((record.tid(), record.pos().try_into()?));
                    Ok(true)
                }
            },
        }
    }
    /// gets contig and start
    pub fn contig_and_start(&self) -> Result<Option<(i32, u64)>, Error> {
        match self.read_state() {
            ReadState::Unknown => Err(Error::UnknownAlignState),
            ReadState::Unmapped => Ok(None),
            _ => Ok(self.contig_and_start),
        }
    }
    /// sets contig name
    pub fn set_contig_name(&mut self, names: &Vec<String>) -> Result<bool, Error> {
        match &self.contig_name {
            Some(_) => Err(Error::InvalidDuplicates(
                "cannot set contig name again!".to_string(),
            )),
            None => match self.contig_and_start() {
                Err(v) => Err(v),
                Ok(None) => Err(Error::InvalidState(
                    "no contig in current read!".to_string(),
                )),
                Ok(Some((v, _))) if 0 <= v && v < names.len() as i32 => {
                    self.contig_name = Some(names[v as usize].clone());
                    Ok(true)
                }
                Ok(Some((_, _))) => Err(Error::InvalidState(
                    "could not assign contig name!".to_string(),
                )),
            },
        }
    }
    /// gets contig name
    pub fn contig_name(&self) -> Result<&str, Error> {
        match (self.read_state(), &self.contig_name) {
            (ReadState::Unknown | ReadState::Unmapped, _) => Err(Error::InvalidState(
                "cannot get contig name from unknown or unmapped read".to_string(),
            )),
            (_, None) => Err(Error::InvalidState("contig name not available".to_string())),
            (_, Some(v)) => Ok(v.as_str()),
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
                    self.read_id()
                }
                Err(_) => Err(Error::InvalidReadID),
            },
        }
    }
    /// gets read id
    pub fn read_id(&self) -> Result<&str, Error> {
        match &self.read_id {
            None => Err(Error::InvalidState("read id not available".to_string())),
            Some(v) => Ok(v.as_str()),
        }
    }
    /// sets modification data using the BAM record
    pub fn set_mod_data(&mut self, record: &Record, mod_thres: ThresholdState) {
        self.mods = Some((
            nanalogue_mm_ml_parser(record, RangeInclusive::from(mod_thres.clone()), None),
            mod_thres,
        ));
    }
    /// sets modification data using BAM record but restricted to the specified
    /// type of modification
    pub fn set_mod_data_one_tag(
        &mut self,
        record: &Record,
        mod_thres: ThresholdState,
        mod_tag: ModChar,
    ) {
        self.mods = Some((
            nanalogue_mm_ml_parser(
                record,
                RangeInclusive::from(mod_thres.clone()),
                Some(mod_tag),
            ),
            mod_thres,
        ));
    }
    /// gets modification data
    pub fn mod_data(&self) -> Result<&Option<(BaseMods, ThresholdState)>, Error> {
        match self.state {
            ReadState::Unknown => Err(Error::UnknownAlignState),
            _ => Ok(&self.mods),
        }
    }
    /// window modification data
    pub fn windowed_mod_data<F>(
        &self,
        window_function: &F,
        win_size: usize,
        slide_size: usize,
        tag: ModChar,
    ) -> Result<Option<Vec<F32Bw0and1>>, Error>
    where
        F: Fn(&[u8]) -> Result<F32Bw0and1, Error>,
    {
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
                            .map(|i| window_function(&data[i..i + win_size]))
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
    fn filter_by_ref_pos(&mut self, start_end: RangeInclusive<i64>) -> Result<bool, Error> {
        macro_rules! subset {
            ( $to_be_subset: expr, $vec_indices:expr ) => {
                $to_be_subset = $vec_indices.iter().map(|&i| $to_be_subset[i]).collect();
            };
        }
        let start = *start_end.start();
        let end = *start_end.end();
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

        if let Some(v) = &self.read_id().ok() {
            output_string = output_string + "\t\"read_id\": \"" + v + "\",\n";
        }

        if let Some(Some(v)) = self.seq_len().ok() {
            output_string = output_string + "\t\"sequence_length\": " + &v.to_string() + ",\n";
        }

        if let Some(Some(v)) = self.align_len().ok() {
            output_string = output_string + "\t\"alignment_length\": " + &v.to_string() + ",\n";
        }

        if let Some(Some((v, w))) = self.contig_and_start().ok() {
            let num_str = &v.to_string();
            output_string = output_string
                + "\t\"contig\": \""
                + if let Some(x) = &self.contig_name().ok() {
                    x
                } else {
                    num_str
                }
                + "\",\n";
            output_string = output_string + "\t\"reference_start\": " + &w.to_string() + ",\n";
        }

        output_string =
            output_string + "\t\"alignment_type\": \"" + &self.state.to_string() + "\",\n";

        let mut mod_count_str = String::from("");
        if let Some((BaseMods { base_mods: v }, w)) = &self.mods {
            if !v.is_empty() {
                for k in v {
                    mod_count_str += format!(
                        "{}{}{}:{};",
                        k.modified_base as char,
                        k.strand,
                        ModChar::new(k.modification_type),
                        k.ranges.qual.len()
                    )
                    .as_str();
                }
                mod_count_str.pop();
            } else {
                mod_count_str += "0";
            }
            mod_count_str += format!(";({})", w).as_str();
        } else {
            mod_count_str += "NA";
        }
        output_string += format!("\t\"mod_count\": \"{}\"\n", mod_count_str).as_str();

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

/// Convert a rust htslib record to our CurrRead struct.
/// NOTE: This operation loads many types of data from the
/// record and you may not be interested in all of them.
/// So, unless you know for sure that you are dealing with
/// a small number of reads, please do not use this function,
/// and call only a subset of the individual invocations below
/// in your program for the sake of speed and/or memory.
impl TryFrom<Record> for CurrRead {
    type Error = crate::Error;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        let mut curr_read_state = CurrRead::default();
        curr_read_state.set_read_state(&record)?;
        curr_read_state.set_read_id(&record)?;
        curr_read_state.set_seq_len(&record)?;
        curr_read_state.set_align_len(&record)?;
        curr_read_state.set_mod_data(&record, ThresholdState::GtEq(128));
        curr_read_state.set_contig_and_start(&record)?;
        Ok(curr_read_state)
    }
}

/// Convert a rust htslib rc record into our struct.
/// I think the rc datatype is just like the normal record,
/// except the record datatype is not destroyed and created
/// every time a new record is read (or something like that).
/// All comments I've made for the `TryFrom<Record>` function
/// apply here as well.
impl TryFrom<Rc<Record>> for CurrRead {
    type Error = crate::Error;

    fn try_from(record: Rc<Record>) -> Result<Self, Self::Error> {
        let mut curr_read_state = CurrRead::default();
        curr_read_state.set_read_state(&record)?;
        curr_read_state.set_read_id(&record)?;
        curr_read_state.set_seq_len(&record)?;
        curr_read_state.set_align_len(&record)?;
        curr_read_state.set_mod_data(&record, ThresholdState::GtEq(128));
        curr_read_state.set_contig_and_start(&record)?;
        Ok(curr_read_state)
    }
}
