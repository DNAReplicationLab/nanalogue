//! # ReadUtils
//!
//! Implements CurrRead Struct for processing information retrieved from BAM files
//! and the mod information in the BAM file using a parser implemented in
//! another module.

use bedrs::prelude::StrandedBed3;
use bedrs::{Coordinates, Strand};
use bio_types::genome::AbstractInterval;
use fibertools_rs::utils::bamranges::Ranges;
use fibertools_rs::utils::basemods::{BaseMod, BaseMods};
use rust_htslib::{bam::ext::BamRecordExtensions, bam::record::Record};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;
use std::num::NonZeroU64;
use std::ops::Range;
use std::rc::Rc;
use std::str::FromStr;

// Import from our crate
use crate::{
    Contains, Error, F32Bw0and1, FilterByRefCoords, InputMods, Intersects, ModChar, OrdPair,
    nanalogue_mm_ml_parser,
};

/// Alignment state of a read; seven possibilities + one unknown state
#[derive(Debug, Clone, Default, Copy, PartialEq, Eq, Serialize, Deserialize)]
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

// Implements conversion of ReadState into the standard BAM flag format
impl TryFrom<ReadState> for u16 {
    type Error = Error;
    /// converts our internal representation to the BAM flag format
    fn try_from(value: ReadState) -> Result<u16, Error> {
        match value {
            ReadState::PrimaryFwd => Ok(0),
            ReadState::Unmapped => Ok(4),
            ReadState::PrimaryRev => Ok(16),
            ReadState::SecondaryFwd => Ok(256),
            ReadState::SecondaryRev => Ok(272),
            ReadState::SupplementaryFwd => Ok(2048),
            ReadState::SupplementaryRev => Ok(2064),
            ReadState::Unknown => Err(Error::UnknownAlignState),
        }
    }
}

// Implements conversion of the standard BAM flag format into ReadState
impl TryFrom<u16> for ReadState {
    type Error = Error;
    /// converts our internal representation to the BAM flag format
    fn try_from(value: u16) -> Result<ReadState, Error> {
        match value {
            0 => Ok(ReadState::PrimaryFwd),
            4 => Ok(ReadState::Unmapped),
            16 => Ok(ReadState::PrimaryRev),
            256 => Ok(ReadState::SecondaryFwd),
            272 => Ok(ReadState::SecondaryRev),
            2048 => Ok(ReadState::SupplementaryFwd),
            2064 => Ok(ReadState::SupplementaryRev),
            _ => Err(Error::UnknownAlignState),
        }
    }
}

/// Implements from string for ReadState
impl FromStr for ReadState {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "primary_forward" => Ok(ReadState::PrimaryFwd),
            "primary_reverse" => Ok(ReadState::PrimaryRev),
            "secondary_forward" => Ok(ReadState::SecondaryFwd),
            "secondary_reverse" => Ok(ReadState::SecondaryRev),
            "supplementary_forward" => Ok(ReadState::SupplementaryFwd),
            "supplementary_reverse" => Ok(ReadState::SupplementaryRev),
            "unmapped" => Ok(ReadState::Unmapped),
            _ => Err(Error::UnknownAlignState),
        }
    }
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
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum ThresholdState {
    /// modification probability >= this value, values are 0 to 255
    GtEq(u8),
    /// modification probability <= this value, values are 0 to 255
    LtEq(u8),
    /// modification probability not within this range.
    /// We expect this to be used to filter out modification calls
    /// around 0.5 i.e. ones with the most uncertainty, although
    /// users of this crate are free to set this to an interval
    /// not including 0.5
    InvertGtEqLtEq(OrdPair<u8>),
}

/// default threshold is >= 0 i.e. all mods are allowed
impl Default for ThresholdState {
    fn default() -> Self {
        ThresholdState::GtEq(0)
    }
}

/// Displays thresholds but using floating point numbers between 0 and 1
///
/// Example 1:
/// ```
/// use nanalogue_core::{ThresholdState, OrdPair};
/// let b = ThresholdState::GtEq(100);
/// assert_eq!("probabilities >= 0.390625", format!("{}", b));
/// ```
/// Example 2:
/// ```
/// # use nanalogue_core::ThresholdState;
/// let b = ThresholdState::LtEq(10);
/// assert_eq!("probabilities <= 0.0390625", format!("{}", b));
/// ```
/// Example 3:
/// ```
/// # use nanalogue_core::{ThresholdState, OrdPair};
/// let b = ThresholdState::InvertGtEqLtEq(OrdPair::new(200, 220).expect("no error"));
/// assert_eq!("probabilities < 0.78125 or > 0.859375", format!("{}", b));
/// ```
impl fmt::Display for ThresholdState {
    /// display the u8 thresholds as a floating point number between 0 and 1
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let printable = match &self {
            ThresholdState::GtEq(v) => format!("probabilities >= {}", F32Bw0and1::from(*v)),
            ThresholdState::LtEq(v) => format!("probabilities <= {}", F32Bw0and1::from(*v)),
            ThresholdState::InvertGtEqLtEq(v) => {
                format!(
                    "probabilities < {} or > {}",
                    F32Bw0and1::from(v.get_low()),
                    F32Bw0and1::from(v.get_high())
                )
            }
        };
        write!(f, "{printable}")
    }
}

/// Check if a given u8 is within the interval covered
///
/// Example 1:
/// ```
/// use nanalogue_core::{Error, OrdPair, ThresholdState, Contains};
/// let b = ThresholdState::GtEq(100);
/// assert!(b.contains(&101));
/// assert!(b.contains(&100));
/// assert!(!b.contains(&99));
/// ```
/// Example 2:
/// ```
/// # use nanalogue_core::{Error, OrdPair, ThresholdState, Contains};
/// let b = ThresholdState::LtEq(10);
/// assert!(!b.contains(&11));
/// assert!(b.contains(&10));
/// assert!(b.contains(&9));
/// ```
/// Example 3:
/// ```
/// # use nanalogue_core::{Error, OrdPair, ThresholdState, Contains};
/// let b = ThresholdState::InvertGtEqLtEq(OrdPair::new(200, 220).expect("no error"));
/// assert!(b.contains(&100));
/// assert!(!b.contains(&200));
/// assert!(!b.contains(&210));
/// assert!(!b.contains(&220));
/// assert!(b.contains(&250));
/// ```
impl Contains<u8> for ThresholdState {
    /// see if value is contained within the interval
    /// specified by the threshold state
    fn contains(&self, val: &u8) -> bool {
        match &self {
            ThresholdState::GtEq(v) => val >= v,
            ThresholdState::LtEq(v) => val <= v,
            ThresholdState::InvertGtEqLtEq(w) => !w.contains(val),
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
    /// Stores alignment type.
    state: ReadState,

    /// Read ID of molecule, also called query name in some contexts.
    read_id: Option<String>,

    /// Length of the stored sequence. This is usually the basecalled
    /// sequence but is not guaranteed to be so.
    seq_len: Option<u64>,

    /// Length of the segment on the reference genome the molecule maps to.
    align_len: Option<u64>,

    /// Stores modification information along with any applied thresholds.
    mods: Option<(BaseMods, ThresholdState)>,

    /// ID of the reference genome contig and the starting position on
    /// contig that the molecule maps or aligns to.
    /// NOTE: the contig here is numeric and refers to an index on the BAM
    /// header. To convert this into an alphanumeric string, you have to
    /// process the header and store it in `contig_name` below.
    /// We have left it this way as it is easier to store and process integers.
    contig_id_and_start: Option<(i32, u64)>,

    /// Contig name.
    contig_name: Option<String>,

    /// Base PHRED-quality threshold (no offset). Mods could have been filtered by this.
    mod_base_qual_thres: u8,
}

impl CurrRead {
    /// sets the alignment of the read using BAM record
    ///
    /// ```
    /// use nanalogue_core::{CurrRead, Error, nanalogue_bam_reader, ReadState};
    /// use rust_htslib::bam::Read;
    /// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    /// let mut count = 0;
    /// for record in reader.records(){
    ///     let r = record?;
    ///     let mut curr_read = CurrRead::default();
    ///     curr_read.set_read_state(&r);
    ///     match count {
    ///         0 => assert_eq!(curr_read.read_state(), ReadState::PrimaryFwd),
    ///         1 => assert_eq!(curr_read.read_state(), ReadState::PrimaryFwd),
    ///         2 => assert_eq!(curr_read.read_state(), ReadState::PrimaryRev),
    ///         3 => assert_eq!(curr_read.read_state(), ReadState::Unmapped),
    ///         _ => unreachable!(),
    ///     }
    ///     count = count + 1;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    pub fn set_read_state(&mut self, record: &Record) -> Result<ReadState, Error> {
        if record.is_paired()
            || record.is_proper_pair()
            || record.is_first_in_template()
            || record.is_last_in_template()
            || record.is_mate_reverse()
            || record.is_mate_unmapped()
            || record.is_duplicate()
        {
            return Err(Error::NotImplementedError(
                "paired-read/mate-read formats not implemented or duplicate reads found!"
                    .to_string(),
            ));
        }
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
                        Ok(self.state)
                    }
                    (true, false, false, false) => {
                        self.state = ReadState::PrimaryRev;
                        Ok(self.state)
                    }
                    (true, false, true, false) => {
                        self.state = ReadState::SecondaryRev;
                        Ok(self.state)
                    }
                    (false, false, true, false) => {
                        self.state = ReadState::SecondaryFwd;
                        Ok(self.state)
                    }
                    (true, false, false, true) => {
                        self.state = ReadState::SupplementaryRev;
                        Ok(self.state)
                    }
                    (false, false, false, true) => {
                        self.state = ReadState::SupplementaryFwd;
                        Ok(self.state)
                    }
                    (false, false, false, false) => {
                        self.state = ReadState::PrimaryFwd;
                        Ok(self.state)
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
    /// resets the read state
    pub fn reset(&mut self) {
        self.state = ReadState::Unknown;
        self.read_id = None;
        self.seq_len = None;
        self.align_len = None;
        self.contig_id_and_start = None;
        self.contig_name = None;
        self.reset_mod_data();
    }
    /// resets the mod data
    pub fn reset_mod_data(&mut self) {
        self.mods = None;
        self.mod_base_qual_thres = 0;
    }
    /// set length of sequence from BAM record
    ///
    /// ```
    /// use nanalogue_core::{CurrRead, Error, nanalogue_bam_reader};
    /// use rust_htslib::bam::Read;
    /// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    /// let mut count = 0;
    /// for record in reader.records(){
    ///     let r = record?;
    ///     let mut curr_read = CurrRead::default();
    ///     curr_read.set_read_state(&r)?;
    ///     let Ok(len) = curr_read.set_seq_len(&r) else { unreachable!() };
    ///     let Ok(len2) = curr_read.seq_len() else { unreachable!() };
    ///     assert_eq!(len, len2);
    ///     match count {
    ///         0 => assert_eq!(len, 8),
    ///         1 => assert_eq!(len, 48),
    ///         2 => assert_eq!(len, 33),
    ///         3 => assert_eq!(len, 48),
    ///         _ => unreachable!(),
    ///     }
    ///     count = count + 1;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    ///
    /// If we call the method twice, we should hit a panic
    /// ```should_panic
    /// # use nanalogue_core::{CurrRead, Error, nanalogue_bam_reader};
    /// # use rust_htslib::bam::Read;
    /// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    /// for record in reader.records(){
    ///     let r = record?;
    ///     let mut curr_read = CurrRead::default();
    ///     curr_read.set_read_state(&r);
    ///     curr_read.set_seq_len(&r)?;
    ///     curr_read.set_seq_len(&r)?;
    ///     break;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    pub fn set_seq_len(&mut self, record: &Record) -> Result<u64, Error> {
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
                self.seq_len.ok_or(Error::UnknownError)
            }
        }
    }
    /// gets length of sequence
    pub fn seq_len(&self) -> Result<u64, Error> {
        match self.read_state() {
            ReadState::Unknown => Err(Error::UnknownAlignState),
            _ => self.seq_len.ok_or(Error::UnavailableData),
        }
    }
    /// set alignment length from BAM record if available
    pub fn set_align_len(&mut self, record: &Record) -> Result<u64, Error> {
        match &self.align_len {
            Some(_) => Err(Error::InvalidDuplicates(
                "cannot set alignment length again!".to_string(),
            )),
            None => match self.read_state() {
                ReadState::Unknown => Err(Error::UnknownAlignState),
                ReadState::Unmapped => Err(Error::Unmapped),
                _ => {
                    let st: i64 = record.pos();
                    let en: i64 = record.reference_end();
                    if en > st && st >= 0 {
                        self.align_len = Some((en - st).try_into()?);
                        self.align_len.ok_or(Error::UnknownError)
                    } else {
                        Err(Error::InvalidAlignLength)
                    }
                }
            },
        }
    }
    /// gets alignment length
    pub fn align_len(&self) -> Result<u64, Error> {
        match self.read_state() {
            ReadState::Unknown => Err(Error::UnknownAlignState),
            ReadState::Unmapped => Err(Error::Unmapped),
            _ => self.align_len.ok_or(Error::UnavailableData),
        }
    }
    /// sets contig ID and start from BAM record if available
    ///
    /// ```
    /// use nanalogue_core::{CurrRead, Error, nanalogue_bam_reader};
    /// use rust_htslib::bam::Read;
    /// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    /// let mut count = 0;
    /// for record in reader.records(){
    ///     let r = record?;
    ///     let mut curr_read = CurrRead::default();
    ///     curr_read.set_read_state(&r)?;
    ///     match (count, curr_read.set_contig_id_and_start(&r)) {
    ///         (0, Ok((0, 9))) |
    ///         (1, Ok((2, 23))) |
    ///         (2, Ok((1, 3))) |
    ///         (3, Err(Error::Unmapped)) => {},
    ///         _ => unreachable!(),
    ///     }
    ///     count = count + 1;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    ///
    /// If we call the method twice, we should hit a panic
    /// ```should_panic
    /// # use nanalogue_core::{CurrRead, Error, nanalogue_bam_reader};
    /// # use rust_htslib::bam::Read;
    /// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    /// for record in reader.records(){
    ///     let r = record?;
    ///     let mut curr_read = CurrRead::default();
    ///     curr_read.set_read_state(&r)?;
    ///     curr_read.set_contig_id_and_start(&r)?;
    ///     curr_read.set_contig_id_and_start(&r)?;
    ///     break;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    pub fn set_contig_id_and_start(&mut self, record: &Record) -> Result<(i32, u64), Error> {
        match &self.contig_id_and_start {
            Some(_) => Err(Error::InvalidDuplicates(
                "cannot set contig and start again!".to_string(),
            )),
            None => match self.read_state() {
                ReadState::Unknown => Err(Error::UnknownAlignState),
                ReadState::Unmapped => Err(Error::Unmapped),
                _ => {
                    self.contig_id_and_start = Some((record.tid(), record.pos().try_into()?));
                    self.contig_id_and_start.ok_or(Error::UnknownError)
                }
            },
        }
    }
    /// gets contig ID and start
    pub fn contig_id_and_start(&self) -> Result<(i32, u64), Error> {
        match self.read_state() {
            ReadState::Unknown => Err(Error::UnknownAlignState),
            ReadState::Unmapped => Err(Error::Unmapped),
            _ => self.contig_id_and_start.ok_or(Error::UnavailableData),
        }
    }
    /// sets contig name
    ///
    /// ```
    /// use nanalogue_core::{CurrRead, Error, nanalogue_bam_reader};
    /// use rust_htslib::bam::Read;
    /// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    /// let mut count = 0;
    /// for record in reader.records(){
    ///     let r = record?;
    ///     let mut curr_read = CurrRead::default();
    ///     curr_read.set_read_state(&r)?;
    ///     curr_read.set_contig_name(&r)?;
    ///     let Ok(contig_name) = curr_read.contig_name() else {unreachable!()};
    ///     match (count, contig_name) {
    ///         (0, "dummyI") |
    ///         (1, "dummyIII") |
    ///         (2, "dummyII") => {},
    ///         _ => unreachable!(),
    ///     }
    ///     count = count + 1;
    ///     if count == 3 { break; } // the fourth entry is unmapped, and will lead to an error.
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    ///
    /// If we try to set contig name on an unmapped read, we will get an error
    ///
    /// ```should_panic
    /// # use nanalogue_core::{CurrRead, Error, nanalogue_bam_reader};
    /// # use rust_htslib::bam::Read;
    /// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    /// let mut count = 0;
    /// for record in reader.records(){
    ///     if count < 3 {
    ///         count = count + 1;
    ///         continue;
    ///     }
    ///     let r = record?;
    ///     let mut curr_read = CurrRead::default();
    ///     curr_read.set_read_state(&r)?;
    ///     curr_read.set_contig_name(&r)?;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    ///
    /// If we call the method twice, we should hit a panic
    /// ```should_panic
    /// # use nanalogue_core::{CurrRead, Error, nanalogue_bam_reader};
    /// # use rust_htslib::bam::Read;
    /// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    /// for record in reader.records(){
    ///     let r = record?;
    ///     let mut curr_read = CurrRead::default();
    ///     curr_read.set_read_state(&r)?;
    ///     curr_read.set_contig_name(&r)?;
    ///     curr_read.set_contig_name(&r)?;
    ///     break;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    pub fn set_contig_name(&mut self, record: &Record) -> Result<&str, Error> {
        match (self.read_state(), &self.contig_name) {
            (ReadState::Unknown, _) => Err(Error::UnknownAlignState),
            (ReadState::Unmapped, _) => Err(Error::Unmapped),
            (_, Some(_)) => Err(Error::InvalidDuplicates(
                "cannot set contig name again!".to_string(),
            )),
            (_, None) => {
                self.contig_name = Some(String::from(record.contig()));
                match &self.contig_name {
                    None => unreachable!(),
                    Some(v) => Ok(v.as_str()),
                }
            }
        }
    }
    /// gets contig name
    pub fn contig_name(&self) -> Result<&str, Error> {
        match (self.read_state(), &self.contig_name) {
            (ReadState::Unknown, _) => Err(Error::UnknownAlignState),
            (ReadState::Unmapped, _) => Err(Error::Unmapped),
            (_, None) => Err(Error::UnavailableData),
            (_, Some(v)) => Ok(v.as_str()),
        }
    }
    /// sets read ID (also called query name) from BAM record
    ///
    /// ```
    /// use nanalogue_core::{CurrRead, Error, nanalogue_bam_reader};
    /// use rust_htslib::bam::Read;
    /// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    /// let mut count = 0;
    /// for record in reader.records(){
    ///     let r = record?;
    ///     let mut curr_read = CurrRead::default();
    ///     curr_read.set_read_state(&r)?;
    ///
    ///     let Ok(read_id) = curr_read.set_read_id(&r) else { unreachable!() };
    ///     let read_id_clone = String::from(read_id);
    ///     drop(read_id);
    ///     let Ok(read_id) = curr_read.read_id() else { unreachable!() };
    ///     assert_eq!(read_id_clone, read_id);
    ///     drop(read_id_clone);
    ///
    ///     match (count, read_id) {
    ///         (0,"5d10eb9a-aae1-4db8-8ec6-7ebb34d32575") |
    ///         (1,"a4f36092-b4d5-47a9-813e-c22c3b477a0c") |
    ///         (2,"fffffff1-10d2-49cb-8ca3-e8d48979001b") |
    ///         (3,"a4f36092-b4d5-47a9-813e-c22c3b477a0c") => {},
    ///         _ => unreachable!(),
    ///     }
    ///     count = count + 1;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    ///
    /// If we call the method twice, we should hit a panic
    /// ```should_panic
    /// # use nanalogue_core::{CurrRead, Error, nanalogue_bam_reader};
    /// # use rust_htslib::bam::Read;
    /// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    /// for record in reader.records(){
    ///     let r = record?;
    ///     let mut curr_read = CurrRead::default();
    ///     curr_read.set_read_state(&r)?;
    ///     curr_read.set_read_id(&r)?;
    ///     curr_read.set_read_id(&r)?;
    ///     break;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
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
            None => Err(Error::UnavailableData),
            Some(v) => Ok(v.as_str()),
        }
    }
    /// sets modification data using the BAM record
    pub fn set_mod_data(
        &mut self,
        record: &Record,
        mod_thres: ThresholdState,
        min_qual: u8,
    ) -> Result<&(BaseMods, ThresholdState), Error> {
        let result = nanalogue_mm_ml_parser(
            record,
            |x| mod_thres.contains(x),
            |_| true,
            |_, _, _| true,
            min_qual,
        )?;
        self.mods = Some((result, mod_thres));
        self.mod_base_qual_thres = min_qual;
        self.mod_data()
    }
    /// sets modification data using BAM record but restricted to the
    /// specified filters
    pub fn set_mod_data_restricted<G, H>(
        &mut self,
        record: &Record,
        mod_thres: ThresholdState,
        mod_fwd_pos_filter: G,
        mod_filter_base_strand_tag: H,
        min_qual: u8,
    ) -> Result<&(BaseMods, ThresholdState), Error>
    where
        G: Fn(&usize) -> bool,
        H: Fn(&u8, &char, &ModChar) -> bool,
    {
        let result = nanalogue_mm_ml_parser(
            record,
            |x| mod_thres.contains(x),
            mod_fwd_pos_filter,
            mod_filter_base_strand_tag,
            min_qual,
        )?;
        self.mods = Some((result, mod_thres));
        self.mod_base_qual_thres = min_qual;
        self.mod_data()
    }
    /// sets modification data using BAM record but with restrictions
    /// applied by the InputMods options
    pub fn set_mod_data_restricted_options(
        &mut self,
        record: &Record,
        mod_options: &InputMods,
    ) -> Result<&(BaseMods, ThresholdState), Error> {
        match (
            &mod_options.trim_read_ends,
            self.seq_len(),
            &mod_options.mod_strand,
        ) {
            (0, _, &Some(v)) => self.set_mod_data_restricted(
                record,
                mod_options.mod_prob_filter,
                |_| true,
                |_, &s, &t| t == mod_options.tag && s == char::from(v),
                mod_options.base_qual_filter,
            ),
            (0, _, None) => self.set_mod_data_restricted(
                record,
                mod_options.mod_prob_filter,
                |_| true,
                |_, _, &t| t == mod_options.tag,
                mod_options.base_qual_filter,
            ),
            (&w, Ok(l), &Some(v)) => self.set_mod_data_restricted(
                record,
                mod_options.mod_prob_filter,
                |x| {
                    (w..usize::try_from(l)
                        .expect("bit conversion error")
                        .checked_sub(w)
                        .unwrap_or_default())
                        .contains(x)
                },
                |_, &s, &t| t == mod_options.tag && s == char::from(v),
                mod_options.base_qual_filter,
            ),
            (&w, Ok(l), None) => self.set_mod_data_restricted(
                record,
                mod_options.mod_prob_filter,
                |x| {
                    (w..usize::try_from(l)
                        .expect("bit conversion error")
                        .checked_sub(w)
                        .unwrap_or_default())
                        .contains(x)
                },
                |_, _, &t| t == mod_options.tag,
                mod_options.base_qual_filter,
            ),
            (_, Err(_), _) => Err(Error::InvalidSeqLength),
        }?;
        self.mod_data()
    }
    /// gets modification data
    pub fn mod_data(&self) -> Result<&(BaseMods, ThresholdState), Error> {
        match self.state {
            ReadState::Unknown => Err(Error::UnknownAlignState),
            _ => match &self.mods {
                Some(v) => Ok(v),
                None => Err(Error::UnavailableData),
            },
        }
    }
    /// window modification data with restrictions.
    /// If a read has the same modification on both the basecalled
    /// strand and its complement, then windows along both are returned.
    /// We make no guarantee about the ordering of the coordinates of the output in this case.
    pub fn windowed_mod_data_restricted<F>(
        &self,
        window_function: &F,
        win_size: usize,
        slide_size: usize,
        tag: ModChar,
    ) -> Result<Vec<F32Bw0and1>, Error>
    where
        F: Fn(&[u8]) -> Result<F32Bw0and1, Error>,
    {
        let mut result = Vec::<F32Bw0and1>::new();
        let mut plus_mod_strand_seen = false;
        let mut minus_mod_strand_seen = false;
        let tag_char = tag.val();
        if let Some((BaseMods { base_mods: v }, _)) = &self.mods {
            for k in v {
                match k {
                    BaseMod {
                        modified_base: _,
                        strand: '+',
                        record_is_reverse: _,
                        modification_type: x,
                        ranges: _,
                    } if *x == tag_char && plus_mod_strand_seen => {
                        return Err(Error::InvalidDuplicates(
                            "mod type has multiple plus tracks!".to_string(),
                        ));
                    }
                    BaseMod {
                        modified_base: _,
                        strand: '-',
                        record_is_reverse: _,
                        modification_type: x,
                        ranges: _,
                    } if *x == tag_char && minus_mod_strand_seen => {
                        return Err(Error::InvalidDuplicates(
                            "mod type has multiple minus tracks!".to_string(),
                        ));
                    }
                    BaseMod {
                        modified_base: _,
                        strand: s,
                        record_is_reverse: _,
                        modification_type: x,
                        ranges: track,
                    } if *x == tag_char => {
                        match s {
                            '+' => plus_mod_strand_seen = true,
                            '-' => minus_mod_strand_seen = true,
                            _ => return Err(Error::InvalidModType),
                        }
                        let mod_data = &track.qual;
                        if win_size > mod_data.len() {
                            continue;
                        }
                        result.extend(
                            (0..=mod_data.len() - win_size)
                                .step_by(slide_size)
                                .map(|i| window_function(&mod_data[i..i + win_size]))
                                .collect::<Result<Vec<F32Bw0and1>, _>>()?,
                        );
                    }
                    _ => {}
                }
            }
        }
        Ok(result)
    }
    /// Performs a count of number of bases per modified type.
    /// Note that this result depends on the type of filtering done
    /// while the struct was created e.g. by modification threshold.
    ///
    /// ```
    /// use nanalogue_core::{CurrRead, Error, ModChar, nanalogue_bam_reader, ThresholdState};
    /// use rust_htslib::bam::Read;
    /// use std::collections::HashMap;
    /// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    /// let mut count = 0;
    /// for record in reader.records(){
    ///     let r = record?;
    ///     let mut curr_read = CurrRead::default();
    ///     curr_read.set_read_state(&r)?;
    ///     curr_read.set_mod_data(&r, ThresholdState::GtEq(180), 0);
    ///
    ///     let mod_count = curr_read.base_count_per_mod();
    ///     let zero_count = Some(HashMap::from([(ModChar::new('T'), 0)]));
    ///     let a = Some(HashMap::from([(ModChar::new('T'), 3)]));
    ///     let b = Some(HashMap::from([(ModChar::new('T'), 1)]));
    ///     let c = Some(HashMap::from([(ModChar::new('T'), 3),(ModChar::new('ᰠ'), 0)]));
    ///     match (count, mod_count) {
    ///         (0, v) => assert_eq!(v, zero_count),
    ///         (1, v) => assert_eq!(v, a),
    ///         (2, v) => assert_eq!(v, b),
    ///         (3, v) => assert_eq!(v, c),
    ///         _ => unreachable!(),
    ///     }
    ///     count = count + 1;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    #[must_use]
    pub fn base_count_per_mod(&self) -> Option<HashMap<ModChar, u32>> {
        let mut output = HashMap::<ModChar, u32>::new();
        match self.mod_data().ok() {
            Some((BaseMods { base_mods: v }, _)) => {
                for k in v {
                    let base_count = k.ranges.qual.len() as u32;
                    output
                        .entry(ModChar::new(k.modification_type))
                        .and_modify(|e| *e += base_count)
                        .or_insert(base_count);
                }
                Some(output)
            }
            None => None,
        }
    }
    /// Uses only alignment information and no modification information to
    /// create the struct. Use this if you want to perform operations that
    /// do not involve reading or manipulating the modification data.
    pub fn try_from_only_alignment(record: &Record) -> Result<Self, Error> {
        let mut curr_read_state = CurrRead::default();
        let read_state: ReadState = curr_read_state.set_read_state(record)?;
        curr_read_state.set_read_id(record)?;
        curr_read_state.set_seq_len(record)?;
        match read_state {
            ReadState::Unknown => unreachable!(),
            ReadState::Unmapped => {}
            _ => {
                curr_read_state.set_align_len(record)?;
                curr_read_state.set_contig_id_and_start(record)?;
                curr_read_state.set_contig_name(record)?;
            }
        }
        Ok(curr_read_state)
    }
    /// Returns the character corresponding to the strand
    ///
    /// ```
    /// use nanalogue_core::{CurrRead, Error, nanalogue_bam_reader};
    /// use rust_htslib::bam::Read;
    /// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    /// let mut count = 0;
    /// for record in reader.records(){
    ///     let r = record?;
    ///     let mut curr_read = CurrRead::default();
    ///     curr_read.set_read_state(&r)?;
    ///     let Ok(strand) = curr_read.strand() else { unreachable!() };
    ///     match (count, strand) {
    ///         (0, '+') | (1, '+') | (2, '-') | (3, '.') => {},
    ///         _ => unreachable!(),
    ///     }
    ///     count = count + 1;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    pub fn strand(&self) -> Result<char, Error> {
        match &self.state {
            ReadState::Unknown => Err(Error::UnknownAlignState),
            ReadState::Unmapped => Ok('.'),
            ReadState::PrimaryFwd | ReadState::SecondaryFwd | ReadState::SupplementaryFwd => {
                Ok('+')
            }
            ReadState::PrimaryRev | ReadState::SecondaryRev | ReadState::SupplementaryRev => {
                Ok('-')
            }
        }
    }
}

impl fmt::Display for CurrRead {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut output_string = String::from("");

        if let Ok(v) = &self.read_id() {
            output_string = output_string + "\t\"read_id\": \"" + v + "\",\n";
        }

        if let Ok(v) = self.seq_len() {
            output_string = output_string + "\t\"sequence_length\": " + &v.to_string() + ",\n";
        }

        if let Ok(v) = self.align_len() {
            output_string = output_string + "\t\"alignment_length\": " + &v.to_string() + ",\n";
        }

        if let Ok((v, w)) = self.contig_id_and_start() {
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
        if let Ok((BaseMods { base_mods: v }, w)) = self.mod_data() {
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

/// Converts CurrRead to StrandedBed3
///
/// ```
/// use bedrs::{Coordinates, Strand};
/// use bedrs::prelude::StrandedBed3;
/// use nanalogue_core::{CurrRead, Error, nanalogue_bam_reader};
/// use rust_htslib::bam::Read;
/// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
/// let mut count = 0;
/// for record in reader.records(){
///     let r = record?;
///     let mut curr_read = CurrRead::default();
///     curr_read.set_read_state(&r);
///     curr_read.set_align_len(&r);
///     curr_read.set_contig_id_and_start(&r);
///     let Ok(bed3_stranded) = StrandedBed3::try_from(curr_read) else {unreachable!()};
///     let exp_bed3_stranded = match count {
///         0 => StrandedBed3::new(0, 9, 17, Strand::Forward),
///         1 => StrandedBed3::new(2, 23, 71, Strand::Forward),
///         2 => StrandedBed3::new(1, 3, 36, Strand::Reverse),
///         3 => StrandedBed3::empty(),
///         _ => unreachable!(),
///     };
///     assert_eq!(*bed3_stranded.chr(), *exp_bed3_stranded.chr());
///     assert_eq!(bed3_stranded.start(), exp_bed3_stranded.start());
///     assert_eq!(bed3_stranded.end(), exp_bed3_stranded.end());
///     assert_eq!(bed3_stranded.strand(), exp_bed3_stranded.strand());
///     count = count + 1;
/// }
/// # Ok::<(), Error>(())
/// ```
impl TryFrom<CurrRead> for StrandedBed3<i32, u64> {
    type Error = crate::Error;

    fn try_from(value: CurrRead) -> Result<Self, Self::Error> {
        match (
            value.read_state(),
            value.align_len().ok(),
            value.contig_id_and_start().ok(),
        ) {
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
        let mut curr_read_state = CurrRead::try_from_only_alignment(&record)?;
        match curr_read_state.set_mod_data(&record, ThresholdState::GtEq(128), 0) {
            Ok(_) | Err(Error::NoModInfo) => {}
            Err(e) => return Err(e),
        };
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
        let mut curr_read_state = CurrRead::try_from_only_alignment(&record)?;
        match curr_read_state.set_mod_data(&record, ThresholdState::GtEq(128), 0) {
            Ok(_) | Err(Error::NoModInfo) => {}
            Err(e) => return Err(e),
        }
        Ok(curr_read_state)
    }
}

/// Implements filter by reference coordinates for the Ranges
/// struct that contains our modification information.
/// NOTE: Ranges does not contain contig information, so we cannot
/// filter by that here.
impl FilterByRefCoords for Ranges {
    /// filters by reference position i.e. all pos such that start <= pos < end
    /// are retained. does not use contig in filtering.
    /// Copied and edited from the fibertools-rs repository.
    fn filter_by_ref_pos(&mut self, start: i64, end: i64) {
        let to_keep = self
            .reference_starts
            .iter()
            .enumerate()
            .filter_map(|(i, &s)| {
                if let Some(s) = s {
                    if s < start || s >= end { None } else { Some(i) }
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();

        self.starts = to_keep.iter().map(|&i| self.starts[i]).collect();
        self.ends = to_keep.iter().map(|&i| self.ends[i]).collect();
        self.lengths = to_keep.iter().map(|&i| self.lengths[i]).collect();
        self.qual = to_keep.iter().map(|&i| self.qual[i]).collect();
        self.reference_starts = to_keep.iter().map(|&i| self.reference_starts[i]).collect();
        self.reference_ends = to_keep.iter().map(|&i| self.reference_ends[i]).collect();
        self.reference_lengths = to_keep.iter().map(|&i| self.reference_lengths[i]).collect();
    }
}

impl Intersects<Range<u64>> for Range<u64> {
    /// Check if a range is within another range
    ///
    /// ```
    /// use nanalogue_core::Intersects;
    /// assert!((0..3).intersects(&(0..1)));
    /// assert!(!(0..3).intersects(&(5..7)));
    /// assert!(!(0..3).intersects(&(1..1)));
    /// assert!((1..3).intersects(&(0..2)));
    /// ```
    fn intersects(&self, val: &Range<u64>) -> bool {
        let end_minus_1 = val.end - 1;
        (self.contains(&val.start) || self.contains(&end_minus_1))
            && !(self.is_empty() || val.is_empty())
    }
}

/// Implements filter by reference coordinates for our CurrRead
impl FilterByRefCoords for CurrRead {
    /// filters by reference position i.e. all pos such that start <= pos < end
    /// are retained. does not use contig in filtering.
    /// First, we check if the read's coordinates are within the interval of interest,
    /// otherwise we just reset the mod data. Then, we call the filter method
    /// in the Ranges struct that contains the mod data.
    fn filter_by_ref_pos(&mut self, start: i64, end: i64) {
        let (start_u64, end_u64) = match (start, end) {
            (st, en) if st >= 0 && en >= 0 && st <= en => (st as u64, en as u64),
            _ => panic!("start and end are invalid!"),
        };
        match (self.contig_id_and_start().ok(), self.align_len().ok()) {
            (Some((_, start_pos)), Some(al))
                if !(start_u64..end_u64).intersects(&(start_pos..start_pos + al)) =>
            {
                self.reset_mod_data();
            }
            _ => match &mut self.mods {
                None => {}
                Some((BaseMods { base_mods: v }, _)) => {
                    for k in v {
                        k.ranges.filter_by_ref_pos(start, end);
                    }
                }
            },
        };
    }
}
