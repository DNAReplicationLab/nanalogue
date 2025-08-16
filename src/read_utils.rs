//! # ReadUtils
//!
//! Implements CurrRead Struct for processing information retrieved from BAM files
//! and the mod information in the BAM file using a parser implemented in
//! another module.

use bedrs::prelude::StrandedBed3;
use bedrs::{Coordinates, Strand};
use fibertools_rs::utils::basemods::{BaseMod, BaseMods};
use rust_htslib::{bam::ext::BamRecordExtensions, bam::record::Record};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;
use std::num::NonZeroU64;
use std::rc::Rc;

// Import from our crate
use crate::{Contains, Error, F32Bw0and1, ModChar, OrdPair, nanalogue_mm_ml_parser};

/// Alignment state of a read; seven possibilities + one unknown state
#[derive(Debug, Clone, Default, Copy, PartialEq, Serialize, Deserialize)]
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

    /// ID of the reference genome contig and the starting position on
    /// contig that the molecule maps or aligns to.
    /// NOTE: the contig here is numeric and refers to an index on the BAM
    /// header. To convert this into an alphanumeric string, you have to
    /// process the header and store it in `contig_name` below.
    /// We have left it this way as it is easier to store and process integers.
    contig_id_and_start: Option<(i32, u64)>,

    /// Contig name.
    contig_name: Option<String>,
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
    ///
    /// ```
    /// use nanalogue_core::{CurrRead, Error, nanalogue_bam_reader};
    /// use rust_htslib::bam::Read;
    /// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    /// let mut count = 0;
    /// for record in reader.records(){
    ///     let r = record?;
    ///     let mut curr_read = CurrRead::default();
    ///     curr_read.set_read_state(&r);
    ///     let Ok(Some(len)) = curr_read.set_seq_len(&r) else { unreachable!() };
    ///     let Ok(Some(len2)) = curr_read.seq_len() else { unreachable!() };
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
    ///     curr_read.set_read_state(&r);
    ///     curr_read.set_contig_id_and_start(&r);
    ///     match (count, curr_read.contig_id_and_start()) {
    ///         (0, Ok(Some((0, 9)))) |
    ///         (1, Ok(Some((2, 23)))) |
    ///         (2, Ok(Some((1, 3)))) |
    ///         (3, Ok(None)) => {},
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
    ///     curr_read.set_contig_id_and_start(&r)?;
    ///     curr_read.set_contig_id_and_start(&r)?;
    ///     break;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    pub fn set_contig_id_and_start(&mut self, record: &Record) -> Result<bool, Error> {
        match &self.contig_id_and_start {
            Some(_) => Err(Error::InvalidDuplicates(
                "cannot set contig and start again!".to_string(),
            )),
            None => match self.read_state() {
                ReadState::Unknown => Err(Error::UnknownAlignState),
                ReadState::Unmapped => Ok(false),
                _ => {
                    self.contig_id_and_start = Some((record.tid(), record.pos().try_into()?));
                    Ok(true)
                }
            },
        }
    }
    /// gets contig ID and start
    pub fn contig_id_and_start(&self) -> Result<Option<(i32, u64)>, Error> {
        match self.read_state() {
            ReadState::Unknown => Err(Error::UnknownAlignState),
            ReadState::Unmapped => Ok(None),
            _ => Ok(self.contig_id_and_start),
        }
    }
    /// sets contig name
    pub fn set_contig_name(&mut self, names: &[String]) -> Result<bool, Error> {
        match &self.contig_name {
            Some(_) => Err(Error::InvalidDuplicates(
                "cannot set contig name again!".to_string(),
            )),
            None => match self.contig_id_and_start() {
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
    ///
    /// ```
    /// use nanalogue_core::{CurrRead, Error, nanalogue_bam_reader};
    /// use rust_htslib::bam::Read;
    /// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    /// let mut count = 0;
    /// for record in reader.records(){
    ///     let r = record?;
    ///     let mut curr_read = CurrRead::default();
    ///     curr_read.set_read_state(&r);
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
    ///     curr_read.set_read_state(&r);
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
            None => Err(Error::InvalidState("read id not available".to_string())),
            Some(v) => Ok(v.as_str()),
        }
    }
    /// sets modification data using the BAM record
    pub fn set_mod_data(&mut self, record: &Record, mod_thres: ThresholdState) {
        self.mods = Some((
            nanalogue_mm_ml_parser(record, |x, _| mod_thres.contains(x), |_, _, _| true),
            mod_thres,
        ));
    }
    /// sets modification data using BAM record but restricted to the
    /// specified filters
    pub fn set_mod_data_restricted<G, H>(
        &mut self,
        record: &Record,
        mod_thres: ThresholdState,
        mod_fwd_pos_filter: G,
        mod_filter_base_strand_tag: H,
    ) where
        G: Fn(&i64) -> bool,
        H: Fn(&u8, &char, &ModChar) -> bool,
    {
        self.mods = Some((
            nanalogue_mm_ml_parser(
                record,
                |x, y| mod_thres.contains(x) && mod_fwd_pos_filter(y),
                mod_filter_base_strand_tag,
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
    /// window modification data with restrictions
    pub fn windowed_mod_data_restricted<F>(
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
        let mut plus_mod_strand_seen = false;
        let mut minus_mod_strand_seen = false;
        let tag_char = tag.get_val();
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
                        result = (0..=mod_data.len() - win_size)
                            .step_by(slide_size)
                            .map(|i| window_function(&mod_data[i..i + win_size]))
                            .collect::<Result<Vec<F32Bw0and1>, _>>()?;
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
    ///     curr_read.set_read_state(&r);
    ///     curr_read.set_mod_data(&r, ThresholdState::GtEq(180));
    ///
    ///     let mod_count = curr_read.base_count_per_mod();
    ///     let a = Some(HashMap::from([(ModChar::new('T'), 3)]));
    ///     let b = Some(HashMap::from([(ModChar::new('T'), 1)]));
    ///     match (count, mod_count) {
    ///         (0, None) => {},
    ///         (1, v) => assert_eq!(v, a),
    ///         (2, v) => assert_eq!(v, b),
    ///         (3, v) => assert_eq!(v, a),
    ///         _ => unreachable!(),
    ///     }
    ///     count = count + 1;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    #[must_use]
    pub fn base_count_per_mod(&self) -> Option<HashMap<ModChar, u32>> {
        let mut output = HashMap::<ModChar, u32>::new();
        match &self.mods {
            Some((BaseMods { base_mods: v }, _)) => {
                if v.is_empty() {
                    None
                } else {
                    for k in v {
                        let base_count = k.ranges.qual.len() as u32;
                        output
                            .entry(ModChar::new(k.modification_type))
                            .and_modify(|e| *e += base_count)
                            .or_insert(base_count);
                    }
                    Some(output)
                }
            }
            None => None,
        }
    }
    /// Uses only alignment information and no modification information to
    /// create the struct. Use this if you want to perform operations that
    /// do not involve reading or manipulating the modification data.
    pub fn try_from_only_alignment(record: Record) -> Result<Self, Error> {
        let mut curr_read_state = CurrRead::default();
        curr_read_state.set_read_state(&record)?;
        curr_read_state.set_align_len(&record)?;
        curr_read_state.set_contig_id_and_start(&record)?;
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
    ///     curr_read.set_read_state(&r);
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

        if let Ok(Some(v)) = self.seq_len() {
            output_string = output_string + "\t\"sequence_length\": " + &v.to_string() + ",\n";
        }

        if let Ok(Some(v)) = self.align_len() {
            output_string = output_string + "\t\"alignment_length\": " + &v.to_string() + ",\n";
        }

        if let Ok(Some((v, w))) = self.contig_id_and_start() {
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
        match (value.state, value.align_len, value.contig_id_and_start) {
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
        curr_read_state.set_contig_id_and_start(&record)?;
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
        curr_read_state.set_contig_id_and_start(&record)?;
        Ok(curr_read_state)
    }
}
