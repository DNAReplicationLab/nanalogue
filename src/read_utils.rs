//! # ReadUtils
//!
//! Implements CurrRead Struct for processing information retrieved from BAM files
//! and the mod information in the BAM file using a parser implemented in
//! another module.

use bedrs::prelude::{Intersect, StrandedBed3};
use bedrs::{Coordinates, Strand};
use bio_types::genome::AbstractInterval;
use fibertools_rs::utils::bamranges::Ranges;
use fibertools_rs::utils::basemods::{BaseMod, BaseMods};
use rust_htslib::{bam::ext::BamRecordExtensions, bam::record::Record};
use std::collections::HashMap;
use std::fmt;
use std::num::NonZeroU64;
use std::ops::Range;
use std::rc::Rc;

// Import from our crate
use crate::{
    Contains, Error, F32Bw0and1, FilterByRefCoords, InputModOptions, Intersects, ModChar,
    ReadState, ThresholdState, nanalogue_mm_ml_parser,
};

/// Shows CurrRead has no data
#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct NoData;

/// Shows CurrRead has only alignment data
#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct OnlyAlignData;

/// Shows CurrRead has only alignment data but with all fields filled
#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct OnlyAlignDataComplete;

/// Shows CurrRead has alignment and modification data
#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct AlignAndModData;

/// Dummy trait
pub trait CurrReadState {}

impl CurrReadState for NoData {}
impl CurrReadState for OnlyAlignData {}
impl CurrReadState for OnlyAlignDataComplete {}
impl CurrReadState for AlignAndModData {}

/// Another dummy trait
pub trait CurrReadStateWithAlign {}

impl CurrReadStateWithAlign for OnlyAlignData {}
impl CurrReadStateWithAlign for OnlyAlignDataComplete {}
impl CurrReadStateWithAlign for AlignAndModData {}

/// Another dummy trait
pub trait CurrReadStateOnlyAlign {}

impl CurrReadStateOnlyAlign for OnlyAlignData {}
impl CurrReadStateOnlyAlign for OnlyAlignDataComplete {}

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
#[derive(Debug, Clone, PartialEq)]
pub struct CurrRead<S: CurrReadState> {
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
    mods: (BaseMods, ThresholdState),

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

    /// PhantomData marker for compiler's sake
    marker: std::marker::PhantomData<S>,
}

/// Implements defaults for CurrRead
impl Default for CurrRead<NoData> {
    fn default() -> Self {
        CurrRead::<NoData> {
            state: ReadState::PrimaryFwd,
            read_id: None,
            seq_len: None,
            align_len: None,
            mods: (BaseMods { base_mods: vec![] }, ThresholdState::default()),
            contig_id_and_start: None,
            contig_name: None,
            mod_base_qual_thres: 0,
            marker: std::marker::PhantomData::<NoData>,
        }
    }
}

impl<S: CurrReadState> CurrRead<S> {
    /// resets the read state
    pub fn reset(self) -> CurrRead<NoData> {
        CurrRead::default()
    }
}

impl CurrRead<NoData> {
    /// sets the alignment of the read using BAM record
    ///
    /// ```
    /// use nanalogue_core::{CurrRead, Error, nanalogue_bam_reader, ReadState};
    /// use rust_htslib::bam::Read;
    /// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    /// let mut count = 0;
    /// for record in reader.records(){
    ///     let r = record?;
    ///     let curr_read = CurrRead::default().set_read_state(&r)?;
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
    pub fn set_read_state(self, record: &Record) -> Result<CurrRead<OnlyAlignData>, Error> {
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
        let state = match (
            record.is_reverse(),
            record.is_unmapped(),
            record.is_secondary(),
            record.is_supplementary(),
        ) {
            (false, true, false, false) => ReadState::Unmapped,
            (true, false, false, false) => ReadState::PrimaryRev,
            (true, false, true, false) => ReadState::SecondaryRev,
            (false, false, true, false) => ReadState::SecondaryFwd,
            (true, false, false, true) => ReadState::SupplementaryRev,
            (false, false, false, true) => ReadState::SupplementaryFwd,
            (false, false, false, false) => ReadState::PrimaryFwd,
            _ => return Err(Error::UnknownAlignState),
        };
        Ok(CurrRead::<OnlyAlignData> {
            state,
            read_id: None,
            seq_len: None,
            align_len: None,
            mods: (BaseMods { base_mods: vec![] }, ThresholdState::default()),
            contig_id_and_start: None,
            contig_name: None,
            mod_base_qual_thres: 0,
            marker: std::marker::PhantomData::<OnlyAlignData>,
        })
    }

    /// Uses only alignment information and no modification information to
    /// create the struct. Use this if you want to perform operations that
    /// do not involve reading or manipulating the modification data.
    pub fn try_from_only_alignment(
        self,
        record: &Record,
    ) -> Result<CurrRead<OnlyAlignDataComplete>, Error> {
        let mut curr_read_state = CurrRead::default()
            .set_read_state(record)?
            .set_seq_len(record)?
            .set_read_id(record)?;
        match curr_read_state.read_state() {
            ReadState::Unmapped => {}
            _ => {
                curr_read_state = curr_read_state
                    .set_align_len(record)?
                    .set_contig_id_and_start(record)?
                    .set_contig_name(record)?;
            }
        }
        let CurrRead::<OnlyAlignData> {
            state,
            read_id,
            seq_len,
            align_len,
            mods: _,
            contig_id_and_start,
            contig_name,
            mod_base_qual_thres: _,
            marker: _,
        } = curr_read_state;
        Ok(CurrRead::<OnlyAlignDataComplete> {
            state,
            read_id,
            seq_len,
            align_len,
            mods: (BaseMods { base_mods: vec![] }, ThresholdState::default()),
            contig_id_and_start,
            contig_name,
            mod_base_qual_thres: 0,
            marker: std::marker::PhantomData::<OnlyAlignDataComplete>,
        })
    }
}

impl<S: CurrReadStateWithAlign + CurrReadState> CurrRead<S> {
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
    ///     let curr_read = CurrRead::default().set_read_state(&r)?.set_seq_len(&r)?;
    ///     let Ok(len) = curr_read.seq_len() else { unreachable!() };
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
    ///     let curr_read = CurrRead::default().set_read_state(&r)?
    ///         .set_seq_len(&r)?.set_seq_len(&r)?;
    ///     break;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    pub fn set_seq_len(mut self, record: &Record) -> Result<Self, Error> {
        self.seq_len = match self.seq_len {
            Some(_) => Err(Error::InvalidDuplicates(
                "cannot set sequence length again!".to_string(),
            )),
            None => Ok(Some(
                NonZeroU64::new(record.seq_len().try_into()?)
                    .ok_or(Error::InvalidSeqLength)?
                    .get(),
            )),
        }?;
        Ok(self)
    }
    /// gets length of sequence
    pub fn seq_len(&self) -> Result<u64, Error> {
        self.seq_len.ok_or(Error::UnavailableData)
    }
    /// set alignment length from BAM record if available
    pub fn set_align_len(mut self, record: &Record) -> Result<Self, Error> {
        self.align_len = match self.align_len {
            Some(_) => Err(Error::InvalidDuplicates(
                "cannot set alignment length again!".to_string(),
            )),
            None => match self.read_state() {
                ReadState::Unmapped => Err(Error::Unmapped),
                _ => {
                    let st: i64 = record.pos();
                    let en: i64 = record.reference_end();
                    if en > st && st >= 0 {
                        Ok(Some((en - st).try_into()?))
                    } else {
                        Err(Error::InvalidAlignLength)
                    }
                }
            },
        }?;
        Ok(self)
    }
    /// gets alignment length
    pub fn align_len(&self) -> Result<u64, Error> {
        match self.read_state() {
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
    ///     let curr_read =
    ///         CurrRead::default().set_read_state(&r)?.set_contig_id_and_start(&r)?;
    ///     match (count, curr_read.contig_id_and_start()) {
    ///         (0, Ok((0, 9))) |
    ///         (1, Ok((2, 23))) |
    ///         (2, Ok((1, 3))) => {},
    ///         _ => unreachable!(),
    ///     }
    ///     count = count + 1;
    ///     if count == 3 { break; } // the fourth entry is unmapped, and will lead to an error.
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    ///
    /// If we call the method on an unmapped read, we should see an error.
    /// ```should_panic
    /// # use nanalogue_core::{CurrRead, Error, nanalogue_bam_reader};
    /// # use rust_htslib::bam::Read;
    /// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    /// let mut count = 0;
    /// for record in reader.records(){
    ///     let r = record?;
    ///     if count < 3 {
    ///         count = count + 1;
    ///         continue;
    ///     }
    ///     // the fourth read is unmapped
    ///     let curr_read =
    ///         CurrRead::default().set_read_state(&r)?.set_contig_id_and_start(&r)?;
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
    ///     let mut curr_read = CurrRead::default().set_read_state(&r)?
    ///         .set_contig_id_and_start(&r)?.set_contig_id_and_start(&r)?;
    ///     break;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    pub fn set_contig_id_and_start(mut self, record: &Record) -> Result<Self, Error> {
        self.contig_id_and_start = match self.contig_id_and_start {
            Some(_) => Err(Error::InvalidDuplicates(
                "cannot set contig and start again!".to_string(),
            )),
            None => match self.read_state() {
                ReadState::Unmapped => Err(Error::Unmapped),
                _ => Ok(Some((record.tid(), record.pos().try_into()?))),
            },
        }?;
        Ok(self)
    }
    /// gets contig ID and start
    pub fn contig_id_and_start(&self) -> Result<(i32, u64), Error> {
        match self.read_state() {
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
    ///     let mut curr_read = CurrRead::default().set_read_state(&r)?.set_contig_name(&r)?;
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
    ///     let mut curr_read = CurrRead::default().set_read_state(&r)?;
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
    ///     let mut curr_read = CurrRead::default().set_read_state(&r)?
    ///         .set_contig_name(&r)?.set_contig_name(&r)?;
    ///     break;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    pub fn set_contig_name(mut self, record: &Record) -> Result<Self, Error> {
        self.contig_name = match (self.read_state(), self.contig_name) {
            (ReadState::Unmapped, _) => Err(Error::Unmapped),
            (_, Some(_)) => Err(Error::InvalidDuplicates(
                "cannot set contig name again!".to_string(),
            )),
            (_, None) => Ok(Some(String::from(record.contig()))),
        }?;
        Ok(self)
    }
    /// gets contig name
    pub fn contig_name(&self) -> Result<&str, Error> {
        match (self.read_state(), &self.contig_name) {
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
    ///     let mut curr_read = CurrRead::default().set_read_state(&r)?.set_read_id(&r)?;
    ///     let Ok(read_id) = curr_read.read_id() else { unreachable!() };
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
    ///     let mut curr_read = CurrRead::default().set_read_state(&r)?
    ///         .set_read_id(&r)?.set_read_id(&r)?;
    ///     break;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    pub fn set_read_id(mut self, record: &Record) -> Result<Self, Error> {
        self.read_id = match self.read_id {
            Some(_) => Err(Error::InvalidDuplicates(
                "cannot set read id again!".to_string(),
            )),
            None => match str::from_utf8(record.qname()) {
                Ok(v) => Ok(Some(v.to_string())),
                Err(_) => Err(Error::InvalidReadID),
            },
        }?;
        Ok(self)
    }
    /// gets read id
    pub fn read_id(&self) -> Result<&str, Error> {
        match &self.read_id {
            None => Err(Error::UnavailableData),
            Some(v) => Ok(v.as_str()),
        }
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
    ///     let mut curr_read = CurrRead::default().set_read_state(&r)?;
    ///     let strand = curr_read.strand() else { unreachable!() };
    ///     match (count, strand) {
    ///         (0, '+') | (1, '+') | (2, '-') | (3, '.') => {},
    ///         _ => unreachable!(),
    ///     }
    ///     count = count + 1;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    pub fn strand(&self) -> char {
        match &self.state {
            ReadState::Unmapped => '.',
            ReadState::PrimaryFwd | ReadState::SecondaryFwd | ReadState::SupplementaryFwd => '+',
            ReadState::PrimaryRev | ReadState::SecondaryRev | ReadState::SupplementaryRev => '-',
        }
    }
    /// sets modification data using the BAM record
    pub fn set_mod_data(
        self,
        record: &Record,
        mod_thres: ThresholdState,
        min_qual: u8,
    ) -> Result<CurrRead<AlignAndModData>, Error> {
        let result = nanalogue_mm_ml_parser(
            record,
            |x| mod_thres.contains(x),
            |_| true,
            |_, _, _| true,
            min_qual,
        )?;
        Ok(CurrRead::<AlignAndModData> {
            state: self.state,
            read_id: self.read_id,
            seq_len: self.seq_len,
            align_len: self.align_len,
            mods: (result, mod_thres),
            contig_id_and_start: self.contig_id_and_start,
            contig_name: self.contig_name,
            mod_base_qual_thres: min_qual,
            marker: std::marker::PhantomData::<AlignAndModData>,
        })
    }
    /// sets modification data using BAM record but restricted to the
    /// specified filters
    pub fn set_mod_data_restricted<G, H>(
        self,
        record: &Record,
        mod_thres: ThresholdState,
        mod_fwd_pos_filter: G,
        mod_filter_base_strand_tag: H,
        min_qual: u8,
    ) -> Result<CurrRead<AlignAndModData>, Error>
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
        Ok(CurrRead::<AlignAndModData> {
            state: self.state,
            read_id: self.read_id,
            seq_len: self.seq_len,
            align_len: self.align_len,
            mods: (result, mod_thres),
            contig_id_and_start: self.contig_id_and_start,
            contig_name: self.contig_name,
            mod_base_qual_thres: min_qual,
            marker: std::marker::PhantomData::<AlignAndModData>,
        })
    }
}

impl CurrRead<OnlyAlignDataComplete> {
    /// sets modification data using BAM record but with restrictions
    /// applied by the InputMods options
    pub fn set_mod_data_restricted_options(
        self,
        record: &Record,
        mod_options: &impl InputModOptions,
    ) -> Result<CurrRead<AlignAndModData>, Error> {
        let l = usize::try_from(self.seq_len().expect("no error")).expect("bit conversion error");
        let w = mod_options.trim_read_ends();
        let interval = if let Some(bed3) = mod_options.region_filter() {
            let stranded_bed3 = StrandedBed3::<i32, u64>::try_from(&self)?;
            if let Some(v) = bed3.intersect(&stranded_bed3) {
                if v.start() == stranded_bed3.start() && v.end() == stranded_bed3.end() {
                    None
                } else {
                    Some(v.start()..v.end())
                }
            } else {
                Some(0..0)
            }
        } else {
            None
        };
        Ok({
            let mut read = self.set_mod_data_restricted(
                record,
                mod_options.mod_prob_filter(),
                |x| w == 0 || (w..l.checked_sub(w).unwrap_or_default()).contains(x),
                |_, &s, &t| {
                    mod_options.tag().is_none_or(|x| x == t)
                        && mod_options.mod_strand().is_none_or(|v| s == char::from(v))
                },
                mod_options.base_qual_filter(),
            )?;
            if let Some(v) = interval {
                read.filter_by_ref_pos(i64::try_from(v.start)?, i64::try_from(v.end)?)
            }
            read
        })
    }
}

impl CurrRead<AlignAndModData> {
    /// gets modification data
    pub fn mod_data(&self) -> &(BaseMods, ThresholdState) {
        &self.mods
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
        let (BaseMods { base_mods: v }, _) = &self.mods;
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
    ///     let curr_read = CurrRead::default().set_read_state(&r)?.set_mod_data(&r, ThresholdState::GtEq(180), 0)?;
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
        let (BaseMods { base_mods: v }, _) = self.mod_data();
        for k in v {
            let base_count = k.ranges.qual.len() as u32;
            let _ = output
                .entry(ModChar::new(k.modification_type))
                .and_modify(|e| *e += base_count)
                .or_insert(base_count);
        }
        Some(output)
    }
}

impl<S: CurrReadStateOnlyAlign + CurrReadState> fmt::Display for CurrRead<S> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut output_string = String::from("");
        let CurrRead::<S> {
            read_id,
            seq_len,
            align_len,
            contig_id_and_start,
            contig_name,
            state,
            mods: _,
            mod_base_qual_thres: _,
            marker: _,
        } = self;

        if let Some(v) = read_id {
            output_string = output_string + "\t\"read_id\": \"" + v + "\",\n";
        }

        if let Some(v) = seq_len {
            output_string = output_string + "\t\"sequence_length\": " + &v.to_string() + ",\n";
        }

        if let Some((v, w)) = contig_id_and_start {
            let num_str = &v.to_string();
            output_string = output_string
                + "\t\"contig\": \""
                + if let Some(x) = contig_name {
                    x
                } else {
                    num_str
                }
                + "\",\n";
            output_string = output_string + "\t\"reference_start\": " + &w.to_string() + ",\n";
            if let Some(x) = align_len {
                output_string =
                    output_string + "\t\"reference_end\": " + &(w + x).to_string() + ",\n";
                output_string = output_string + "\t\"alignment_length\": " + &x.to_string() + ",\n";
            }
        }

        output_string = output_string + "\t\"alignment_type\": \"" + &state.to_string() + "\",\n";

        write!(f, "{{\n{output_string}}}")
    }
}

impl fmt::Display for CurrRead<AlignAndModData> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut output_string = String::from("");
        let CurrRead::<AlignAndModData> {
            read_id,
            seq_len,
            align_len,
            contig_id_and_start,
            contig_name,
            state,
            mods,
            mod_base_qual_thres,
            marker: _,
        } = self;

        if let Some(v) = read_id {
            output_string = output_string + "\t\"read_id\": \"" + v + "\",\n";
        }

        if let Some(v) = seq_len {
            output_string = output_string + "\t\"sequence_length\": " + &v.to_string() + ",\n";
        }

        if let Some((v, w)) = contig_id_and_start {
            let num_str = &v.to_string();
            output_string = output_string
                + "\t\"contig\": \""
                + if let Some(x) = contig_name {
                    x
                } else {
                    num_str
                }
                + "\",\n";
            output_string = output_string + "\t\"reference_start\": " + &w.to_string() + ",\n";
            if let Some(x) = align_len {
                output_string =
                    output_string + "\t\"reference_end\": " + &(w + x).to_string() + ",\n";
                output_string = output_string + "\t\"alignment_length\": " + &x.to_string() + ",\n";
            }
        }

        output_string = output_string + "\t\"alignment_type\": \"" + &state.to_string() + "\",\n";

        let mut mod_count_str = String::from("");
        let (BaseMods { base_mods: v }, w) = mods;
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
        if mod_count_str.is_empty() {
            mod_count_str += "NA";
        } else {
            mod_count_str +=
                format!("({}, PHRED base qual >= {})", w, mod_base_qual_thres).as_str();
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
///     let mut curr_read = CurrRead::default().try_from_only_alignment(&r)?;
///     let Ok(bed3_stranded) = StrandedBed3::try_from(&curr_read) else {unreachable!()};
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
impl<S: CurrReadStateWithAlign + CurrReadState> TryFrom<&CurrRead<S>> for StrandedBed3<i32, u64> {
    type Error = Error;

    fn try_from(value: &CurrRead<S>) -> Result<Self, Self::Error> {
        match (
            value.read_state(),
            value.align_len().ok(),
            value.contig_id_and_start().ok(),
        ) {
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
impl TryFrom<Record> for CurrRead<AlignAndModData> {
    type Error = Error;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        let curr_read_state = CurrRead::default()
            .try_from_only_alignment(&record)?
            .set_mod_data(&record, ThresholdState::GtEq(128), 0)?;
        Ok(curr_read_state)
    }
}

/// Convert a rust htslib rc record into our struct.
/// I think the rc datatype is just like the normal record,
/// except the record datatype is not destroyed and created
/// every time a new record is read (or something like that).
/// All comments I've made for the `TryFrom<Record>` function
/// apply here as well.
impl TryFrom<Rc<Record>> for CurrRead<AlignAndModData> {
    type Error = Error;

    fn try_from(record: Rc<Record>) -> Result<Self, Self::Error> {
        let curr_read_state = CurrRead::default()
            .try_from_only_alignment(&record)?
            .set_mod_data(&record, ThresholdState::GtEq(128), 0)?;
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
        let (start_idx, end_idx) = {
            let mut start_idx = 0;
            let mut end_idx = 0;
            for k in self.reference_starts.iter().enumerate() {
                if (*k.1).is_some_and(|x| x < start) {
                    start_idx = k.0;
                }
                if (*k.1).is_some_and(|x| x < end) {
                    end_idx = k.0;
                }
            }
            (start_idx, end_idx)
        };

        for k in [
            &mut self.starts,
            &mut self.ends,
            &mut self.lengths,
            &mut self.reference_starts,
            &mut self.reference_ends,
            &mut self.reference_lengths,
        ] {
            k.extract_if(end_idx.., |_| true).for_each(drop);
            k.extract_if(0..start_idx, |_| true).for_each(drop);
        }
        self.qual.extract_if(end_idx.., |_| true).for_each(drop);
        self.qual.extract_if(0..start_idx, |_| true).for_each(drop);
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
impl FilterByRefCoords for CurrRead<AlignAndModData> {
    /// filters by reference position i.e. all pos such that start <= pos < end
    /// are retained. does not use contig in filtering.
    fn filter_by_ref_pos(&mut self, start: i64, end: i64) {
        let (BaseMods { base_mods: v }, _) = &mut self.mods;
        for k in v {
            k.ranges.filter_by_ref_pos(start, end);
        }
    }
}
