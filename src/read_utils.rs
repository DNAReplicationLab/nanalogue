//! # `ReadUtils`
//!
//! Implements `CurrRead` Struct for processing information retrieved from BAM files
//! and the mod information in the BAM file using a parser implemented in
//! another module.

use bedrs::prelude::{Intersect, StrandedBed3};
use bedrs::{Bed3, Coordinates, Strand};
use bio_types::genome::AbstractInterval;
use fibertools_rs::utils::bamranges::Ranges;
use fibertools_rs::utils::basemods::{BaseMod, BaseMods};
use rust_htslib::{bam::ext::BamRecordExtensions, bam::record::Record};
use std::collections::{HashMap, HashSet};
use std::fmt::{self, Write};
use std::num::NonZeroU64;
use std::rc::Rc;

// Import from our crate
use crate::{
    Contains, Error, F32Bw0and1, FilterByRefCoords, InputModOptions, InputRegionOptions,
    InputWindowing, ModChar, ReadState, ThresholdState, nanalogue_mm_ml_parser,
};
use serde::{Deserialize, Serialize};

/// Shows `CurrRead` has no data
#[derive(Debug, Default, Copy, Clone, PartialEq)]
#[non_exhaustive]
pub struct NoData;

/// Shows `CurrRead` has only alignment data
#[derive(Debug, Default, Copy, Clone, PartialEq)]
#[non_exhaustive]
pub struct OnlyAlignData;

/// Shows `CurrRead` has only alignment data but with all fields filled
#[derive(Debug, Default, Copy, Clone, PartialEq)]
#[non_exhaustive]
pub struct OnlyAlignDataComplete;

/// Shows `CurrRead` has alignment and modification data
#[derive(Debug, Default, Copy, Clone, PartialEq)]
#[non_exhaustive]
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
/// The information within the struct is hard to access without
/// the methods defined here. This is to ensure the struct
/// doesn't fall into an invalid state, which could cause mistakes
/// in calculations associated with the struct. For example:
/// if I want to measure mean modification density along windows
/// of the raw modification data, I need a guarantee that the
/// modification data is sorted by position. We can guarantee
/// this when the modification data is parsed, but we cannot
/// guarantee this if we allow free access to the struct.
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

    /// `PhantomData` marker for compiler's sake
    marker: std::marker::PhantomData<S>,
}

/// Implements defaults for `CurrRead`
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

impl CurrRead<NoData> {
    /// sets the alignment of the read using BAM record
    ///
    /// # Errors
    /// While we support normal BAM reads from `ONT`, `PacBio` etc. that contain modifications,
    /// we do not support some BAM flags like paired, duplicate, quality check failed etc.
    /// This is because of our design choices e.g. if mods are called on paired reads,
    /// then we'll have to include both records as one read in our statistics
    /// and we do not have functionality in place to do this.
    /// So, we return errors if such flags or an invalid combination of flags (e.g.
    /// secondary and supplementary bits are set) are encountered.
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
            || record.is_quality_check_failed()
        {
            return Err(Error::NotImplementedError(
                "paired-read/mate-read/duplicate/qual-check-failed flags not supported!"
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
    ///
    /// # Errors
    /// Errors are returned if getting record information fails e.g. read id
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
    /// # Errors
    /// Errors are returned if sequence length is already set or
    /// sequence length is not non-zero.
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
    ///
    /// # Errors
    /// Error if sequence length is not set
    pub fn seq_len(&self) -> Result<u64, Error> {
        self.seq_len.ok_or(Error::UnavailableData)
    }
    /// set alignment length from BAM record if available
    ///
    /// # Errors
    /// Returns errors if alignment len is already set, instance is
    /// unmapped, or if alignment coordinates are malformed
    /// (e.g. end < start).
    pub fn set_align_len(mut self, record: &Record) -> Result<Self, Error> {
        self.align_len = match self.align_len {
            Some(_) => Err(Error::InvalidDuplicates(
                "cannot set alignment length again!".to_string(),
            )),
            None => {
                if self.read_state() == ReadState::Unmapped {
                    Err(Error::Unmapped)
                } else {
                    // NOTE: right now, I don't know of a way to test the error below
                    // as rust htslib initializes an empty record with an alignment
                    // length of 1 (see the code below). This is only a note about
                    // the error variant, not the normal function of this code block
                    // which is fine.
                    // ```
                    // use rust_htslib::bam::ext::BamRecordExtensions;
                    // let r = Record::new();
                    // assert_eq!(r.seq_len(), 0);
                    // assert_eq!(r.pos(), 0);
                    // assert_eq!(r.reference_end(), 1);
                    // ```
                    let st: i64 = record.pos();
                    let en: i64 = record.reference_end();
                    if en > st && st >= 0 {
                        #[expect(
                            clippy::arithmetic_side_effects,
                            reason = "en > st && st >= 0 guarantee no i64 overflows"
                        )]
                        #[expect(
                            clippy::missing_panics_doc,
                            reason = "en > st && st >= 0 guarantee no panic"
                        )]
                        Ok(Some(
                            (en - st)
                                .try_into()
                                .expect("en>st && st>=0 guarantee no problems i64->u64"),
                        ))
                    } else {
                        Err(Error::InvalidAlignLength)
                    }
                }
            }
        }?;
        Ok(self)
    }
    /// gets alignment length
    ///
    /// # Errors
    /// if instance is unmapped or alignment length is not set
    pub fn align_len(&self) -> Result<u64, Error> {
        match self.read_state() {
            ReadState::Unmapped => Err(Error::Unmapped),
            _ => self.align_len.ok_or(Error::UnavailableData),
        }
    }
    /// sets contig ID and start from BAM record if available
    ///
    /// # Errors
    /// if instance is unmapped, if these data are already set and
    /// the user is trying to set them again, or if coordinates
    /// are malformed (start position is negative)
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
    ///
    /// # Errors
    /// If instance is unmapped or if data (contig id and start) are not set
    pub fn contig_id_and_start(&self) -> Result<(i32, u64), Error> {
        match self.read_state() {
            ReadState::Unmapped => Err(Error::Unmapped),
            _ => self.contig_id_and_start.ok_or(Error::UnavailableData),
        }
    }
    /// sets contig name
    ///
    /// # Errors
    /// Returns error if instance is unmapped or contig name has already been set
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
    ///
    /// # Errors
    /// If instance is unmapped or contig name has not been set
    pub fn contig_name(&self) -> Result<&str, Error> {
        match (self.read_state(), &self.contig_name) {
            (ReadState::Unmapped, _) => Err(Error::Unmapped),
            (_, None) => Err(Error::UnavailableData),
            (_, Some(v)) => Ok(v.as_str()),
        }
    }
    /// sets read ID (also called query name) from BAM record
    ///
    /// # Errors
    /// If read id has already been set for the instance, or
    /// if read id is not valid UTF-8
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
    ///
    /// # Errors
    /// If read id has not been set
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
    #[must_use]
    pub fn strand(&self) -> char {
        match &self.state {
            ReadState::Unmapped => '.',
            ReadState::PrimaryFwd | ReadState::SecondaryFwd | ReadState::SupplementaryFwd => '+',
            ReadState::PrimaryRev | ReadState::SecondaryRev | ReadState::SupplementaryRev => '-',
        }
    }
    /// Returns subset of sequence using reference coordinates
    ///
    /// # Errors
    /// If conversion to Bed3 fails, malformed coordinates,
    /// or no intersection with given region, or conversion to usize errors.
    /// Absence of a sequence due to most of the above reasons is reported
    /// with `Error::UnavailableData`.
    ///
    /// ```
    /// use bedrs::Bed3;
    /// use nanalogue_core::{CurrRead, Error, nanalogue_bam_reader};
    /// use rust_htslib::bam::Read;
    ///
    /// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
    /// let mut count = 0;
    /// for record in reader.records() {
    ///     let r = record?;
    ///     let curr_read = CurrRead::default().try_from_only_alignment(&r)?;
    ///
    ///     // Skip unmapped reads
    ///     if curr_read.read_state().to_string() != "unmapped" {
    ///         let (contig_id, start) = curr_read.contig_id_and_start()?;
    ///         let align_len = curr_read.align_len()?;
    ///
    ///         // Create a region that overlaps with the read but is short of one bp.
    ///         // Note that this BAM file has reads with all bases matching perfectly
    ///         // with the reference.
    ///         let region = Bed3::new(contig_id, start, start + align_len - 1);
    ///         let seq_subset = curr_read.seq_on_ref_coords(&r, &region)?;
    ///
    ///         // Check for sequence length match
    ///         assert_eq!(curr_read.seq_len()? - 1, u64::try_from(seq_subset.len())?);
    ///
    ///         // Create a region with no overlap at all and check we get no data
    ///         let region = Bed3::new(contig_id, start + align_len, start + align_len + 2);
    ///         match curr_read.seq_on_ref_coords(&r, &region){
    ///             Err(Error::UnavailableData) => Ok(()),
    ///             _ => Err(Error::UnknownError),
    ///         }?;
    ///
    ///     }
    ///     count += 1;
    /// }
    /// # Ok::<(), Error>(())
    /// ```
    pub fn seq_on_ref_coords(
        &self,
        record: &Record,
        region: &Bed3<i32, u64>,
    ) -> Result<Vec<u8>, Error> {
        #[expect(
            clippy::missing_panics_doc,
            reason = "genomic coordinates are far less than (2^64-1)/2 i.e. u64->i64 should be ok"
        )]
        let interval = {
            let stranded_bed3 = StrandedBed3::<i32, u64>::try_from(self)?;
            if let Some(v) = region.intersect(&stranded_bed3) {
                let start = i64::try_from(v.start())
                    .expect("genomic coordinates are far less than (2^64 - 1)/2");
                let end = i64::try_from(v.end())
                    .expect("genomic coordinates are far less than (2^64 - 1)/2");
                if start < end {
                    Ok(start..end)
                } else {
                    Err(Error::UnavailableData)
                }
            } else {
                Err(Error::UnavailableData)
            }
        }?;

        // Get sequence and initialize subset.
        // We don't know how long the subset will be, we initialize with a guess
        // of 2 * interval size
        let seq = record.seq();
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "genomic coordinates far less than i64::MAX (approx (2^64-1)/2)"
        )]
        let mut s: Vec<u8> =
            Vec::with_capacity(usize::try_from(2 * (interval.end - interval.start))?);

        // we may have to trim the sequence if we hit a bunch of unaligned base
        // pairs right at the end e.g. a softclip.
        let mut trim_end_bp: u64 = 0;

        for w in record
            .aligned_pairs_full()
            .skip_while(|x| x[1].is_none_or(|y| !interval.contains(&y)))
            .take_while(|x| x[1].is_none_or(|y| interval.contains(&y)))
        {
            // the logic below is as follows:
            // * matches or mismatches, we show the base on the sequence.
            //   so SNPs for example (i.e. a 1 bp difference from the ref) will show up
            //   as a base different from the reference.
            // * a deletion or a ref skip ("N" in cigar) will show up as dot(s) i.e. ".".
            // * insertions are displayed i.e. bases in the middle of a read present
            //   on the read but not on the reference
            // * clipped bases at the end of the read are not displayed. these are bp
            //   on the read but not on the reference and are denoted as soft or hard
            //   clips on the CIGAR string e.g. barcodes from sequencing
            // * some CIGAR combinations are illogical and we are assuming they do not happen
            //   e.g. a read's CIGAR can end with, say 10D20S, this means last 10 bp
            //   are in a deletion and the next 20 are a softclip. This is illogical,
            //   as they must be combined into a 30-bp softclip i.e. 30S. So if the
            //   aligner produces such illogical states, then the sequences reported
            //   here may be erroneous.
            #[expect(
                clippy::arithmetic_side_effects,
                reason = "coordinates far less than u64::MAX (2^64-1) so no chance of counter overflow"
            )]
            match w {
                [Some(x), Some(_)] => {
                    s.push(seq[usize::try_from(x)?]);
                    trim_end_bp = 0;
                }
                [Some(x), None] => {
                    s.push(seq[usize::try_from(x)?]);
                    trim_end_bp += 1;
                }
                [None, Some(_)] => {
                    s.push(b'.');
                    trim_end_bp = 0;
                }
                _ => {}
            }
        }

        // if last few bp in sequence are all unmapped, we remove them here.
        for _ in 0..trim_end_bp {
            let Some(_) = s.pop() else { unreachable!() };
        }

        if s.is_empty() {
            Err(Error::UnavailableData)
        } else {
            Ok(s)
        }
    }
    /// sets modification data using the BAM record
    ///
    /// # Errors
    /// If tags in the BAM record containing the modification information (MM, ML)
    /// contain mistakes.
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
    /// sets modification data using BAM record but restricted to the specified filters
    ///
    /// # Errors
    /// If tags in the BAM record containing the modification information (MM, ML)
    /// contain mistakes.
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
    /// applied by the `InputMods` options
    ///
    /// # Errors
    /// If a region filter is specified and we fail to convert current instance to Bed,
    /// and if parsing the MM/ML BAM tags fails (presumably because they are malformed).
    #[expect(
        clippy::missing_panics_doc,
        reason = "integer conversions (u64->usize, u64->i64) are expected to not fail as \
genomic coordinates are far smaller than ~2^63"
    )]
    pub fn set_mod_data_restricted_options<S: InputModOptions + InputRegionOptions>(
        self,
        record: &Record,
        mod_options: &S,
    ) -> Result<CurrRead<AlignAndModData>, Error> {
        let l = usize::try_from(self.seq_len().expect("no error"))
            .expect("bit conversion errors unlikely");
        let w = mod_options.trim_read_ends_mod();
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
                mod_options.base_qual_filter_mod(),
            )?;
            if let Some(v) = interval {
                read.filter_by_ref_pos(
                    i64::try_from(v.start)
                        .expect("no error as genomic coordinates far less than ~2^63"),
                    i64::try_from(v.end)
                        .expect("no error as genomic coordinates far less than ~2^63"),
                );
            }
            read
        })
    }
}

impl CurrRead<AlignAndModData> {
    /// gets modification data
    #[must_use]
    pub fn mod_data(&self) -> &(BaseMods, ThresholdState) {
        &self.mods
    }
    /// window modification data with restrictions.
    /// If a read has the same modification on both the basecalled
    /// strand and its complement, then windows along both are returned.
    ///
    /// If `win_size` exceeds the modification data length, no windows are produced.
    ///
    /// # Errors
    /// Returns an error if the window function returns an error.
    pub fn windowed_mod_data_restricted<F>(
        &self,
        window_function: &F,
        win_options: InputWindowing,
        tag: ModChar,
    ) -> Result<Vec<F32Bw0and1>, Error>
    where
        F: Fn(&[u8]) -> Result<F32Bw0and1, Error>,
    {
        let win_size = win_options.win.get();
        let slide_size = win_options.step.get();
        let mut result = Vec::<F32Bw0and1>::new();
        let tag_char = tag.val();
        let (BaseMods { base_mods: v }, _) = &self.mods;

        // we make a few assumptions below:
        // * data is sorted by coordinate along sequence
        // * illegal types like strand not '+' or '-', or multiple entries
        //   corresponding to the same modification strand combination
        //   e.g. C+m occuring twice.
        // we control data flow into CurrRead, checking these do not happen
        // during ingress. To be future-proof etc., we should check these things
        // here but we do not as there is no way right now to test error checking
        // as there is no way to make CurrRead fall into these illegal states.
        #[expect(
            clippy::missing_panics_doc,
            reason = "checked_sub ensures win_size <= mod_data.len() before windowing"
        )]
        for k in v {
            match k {
                BaseMod {
                    modified_base: _,
                    strand: _,
                    record_is_reverse: _,
                    modification_type: x,
                    ranges: track,
                } if *x == tag_char => {
                    let mod_data = &track.qual;
                    if let Some(v) = mod_data.len().checked_sub(win_size) {
                        result.extend(
                            (0..=v)
                                .step_by(slide_size)
                                .map(|i| {
                                    window_function(
                                        mod_data[i..]
                                            .get(0..win_size)
                                            .expect("checked len>=win_size so no error"),
                                    )
                                })
                                .collect::<Result<Vec<F32Bw0and1>, _>>()?,
                        );
                    }
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
    /// # Panics
    /// Panics if the number of modifications exceeds `u32::MAX` (approximately 4.2 billion).
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
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "u32::MAX approx 4.2 Gb, v unlikely 1 molecule is this modified"
        )]
        for k in v {
            let base_count = u32::try_from(k.ranges.qual.len()).expect("number conversion error");
            let _: &mut u32 = output
                .entry(ModChar::new(k.modification_type))
                .and_modify(|e| *e += base_count)
                .or_insert(base_count);
        }
        Some(output)
    }
}

/// To format and display modification data in a condensed manner.
trait DisplayCondensedModData {
    fn mod_data_section(&self) -> Result<String, fmt::Error>;
}

/// No mod data means no display is produced
impl<S> DisplayCondensedModData for CurrRead<S>
where
    S: CurrReadStateOnlyAlign + CurrReadState,
{
    fn mod_data_section(&self) -> Result<String, fmt::Error> {
        Ok(String::new())
    }
}

/// Implements display when mod data is available
impl DisplayCondensedModData for CurrRead<AlignAndModData> {
    fn mod_data_section(&self) -> Result<String, fmt::Error> {
        let mut mod_count_str = String::new();
        let (BaseMods { base_mods: v }, w) = &self.mods;
        for k in v {
            write!(
                mod_count_str,
                "{}{}{}:{};",
                k.modified_base as char,
                k.strand,
                ModChar::new(k.modification_type),
                k.ranges.qual.len()
            )?;
        }
        if mod_count_str.is_empty() {
            write!(mod_count_str, "NA")?;
        } else {
            write!(
                mod_count_str,
                "({}, PHRED base qual >= {})",
                w, self.mod_base_qual_thres
            )?;
        }
        Ok(format!(",\n\t\"mod_count\": \"{mod_count_str}\""))
    }
}

impl<S> fmt::Display for CurrRead<S>
where
    S: CurrReadState,
    CurrRead<S>: DisplayCondensedModData,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut output_string = String::from("{\n");

        if let Some(v) = &self.read_id {
            writeln!(output_string, "\t\"read_id\": \"{v}\",")?;
        }

        if let Some(v) = self.seq_len {
            writeln!(output_string, "\t\"sequence_length\": {v},")?;
        }

        if let Some((v, w)) = self.contig_id_and_start {
            let num_str = &v.to_string();
            writeln!(
                output_string,
                "\t\"contig\": \"{}\",",
                if let Some(x) = &self.contig_name {
                    x
                } else {
                    num_str
                }
            )?;
            writeln!(output_string, "\t\"reference_start\": {w},")?;
            if let Some(x) = self.align_len {
                writeln!(
                    output_string,
                    "\t\"reference_end\": {},",
                    w.checked_add(x)
                        .expect("numeric overflow in calculating reference_end")
                )?;
                writeln!(output_string, "\t\"alignment_length\": {x},")?;
            }
        }

        write!(output_string, "\t\"alignment_type\": \"{}\"", self.state)?;
        writeln!(output_string, "{}", &self.mod_data_section()?)?;
        output_string.push('}');
        output_string.fmt(f)
    }
}

/// Converts `CurrRead` to `StrandedBed3`
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

    #[expect(
        clippy::arithmetic_side_effects,
        reason = "u64 variables won't overflow with genomic coords (<2^64-1)"
    )]
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

/// Convert a rust htslib record to our `CurrRead` struct.
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

/// Implements filter by reference coordinates for our `CurrRead`
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

/// Serialized representation of `CurrRead` with condensed JSON format
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
struct SerializedCurrRead {
    /// The type of alignment (primary, secondary, supplementary, unmapped)
    alignment_type: ReadState,
    /// Alignment information, None for unmapped reads
    #[serde(skip_serializing_if = "Option::is_none")]
    alignment: Option<AlignmentInfo>,
    /// Condensed modification data table
    mod_table: Vec<ModTableEntry>,
    /// Read identifier
    read_id: String,
    /// Sequence length
    seq_len: u64,
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
struct AlignmentInfo {
    /// Start position on reference
    start: u64,
    /// End position on reference
    end: u64,
    /// Contig/chromosome name
    contig: String,
    /// Contig/chromosome ID
    contig_id: i32,
}

/// Individual modification table entry
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(default)]
struct ModTableEntry {
    /// Base that is modified (A, C, G, T, etc.)
    base: char,
    /// Whether this is on the plus strand
    is_strand_plus: bool,
    /// Modification code (character or numeric)
    mod_code: ModChar,
    /// Whether this modification data is implicit
    implicit: bool,
    /// Modification data as [start, `ref_start`, qual] tuples
    data: Vec<(u64, i64, u8)>,
}

impl SerializedCurrRead {
    /// Validates that coordinates in modification data are within expected ranges
    fn check_coordinates(&self) -> Result<(), Error> {
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
                        let align_range =
                            i64::try_from(alignment.start)?..i64::try_from(alignment.end)?;
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
        let serialized_curr_read: SerializedCurrRead =
            self.clone().try_into().map_err(serde::ser::Error::custom)?;
        serialized_curr_read.serialize(serializer)
    }
}

impl TryFrom<CurrRead<AlignAndModData>> for SerializedCurrRead {
    type Error = Error;

    fn try_from(curr_read: CurrRead<AlignAndModData>) -> Result<Self, Self::Error> {
        let alignment_type = curr_read.read_state();

        #[expect(
            clippy::arithmetic_side_effects,
            reason = "u64 variables won't overflow with genomic coords (<2^64-1)"
        )]
        let alignment = if curr_read.read_state() == ReadState::Unmapped {
            None
        } else {
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
                let align_len = {
                    if let Some(v) = alignment.end.checked_sub(alignment.start) {
                        Ok(Some(v))
                    } else {
                        Err(Error::InvalidAlignCoords)
                    }
                }?;
                let contig_id_and_start = Some((alignment.contig_id, alignment.start));
                let contig_name = Some(alignment.contig.clone());
                (align_len, contig_id_and_start, contig_name)
            }
            None => (None, None, None),
        };

        // Create CurrRead directly
        Ok(CurrRead {
            state: serialized.alignment_type,
            read_id: Some(serialized.read_id),
            seq_len: Some(serialized.seq_len),
            align_len,
            mods: (base_mods, ThresholdState::GtEq(0)), // Default threshold
            contig_id_and_start,
            contig_name,
            mod_base_qual_thres: 0, // Default value
            marker: std::marker::PhantomData,
        })
    }
}

/// Convert `BaseMods` to condensed `mod_table` format
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

/// Helper function to check if starts vector is sorted in ascending order
fn validate_starts_sorting(starts: &[Option<i64>]) -> Result<(), Error> {
    // Extract non-None values for sorting validation
    let positions: Vec<i64> = starts.iter().filter_map(|&x| x).collect();

    // Skip validation for empty or single-element vectors
    if positions.len() <= 1 {
        return Ok(());
    }

    // Always check ascending order regardless of alignment type
    let is_sorted_ascending = positions.windows(2).all(|w| w[0] <= w[1]);
    if !is_sorted_ascending {
        let sample_positions: Vec<i64> = positions.iter().take(5).copied().collect();
        return Err(Error::InvalidSorting(format!(
            "Expected ascending order, but got first {} positions: {:?}{}",
            sample_positions.len(),
            sample_positions,
            if positions.len() > 5 { "..." } else { "" }
        )));
    }

    Ok(())
}

/// Reconstruct `BaseMods` from condensed `mod_table` format
fn reconstruct_base_mods(
    mod_table: &[ModTableEntry],
    alignment_type: ReadState,
    seq_len: u64,
) -> Result<BaseMods, Error> {
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

        // Validate sorting after starts vector is populated
        validate_starts_sorting(&starts)?;

        // Calculate ends: starts + 1 where available, None otherwise
        #[expect(
            clippy::arithmetic_side_effects,
            reason = "overflow error impossible as genomic coords do not exceed (2^64-1)/2"
        )]
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

        // Marks if a read is reverse aligned
        let is_reverse = matches!(
            alignment_type,
            ReadState::PrimaryRev | ReadState::SecondaryRev | ReadState::SupplementaryRev
        );

        let ranges = Ranges {
            starts,
            ends,
            lengths,
            reference_starts,
            reference_ends,
            reference_lengths,
            qual,
            seq_len: seq_len.try_into()?,
            reverse: is_reverse,
        };

        let strand = if entry.is_strand_plus { '+' } else { '-' };

        base_mods.push(BaseMod {
            modified_base: entry.base as u8,
            strand,
            modification_type: entry.mod_code.val(),
            ranges,
            record_is_reverse: is_reverse,
        });
    }

    base_mods.sort();

    // Check for duplicate strand, modification_type combinations
    let mut seen_combinations = HashSet::new();
    for base_mod in &base_mods {
        let combination = (base_mod.strand, base_mod.modification_type);
        if seen_combinations.contains(&combination) {
            return Err(Error::InvalidDuplicates(format!(
                "Duplicate strand '{}' and modification_type '{}' combination found",
                base_mod.strand, base_mod.modification_type
            )));
        }
        let _: bool = seen_combinations.insert(combination);
    }

    Ok(BaseMods { base_mods })
}

#[cfg(test)]
mod test_error_handling {
    use super::*;

    #[test]
    fn test_set_read_state_not_implemented_error() {
        for flag_value in 0..4096u16 {
            let mut record = Record::new();
            record.set_flags(flag_value);
            let curr_read = CurrRead::default();
            // * first line below is usual primary, secondary, supplementary,
            //   unmapped with reversed set or not with the first 3 categories
            // * second line is if unmapped is set with a mapped state, or
            //   secondary and supplementary are both set
            // * third line are flags we have chosen to exclude from our program.
            match (flag_value, curr_read.set_read_state(&record)) {
                (0 | 4 | 16 | 256 | 272 | 2048 | 2064, Ok(_))
                | (
                    20 | 260 | 276 | 2052 | 2068 | 2304 | 2308 | 2320 | 2324,
                    Err(Error::UnknownAlignState),
                )
                | (_, Err(Error::NotImplementedError(_))) => {}
                (_, _) => unreachable!(),
            }
        }
    }

    #[test]
    #[should_panic(expected = "InvalidDuplicates")]
    fn test_reconstruct_base_mods_detects_duplicates() {
        // Create a test case with duplicate strand and modification_type combinations
        // i.e. we have two entries for T+T below with different data
        let mod_entries = vec![
            ModTableEntry {
                base: 'T',
                is_strand_plus: true,
                mod_code: ModChar::new('T'),
                implicit: false,
                data: vec![(0, 0, 10)],
            },
            ModTableEntry {
                base: 'T',
                is_strand_plus: true,
                mod_code: ModChar::new('T'),
                implicit: false,
                data: vec![(1, 1, 20)],
            },
        ];

        // Test reconstruct_base_mods with duplicate entries
        let _: BaseMods = reconstruct_base_mods(&mod_entries, ReadState::PrimaryFwd, 10).unwrap();
    }

    #[test]
    fn test_reconstruct_base_mods_accepts_unique_combinations() {
        // Create a valid case with different combinations
        // First entry: T+ modification
        // Second entry: C+ modification (different base, same strand)
        // Third entry: T- modification (same base, different strand)
        let mod_entries = vec![
            ModTableEntry {
                base: 'T',
                is_strand_plus: true,
                mod_code: ModChar::new('T'),
                implicit: false,
                data: vec![(0, 0, 10)],
            },
            ModTableEntry {
                base: 'C',
                is_strand_plus: true,
                mod_code: ModChar::new('m'),
                implicit: false,
                data: vec![(1, 1, 20)],
            },
            ModTableEntry {
                base: 'T',
                is_strand_plus: false,
                mod_code: ModChar::new('T'),
                implicit: false,
                data: vec![(2, 2, 30)],
            },
        ];

        // This should succeed since all combinations are unique
        let _: BaseMods = reconstruct_base_mods(&mod_entries, ReadState::PrimaryFwd, 10).unwrap();
    }
}

#[cfg(test)]
mod test_serde {
    use super::*;
    use crate::nanalogue_bam_reader;
    use indoc::indoc;
    use rust_htslib::bam::Read;

    #[test]
    fn test_first_record_serde() -> Result<(), Error> {
        // Read the first record from the example BAM file
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
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
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
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
        let json_str = indoc! {r"
            {
            }"};

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

    #[test]
    fn test_example_1_unmapped() -> Result<(), Error> {
        // Read the fourth record from the example BAM file (this is the unmapped read)
        let fourth_record = {
            let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
            let mut records = reader.records();
            for _ in 0..3 {
                let _: Record = records.next().unwrap()?;
            }
            records.next().unwrap()?
        };

        // Create CurrRead with alignment and modification data
        let curr_read = CurrRead::default()
            .try_from_only_alignment(&fourth_record)?
            .set_mod_data(&fourth_record, ThresholdState::GtEq(0), 0)?;

        let actual_json: serde_json::Value = serde_json::to_value(&curr_read)?;

        let expected_json_str = indoc! {r#"
            {
              "alignment_type": "unmapped",
              "mod_table": [
                {
                  "base": "G",
                  "is_strand_plus": false,
                  "mod_code": "7200",
                  "implicit": false,
                  "data": [
                    [28, -1, 0],
                    [29, -1, 0],
                    [30, -1, 0],
                    [32, -1, 0],
                    [43, -1, 77],
                    [44, -1, 0]
                  ]
                },
                {
                  "base": "T",
                  "is_strand_plus": true,
                  "mod_code": "T",
                  "implicit": false,
                  "data": [
                    [3, -1, 221],
                    [8, -1, 242],
                    [27, -1, 0],
                    [39, -1, 47],
                    [47, -1, 239]
                  ]
                }
              ],
              "read_id": "a4f36092-b4d5-47a9-813e-c22c3b477a0c",
              "seq_len": 48
            }"#};

        let expected_json: serde_json::Value = serde_json::from_str(expected_json_str)?;

        // Compare expected and actual outputs
        assert_eq!(actual_json, expected_json);

        Ok(())
    }

    #[test]
    #[should_panic(expected = "invalid alignment coordinates")]
    fn test_invalid_align_coords_unmapped_with_reference_positions() {
        let invalid_json = r#"{
            "mod_table": [
                {
                    "data": [[2, 3, 200]]
                }
            ],
            "seq_len": 10
        }"#;

        // Deserialize JSON to CurrRead - this should panic with InvalidAlignCoords
        let _: CurrRead<AlignAndModData> = serde_json::from_str(invalid_json).unwrap();
    }

    #[test]
    #[should_panic(expected = "invalid sequence length")]
    fn test_invalid_sequence_length() {
        let invalid_json = r#"{
            "mod_table": [
                {
                    "data": [[20, -1, 200]]
                }
            ],
            "seq_len": 10
        }"#;

        // Deserialize JSON to CurrRead - this should panic with InvalidSeqLength
        let _: CurrRead<AlignAndModData> = serde_json::from_str(invalid_json).unwrap();
    }

    #[test]
    #[should_panic(expected = "cannot deserialize when implicit")]
    fn test_invalid_implicit() {
        let invalid_json = r#"{
            "mod_table": [
                {
                    "implicit": true
                }
            ]
        }"#;

        // Deserialize JSON to CurrRead - this should panic with InvalidSeqLength
        let _: CurrRead<AlignAndModData> = serde_json::from_str(invalid_json).unwrap();
    }

    #[test]
    #[should_panic(expected = "invalid value")]
    fn test_invalid_sequence_coordinate() {
        let invalid_json = r#"{
            "alignment_type": "primary_forward",
            "alignment": {
                "start": 0,
                "end": 30
            },
            "mod_table": [
                {
                    "data": [[-1, 1, 200]]
                }
            ],
            "seq_len": 30
        }"#;

        // Deserialize JSON to CurrRead - this should panic with InvalidSeqLength
        let _: CurrRead<AlignAndModData> = serde_json::from_str(invalid_json).unwrap();
    }

    #[test]
    #[should_panic(expected = "invalid alignment coordinates")]
    fn test_invalid_alignment_coordinates() {
        let invalid_json = r#"{
            "alignment_type": "primary_forward",
            "alignment": {
                "start": 10,
                "end": 25
            },
            "mod_table": [
                {
                    "data": [[0, 10, 200], [1, 20, 180], [2, 30, 220]]
                }
            ],
            "seq_len": 3
        }"#;

        // Deserialize JSON to CurrRead - this should panic with InvalidAlignCoords
        let _: CurrRead<AlignAndModData> = serde_json::from_str(invalid_json).unwrap();
    }

    #[test]
    #[should_panic(expected = "invalid sorting")]
    fn test_invalid_sorting_forward_alignment() {
        let invalid_json = r#"{
            "alignment_type": "primary_forward",
            "alignment": {
                "start": 10,
                "end": 40
            },
            "mod_table": [
                {
                    "data": [[0, 10, 200], [2, 30, 180], [1, 20, 220]]
                }
            ],
            "seq_len": 3
        }"#;

        // Deserialize JSON to CurrRead - this should panic with InvalidSorting
        let _: CurrRead<AlignAndModData> = serde_json::from_str(invalid_json).unwrap();
    }

    #[test]
    #[should_panic(expected = "invalid sorting")]
    fn test_invalid_sorting_reverse_alignment() {
        let invalid_json = r#"{
            "alignment_type": "primary_reverse",
            "alignment": {
                "start": 10,
                "end": 40
            },
            "mod_table": [
                {
                    "data": [[2, 30, 180], [1, 20, 220], [0, 10, 200]]
                }
            ],
            "seq_len": 3
        }"#;

        // Deserialize JSON to CurrRead - this should panic with InvalidSorting
        let _: CurrRead<AlignAndModData> = serde_json::from_str(invalid_json).unwrap();
    }
}
