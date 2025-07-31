use std::collections::HashMap;
use std::num::NonZeroU64;
use std::fmt;
use rust_htslib::{bam::record::Record, bam::ext::BamRecordExtensions};
use fibertools_rs::utils::basemods::{BaseMod, BaseMods};
use bedrs::prelude::*;

// Import from our crate
use crate::nanalogue_mm_ml_parser;
use crate::Error;

// A read can exist in seven states + one unknown state
#[derive(Debug, Clone, Default, Copy, PartialEq)]
pub enum ReadState {
    #[default]
    Unknown,
    PrimaryFwd,
    PrimaryRev,
    SecondaryFwd,
    SecondaryRev,
    SupplementaryFwd,
    SupplementaryRev,
    Unmapped,
}

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

// Three types of thresholds can be applied to modification data
#[derive(Debug, Clone, PartialEq)]
pub enum ThresholdState{
    Above(u8),
    Below(u8),
    Between(u8, u8),
}

impl Default for ThresholdState {
    fn default () -> Self {
        ThresholdState::Above(0)
    }
}

#[derive(Debug, Default, Clone, PartialEq)]
pub struct CurrRead {
    state: ReadState,
    read_id: Option<String>,
    seq_len: Option<u64>,
    align_len: Option<u64>,
    mods: Option<(Vec<BaseMod>, ThresholdState)>,
    contig_and_start: Option<(i32, u64)>,
}

impl CurrRead {
    pub fn new() -> Self {
        Default::default()
    }
    pub fn set_read_state(&mut self, record: &Record) -> Result<bool, Error> {
        // set read state
        match &self.state {
            ReadState::Unknown => {
                match (record.is_reverse(), record.is_unmapped(), 
                        record.is_secondary(), record.is_supplementary()) {
                    (false, true, false, false) => { self.state = ReadState::Unmapped; Ok(true) },
                    (true, false, false, false) => { self.state = ReadState::PrimaryRev; Ok(true) },
                    (true, false, true, false) => { self.state = ReadState::SecondaryRev; Ok(true) },
                    (false, false, true, false) => { self.state = ReadState::SecondaryFwd; Ok(true) },
                    (true, false, false, true) => { self.state = ReadState::SupplementaryRev; Ok(true) },
                    (false, false, false, true) => { self.state = ReadState::SupplementaryFwd; Ok(true) },
                    (false, false, false, false) => { self.state = ReadState::PrimaryFwd; Ok(true) },
                    _ => Err(Error::UnknownAlignState),
                }
            },
            _ => Err(Error::InvalidDuplicates("cannot set align state again!".to_string())),
        }
    }
    pub fn get_read_state(&self) -> ReadState {
        // get read state
        self.state
    }
    pub fn set_seq_len(&mut self, record: &Record) -> Result<Option<u64>, Error>{
        // set length of sequence from BAM record
        match &self.seq_len {
            Some(_) => Err(Error::InvalidDuplicates("cannot set sequence length again!".to_string())),
            None => {
                self.seq_len = Some(NonZeroU64::new(record.seq_len().try_into()?)
                    .ok_or(Error::InvalidSeqLength)?.get());
                Ok(self.seq_len)
            }
        }
    }
    pub fn set_align_len(&mut self, record: &Record) -> Result<Option<u64>, Error>{
        // set alignment length from BAM record if available
        match &self.align_len {
            Some(_) => Err(Error::InvalidDuplicates("cannot set alignment length again!".to_string())),
            None => {
                match self.state {
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
                    },
                }
            },
        }
    }
    pub fn set_contig_and_start(&mut self, record: &Record) -> Result<bool, Error> {
        match &self.contig_and_start{
            Some(_) => Err(Error::InvalidDuplicates("cannot set contig and start again!".to_string())),
            None => {
                match self.state {
                    ReadState::Unknown => Err(Error::UnknownAlignState),
                    ReadState::Unmapped => Ok(false),
                    _ => {
                        self.contig_and_start = Some((record.tid(), record.pos().try_into().unwrap()));
                        Ok(true)
                    },
                }
            }
        }
    }
    pub fn set_read_id(&mut self, record: &Record) -> Result<&str, Error>{
        // sets read id
        match &self.read_id {
            Some(_) => Err(Error::InvalidDuplicates("cannot set read id again!".to_string())),
            None => {
                match str::from_utf8(record.qname()) {
                    Ok(v) => {
                        self.read_id = Some(v.to_string());
                        self.get_read_id()
                    },
                    Err(_) => Err(Error::InvalidReadID),
                }
            },
        }
    }
    pub fn get_read_id(&self) -> Result<&str, Error> {
        // get read id
        match &self.read_id {
            None => Err(Error::InvalidState("read id not available".to_string())),
            Some(v) => Ok(v.as_str()),
        }
    }
    pub fn set_mod_data(&mut self, record: &Record, mod_thres: u8){
        let BaseMods { base_mods: v } = nanalogue_mm_ml_parser(record, mod_thres, None);
        self.mods = Some((v, ThresholdState::Above(mod_thres)));
    }
    pub fn set_mod_data_one_tag(&mut self, record: &Record, mod_thres: u8, mod_tag: char){
        let BaseMods { base_mods: v } = nanalogue_mm_ml_parser(record, mod_thres, Some(mod_tag));
        self.mods = Some((v, ThresholdState::Above(mod_thres)));
    }
    pub fn windowed_mod_data(&self, win_size: usize, slide_size: usize, tag_char: char) 
        -> Result<Option<Vec<f32>>, Error>{
        let mut result = Vec::<f32>::new();
        let mut is_track_seen: bool = false;
        if let Some((v, _)) = &self.mods {
            for k in v {
                match k {
                    BaseMod {
                        modified_base: _, strand: _, record_is_reverse: _,
                        modification_type: x, ranges: track,
                    } if *x == tag_char && !is_track_seen => {
                        is_track_seen = true;
                        let data = &track.qual;
                        if win_size > data.len(){
                            continue;
                        }
                        result = (0..=data.len() - win_size).step_by(slide_size)
                            .map(|i| {
                                let window_slice = &data[i..i + win_size];
                                let sum: f32 = window_slice.iter().map(|&val| val as f32).sum();
                                sum / (256.0 * win_size as f32)
                            }).collect();
                    },
                    BaseMod {
                        modified_base: _, strand: _, record_is_reverse: _,
                        modification_type: x, ranges: _,
                    } if *x == tag_char && is_track_seen => {
                        return Err(Error::NotImplementedError("cannot window on data on multiple tracks".to_string()));
                    },
                    _ => {},
                }
            }
        }
        if ! result.is_empty(){
            Ok(Some(result))
        } else {
            Ok(None)
        }
    }

    pub fn mod_count_per_mod(&self) -> Option<HashMap<char, u32>> {
        let mut output = HashMap::<char, u32>::new();
        match &self.mods {
            Some((v, _)) => {
                if v.is_empty() {
                    None
                } else {
                    for k in v {
                        let mod_count = k.ranges.qual.len() as u32;
                        output.entry(k.modification_type).and_modify(|e| *e += mod_count ).or_insert(mod_count);
                    }
                    Some(output)
                }
            }
            None => None,
        }
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

        output_string = output_string + "\t\"alignment_type\": \"" + &self.state.to_string() + "\",\n";

        if let Some((v, _)) = &self.mods {
            if !v.is_empty() {
                output_string += "\t\"mod_count\": \"";
                for k in v {
                    output_string += format!("{}{}{}:{};",
                        k.modified_base as char,
                        k.strand,
                        match k.modification_type {
                            w @ ('A'..='Z' | 'a'..='z') => w.to_string(),
                            w => format!("{}", w as u32),
                        },
                        k.ranges.qual.len()
                    ).as_str();
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

impl TryFrom<CurrRead> for StrandedBed3<i32, u64>{
    type Error = crate::Error;

    fn try_from(value: CurrRead) -> Result<Self, Self::Error> {
        match (value.state, value.align_len, value.contig_and_start) {
            (ReadState::Unknown | ReadState::Unmapped, _, _) => Err(Error::UnknownAlignState),
            (_, None, _)  => Err(Error::InvalidAlignLength),
            (_, _, None) => Err(Error::InvalidContigAndStart),
            (ReadState::PrimaryFwd | ReadState::SecondaryFwd | ReadState::SupplementaryFwd, Some(al), Some((cg, st))) => {
                Ok(StrandedBed3::new(cg, st, st + al, Strand::Forward))
            },
            (ReadState::PrimaryRev | ReadState::SecondaryRev | ReadState::SupplementaryRev, Some(al), Some((cg, st))) => {
                Ok(StrandedBed3::new(cg, st, st + al, Strand::Reverse))
            },
        }
    }
}
