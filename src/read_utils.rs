use std::collections::HashMap;
use std::num::NonZeroU64;
use std::fmt;
use rust_htslib::{bam::record::Record, bam::ext::BamRecordExtensions};
use fibertools_rs::utils::basemods::{BaseMod, BaseMods};
use bio_types::genome::AbstractInterval;

// import from our crate
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

#[readonly::make]
#[derive(Debug, Default, Clone, PartialEq)]
pub struct CurrRead {
    state: ReadState,
    read_id: Option<String>,
    seq_len: Option<u64>,
    align_len: Option<u64>,
    mods: Option<Vec<BaseMod>>,
    mod_thres: Option<u8>,
    contig_and_start: Option<(String, u64)>,
}

impl CurrRead {
    pub fn new() -> Self {
        Default::default()
    }
    pub fn set_read_state(&mut self, record: &Record) -> Result<bool, Error> {
        match self.state {
            ReadState::Unknown => {},
            _ => Err(Error::InvalidDuplicates("cannot set align state again!".to_string()))?,
        };
        match (record.is_reverse(), record.is_unmapped(), 
                record.is_secondary(), record.is_supplementary()) {
            (false, true, false, false) => self.state = ReadState::Unmapped,
            (true, false, false, false) => self.state = ReadState::PrimaryRev,
            (true, false, true, false) => self.state = ReadState::SecondaryRev,
            (false, false, true, false) => self.state = ReadState::SecondaryFwd,
            (true, false, false, true) => self.state = ReadState::SupplementaryRev,
            (false, false, false, true) => self.state = ReadState::SupplementaryFwd,
            (false, false, false, false) => self.state = ReadState::PrimaryFwd,
            _ => self.state = ReadState::Unknown,
        };
        match self.state {
            ReadState::Unknown => Err(Error::UnknownAlignState),
            _ => Ok(true),
        }
    }
    pub fn get_read_state(&self) -> ReadState {
        self.state
    }
    pub fn set_seq_len(&mut self, record: &Record) -> Result<Option<u64>, Error>{
        // get length of sequence
        self.seq_len = Some(NonZeroU64::new(record.seq_len().try_into()?)
            .ok_or(Error::InvalidSeqLength)?.get());
        Ok(self.seq_len)
    }
    pub fn set_align_len(&mut self, record: &Record) -> Result<Option<u64>, Error>{
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
    }
    pub fn set_contig_and_start(&mut self, record: &Record) -> Result<bool, Error> {
        match self.state {
            ReadState::Unknown => Err(Error::UnknownAlignState),
            ReadState::Unmapped => Ok(false),
            _ => {
                self.contig_and_start = Some((record.contig().to_string(),
                    record.pos().try_into().unwrap()));
                Ok(true)
            },
        }
    }
    pub fn set_read_id(&mut self, record: &Record) -> Result<bool, Error>{
        // get read id from BAM and set it
        let qname: String = match str::from_utf8(record.qname()) {
            Ok(v) => v.to_string(),
            Err(_) => String::from(""),
        };
        if ! qname.is_empty() {
            self.read_id = Some(qname.clone());
            Ok(true)
        } else {
            Err(Error::InvalidReadID)
        }
    }
    pub fn get_read_id(&self) -> Result<&str, Error> {
        match &self.read_id {
            None => Err(Error::InvalidState("read id not available".to_string())),
            Some(v) => Ok(v.as_str()),
        }
    }
    pub fn set_mod_data(&mut self, record: &Record, mod_thres: u8){
        let BaseMods { base_mods: v } = nanalogue_mm_ml_parser(record, mod_thres, None);
        self.mods = Some(v);
        self.mod_thres = Some(mod_thres);
    }
    pub fn set_mod_data_one_tag(&mut self, record: &Record, mod_thres: u8, mod_tag: char){
        let BaseMods { base_mods: v } = nanalogue_mm_ml_parser(record, mod_thres, Some(mod_tag));
        self.mods = Some(v);
        self.mod_thres = Some(mod_thres);
    }
    pub fn windowed_mod_data(&self, win_size: usize, slide_size: usize) -> Result<Vec<f32>, Error>{
        if let Some(v) = &self.mods {
            if v.len() > 1 {
                Err(Error::NotImplementedError("Cannot window on mod data on multiple tracks".to_string()))
            } else if let Some(track) = v.first() {
                let data = &track.ranges.qual;
                if win_size > data.len(){
                    return Ok(Vec::new());
                }
                let result = (0..=data.len() - win_size).step_by(slide_size)
                    .map(|i| {
                        let window_slice = &data[i..i + win_size];
                        let sum: f32 = window_slice.iter().map(|&val| val as f32).sum();
                        sum / (256.0 * win_size as f32)
                    }).collect();
                Ok(result)
            } else {
                Ok(Vec::new())
            }
        } else {
            Ok(Vec::new())
        }
    }
    pub fn mod_count_per_mod(&self) -> Option<HashMap<char, u32>> {
        let mut output = HashMap::<char, u32>::new();
        match &self.mods {
            Some(v) => {
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

        if let Some(v) = &self.mods {
            if !v.is_empty() {
                output_string += "\t\"mod_count\": \"";
                for k in v {
                    output_string += format!("{}{}{}:{};",
                        k.modified_base as char,
                        k.strand,
                        match k.modification_type {
                            'A'..='Z' | 'a'..='z' => k.modification_type.to_string(),
                            _ => format!("{}", k.modification_type as u32),
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
