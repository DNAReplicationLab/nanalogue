//! # Error
//!
//! Covers all errors in our module. These errors arise from us processing
//! and calculating data associated with DNA molecules, their alignments to
//! reference genomes, modification information on them, and other miscellaneous
//! information. We convert errors from other packages to this error type so that
//! error handling in our package becomes easier.

use crate::F32Bw0and1;
use derive_builder::UninitializedFieldError;
use std::char::TryFromCharError;
use std::fmt;
use std::io;
use std::num::{ParseFloatError, ParseIntError, TryFromIntError};
use std::str::Utf8Error;
use std::string::FromUtf8Error;

/// Shared CLI guidance for contig/region inspection failures.
const PEEK_HELP: &str =
    "In command line tool, use `nanalogue peek` to check contig names and lengths.";

/// Shared CLI guidance for alignment/header-related failures.
const PEEK_AND_HEADER_HELP: &str = "In command line tool, use `nanalogue peek` to check contig names and lengths.\n\
In command line tool, if piping in a samtools view command, please include header with -h in samtools.";

/// Truncate display strings to 100 characters.
fn trunc_100(s: &str) -> String {
    let truncated: String = s.chars().take(100).collect();
    if s.len() > truncated.len() {
        format!("{truncated}...")
    } else {
        truncated
    }
}

/// Enum that covers errors in our module.
///
/// Any error arising from our crate does not have the
/// suffix 'Error'. If we are deriving an error from
/// an error from another crate, and that has a suffix
/// 'Error', we have let it be in our naming of the error.
#[derive(Debug)]
#[non_exhaustive]
pub enum Error {
    /// Alignment of sequence is not known
    UnknownAlignState(String),
    /// Failure upon extracting sequence length
    InvalidSeqLength(String),
    /// Failure upon extracting or calculating alignment length of molecule
    InvalidAlignLength(String),
    /// Contig and start of alignment of molecule are invalid
    InvalidContigAndStart(String),
    /// Alignment coordinates (contig/start/end) are invalid.
    InvalidAlignCoords(String),
    /// Modification coordinates are invalid
    InvalidModCoords(String),
    /// Modification probabilities are invalid
    InvalidModProbs(String),
    /// Sequence is invalid
    InvalidSeq(String),
    /// Base is invalid (not A, G, C, T, or N)
    InvalidBase(String),
    /// Read id of molecule is invalid
    InvalidReadID(String),
    /// Modification type is invalid.
    InvalidModType(String),
    /// Modification type is empty
    EmptyModType(String),
    /// Some error from the rust htslib library we use to read BAM files
    RustHtslibError(Box<rust_htslib::errors::Error>),
    /// Error upon conversion from integer
    IntConversionError(Box<TryFromIntError>),
    /// Error involving string conversion
    StringConversionError(Box<Utf8Error>),
    /// Error converting from UTF-8 bytes to string
    Utf8ConversionError(Box<FromUtf8Error>),
    /// Error parsing JSON
    JsonParseError(Box<serde_json::Error>),
    /// `OrdPair` string conversion failed.
    OrdPairConversion(String),
    /// Problem parsing integers
    IntParseError(Box<ParseIntError>),
    /// Problem parsing floats
    FloatParseError(Box<ParseFloatError>),
    /// Generic Input-Output error
    InputOutputError(Box<io::Error>),
    /// Generic formatting error
    FormattingError(Box<fmt::Error>),
    /// Error when unexpected duplicates are seen
    InvalidDuplicates(String),
    /// Generic error used when program hits an invalid state
    InvalidState(String),
    /// Error while writing output
    WriteOutput(String),
    /// Generic not implemented error
    NotImplemented(String),
    /// General error when ordering of items in some context is wrong.
    WrongOrder(String),
    /// Data not available
    UnavailableData(String),
    /// Read is unmapped
    Unmapped(String),
    /// Zero values used where they should not be
    Zero(String),
    /// Zero sequence length
    ZeroSeqLen(String),
    /// Genomic region coordinates exceed contig boundaries
    InvalidRegion {
        /// The original region string provided by the user
        region: String,
        /// The position that exceeds the contig boundary
        pos: u64,
        /// The actual length of the contig
        contig_length: u64,
    },
    /// Sorting validation failure
    InvalidSorting(String),
    /// Window density is below threshold
    WindowDensBelowThres {
        /// The density value that was below threshold
        density: F32Bw0and1,
        /// The threshold value
        threshold: F32Bw0and1,
    },
    /// Window does not contain any data
    EmptyWindow(String),
    /// Data is not of sufficient size (e.g. in a window)
    InsufficientDataSize(String),
    /// Arithmetic error
    Arithmetic(String),
    /// Problem parsing items while building structs with Builder methods
    BuilderError(Box<UninitializedFieldError>),
    /// Problem validating items with a Builder method
    BuilderValidation(String),
    /// Problem parsing items while converting between DNA base representations
    FromCharError(Box<TryFromCharError>),
    /// Error from Polars during `DataFrame` construction or manipulation
    PolarsError(Box<polars::error::PolarsError>),
    /// Error simulating DNA sequences
    SimulateDNASeqCIGAREndProblem(String),
}

impl fmt::Display for Error {
    #[expect(
        clippy::pattern_type_mismatch,
        clippy::too_many_lines,
        reason = "many variants, so too many lines is fine. And matching patterns written better for readability"
    )]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownAlignState(v) => {
                write!(f, "unknown alignment state: `{}`", trunc_100(v))
            }
            Self::InvalidSeqLength(v) => {
                write!(f, "invalid sequence length: `{}`", trunc_100(v))
            }
            Self::InvalidAlignLength(v) => {
                write!(f, "invalid alignment length: `{}`", trunc_100(v))
            }
            Self::InvalidContigAndStart(v) => {
                write!(f, "invalid contig and/or start: `{}`", trunc_100(v))
            }
            Self::InvalidAlignCoords(v) => {
                write!(
                    f,
                    "invalid alignment coordinates (contig/start/end): `{}`.\n{PEEK_AND_HEADER_HELP}",
                    trunc_100(v)
                )
            }
            Self::InvalidModCoords(v) => write!(f, "invalid mod coordinates: `{}`", trunc_100(v)),
            Self::InvalidModProbs(v) => write!(f, "invalid mod probabilities: `{}`", trunc_100(v)),
            Self::InvalidSeq(v) => write!(f, "invalid sequence: `{}`", trunc_100(v)),
            Self::InvalidBase(v) => write!(f, "invalid base: `{}`", trunc_100(v)),
            Self::InvalidReadID(v) => write!(f, "invalid read id: `{}`", trunc_100(v)),
            Self::InvalidModType(v) => write!(f, "invalid mod type: `{}`", trunc_100(v)),
            Self::EmptyModType(v) => write!(f, "empty mod type: `{}`", trunc_100(v)),
            Self::RustHtslibError(v) => {
                write!(f, "rust_htslib error: `{v}`\n{PEEK_AND_HEADER_HELP}")
            }
            Self::IntConversionError(v) => write!(f, "integer conversion error: `{v}`"),
            Self::StringConversionError(v) => {
                write!(f, "error involving string conversion: `{v}`")
            }
            Self::Utf8ConversionError(v) => write!(f, "UTF-8 conversion error: `{v}`"),
            Self::JsonParseError(v) => write!(f, "JSON parsing error: `{v}`"),
            Self::OrdPairConversion(v) => {
                write!(f, "ordered pair conversion error: `{}`", trunc_100(v))
            }
            Self::IntParseError(v) => write!(f, "integer parsing error: `{v}`"),
            Self::FloatParseError(v) => write!(f, "float parsing error: `{v}`"),
            Self::InputOutputError(v) => write!(f, "input output error: `{v}`"),
            Self::FormattingError(v) => write!(f, "formatting error: `{v}`"),
            Self::InvalidDuplicates(v) => write!(f, "duplicates detected: `{}`", trunc_100(v)),
            Self::InvalidState(v) => write!(f, "`{}`", trunc_100(v)),
            Self::WriteOutput(v) => write!(f, "error while writing output: `{}`", trunc_100(v)),
            Self::NotImplemented(v) => write!(f, "not implemented: `{}`", trunc_100(v)),
            Self::WrongOrder(v) => write!(f, "items in wrong order: `{}`", trunc_100(v)),
            Self::UnavailableData(v) => write!(f, "data not available: `{}`", trunc_100(v)),
            Self::Unmapped(v) => write!(f, "read is unmapped: `{}`", trunc_100(v)),
            Self::Zero(v) => write!(f, "zero values not allowed: `{}`", trunc_100(v)),
            Self::ZeroSeqLen(v) => write!(f, "zero sequence length: `{}`", trunc_100(v)),
            Self::InvalidRegion {
                region,
                pos,
                contig_length,
            } => {
                write!(
                    f,
                    "invalid region '{region}': position {pos} exceeds contig length {contig_length}\n{PEEK_HELP}"
                )
            }
            Self::InvalidSorting(v) => write!(f, "invalid sorting: {}", trunc_100(v)),
            Self::WindowDensBelowThres { density, threshold } => {
                write!(f, "window density {density} below threshold {threshold}")
            }
            Self::EmptyWindow(v) => {
                write!(f, "window does not contain any data: `{}`", trunc_100(v))
            }
            Self::InsufficientDataSize(v) => {
                write!(
                    f,
                    "data is not of sufficient size (e.g. in a window): `{}`",
                    trunc_100(v)
                )
            }
            Self::Arithmetic(v) => {
                write!(
                    f,
                    "unanticipated arithmetic error e.g. overflow: `{}`",
                    trunc_100(v)
                )
            }
            Self::BuilderError(v) => write!(f, "building error, are you missing inputs?: `{v}`"),
            Self::BuilderValidation(v) => {
                write!(
                    f,
                    "building error, input validation failed: `{}`",
                    trunc_100(v)
                )
            }
            Self::FromCharError(v) => write!(f, "error converting between DNA bases: `{v}`"),
            Self::PolarsError(v) => write!(f, "Polars error: `{v}`"),
            Self::SimulateDNASeqCIGAREndProblem(v) => {
                write!(
                    f,
                    "Simulate DNA sequence error; problem at end of CIGAR: {}",
                    trunc_100(v)
                )
            }
        }
    }
}

impl std::error::Error for Error {
    #[expect(
        clippy::pattern_type_mismatch,
        reason = "matching borrowed enum variants keeps the implementation concise"
    )]
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::RustHtslibError(err) => Some(err.as_ref()),
            Self::IntConversionError(err) => Some(err.as_ref()),
            Self::StringConversionError(err) => Some(err.as_ref()),
            Self::Utf8ConversionError(err) => Some(err.as_ref()),
            Self::JsonParseError(err) => Some(err.as_ref()),
            Self::IntParseError(err) => Some(err.as_ref()),
            Self::FloatParseError(err) => Some(err.as_ref()),
            Self::InputOutputError(err) => Some(err.as_ref()),
            Self::FormattingError(err) => Some(err.as_ref()),
            Self::BuilderError(err) => Some(err.as_ref()),
            Self::FromCharError(err) => Some(err.as_ref()),
            Self::PolarsError(err) => Some(err.as_ref()),
            Self::UnknownAlignState(_)
            | Self::InvalidSeqLength(_)
            | Self::InvalidAlignLength(_)
            | Self::InvalidContigAndStart(_)
            | Self::InvalidAlignCoords(_)
            | Self::InvalidModCoords(_)
            | Self::InvalidModProbs(_)
            | Self::InvalidSeq(_)
            | Self::InvalidBase(_)
            | Self::InvalidReadID(_)
            | Self::InvalidModType(_)
            | Self::EmptyModType(_)
            | Self::OrdPairConversion(_)
            | Self::InvalidDuplicates(_)
            | Self::InvalidState(_)
            | Self::WriteOutput(_)
            | Self::NotImplemented(_)
            | Self::WrongOrder(_)
            | Self::UnavailableData(_)
            | Self::Unmapped(_)
            | Self::Zero(_)
            | Self::ZeroSeqLen(_)
            | Self::InvalidRegion { .. }
            | Self::InvalidSorting(_)
            | Self::WindowDensBelowThres { .. }
            | Self::EmptyWindow(_)
            | Self::InsufficientDataSize(_)
            | Self::Arithmetic(_)
            | Self::BuilderValidation(_)
            | Self::SimulateDNASeqCIGAREndProblem(_) => None,
        }
    }
}

impl From<rust_htslib::errors::Error> for Error {
    fn from(value: rust_htslib::errors::Error) -> Self {
        Self::RustHtslibError(Box::new(value))
    }
}

impl From<TryFromIntError> for Error {
    fn from(value: TryFromIntError) -> Self {
        Self::IntConversionError(Box::new(value))
    }
}

impl From<Utf8Error> for Error {
    fn from(value: Utf8Error) -> Self {
        Self::StringConversionError(Box::new(value))
    }
}

impl From<FromUtf8Error> for Error {
    fn from(value: FromUtf8Error) -> Self {
        Self::Utf8ConversionError(Box::new(value))
    }
}

impl From<serde_json::Error> for Error {
    fn from(value: serde_json::Error) -> Self {
        Self::JsonParseError(Box::new(value))
    }
}

impl From<ParseIntError> for Error {
    fn from(value: ParseIntError) -> Self {
        Self::IntParseError(Box::new(value))
    }
}

impl From<ParseFloatError> for Error {
    fn from(value: ParseFloatError) -> Self {
        Self::FloatParseError(Box::new(value))
    }
}

impl From<io::Error> for Error {
    fn from(value: io::Error) -> Self {
        Self::InputOutputError(Box::new(value))
    }
}

impl From<fmt::Error> for Error {
    fn from(value: fmt::Error) -> Self {
        Self::FormattingError(Box::new(value))
    }
}

impl From<UninitializedFieldError> for Error {
    fn from(value: UninitializedFieldError) -> Self {
        Self::BuilderError(Box::new(value))
    }
}

impl From<TryFromCharError> for Error {
    fn from(value: TryFromCharError) -> Self {
        Self::FromCharError(Box::new(value))
    }
}

impl From<polars::error::PolarsError> for Error {
    fn from(value: polars::error::PolarsError) -> Self {
        Self::PolarsError(Box::new(value))
    }
}
