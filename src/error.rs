//! # Error
//!
//! Covers all errors in our module. These errors arise from us processing
//! and calculating data associated with DNA molecules, their alignments to
//! reference genomes, modification information on them, and other miscellaneous
//! information. We convert errors from other packages to this error type so that
//! error handling in our package becomes easier.

use crate::F32Bw0and1;
use std::fmt;
use std::io;
use std::num::{ParseFloatError, ParseIntError, TryFromIntError};
use std::str::Utf8Error;
use std::string::FromUtf8Error;
use thiserror::Error;

/// Enum that covers errors in our module.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Error {
    /// Alignment of sequence is not known
    #[error("unknown alignment state")]
    UnknownAlignState,

    /// Failure upon extracting sequence length
    #[error("invalid sequence length")]
    InvalidSeqLength,

    /// Failure upon extracting or calculating alignment length of molecule
    #[error("invalid alignment length")]
    InvalidAlignLength,

    /// Contig and start of alignment of molecule are invalid
    #[error("invalid contig and start")]
    InvalidContigAndStart,

    /// Alignment coordinates (contig/start/end) are invalid.
    #[error("invalid alignment coordinates (contig/start/end)")]
    InvalidAlignCoords,

    /// Modification coordinates are invalid
    #[error("invalid mod coordinates")]
    InvalidModCoords,

    /// Modification probabilities are invalid
    #[error("invalid mod probabilities")]
    InvalidModProbs,

    /// Sequence is invalid
    #[error("invalid sequence")]
    InvalidSeq,

    /// Read id of molecule is invalid
    #[error("invalid read id")]
    InvalidReadID,

    /// Modification type is invalid. Mod types are indicated in
    /// mod BAM files like so: ...C+m... where C is the base and
    /// m is the modification type, in this case methylation.
    #[error("invalid mod type")]
    InvalidModType,

    /// Modification type is empty
    #[error("empty mod type")]
    EmptyModType,

    /// Some error from the rust htslib library we use to read BAM files
    #[error(
        "rust_htslib error: `{0}` \nIf piping in a samtools view command, please include header with -h in samtools. "
    )]
    RustHtslibError(#[from] rust_htslib::errors::Error),

    /// Error upon conversion from integer
    #[error("integer conversion error: `{0}`")]
    IntConversionError(#[from] TryFromIntError),

    /// Error involving string conversion
    #[error("error involving string conversion: `{0}`")]
    StringConversionError(#[from] Utf8Error),

    /// Error converting from UTF-8 bytes to string
    #[error("UTF-8 conversion error: `{0}`")]
    Utf8ConversionError(#[from] FromUtf8Error),

    /// Error parsing JSON
    #[error("JSON parsing error: `{0}`")]
    JsonParseError(#[from] serde_json::Error),

    /// `OrdPair` is an ordered pair, which can be obtained from
    /// a string of the correct format. This error says string
    /// conversion failed.
    #[error("ordered pair conversion error: `{0}`")]
    OrdPairConversionError(String),

    /// Problem parsing integers
    #[error("integer parsing error: `{0}`")]
    IntParseError(#[from] ParseIntError),

    /// Problem parsing floats
    #[error("float parsing error: `{0}`")]
    FloatParseError(#[from] ParseFloatError),

    /// Generic Input-Output error
    #[error("input output error: `{0}`")]
    InputOutputError(#[from] io::Error),

    /// Generic formatting error
    #[error("formatting error: `{0}`")]
    FormattingError(#[from] fmt::Error),

    /// Problem reading or parsing CSV files
    #[error("error parsing csv: `{0}`")]
    CsvError(#[from] csv::Error),

    /// Error when unexpected duplicates are seen
    #[error("duplicates detected: `{0}`")]
    InvalidDuplicates(String),

    /// Generic error used when program hits an invalid state
    #[error("`{0}`")]
    InvalidState(String),

    /// Error while writing output
    #[error("error while writing output")]
    WriteOutputError,

    /// Generic not implemented error
    #[error("not implemented: `{0}`")]
    NotImplementedError(String),

    /// General error when ordering of items in some context is wrong.
    #[error("items in wrong order")]
    WrongOrder,

    /// Data not available
    #[error("data not available")]
    UnavailableData,

    /// Unknown error
    #[error("unknown error")]
    UnknownError,

    /// Read is unmapped, use this whenever some function
    /// meant for a mapped read is called on an unmapped read
    #[error("read is unmapped")]
    Unmapped,

    /// Zero values used where they should not be
    #[error("zero values not allowed")]
    Zero,

    /// Genomic region coordinates exceed contig boundaries
    #[error("invalid region '{region}': position {pos} exceeds contig length {contig_length}")]
    InvalidRegionError {
        /// The original region string provided by the user
        region: String,
        /// The position that exceeds the contig boundary
        pos: u64,
        /// The actual length of the contig
        contig_length: u64,
    },

    /// Cannot de-serialize implicit modBAM format
    #[error("cannot deserialize when implicit is set to true")]
    DeSerializeImplicit,

    /// Sorting validation failure
    #[error("invalid sorting: {0}")]
    InvalidSorting(String),

    /// Window density is below threshold
    #[error("window density {density} below threshold {threshold}")]
    WindowDensBelowThres {
        /// The density value that was below threshold
        density: F32Bw0and1,
        /// The threshold value
        threshold: F32Bw0and1,
    },

    /// Window does not contain any data
    #[error("window does not contain any data")]
    EmptyWindow,

    /// Data is not of sufficient size (e.g. in a window)
    #[error("data is not of sufficient size (e.g. in a window)")]
    InsufficientDataSize,

    /// Arithmetic error
    #[error("unanticipated arithmetic error e.g. overflow")]
    Arithmetic,
}
