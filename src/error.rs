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
    let mut chars = s.chars();
    let truncated: String = chars.by_ref().take(100).collect();
    if chars.next().is_some() {
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
        pos: u32,
        /// The actual length of the contig
        contig_length: u32,
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

#[cfg(test)]
mod tests {
    use super::*;
    use polars::error::PolarsError;
    use rust_htslib::bam;

    /// Helper to create a `rust_htslib` error suitable for conversion tests.
    fn sample_rust_htslib_error() -> rust_htslib::errors::Error {
        let missing_bam = std::env::temp_dir()
            .join(format!("nanalogue-error-test-{}", std::process::id()))
            .join("missing.bam");
        match bam::Reader::from_path(&missing_bam) {
            Ok(_) => unreachable!("expected missing BAM path to produce a rust_htslib error"),
            Err(error) => error,
        }
    }

    /// Helper to create a `TryFromIntError` suitable for conversion tests.
    fn sample_try_from_int_error() -> TryFromIntError {
        match u8::try_from(u16::MAX) {
            Ok(_) => unreachable!("expected integer conversion to fail"),
            Err(error) => error,
        }
    }

    /// Helper to create a `Utf8Error` suitable for conversion tests.
    fn sample_utf8_error() -> Utf8Error {
        let bytes = [u8::MAX];
        match std::str::from_utf8(&bytes) {
            Ok(_) => unreachable!("expected invalid UTF-8 bytes to fail"),
            Err(error) => error,
        }
    }

    /// Helper to create a `FromUtf8Error` suitable for conversion tests.
    fn sample_from_utf8_error() -> FromUtf8Error {
        match String::from_utf8(vec![0xFF]) {
            Ok(_) => unreachable!("expected invalid UTF-8 bytes to fail"),
            Err(error) => error,
        }
    }

    /// Helper to create a JSON parse error suitable for conversion tests.
    fn sample_json_parse_error() -> serde_json::Error {
        match serde_json::from_str::<serde_json::Value>("{") {
            Ok(_) => unreachable!("expected invalid JSON to fail"),
            Err(error) => error,
        }
    }

    /// Helper to create an integer parse error suitable for conversion tests.
    fn sample_int_parse_error() -> ParseIntError {
        match "abc".parse::<u32>() {
            Ok(_) => unreachable!("expected invalid integer parsing to fail"),
            Err(error) => error,
        }
    }

    /// Helper to create a float parse error suitable for conversion tests.
    fn sample_float_parse_error() -> ParseFloatError {
        match "abc".parse::<f32>() {
            Ok(_) => unreachable!("expected invalid float parsing to fail"),
            Err(error) => error,
        }
    }

    /// Helper to create a `TryFromCharError` suitable for conversion tests.
    fn sample_try_from_char_error() -> TryFromCharError {
        match u8::try_from('💥') {
            Ok(_) => unreachable!("expected out-of-range char conversion to fail"),
            Err(error) => error,
        }
    }

    fn assert_string_variant_display(error: &Error, expected_prefix: &str, expected_payload: &str) {
        let rendered = error.to_string();
        assert!(
            rendered.starts_with(expected_prefix),
            "expected prefix `{expected_prefix}` in `{rendered}`"
        );
        assert!(
            rendered.contains(expected_payload),
            "expected payload `{expected_payload}` in `{rendered}`"
        );
    }

    fn assert_backticked(error: &Error, expected_prefix: &str) {
        assert_string_variant_display(error, expected_prefix, "`payload`");
    }

    fn assert_plain(error: &Error, expected_prefix: &str, expected_payload: &str) {
        assert_string_variant_display(error, expected_prefix, expected_payload);
        let rendered = error.to_string();
        assert!(
            !rendered.contains(&format!("`{expected_payload}`")),
            "did not expect backticks around `{expected_payload}` in `{rendered}`"
        );
    }

    #[test]
    fn string_backed_variants_display_expected_prefixes() {
        macro_rules! assert_backticked_cases {
            ($( $error:expr => $prefix:literal; )+) => {
                $(assert_backticked(&$error, $prefix);)+
            };
        }

        assert_backticked_cases!(
            Error::UnknownAlignState(String::from("payload")) => "unknown alignment state: ";
            Error::InvalidSeqLength(String::from("payload")) => "invalid sequence length: ";
            Error::InvalidAlignLength(String::from("payload")) => "invalid alignment length: ";
            Error::InvalidContigAndStart(String::from("payload")) => "invalid contig and/or start: ";
            Error::InvalidModCoords(String::from("payload")) => "invalid mod coordinates: ";
            Error::InvalidModProbs(String::from("payload")) => "invalid mod probabilities: ";
            Error::InvalidSeq(String::from("payload")) => "invalid sequence: ";
            Error::InvalidBase(String::from("payload")) => "invalid base: ";
            Error::InvalidReadID(String::from("payload")) => "invalid read id: ";
            Error::InvalidModType(String::from("payload")) => "invalid mod type: ";
            Error::EmptyModType(String::from("payload")) => "empty mod type: ";
            Error::OrdPairConversion(String::from("payload")) => "ordered pair conversion error: ";
            Error::InvalidDuplicates(String::from("payload")) => "duplicates detected: ";
            Error::WriteOutput(String::from("payload")) => "error while writing output: ";
            Error::NotImplemented(String::from("payload")) => "not implemented: ";
            Error::WrongOrder(String::from("payload")) => "items in wrong order: ";
            Error::UnavailableData(String::from("payload")) => "data not available: ";
            Error::Unmapped(String::from("payload")) => "read is unmapped: ";
            Error::Zero(String::from("payload")) => "zero values not allowed: ";
            Error::ZeroSeqLen(String::from("payload")) => "zero sequence length: ";
            Error::EmptyWindow(String::from("payload")) => "window does not contain any data: ";
            Error::InsufficientDataSize(String::from("payload")) => "data is not of sufficient size (e.g. in a window): ";
            Error::Arithmetic(String::from("payload")) => "unanticipated arithmetic error e.g. overflow: ";
            Error::BuilderValidation(String::from("payload")) => "building error, input validation failed: ";
        );

        assert_eq!(
            Error::InvalidState(String::from("payload")).to_string(),
            "`payload`"
        );
        assert_plain(
            &Error::InvalidSorting(String::from("payload")),
            "invalid sorting: payload",
            "payload",
        );
        assert_plain(
            &Error::SimulateDNASeqCIGAREndProblem(String::from("payload")),
            "Simulate DNA sequence error; problem at end of CIGAR: payload",
            "payload",
        );
    }

    #[test]
    fn trunc_100_preserves_boundary_and_truncates_unicode() {
        assert_eq!(trunc_100(""), "");
        assert_eq!(trunc_100("x"), "x");

        let exact = "é".repeat(100);
        assert_eq!(trunc_100(&exact), exact);

        let longer = "界".repeat(101);
        let expected = format!("{}...", "界".repeat(100));
        assert_eq!(trunc_100(&longer), expected);
        assert_eq!(
            Error::UnknownAlignState(longer).to_string(),
            format!("unknown alignment state: `{expected}`")
        );
    }

    #[test]
    fn display_includes_peek_guidance_for_alignment_and_region_errors() {
        let align = Error::InvalidAlignCoords(String::from("bad region")).to_string();
        assert_eq!(
            align,
            format!(
                "invalid alignment coordinates (contig/start/end): `bad region`.\n{PEEK_AND_HEADER_HELP}"
            )
        );

        let region = Error::InvalidRegion {
            region: String::from("chr1:1-1000"),
            pos: 1000,
            contig_length: 100,
        }
        .to_string();
        assert_eq!(
            region,
            format!(
                "invalid region 'chr1:1-1000': position 1000 exceeds contig length 100\n{PEEK_HELP}"
            )
        );
    }

    #[test]
    fn display_formats_rust_htslib_error() {
        let rust_htslib = Error::from(sample_rust_htslib_error()).to_string();
        assert!(
            rust_htslib.starts_with("rust_htslib error: `"),
            "expected rust_htslib prefix in `{rust_htslib}`"
        );
        assert!(
            rust_htslib.ends_with(PEEK_AND_HEADER_HELP),
            "expected shared help text in `{rust_htslib}`"
        );
    }

    #[test]
    fn display_formats_conversion_and_parse_errors() {
        let int_conversion = Error::from(sample_try_from_int_error()).to_string();
        assert!(
            int_conversion.starts_with("integer conversion error: `"),
            "expected integer conversion prefix in `{int_conversion}`"
        );

        let string_conversion = Error::from(sample_utf8_error()).to_string();
        assert!(
            string_conversion.starts_with("error involving string conversion: `"),
            "expected string conversion prefix in `{string_conversion}`"
        );

        let utf8_conversion = Error::from(sample_from_utf8_error()).to_string();
        assert!(
            utf8_conversion.starts_with("UTF-8 conversion error: `"),
            "expected UTF-8 conversion prefix in `{utf8_conversion}`"
        );

        let json_parse = Error::from(sample_json_parse_error()).to_string();
        assert!(
            json_parse.starts_with("JSON parsing error: `"),
            "expected JSON parsing prefix in `{json_parse}`"
        );

        let int_parse = Error::from(sample_int_parse_error()).to_string();
        assert!(
            int_parse.starts_with("integer parsing error: `"),
            "expected integer parsing prefix in `{int_parse}`"
        );

        let float_parse = Error::from(sample_float_parse_error()).to_string();
        assert!(
            float_parse.starts_with("float parsing error: `"),
            "expected float parsing prefix in `{float_parse}`"
        );
    }

    #[test]
    fn display_formats_io_and_formatting_errors() {
        let io_error = Error::from(io::Error::other("boom")).to_string();
        assert_eq!(io_error, "input output error: `boom`");

        let formatting_error = Error::from(fmt::Error).to_string();
        assert_eq!(
            formatting_error,
            "formatting error: `an error occurred when formatting an argument`"
        );
    }

    #[test]
    fn display_formats_builder_polars_and_density_errors() {
        let builder_error = Error::from(UninitializedFieldError::new("bam_path")).to_string();
        assert_eq!(
            builder_error,
            "building error, are you missing inputs?: `Field not initialized: bam_path`"
        );

        let from_char = Error::from(sample_try_from_char_error()).to_string();
        assert!(
            from_char.starts_with("error converting between DNA bases: `"),
            "expected char conversion prefix in `{from_char}`"
        );

        let polars_error =
            Error::from(PolarsError::ComputeError(String::from("boom").into())).to_string();
        assert!(
            polars_error.starts_with("Polars error: `"),
            "expected Polars prefix in `{polars_error}`"
        );

        let density = Error::WindowDensBelowThres {
            density: F32Bw0and1::new(0.25).expect("density should be valid"),
            threshold: F32Bw0and1::new(0.75).expect("threshold should be valid"),
        }
        .to_string();
        assert_eq!(density, "window density 0.25 below threshold 0.75");
    }

    fn assert_sources_present(cases: &[Error]) {
        for error in cases {
            assert!(
                std::error::Error::source(error).is_some(),
                "expected source for `{error}`"
            );
        }
    }

    #[test]
    fn wrapped_error_variants_expose_sources() {
        let cases = [
            Error::from(sample_rust_htslib_error()),
            Error::from(sample_try_from_int_error()),
            Error::from(sample_utf8_error()),
            Error::from(sample_from_utf8_error()),
            Error::from(sample_json_parse_error()),
            Error::from(sample_int_parse_error()),
            Error::from(sample_float_parse_error()),
            Error::from(io::Error::other("boom")),
            Error::from(fmt::Error),
            Error::from(UninitializedFieldError::new("bam_path")),
            Error::from(sample_try_from_char_error()),
            Error::from(PolarsError::ComputeError(String::from("boom").into())),
        ];

        assert_sources_present(&cases);
    }

    fn assert_sources_absent(cases: &[Error]) {
        for error in cases {
            assert!(
                std::error::Error::source(error).is_none(),
                "did not expect source for `{error}`"
            );
        }
    }

    #[test]
    fn native_error_variants_have_no_source() {
        // Keep this list in sync with `source()` and the `Error` enum so new native
        // variants explicitly decide whether they should expose an underlying source.
        let cases = [
            Error::UnknownAlignState(String::from("payload")),
            Error::InvalidSeqLength(String::from("payload")),
            Error::InvalidAlignLength(String::from("payload")),
            Error::InvalidContigAndStart(String::from("payload")),
            Error::InvalidAlignCoords(String::from("payload")),
            Error::InvalidModCoords(String::from("payload")),
            Error::InvalidModProbs(String::from("payload")),
            Error::InvalidSeq(String::from("payload")),
            Error::InvalidBase(String::from("payload")),
            Error::InvalidReadID(String::from("payload")),
            Error::InvalidModType(String::from("payload")),
            Error::EmptyModType(String::from("payload")),
            Error::OrdPairConversion(String::from("payload")),
            Error::InvalidDuplicates(String::from("payload")),
            Error::InvalidState(String::from("payload")),
            Error::WriteOutput(String::from("payload")),
            Error::NotImplemented(String::from("payload")),
            Error::WrongOrder(String::from("payload")),
            Error::UnavailableData(String::from("payload")),
            Error::Unmapped(String::from("payload")),
            Error::Zero(String::from("payload")),
            Error::ZeroSeqLen(String::from("payload")),
            Error::InvalidRegion {
                region: String::from("chr1:1-1000"),
                pos: 1000,
                contig_length: 100,
            },
            Error::InvalidSorting(String::from("payload")),
            Error::WindowDensBelowThres {
                density: F32Bw0and1::new(0.25).expect("density should be valid"),
                threshold: F32Bw0and1::new(0.75).expect("threshold should be valid"),
            },
            Error::EmptyWindow(String::from("payload")),
            Error::InsufficientDataSize(String::from("payload")),
            Error::Arithmetic(String::from("payload")),
            Error::BuilderValidation(String::from("payload")),
            Error::SimulateDNASeqCIGAREndProblem(String::from("payload")),
        ];

        assert_sources_absent(&cases);
    }

    #[test]
    fn from_conversions_map_to_expected_variants() {
        assert!(matches!(
            Error::from(sample_rust_htslib_error()),
            Error::RustHtslibError(_)
        ));
        assert!(matches!(
            Error::from(sample_try_from_int_error()),
            Error::IntConversionError(_)
        ));
        assert!(matches!(
            Error::from(sample_utf8_error()),
            Error::StringConversionError(_)
        ));
        assert!(matches!(
            Error::from(sample_from_utf8_error()),
            Error::Utf8ConversionError(_)
        ));
        assert!(matches!(
            Error::from(sample_json_parse_error()),
            Error::JsonParseError(_)
        ));
        assert!(matches!(
            Error::from(sample_int_parse_error()),
            Error::IntParseError(_)
        ));
        assert!(matches!(
            Error::from(sample_float_parse_error()),
            Error::FloatParseError(_)
        ));
        assert!(matches!(
            Error::from(io::Error::other("boom")),
            Error::InputOutputError(_)
        ));
        assert!(matches!(Error::from(fmt::Error), Error::FormattingError(_)));
        assert!(matches!(
            Error::from(UninitializedFieldError::new("bam_path")),
            Error::BuilderError(_)
        ));
        assert!(matches!(
            Error::from(sample_try_from_char_error()),
            Error::FromCharError(_)
        ));
        assert!(matches!(
            Error::from(PolarsError::ComputeError(String::from("boom").into())),
            Error::PolarsError(_)
        ));
    }
}
