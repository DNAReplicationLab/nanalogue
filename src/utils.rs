//! # Utils implementing simple, shared datatypes
//!
//! Implements simple datatypes like an ordered pair or mod tag
//! used by all modules in our crate.

use crate::Error;
use rust_htslib::bam::FetchDefinition;
use serde::{Deserialize, Serialize};
use std::convert::From;
use std::fmt;
use std::fmt::Debug;
use std::ops::RangeInclusive;
use std::str::FromStr;

/// Implements test if a value is within some interval
pub trait Contains<T> {
    /// see if value is contained within
    fn contains(&self, val: &T) -> bool;
}

/// Implements test if an interval intersects with self
pub trait Intersects<T> {
    /// see if interval intersects with self
    fn intersects(&self, val: &T) -> bool;
}

/// Implements filter by coordinates on the reference genome.
pub trait FilterByRefCoords {
    /// filters by reference position i.e. all pos such that start <= pos < end
    /// are retained. does not use contig in filtering.
    fn filter_by_ref_pos(&mut self, _: i64, _: i64) {
        todo!()
    }
    /// filters by coordinates contained in FetchDefinition
    fn filter_by_ref_pos_fd(&mut self, _: FetchDefinition) {
        todo!()
    }
}

/// Datatype holding two values low, high such that low <= high is guaranteed at creation.
#[derive(Debug, Clone, Default, Copy, PartialOrd, PartialEq, Serialize, Deserialize)]
pub struct OrdPair<T: Clone + Copy + Debug + Default> {
    low: T,
    high: T,
}

impl<T: Clone + Copy + Debug + Default + PartialEq + PartialOrd> OrdPair<T> {
    /// Constructor with two values, will fail if ordering in input is not respected.
    ///
    /// ```should_panic
    /// use nanalogue_core::OrdPair;
    /// let x = OrdPair::<f32>::new(1.0,0.0)?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    /// ```
    /// # use nanalogue_core::OrdPair;
    /// let x = OrdPair::<f32>::new(0.0, 1.0)?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    pub fn new(low: T, high: T) -> Result<Self, Error> {
        if low <= high {
            Ok(OrdPair { low, high })
        } else {
            Err(Error::WrongOrder)
        }
    }
    /// Gets the low value
    ///
    /// ```
    /// use nanalogue_core::OrdPair;
    /// let x = OrdPair::<u8>::new(10,11).expect("no failure");
    /// assert_eq!(x.get_low(),10);
    /// ```
    pub fn get_low(&self) -> T {
        self.low
    }
    /// Gets the high value
    ///
    /// ```
    /// use nanalogue_core::OrdPair;
    /// let x = OrdPair::<u8>::new(10,11).expect("no failure");
    /// assert_eq!(x.get_high(),11);
    /// ```
    pub fn get_high(&self) -> T {
        self.high
    }
}

impl<T: Clone + Copy + Debug + Default + PartialEq + PartialOrd + FromStr> FromStr for OrdPair<T> {
    type Err = Error;

    /// Parse a string to obtain an Ordered Pair, return Error if cannot be done.
    fn from_str(val_str: &str) -> Result<Self, Self::Err> {
        macro_rules! parse_error {
            () => {
                Err(Error::OrdPairConversionError(
                    "Bad ordered pair inputs!".to_string(),
                ))
            };
        }
        let v: Vec<&str> = val_str.split(",").map(|s| s.trim()).collect();
        match v.len() {
            2 => {
                let Ok(low) = T::from_str(v[0]) else {
                    parse_error!()?
                };
                let Ok(high) = T::from_str(v[1]) else {
                    parse_error!()?
                };
                OrdPair::<T>::new(low, high)
            }
            _ => parse_error!(),
        }
    }
}

impl<T: Clone + Copy + Debug + Default + PartialEq + PartialOrd + FromStr> From<OrdPair<T>>
    for RangeInclusive<T>
{
    /// Convert the OrdPair into a RangeInclusive i.e. (start..=end)
    fn from(value: OrdPair<T>) -> Self {
        RangeInclusive::<T>::new(value.get_low(), value.get_high())
    }
}

impl<T: Clone + Copy + Debug + Default + PartialEq + PartialOrd + FromStr> Contains<T>
    for OrdPair<T>
{
    /// Check if the provided value is within the Range of the OrdPair
    fn contains(&self, val: &T) -> bool {
        RangeInclusive::<T>::from(*self).contains(val)
    }
}

impl<T: Clone + Copy + Debug + Default + fmt::Display + PartialEq + PartialOrd> fmt::Display
    for OrdPair<T>
{
    /// converts to string for display i.e. "low, high"
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}, {}", self.get_low(), self.get_high())
    }
}

/// Datatype holding a float (f32) between 0 and 1 (both inclusive) guaranteed at creation.
#[derive(Debug, Clone, Default, Copy, PartialOrd, PartialEq, Serialize, Deserialize)]
pub struct F32Bw0and1(f32);

impl F32Bw0and1 {
    /// Constructor, will fail if float is not between 0 and 1
    ///
    /// ```should_panic
    /// use nanalogue_core::Error;
    /// use nanalogue_core::F32Bw0and1;
    /// let x = F32Bw0and1::new(-0.1)?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    /// ```
    /// # use nanalogue_core::Error;
    /// # use nanalogue_core::F32Bw0and1;
    /// let x = F32Bw0and1::new(0.1)?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    pub fn new(val: f32) -> Result<Self, Error> {
        if (0.0..=1.0).contains(&val) {
            Ok(F32Bw0and1(val))
        } else {
            Err(Error::InvalidState("Num not b/w 0 and 1!".to_string()))
        }
    }
    /// Returns the value of the float.
    ///
    /// ```
    /// use nanalogue_core::F32Bw0and1;
    /// for y in vec![0.0,0.1,0.7,1.0]{
    ///     let x = F32Bw0and1::new(y.clone())?;
    ///     assert_eq!(x.val(), y);
    /// }
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    pub fn val(&self) -> f32 {
        self.0
    }
    /// Shortcut for 1.0
    pub fn one() -> Self {
        F32Bw0and1::new(1.0).expect("no error")
    }
    /// Shortcut for 0.0
    pub fn zero() -> Self {
        F32Bw0and1::new(0.0).expect("no error")
    }
    /// Converts from F32AbsValBelow1 using the absolute value
    pub fn abs_f32_abs_val_below_1(val: F32AbsValBelow1) -> Self {
        F32Bw0and1::new(f32::abs(val.val())).expect("no error")
    }
}

impl FromStr for F32Bw0and1 {
    type Err = Error;

    /// Parse a string to obtain float and then convert if b/w 0 and 1
    fn from_str(val_str: &str) -> Result<Self, Self::Err> {
        Self::new(f32::from_str(val_str)?)
    }
}

impl From<u8> for F32Bw0and1 {
    /// Convert from a u8 i.e. a number >= 0 and <= 255
    fn from(value: u8) -> Self {
        F32Bw0and1::new((value as f32) / (u8::MAX as f32 + 1.0)).expect("no F32 conversion error")
    }
}

impl fmt::Display for F32Bw0and1 {
    /// converts to string for display.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.val())
    }
}

/// Datatype holding a float (f32) between -1 and 1 (both inclusive) guaranteed at creation.
#[derive(Debug, Clone, Default, Copy, PartialOrd, PartialEq)]
pub struct F32AbsValBelow1(f32);

impl F32AbsValBelow1 {
    /// Constructor, will fail if float is not between -1 and 1
    ///
    /// ```should_panic
    /// use nanalogue_core::Error;
    /// use nanalogue_core::F32AbsValBelow1;
    /// let x = F32AbsValBelow1::new(-1.1)?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    /// ```should_panic
    /// # use nanalogue_core::Error;
    /// # use nanalogue_core::F32AbsValBelow1;
    /// let x = F32AbsValBelow1::new(1.1)?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    /// ```
    /// # use nanalogue_core::Error;
    /// # use nanalogue_core::F32AbsValBelow1;
    /// let x = F32AbsValBelow1::new(0.1)?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    /// ```
    /// # use nanalogue_core::Error;
    /// # use nanalogue_core::F32AbsValBelow1;
    /// let x = F32AbsValBelow1::new(-0.5)?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    pub fn new(val: f32) -> Result<Self, Error> {
        if (-1.0..=1.0).contains(&val) {
            Ok(F32AbsValBelow1(val))
        } else {
            Err(Error::InvalidState("Num not b/w -1 and 1!".to_string()))
        }
    }
    /// Returns the value of the float.
    ///
    /// ```
    /// use nanalogue_core::F32AbsValBelow1;
    /// for y in vec![0.0,0.1,-0.7,1.0,-1.0]{
    ///     let x = F32AbsValBelow1::new(y.clone())?;
    ///     assert_eq!(x.val(), y);
    /// }
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    pub fn val(&self) -> f32 {
        self.0
    }
}

impl FromStr for F32AbsValBelow1 {
    type Err = Error;

    /// Parse a string to obtain float and then convert if b/w -1 and 1
    fn from_str(val_str: &str) -> Result<Self, Self::Err> {
        Self::new(f32::from_str(val_str)?)
    }
}

impl From<F32Bw0and1> for F32AbsValBelow1 {
    /// Convert between the two types of floats
    fn from(value: F32Bw0and1) -> Self {
        F32AbsValBelow1::new(value.val()).expect("no F32 conversion error")
    }
}

impl fmt::Display for F32AbsValBelow1 {
    /// converts to string for display.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.val())
    }
}

/// Our struct to hold a modification tag.
/// The BAM file format uses the syntax base+mod_code in its ML tag
/// to show which modification is represented e.g. C+m, A+a, T+T, ...
/// This can be a letter or a number e.g. T+472232 represents BrdU as that is its CheBI code.
/// As we rely on a fibertools-rs data structure to store mod information (BaseMod),
/// which uses a char datatype to represent the mod code, representing a one letter
/// mod code is easy, but numbers need to be converted to char first.
/// Fortunately, rust's char datatype is almost equivalent to u32, so we just
/// convert numbers to char before storing.
/// NOTE: the above conversion has some problems e.g. A+a is equivalent to A+97 etc.
/// (as 97 is the ascii code of a),
/// and the rust char datatype does not allow a set of u32s somewhere between 55000
/// and 59000. We have chosen to live with this problem. I think the probability of
/// having a DNA modification with a CheBI code overlapping with ASCII values or
/// within this narrow range of values near 59000 is very small.
#[derive(Debug, Clone, Copy, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct ModChar {
    val: char,
}

/// Defaults to mod tag 'N'
impl Default for ModChar {
    fn default() -> Self {
        ModChar::new('N')
    }
}

impl ModChar {
    /// We initialize with a character
    pub fn new(val: char) -> Self {
        ModChar { val }
    }
    /// Returns the character
    ///
    /// ```
    /// use nanalogue_core::ModChar;
    /// for y in vec!['a','b','c','\u{D000}']{
    ///     let x = ModChar::new(y.clone());
    ///     assert_eq!(x.val(), y);
    /// }
    /// ```
    pub fn val(&self) -> char {
        self.val
    }
}

impl FromStr for ModChar {
    type Err = Error;

    /// process the modification type from a string,
    /// returning the first character if it is a letter,
    /// or converting it to a character if the first character is a number
    fn from_str(mod_type: &str) -> Result<Self, Self::Err> {
        let first_char = mod_type.chars().next().ok_or(Error::EmptyModType)?;
        match first_char {
            'A'..='Z' | 'a'..='z' => Ok(ModChar { val: first_char }),
            '0'..='9' => {
                let val = char::from_u32(mod_type.parse()?).ok_or(Error::InvalidModType)?;
                Ok(ModChar { val })
            }
            _ => Err(Error::InvalidModType),
        }
    }
}

impl fmt::Display for ModChar {
    /// converts to string for display. If the value is in the alphabet,
    /// display it. Otherwise, display the equivalent u8 number.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}",
            match self.val() {
                w @ ('A'..='Z' | 'a'..='z') => w.to_string(),
                w => format!("{}", w as u32),
            }
        )
    }
}

/// Struct to mark if we want data only from the "+"
/// or the "-" mod-called strand.
/// NOTE: this doesn't mean the strand from the alignment,
/// but the strand from the modification data i.e. the strand
/// in the MM tag e.g. T+T, T-a where the strand is "+" and "-"
/// respectively. This denotes if mod data is from the basecalled
/// strand or the complementary strand. To not confuse users of the
/// command line interface, we allow creation of this from a string
/// "bc" or "bc_comp" i.e. basecalled or basecalled complement instead
/// of "+" and "-", which may be mistaken for the alignment strand.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct RestrictModCalledStrand(bool);

/// default mod strand restriction is +
impl Default for RestrictModCalledStrand {
    fn default() -> Self {
        RestrictModCalledStrand(true)
    }
}

/// Set states according to if we receive "bc" or "bc_comp"
///
/// ```
/// use nanalogue_core::RestrictModCalledStrand;
/// use std::str::FromStr;
/// assert_eq!('+', char::from(RestrictModCalledStrand::from_str("bc")?));
/// assert_eq!('-', char::from(RestrictModCalledStrand::from_str("bc_comp")?));
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
impl FromStr for RestrictModCalledStrand {
    type Err = Error;

    /// Parse a string to obtain float and then convert if b/w 0 and 1
    fn from_str(val_str: &str) -> Result<Self, Self::Err> {
        match val_str {
            "bc" => Ok(RestrictModCalledStrand(true)),
            "bc_comp" => Ok(RestrictModCalledStrand(false)),
            _ => Err(Error::InvalidState(
                "Please specify bc or bc_comp for mod-called strand!".to_string(),
            )),
        }
    }
}

/// Prints the true state as "+" and the false state as "-"
///
/// ```
/// use nanalogue_core::RestrictModCalledStrand;
/// use std::fmt::Display;
/// use std::str::FromStr;
/// assert_eq!("+", format!("{}",RestrictModCalledStrand::from_str("bc")?));
/// assert_eq!("-", format!("{}",RestrictModCalledStrand::from_str("bc_comp")?));
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
impl fmt::Display for RestrictModCalledStrand {
    /// converts to string for display i.e. "+" or "-"
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", char::from(*self))
    }
}

/// Converts the true state to the '+' character, and false to '-'
///
/// ```
/// use nanalogue_core::RestrictModCalledStrand;
/// use std::str::FromStr;
/// assert_eq!('+', char::from(RestrictModCalledStrand::from_str("bc")?));
/// assert_eq!('-', char::from(RestrictModCalledStrand::from_str("bc_comp")?));
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
impl From<RestrictModCalledStrand> for char {
    /// converts to char
    fn from(val: RestrictModCalledStrand) -> char {
        match val {
            RestrictModCalledStrand(true) => '+',
            RestrictModCalledStrand(false) => '-',
        }
    }
}

/// Alignment state of a read; seven possibilities + one unknown state
#[derive(Debug, Clone, Default, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReadState {
    #[default]
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
    /// modification probability not within this range.
    /// We expect this to be used to filter out modification calls
    /// around 0.5 i.e. ones with the most uncertainty, although
    /// users of this crate are free to set this to an interval
    /// not including 0.5
    InvertGtEqLtEq(OrdPair<u8>),
    /// modification probability >= first value, and mod prob
    /// not within the second range i.e. the 'and' combination
    /// of the two possibilities above
    Both((u8, OrdPair<u8>)),
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
/// # use nanalogue_core::{ThresholdState, OrdPair};
/// let b = ThresholdState::InvertGtEqLtEq(OrdPair::new(200, 220).expect("no error"));
/// assert_eq!("probabilities < 0.78125 or > 0.859375", format!("{}", b));
/// ```
///
/// Example 3:
/// ```
/// # use nanalogue_core::{ThresholdState, OrdPair};
/// let b = ThresholdState::Both((100, OrdPair::new(200, 220).expect("no error")));
/// assert_eq!("probabilities >= 0.390625 and (probabilities < 0.78125 or > 0.859375)", format!("{}", b));
/// ```
impl fmt::Display for ThresholdState {
    /// display the u8 thresholds as a floating point number between 0 and 1
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let printable = match &self {
            ThresholdState::GtEq(v) => format!("probabilities >= {}", F32Bw0and1::from(*v)),
            ThresholdState::InvertGtEqLtEq(v) => {
                format!(
                    "probabilities < {} or > {}",
                    F32Bw0and1::from(v.get_low()),
                    F32Bw0and1::from(v.get_high())
                )
            }
            ThresholdState::Both((a, b)) => {
                format!(
                    "{} and ({})",
                    ThresholdState::GtEq(*a),
                    ThresholdState::InvertGtEqLtEq(*b)
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
/// let b = ThresholdState::InvertGtEqLtEq(OrdPair::new(200, 220).expect("no error"));
/// assert!(b.contains(&100));
/// assert!(!b.contains(&200));
/// assert!(!b.contains(&210));
/// assert!(!b.contains(&220));
/// assert!(b.contains(&250));
/// ```
/// Example 3:
/// ```
/// # use nanalogue_core::{Error, OrdPair, ThresholdState, Contains};
/// let b = ThresholdState::Both((100, OrdPair::new(200, 220).expect("no error")));
/// assert!(!b.contains(&99));
/// assert!(b.contains(&100));
/// assert!(b.contains(&101));
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
            ThresholdState::InvertGtEqLtEq(w) => !w.contains(val),
            ThresholdState::Both((a, b)) => {
                ThresholdState::GtEq(*a).contains(val)
                    && ThresholdState::InvertGtEqLtEq(*b).contains(val)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tests if our Ordered Pair struct can be correctly obtained from strings
    #[test]
    fn test_ord_pair_from_str() {
        assert_eq!(
            OrdPair::<f32>::from_str("1.0,2.0")
                .expect("no failure")
                .get_low(),
            1.0
        );
        assert_eq!(
            OrdPair::<f32>::from_str("1.0,2.0")
                .expect("no failure")
                .get_high(),
            2.0
        );
        assert_eq!(
            OrdPair::<u8>::from_str("1, 2")
                .expect("no failure")
                .get_low(),
            1
        );
        assert_eq!(
            OrdPair::<u8>::from_str("1, 2")
                .expect("no failure")
                .get_high(),
            2
        );
        assert!(matches!(
            OrdPair::<u8>::from_str(",2"),
            Err(Error::OrdPairConversionError(_))
        ));
        assert!(matches!(
            OrdPair::<u8>::from_str("2,1"),
            Err(Error::WrongOrder)
        ));
    }

    /// Tests if OrdPair can be converted into a range
    #[test]
    fn test_ord_pair_to_range() {
        assert_eq!(
            (3..=5),
            RangeInclusive::from(OrdPair::new(3, 5).expect("no failure"))
        )
    }

    /// Tests if ModChar is displayed correctly
    #[test]
    fn test_display_mod_char() {
        assert_eq!(
            format!("{}", ModChar::from_str("a").expect("no failure")),
            "a"
        );
        assert_eq!(
            format!("{}", ModChar::from_str("T").expect("no failure")),
            "T"
        );
        assert_eq!(
            format!("{}", ModChar::from_str("77000").expect("no failure")),
            "77000"
        );
    }
}
