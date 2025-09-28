//! # Utils implementing simple, shared datatypes
//!
//! Implements simple datatypes like an ordered pair or mod tag
//! used by all modules in our crate.

use crate::Error;
use bedrs::prelude::Bed3;
use rust_htslib::bam;
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

impl OrdPair<u64> {
    /// Parse an interval string specifically for genomic regions.
    /// Supports formats like "1000-2000" and "1000-" (where end defaults to u64::MAX).
    /// Enforces strict inequality (start < end).
    ///
    /// ```
    /// use nanalogue_core::OrdPair;
    ///
    /// // Standard interval
    /// let interval = OrdPair::<u64>::from_interval("1000-2000")?;
    /// assert_eq!(interval.get_low(), 1000);
    /// assert_eq!(interval.get_high(), 2000);
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```
    /// # use nanalogue_core::OrdPair;
    /// #
    /// // Open-ended interval (end defaults to u64::MAX)
    /// let interval = OrdPair::<u64>::from_interval("1000-")?;
    /// assert_eq!(interval.get_low(), 1000);
    /// assert_eq!(interval.get_high(), u64::MAX);
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```should_panic
    /// # use nanalogue_core::OrdPair;
    /// #
    /// // Equal start and end should fail
    /// let interval = OrdPair::<u64>::from_interval("1000-1000")?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    pub fn from_interval(interval_str: &str) -> Result<Self, Error> {
        let parts: Vec<&str> = interval_str.split('-').collect();

        match parts.len() {
            2 => {
                let start = parts[0].trim().parse::<u64>().map_err(|_| {
                    Error::OrdPairConversionError(
                        "Invalid start coordinate in interval!".to_string(),
                    )
                })?;

                let end = if parts[1].trim().is_empty() {
                    // Open-ended interval: "1000-"
                    u64::MAX
                } else {
                    // Closed interval: "1000-2000"
                    parts[1].trim().parse::<u64>().map_err(|_| {
                        Error::OrdPairConversionError(
                            "Invalid end coordinate in interval!".to_string(),
                        )
                    })?
                };

                // Enforce strict inequality (start < end)
                if start < end {
                    Ok(OrdPair {
                        low: start,
                        high: end,
                    })
                } else {
                    Err(Error::WrongOrder)
                }
            }
            _ => Err(Error::OrdPairConversionError(
                "Invalid interval format! Expected 'start-end' or 'start-'".to_string(),
            )),
        }
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

/// Datatype holding a genomic region
#[derive(Debug, Default, Clone, PartialOrd, PartialEq, Serialize, Deserialize)]
pub struct GenomicRegion(pub (String, Option<OrdPair<u64>>));

/// Obtains genomic region from a string with the standard region format of name[:begin[-end]].
///
/// ```
/// use nanalogue_core::GenomicRegion;
/// use std::str::FromStr;
///
/// // Simple contig name only
/// let region = GenomicRegion::from_str("chr1")?;
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
///
/// ```
/// # use nanalogue_core::GenomicRegion;
/// # use std::str::FromStr;
/// #
/// // Contig with coordinates
/// let region = GenomicRegion::from_str("chr1:1000-2000")?;
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
///
/// ```
/// # use nanalogue_core::GenomicRegion;
/// # use std::str::FromStr;
/// #
/// // Contig name with colons (e.g., from some assemblies)
/// let region = GenomicRegion::from_str("chr1:alternate:1000-2000")?;
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
impl FromStr for GenomicRegion {
    type Err = Error;

    fn from_str(val_str: &str) -> Result<Self, Self::Err> {
        let mut colon_split: Vec<&str> = val_str.split(":").collect();
        match colon_split.len() {
            0 => Err(Error::InvalidAlignCoords),
            1 => Ok(GenomicRegion((val_str.to_string(), None))),
            _ => {
                let interval_str = colon_split.pop().ok_or(Error::UnknownError)?;
                Ok(GenomicRegion((
                    colon_split.join(":").to_string(),
                    Some(OrdPair::<u64>::from_interval(interval_str)?),
                )))
            }
        }
    }
}

/// Converts genomic region from genomic string representation to bed3 representation
impl GenomicRegion {
    /// converts genomic region from genomic string representation to bed3 representation
    pub fn try_to_bed3(self, header: bam::HeaderView) -> Result<Bed3<i32, u64>, Error> {
        let region_bed = {
            let GenomicRegion((contig_name, coords)) = &self;
            let numeric_contig: i32 = header
                .tid(contig_name.as_bytes())
                .ok_or(Error::InvalidAlignCoords)?
                .try_into()?;

            let (start, end) = if let Some(c) = coords {
                let start = c.get_low();
                let end = c.get_high();

                // Check if start position exceeds contig length
                if let Some(contig_length) = header.target_len(u32::try_from(numeric_contig)?) {
                    if start >= contig_length {
                        let region_str = if end == u64::MAX {
                            format!("{}:{}-", contig_name, start)
                        } else {
                            format!("{}:{}-{}", contig_name, start, end)
                        };

                        return Err(Error::InvalidRegionError {
                            region: region_str,
                            start,
                            contig_length,
                        });
                    }
                }

                (start, end)
            } else {
                (u64::MIN, u64::MAX)
            };

            Bed3::<i32, u64>::new(numeric_contig, start, end)
        };
        Ok(region_bed)
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
    ///
    /// ```
    /// use nanalogue_core::F32Bw0and1;
    /// use std::str::FromStr;
    ///
    /// // Boundary values - exactly 0.0 and 1.0 should work
    /// let zero = F32Bw0and1::from_str("0.0")?;
    /// assert_eq!(zero.val(), 0.0);
    /// let one = F32Bw0and1::from_str("1.0")?;
    /// assert_eq!(one.val(), 1.0);
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```
    /// # use nanalogue_core::F32Bw0and1;
    /// # use std::str::FromStr;
    /// #
    /// // Near-boundary values
    /// let near_zero = F32Bw0and1::from_str("0.000001")?;
    /// let near_one = F32Bw0and1::from_str("0.999999")?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```should_panic
    /// # use nanalogue_core::F32Bw0and1;
    /// # use std::str::FromStr;
    /// #
    /// // Just outside boundaries should fail
    /// let outside = F32Bw0and1::from_str("1.000001")?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
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
    ///
    /// ```
    /// use nanalogue_core::F32AbsValBelow1;
    /// use std::str::FromStr;
    ///
    /// // Boundary values - exactly -1.0, 0.0, and 1.0 should work
    /// let neg_one = F32AbsValBelow1::from_str("-1.0")?;
    /// assert_eq!(neg_one.val(), -1.0);
    /// let zero = F32AbsValBelow1::from_str("0.0")?;
    /// assert_eq!(zero.val(), 0.0);
    /// let one = F32AbsValBelow1::from_str("1.0")?;
    /// assert_eq!(one.val(), 1.0);
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```
    /// # use nanalogue_core::F32AbsValBelow1;
    /// # use std::str::FromStr;
    /// #
    /// // Near-boundary values
    /// let near_neg_one = F32AbsValBelow1::from_str("-0.999999")?;
    /// let near_one = F32AbsValBelow1::from_str("0.999999")?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```should_panic
    /// # use nanalogue_core::F32AbsValBelow1;
    /// # use std::str::FromStr;
    /// #
    /// // Just outside boundaries should fail
    /// let outside = F32AbsValBelow1::from_str("-1.000001")?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
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
#[derive(Debug, Clone, Copy, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct ModChar(char);

/// Defaults to mod tag 'N'
impl Default for ModChar {
    fn default() -> Self {
        ModChar::new('N')
    }
}

impl ModChar {
    /// We initialize with a character
    pub fn new(val: char) -> Self {
        ModChar(val)
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
        self.0
    }
}

impl FromStr for ModChar {
    type Err = Error;

    /// process the modification type from a string,
    /// returning the first character if it is a letter,
    /// or converting it to a character if the first character is a number
    ///
    /// ```
    /// use nanalogue_core::ModChar;
    /// use std::str::FromStr;
    ///
    /// // Single letter modification codes
    /// let mod_char = ModChar::from_str("m")?;
    /// assert_eq!(mod_char.val(), 'm');
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```
    /// # use nanalogue_core::ModChar;
    /// # use std::str::FromStr;
    /// #
    /// // CheBI code for BrdU (5-bromo-2'-deoxyuridine)
    /// let mod_char = ModChar::from_str("472232")?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```
    /// # use nanalogue_core::ModChar;
    /// # use std::str::FromStr;
    /// #
    /// // Small numeric codes
    /// let mod_char = ModChar::from_str("123")?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    ///
    /// ```should_panic
    /// # use nanalogue_core::ModChar;
    /// # use std::str::FromStr;
    /// #
    /// // Invalid: starts with special character
    /// let mod_char = ModChar::from_str("@123")?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    fn from_str(mod_type: &str) -> Result<Self, Self::Err> {
        let first_char = mod_type.chars().next().ok_or(Error::EmptyModType)?;
        match first_char {
            'A'..='Z' | 'a'..='z' => Ok(ModChar(first_char)),
            '0'..='9' => {
                let val = char::from_u32(mod_type.parse()?).ok_or(Error::InvalidModType)?;
                Ok(ModChar(val))
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

impl Serialize for ModChar {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        // Use the Display implementation to serialize as a string
        serializer.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for ModChar {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        ModChar::from_str(&s).map_err(serde::de::Error::custom)
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
    #[serde(rename = "primary_forward")]
    PrimaryFwd,
    /// Primary alignment opposite the reference strand
    #[serde(rename = "primary_reverse")]
    PrimaryRev,
    /// Secondary alignment to the reference strand
    #[serde(rename = "secondary_forward")]
    SecondaryFwd,
    /// Secondary alignment opposite the reference strand
    #[serde(rename = "secondary_reverse")]
    SecondaryRev,
    /// Supplementary alignment to the reference strand
    #[serde(rename = "supplementary_forward")]
    SupplementaryFwd,
    /// Supplementary alignment opposite the reference strand
    #[serde(rename = "supplementary_reverse")]
    SupplementaryRev,
    /// Marked as unmapped in the BAM file. We are assuming
    /// that unmapped sequences will not be stored as reversed
    /// complements, as what would be the point of that?
    #[serde(rename = "unmapped")]
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
///
/// ```
/// use nanalogue_core::ReadState;
/// use std::str::FromStr;
///
/// // Primary alignments
/// let state = ReadState::from_str("primary_forward")?;
/// assert_eq!(state, ReadState::PrimaryFwd);
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
///
/// ```
/// # use nanalogue_core::ReadState;
/// # use std::str::FromStr;
/// #
/// // Secondary alignments
/// let state = ReadState::from_str("secondary_reverse")?;
/// assert_eq!(state, ReadState::SecondaryRev);
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
///
/// ```
/// # use nanalogue_core::ReadState;
/// # use std::str::FromStr;
/// #
/// // Supplementary alignments
/// let state = ReadState::from_str("supplementary_forward")?;
/// assert_eq!(state, ReadState::SupplementaryFwd);
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
///
/// ```
/// # use nanalogue_core::ReadState;
/// # use std::str::FromStr;
/// #
/// // Unmapped reads
/// let state = ReadState::from_str("unmapped")?;
/// assert_eq!(state, ReadState::Unmapped);
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
///
/// ```should_panic
/// # use nanalogue_core::ReadState;
/// # use std::str::FromStr;
/// #
/// // Invalid string should error
/// let state = ReadState::from_str("invalid_state")?;
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
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

    /// Tests OrdPair::from_interval method for genomic intervals
    #[test]
    fn test_ord_pair_from_interval() {
        // Standard interval
        let interval = OrdPair::<u64>::from_interval("1000-2000").expect("should parse");
        assert_eq!(interval.get_low(), 1000);
        assert_eq!(interval.get_high(), 2000);

        // Open-ended interval
        let interval = OrdPair::<u64>::from_interval("1000-").expect("should parse");
        assert_eq!(interval.get_low(), 1000);
        assert_eq!(interval.get_high(), u64::MAX);
    }

    /// Tests OrdPair::from_interval error cases
    #[test]
    fn test_ord_pair_from_interval_errors() {
        // Equal start and end should fail (strict inequality enforced)
        assert!(matches!(
            OrdPair::<u64>::from_interval("1000-1000"),
            Err(Error::WrongOrder)
        ));

        // Start greater than end should fail
        assert!(matches!(
            OrdPair::<u64>::from_interval("2000-1000"),
            Err(Error::WrongOrder)
        ));

        // Invalid format - no dash
        assert!(matches!(
            OrdPair::<u64>::from_interval("1000"),
            Err(Error::OrdPairConversionError(_))
        ));

        // Invalid format - too many dashes
        assert!(matches!(
            OrdPair::<u64>::from_interval("1000-2000-3000"),
            Err(Error::OrdPairConversionError(_))
        ));

        // Invalid start coordinate
        assert!(matches!(
            OrdPair::<u64>::from_interval("abc-2000"),
            Err(Error::OrdPairConversionError(_))
        ));

        // Invalid end coordinate
        assert!(matches!(
            OrdPair::<u64>::from_interval("1000-xyz"),
            Err(Error::OrdPairConversionError(_))
        ));

        // Empty string
        assert!(matches!(
            OrdPair::<u64>::from_interval(""),
            Err(Error::OrdPairConversionError(_))
        ));

        // Just a dash
        assert!(matches!(
            OrdPair::<u64>::from_interval("-"),
            Err(Error::OrdPairConversionError(_))
        ));

        // Negative numbers (should fail for u64)
        assert!(matches!(
            OrdPair::<u64>::from_interval("-100-200"),
            Err(Error::OrdPairConversionError(_))
        ));
    }

    /// Tests comprehensive GenomicRegion parsing
    #[test]
    fn test_genomic_region_parsing() {
        // Simple contig name only
        let region = GenomicRegion::from_str("chr1").expect("should parse");
        assert_eq!(region.0.0, "chr1");
        assert_eq!(region.0.1, None);

        // Contig with coordinates
        let region = GenomicRegion::from_str("chr1:1000-2000").expect("should parse");
        assert_eq!(region.0.0, "chr1");
        assert!(region.0.1.is_some());
        let coords = region.0.1.unwrap();
        assert_eq!(coords.get_low(), 1000);
        assert_eq!(coords.get_high(), 2000);

        // Contig name with colons (e.g., from some assemblies)
        let region = GenomicRegion::from_str("chr1:alternate:1000-2000").expect("should parse");
        assert_eq!(region.0.0, "chr1:alternate");
        assert!(region.0.1.is_some());
        let coords = region.0.1.unwrap();
        assert_eq!(coords.get_low(), 1000);
        assert_eq!(coords.get_high(), 2000);

        // Complex contig names
        let region = GenomicRegion::from_str("scaffold_123:456-789").expect("should parse");
        assert_eq!(region.0.0, "scaffold_123");
        let coords = region.0.1.unwrap();
        assert_eq!(coords.get_low(), 456);
        assert_eq!(coords.get_high(), 789);

        // Mitochondrial chromosome
        let region = GenomicRegion::from_str("chrM:1-16569").expect("should parse");
        assert_eq!(region.0.0, "chrM");
        let coords = region.0.1.unwrap();
        assert_eq!(coords.get_low(), 1);
        assert_eq!(coords.get_high(), 16569);
    }

    /// Tests GenomicRegion parsing with open-ended support
    #[test]
    fn test_genomic_region_open_ended() {
        // Open-ended interval support
        let region = GenomicRegion::from_str("chr1:1000-").expect("should parse");
        assert_eq!(region.0.0, "chr1");
        assert!(region.0.1.is_some());
        let coords = region.0.1.unwrap();
        assert_eq!(coords.get_low(), 1000);
        assert_eq!(coords.get_high(), u64::MAX);
    }

    /// Tests GenomicRegion parsing error cases
    #[test]
    fn test_genomic_region_parsing_errors() {
        // Wrong order coordinates should fail
        assert!(matches!(
            GenomicRegion::from_str("chr1:2000-1000"),
            Err(Error::WrongOrder)
        ));

        // Equal start and end should now fail (strict inequality)
        assert!(matches!(
            GenomicRegion::from_str("chr1:1000-1000"),
            Err(Error::WrongOrder)
        ));

        // Invalid coordinate format should fail
        assert!(GenomicRegion::from_str("chr1:abc-def").is_err());

        // Invalid format should fail
        assert!(GenomicRegion::from_str("chr1:1000-2000-3000").is_err());
    }

    /// Tests ReadState u16 conversion round-trip
    #[test]
    fn test_readstate_u16_conversion_roundtrip() {
        let states = vec![
            ReadState::PrimaryFwd,
            ReadState::PrimaryRev,
            ReadState::SecondaryFwd,
            ReadState::SecondaryRev,
            ReadState::SupplementaryFwd,
            ReadState::SupplementaryRev,
            ReadState::Unmapped,
        ];

        for state in states {
            // Convert to u16 and back
            let flag: u16 = state.try_into().expect("conversion to u16 should work");
            let recovered_state: ReadState =
                flag.try_into().expect("conversion from u16 should work");
            assert_eq!(state, recovered_state);
        }
    }

    /// Tests specific ReadState u16 flag values
    #[test]
    fn test_readstate_specific_flag_values() {
        assert_eq!(u16::try_from(ReadState::PrimaryFwd).unwrap(), 0);
        assert_eq!(u16::try_from(ReadState::Unmapped).unwrap(), 4);
        assert_eq!(u16::try_from(ReadState::PrimaryRev).unwrap(), 16);
        assert_eq!(u16::try_from(ReadState::SecondaryFwd).unwrap(), 256);
        assert_eq!(u16::try_from(ReadState::SecondaryRev).unwrap(), 272);
        assert_eq!(u16::try_from(ReadState::SupplementaryFwd).unwrap(), 2048);
        assert_eq!(u16::try_from(ReadState::SupplementaryRev).unwrap(), 2064);
    }

    /// Tests ReadState from invalid u16 values
    #[test]
    fn test_readstate_invalid_u16_values() {
        // Test various invalid flag combinations
        let invalid_flags = vec![1, 2, 8, 32, 64, 128, 512, 1024, 4096, 8192];
        for flag in invalid_flags {
            assert!(matches!(
                ReadState::try_from(flag),
                Err(Error::UnknownAlignState)
            ));
        }
    }

    /// Tests ModChar numeric conversion and edge cases
    #[test]
    fn test_modchar_numeric_conversion() {
        // Test letter codes
        let mod_char = ModChar::from_str("m").expect("should parse");
        assert_eq!(mod_char.val(), 'm');
        assert_eq!(format!("{}", mod_char), "m");

        let mod_char = ModChar::from_str("T").expect("should parse");
        assert_eq!(mod_char.val(), 'T');
        assert_eq!(format!("{}", mod_char), "T");

        // Test small numeric codes
        let mod_char = ModChar::from_str("123").expect("should parse");
        assert_eq!(format!("{}", mod_char), "123");

        // Test CheBI code for BrdU
        let mod_char = ModChar::from_str("472232").expect("should parse");
        assert_eq!(format!("{}", mod_char), "472232");

        // Test ASCII boundary - 97 is 'a'
        let mod_char = ModChar::from_str("97").expect("should parse");
        assert_eq!(mod_char.val(), 'a');
        // When the char value is in alphabet range, it displays as the letter, not the number
        assert_eq!(format!("{}", mod_char), "a");

        // Test very large numbers that are valid unicode
        let mod_char = ModChar::from_str("65536").expect("should parse");
        assert_eq!(format!("{}", mod_char), "65536");
    }

    /// Tests ModChar error cases
    #[test]
    fn test_modchar_error_cases() {
        // Empty string should fail
        assert!(matches!(ModChar::from_str(""), Err(Error::EmptyModType)));

        // Starting with special characters should fail
        assert!(matches!(
            ModChar::from_str("@123"),
            Err(Error::InvalidModType)
        ));
        assert!(matches!(
            ModChar::from_str("#abc"),
            Err(Error::InvalidModType)
        ));

        // Invalid numeric values that can't convert to char should fail
        // (Note: this would be numbers > u32::MAX or invalid unicode ranges)
        // Most very large numbers should still work due to char's range
    }

    /// Tests ModChar display format consistency
    #[test]
    fn test_modchar_display_consistency() {
        // Letters should display as letters
        for letter in ['a', 'b', 'z', 'A', 'B', 'Z'] {
            let mod_char = ModChar::new(letter);
            assert_eq!(format!("{}", mod_char), letter.to_string());
        }

        // Numbers converted to char should display as their numeric value
        let test_numbers = vec![123, 456, 789, 472232];
        for num in test_numbers {
            let mod_char = ModChar::from_str(&num.to_string()).expect("should parse");
            assert_eq!(format!("{}", mod_char), num.to_string());
        }
    }

    /// Tests integration between F32Bw0and1 and F32AbsValBelow1 types
    #[test]
    fn test_f32_types_integration() {
        // Test conversion from F32Bw0and1 to F32AbsValBelow1
        let pos_values = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        for val in pos_values {
            let bw_val = F32Bw0and1::new(val).expect("should create");
            let abs_val: F32AbsValBelow1 = bw_val.into();
            assert_eq!(abs_val.val(), val);
        }

        // Test conversion via absolute value function
        let neg_val = F32AbsValBelow1::new(-0.5).expect("should create");
        let abs_converted = F32Bw0and1::abs_f32_abs_val_below_1(neg_val);
        assert_eq!(abs_converted.val(), 0.5);

        let pos_val = F32AbsValBelow1::new(0.7).expect("should create");
        let abs_converted = F32Bw0and1::abs_f32_abs_val_below_1(pos_val);
        assert_eq!(abs_converted.val(), 0.7);
    }

    /// Tests u8 to F32Bw0and1 conversion
    #[test]
    fn test_u8_to_f32bw0and1_conversion() {
        // Test boundary values
        let zero = F32Bw0and1::from(0u8);
        assert_eq!(zero.val(), 0.0);

        let max = F32Bw0and1::from(255u8);
        assert!(max.val() < 1.0); // Should be 255/256 = 0.99609375

        // Test some intermediate values
        let half = F32Bw0and1::from(128u8);
        assert!(half.val() > 0.49 && half.val() < 0.51);

        // Test exact calculation
        let test_val = 100u8;
        let converted = F32Bw0and1::from(test_val);
        let expected = (test_val as f32) / (u8::MAX as f32 + 1.0);
        assert_eq!(converted.val(), expected);
    }

    /// Tests shortcut constructors for F32Bw0and1
    #[test]
    fn test_f32bw0and1_shortcuts() {
        let zero = F32Bw0and1::zero();
        assert_eq!(zero.val(), 0.0);

        let one = F32Bw0and1::one();
        assert_eq!(one.val(), 1.0);
    }

    /// Tests error handling for all FromStr implementations
    #[test]
    fn test_fromstr_error_handling() {
        // OrdPair error cases
        assert!(OrdPair::<i32>::from_str("").is_err());
        assert!(OrdPair::<i32>::from_str("1").is_err()); // Single value
        assert!(OrdPair::<i32>::from_str("1,2,3").is_err()); // Too many values
        assert!(OrdPair::<i32>::from_str("abc,def").is_err()); // Non-numeric
        assert!(matches!(
            OrdPair::<i32>::from_str("2,1"),
            Err(Error::WrongOrder)
        )); // Wrong order

        // F32Bw0and1 error cases
        assert!(F32Bw0and1::from_str("").is_err());
        assert!(F32Bw0and1::from_str("abc").is_err());
        assert!(matches!(
            F32Bw0and1::from_str("-0.1"),
            Err(Error::InvalidState(_))
        ));
        assert!(matches!(
            F32Bw0and1::from_str("1.1"),
            Err(Error::InvalidState(_))
        ));

        // F32AbsValBelow1 error cases
        assert!(F32AbsValBelow1::from_str("").is_err());
        assert!(F32AbsValBelow1::from_str("xyz").is_err());
        assert!(matches!(
            F32AbsValBelow1::from_str("-1.1"),
            Err(Error::InvalidState(_))
        ));
        assert!(matches!(
            F32AbsValBelow1::from_str("1.1"),
            Err(Error::InvalidState(_))
        ));

        // RestrictModCalledStrand error cases
        assert!(matches!(
            RestrictModCalledStrand::from_str("invalid"),
            Err(Error::InvalidState(_))
        ));
        assert!(matches!(
            RestrictModCalledStrand::from_str(""),
            Err(Error::InvalidState(_))
        ));
        assert!(matches!(
            RestrictModCalledStrand::from_str("+"),
            Err(Error::InvalidState(_))
        )); // Should be "bc", not "+"

        // ReadState error cases
        assert!(matches!(
            ReadState::from_str("invalid_state"),
            Err(Error::UnknownAlignState)
        ));
        assert!(matches!(
            ReadState::from_str(""),
            Err(Error::UnknownAlignState)
        ));
        assert!(matches!(
            ReadState::from_str("primary"), // Incomplete
            Err(Error::UnknownAlignState)
        ));
    }

    /// Tests display format consistency across all types
    #[test]
    fn test_display_format_consistency() {
        // OrdPair display format
        let pair = OrdPair::new(10, 20).expect("should create");
        assert_eq!(format!("{}", pair), "10, 20");

        let float_pair = OrdPair::new(1.5, 2.5).expect("should create");
        assert_eq!(format!("{}", float_pair), "1.5, 2.5");

        // F32Bw0and1 display format
        let f32_val = F32Bw0and1::new(0.5).expect("should create");
        assert_eq!(format!("{}", f32_val), "0.5");

        // F32AbsValBelow1 display format
        let abs_val = F32AbsValBelow1::new(-0.5).expect("should create");
        assert_eq!(format!("{}", abs_val), "-0.5");

        // RestrictModCalledStrand display format
        let bc = RestrictModCalledStrand::from_str("bc").expect("should parse");
        assert_eq!(format!("{}", bc), "+");

        let bc_comp = RestrictModCalledStrand::from_str("bc_comp").expect("should parse");
        assert_eq!(format!("{}", bc_comp), "-");

        // ReadState display format
        assert_eq!(format!("{}", ReadState::PrimaryFwd), "primary_forward");
        assert_eq!(format!("{}", ReadState::PrimaryRev), "primary_reverse");
        assert_eq!(format!("{}", ReadState::SecondaryFwd), "secondary_forward");
        assert_eq!(format!("{}", ReadState::SecondaryRev), "secondary_reverse");
        assert_eq!(
            format!("{}", ReadState::SupplementaryFwd),
            "supplementary_forward"
        );
        assert_eq!(
            format!("{}", ReadState::SupplementaryRev),
            "supplementary_reverse"
        );
        assert_eq!(format!("{}", ReadState::Unmapped), "unmapped");

        // ThresholdState display format
        let threshold = ThresholdState::GtEq(128);
        assert!(format!("{}", threshold).contains("probabilities >="));
        assert!(format!("{}", threshold).contains("0.5")); // 128/256 = 0.5

        let range_threshold =
            ThresholdState::InvertGtEqLtEq(OrdPair::new(100, 150).expect("should create"));
        let display_str = format!("{}", range_threshold);
        assert!(display_str.contains("probabilities <"));
        assert!(display_str.contains("or >"));
    }
}
