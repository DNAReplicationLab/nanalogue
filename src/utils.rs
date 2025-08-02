//! # Utils implementing simple, shared datatypes
//!
//! Implements simple datatypes like an ordered pair or mod tag
//! used by all modules in our crate.

use crate::Error;
use std::fmt;
use std::fmt::Debug;
use std::str::FromStr;

/// Datatype holding two values low, high such that low <= high is guaranteed at creation.
#[derive(Debug, Clone, Default, Copy, PartialOrd, PartialEq)]
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
    /// use nanalogue_core::OrdPair;
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

/// Datatype holding a float (f32) between 0 and 1 (both inclusive) guaranteed at creation.
#[derive(Debug, Clone, Default, Copy, PartialOrd, PartialEq)]
pub struct F32Bw0and1 {
    val: f32,
}

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
    /// use nanalogue_core::Error;
    /// use nanalogue_core::F32Bw0and1;
    /// let x = F32Bw0and1::new(0.1)?;
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    pub fn new(val: f32) -> Result<Self, Error> {
        if (0.0..=1.0).contains(&val) {
            Ok(F32Bw0and1 { val })
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
    ///     assert_eq!(x.get_val(), y);
    /// }
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    pub fn get_val(&self) -> f32 {
        self.val
    }
}

impl FromStr for F32Bw0and1 {
    type Err = Error;

    /// Parse a string to obtain float and then convert if b/w 0 and 1
    fn from_str(val_str: &str) -> Result<Self, Self::Err> {
        Self::new(f32::from_str(val_str)?)
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
#[derive(Debug, Clone, Default, Copy, Eq, Hash, PartialEq, PartialOrd)]
pub struct ModChar {
    val: char,
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
    ///     assert_eq!(x.get_val(), y);
    /// }
    /// ```
    pub fn get_val(&self) -> char {
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
            match self.get_val() {
                w @ ('A'..='Z' | 'a'..='z') => w.to_string(),
                w => format!("{}", w as u32),
            }
        )
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
