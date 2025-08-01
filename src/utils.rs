use crate::Error;
use std::fmt::Debug;
use std::fmt;
use std::str::FromStr;

#[derive(Debug, Clone, Default, Copy, PartialOrd, PartialEq)]
pub struct OrdPair<T: Clone + Copy + Debug + Default> {
    low: T,
    high: T,
}

impl<T: Clone + Copy + Debug + Default + PartialEq + PartialOrd> OrdPair<T>{
    pub fn new(low: T, high: T) -> Result<Self, Error> {
        if low <= high {
            Ok(OrdPair { low, high, })
        } else {
            Err(Error::WrongOrder)
        }
    }
    pub fn get_low(&self) -> T {
        self.low
    }
    pub fn get_high(&self) -> T {
        self.high
    }
}

impl<T: Clone + Copy + Debug + Default + PartialEq + PartialOrd + FromStr> FromStr for OrdPair<T>{
    type Err = Error;

    fn from_str(val_str: &str) -> Result<Self, Self::Err>{
        macro_rules! parse_error {
            () => {
                Err(Error::OrdPairConversionError("Bad ordered pair inputs!".to_string()))
            }
        }
        let v: Vec<&str> = val_str.split(",").collect();
        match v.len() {
            2 => {
                let Ok(low) = T::from_str(v[0]) else { parse_error!()? };
                let Ok(high) = T::from_str(v[1]) else { parse_error!()? };
                OrdPair::<T>::new(low, high)
            }
            _ => parse_error!(),
        }
    }
}

#[derive(Debug, Clone, Default, Copy, PartialOrd, PartialEq)]
pub struct F32Bw0and1{
    val: f32
}

impl F32Bw0and1{
    pub fn new(val: f32) -> Result<Self, Error>{
        if (0.0..=1.0).contains(&val) {
            Ok(F32Bw0and1{ val })
        } else {
            Err(Error::InvalidState("Num not b/w 0 and 1!".to_string()))
        }
    }
    pub fn get_val(&self) -> f32 {
        self.val
    }
}

impl FromStr for F32Bw0and1 {
    type Err = Error;

    fn from_str(val_str: &str) -> Result<Self, Self::Err> {
        Self::new(f32::from_str(val_str)?)
    }
}

#[derive(Debug, Clone, Default, Copy, PartialEq, PartialOrd)]
pub struct ModChar{
    val: char
}

impl ModChar{
    pub fn new(val: char) -> Self{
        ModChar { val }
    }
    pub fn get_val(&self) -> char {
        self.val
    }
}

impl FromStr for ModChar{
    type Err = Error;

    fn from_str(mod_type: &str) -> Result<Self, Self::Err> {
        // process the modification type, returning the first character if it is a letter,
        // or converting it to a character if it is a number
        let first_char = mod_type.chars().next().ok_or(Error::EmptyModType)?; 
        match first_char {
            'A' ..= 'Z' | 'a' ..= 'z' => Ok(ModChar{ val: first_char}),
            '0' ..= '9' => {
                let val = char::from_u32(mod_type.parse()?).ok_or(Error::InvalidModType)?;
                Ok(ModChar{ val })
            }, 
            _ => Err(Error::InvalidModType),
        }
    }
}

impl fmt::Display for ModChar{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.get_val())
    }
}
