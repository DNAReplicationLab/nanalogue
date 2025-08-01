use crate::Error;
use std::fmt::Debug;
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
