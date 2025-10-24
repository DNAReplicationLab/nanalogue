//! `RestrictModCalledStrand` struct for strand-specific modification filtering
//! Controls whether to use modifications from basecalled or complement strand

use crate::Error;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

/// Struct to mark if we want data only from the "+"
/// or the "-" mod-called strand.
/// NOTE: this doesn't mean the strand from the alignment,
/// but the strand from the modification data i.e. the strand
/// in the MM tag e.g. T+T, T-a where the strand is "+" and "-"
/// respectively. This denotes if mod data is from the basecalled
/// strand or the complementary strand. To not confuse users of the
/// command line interface, we allow creation of this from a string
/// "bc" or "`bc_comp`" i.e. basecalled or basecalled complement instead
/// of "+" and "-", which may be mistaken for the alignment strand.
#[derive(Debug, Clone, Copy, Eq, PartialEq, Serialize, Deserialize)]
pub struct RestrictModCalledStrand(bool);

/// default mod strand restriction is +
impl Default for RestrictModCalledStrand {
    fn default() -> Self {
        RestrictModCalledStrand(true)
    }
}

/// Set states according to if we receive "bc" or "`bc_comp`"
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

    /// Parse a string to create a `RestrictModCalledStrand` ("bc or "`bc_comp`")
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_restrict_mod_called_strand_basic() {
        let bc = RestrictModCalledStrand::from_str("bc").expect("should parse");
        assert_eq!(format!("{bc}"), "+");
        assert_eq!(char::from(bc), '+');

        let bc_comp = RestrictModCalledStrand::from_str("bc_comp").expect("should parse");
        assert_eq!(format!("{bc_comp}"), "-");
        assert_eq!(char::from(bc_comp), '-');
    }

    #[test]
    fn test_restrict_mod_called_strand_errors() {
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
    }

    #[test]
    fn test_restrict_mod_called_strand_default() {
        let default_strand = RestrictModCalledStrand::default();
        assert_eq!(char::from(default_strand), '+');
    }
}
