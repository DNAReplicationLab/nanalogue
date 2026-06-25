//! Helpers for protecting loops against runaway iteration.

use crate::Error;

/// Increment a counter and assert that it stays within a configured bound.
///
/// # Errors
/// Returns an error if the counter overflows or exceeds `max`.
pub fn assert_bounded_counter(idx: &mut u32, max: u32, what: &str) -> Result<(), Error> {
    *idx = idx.checked_add(1).ok_or_else(|| {
        Error::InvalidState(format!("{what} record counter overflowed at u32::MAX"))
    })?;

    if *idx > max {
        return Err(Error::InvalidState(format!(
            "{what} record limit exceeded: {max}"
        )));
    }

    Ok(())
}

/// Assert that a counter is non-zero.
///
/// # Errors
/// Returns an error if `idx` is zero.
pub fn assert_nonzero_counter(idx: u32, what: &str) -> Result<(), Error> {
    if idx == 0 {
        return Err(Error::InvalidState(format!("no {what} found!")));
    }

    Ok(())
}

/// Assert that a BAM record's internal data capacity stays within a configured bound.
///
/// # Errors
/// Returns an error if `m_data` exceeds `max`.
pub fn assert_record_data_capacity(m_data: u32, max: u32, what: &str) -> Result<(), Error> {
    if m_data > max {
        return Err(Error::InvalidState(format!(
            "{what} record capacity limit exceeded: {max}"
        )));
    }

    Ok(())
}

/// Assert that a flag is true.
///
/// # Errors
/// Returns an error if `flag` is false.
pub fn assert_flag(flag: bool, msg: &str) -> Result<(), Error> {
    if !flag {
        return Err(Error::InvalidState(msg.to_owned()));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        assert_bounded_counter, assert_flag, assert_nonzero_counter, assert_record_data_capacity,
    };
    use crate::Error;
    use crate::constants::shared::MAX_RECORD_CAPACITY_BYTES;

    #[test]
    fn assert_bounded_counter_increments_within_bound() -> Result<(), Error> {
        let mut idx = 0;

        assert_bounded_counter(&mut idx, 2, "peek")?;
        assert_eq!(idx, 1);

        assert_bounded_counter(&mut idx, 2, "peek")?;
        assert_eq!(idx, 2);

        Ok(())
    }

    #[test]
    fn assert_bounded_counter_accepts_exact_limit() -> Result<(), Error> {
        let mut idx = 0;

        assert_bounded_counter(&mut idx, 1, "peek")?;
        assert_eq!(idx, 1);

        Ok(())
    }

    #[test]
    #[should_panic(expected = "peek record limit exceeded: 1")]
    fn assert_bounded_counter_errors_when_limit_exceeded() {
        let mut idx = 1;

        assert_bounded_counter(&mut idx, 1, "peek").unwrap();
    }

    #[test]
    #[should_panic(expected = "peek record counter overflowed at u32::MAX")]
    fn assert_bounded_counter_errors_on_overflow() {
        let mut idx = u32::MAX;

        assert_bounded_counter(&mut idx, u32::MAX, "peek").unwrap();
    }

    #[test]
    fn assert_nonzero_counter_accepts_nonzero() -> Result<(), Error> {
        assert_nonzero_counter(1, "records")
    }

    #[test]
    #[should_panic(expected = "no records found!")]
    fn assert_nonzero_counter_errors_on_zero() {
        assert_nonzero_counter(0, "records").unwrap();
    }

    #[test]
    fn assert_record_data_capacity_accepts_exact_limit() -> Result<(), Error> {
        assert_record_data_capacity(MAX_RECORD_CAPACITY_BYTES, MAX_RECORD_CAPACITY_BYTES, "peek")
    }

    #[test]
    #[should_panic(expected = "peek record capacity limit exceeded")]
    fn assert_record_data_capacity_errors_when_limit_exceeded() {
        assert_record_data_capacity(
            MAX_RECORD_CAPACITY_BYTES + 1,
            MAX_RECORD_CAPACITY_BYTES,
            "peek",
        )
        .unwrap();
    }

    #[test]
    fn assert_flag_accepts_true() -> Result<(), Error> {
        assert_flag(true, "should not fail")
    }

    #[test]
    #[should_panic(expected = "flag failed")]
    fn assert_flag_errors_on_false() {
        assert_flag(false, "flag failed").unwrap();
    }
}
