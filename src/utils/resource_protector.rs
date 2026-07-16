//! Helpers for (1) protecting loops against runaway iteration.
//! and for (2) protecting against memory exhaustion.

use crate::Error;

/// Increment a counter and assert that it stays within a configured bound.
///
/// # Errors
/// Returns an error if the counter exceeds `max`.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "no overflow as bounded by max"
)]
pub fn assert_bounded_counter(idx: &mut u32, max: u32, what: &str) -> Result<(), Error> {
    if *idx < max {
        *idx += 1;
    } else {
        return Err(Error::InvalidState(format!("{what} limit exceeded: {max}")));
    }
    Ok(())
}

/// Assert that a counter is non-zero.
///
/// # Errors
/// Returns an error if `idx` is zero.
pub fn assert_nonzero_counter(idx: u32, msg: &str) -> Result<(), Error> {
    if idx == 0 {
        return Err(Error::InvalidState(msg.to_owned()));
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

/// Verify that a read id is valid.
/// We apply stricter conditions than the BAM standard i.e. the BAM standard allows
/// read ids up to 255 characters I believe, and allows the various quotes shown
/// below, and allows a read id starting with a #. We deliberately apply stricter
/// standards because we don't want to run into problems downstream
/// e.g. in a table with the first column as read id, a read starting with a '#'
/// could be interpreted as a comment.
///
/// # Errors
/// - if empty
/// - if above a max read id length
/// - if non ASCII characters or strange characters
/// - if starts with a reserved leading character
/// - if contains a quote-like character
#[expect(clippy::else_if_without_else, reason = "simple enough structure")]
pub fn assert_valid_read_id(qname: &[u8], max_len: u8) -> Result<(), Error> {
    #[expect(
        clippy::indexing_slicing,
        reason = "the first branch returns on empty input, so later `qname[0]` is guarded"
    )]
    if qname.is_empty() {
        return Err(Error::InvalidReadID("read id is blank".to_owned()));
    } else if qname.len() > usize::from(max_len) {
        return Err(Error::InvalidState(format!(
            "error in setting read id, length > {max_len}"
        )));
    } else if qname.starts_with(b"#") {
        return Err(Error::InvalidState(
            "we do not accept read ids starting with a # symbol".to_owned(),
        ));
    } else if matches!(qname[0], b'=' | b'+' | b'-' | b'@') {
        // These are the classic CSV/spreadsheet formula-injection trigger characters (Excel/Sheets).
        // We don't want to deal with these issues.
        return Err(Error::InvalidState(
            "we do not accept read ids starting with reserved leading characters".to_owned(),
        ));
    }
    for k in qname {
        if (0..33).contains(k) || (127..).contains(k) || *k == b'`' || *k == b'"' || *k == b'\'' {
            return Err(Error::InvalidReadID(
                "read_id contains strange characters and/or quotes!".to_owned(),
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        assert_bounded_counter, assert_flag, assert_nonzero_counter, assert_record_data_capacity,
        assert_valid_read_id,
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

        assert_bounded_counter(&mut idx, 1, "peek record").unwrap();
    }

    #[test]
    #[should_panic(expected = "peek record limit exceeded: 4294967295")]
    fn assert_bounded_counter_errors_on_overflow() {
        let mut idx = u32::MAX;
        assert_bounded_counter(&mut idx, u32::MAX, "peek record").unwrap();
    }

    #[test]
    fn assert_nonzero_counter_accepts_nonzero() -> Result<(), Error> {
        assert_nonzero_counter(1, "records")
    }

    #[test]
    #[should_panic(expected = "no records found!")]
    fn assert_nonzero_counter_errors_on_zero() {
        assert_nonzero_counter(0, "no records found!").unwrap();
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

    #[test]
    fn assert_valid_read_id_accepts_valid_read_id() -> Result<(), Error> {
        assert_valid_read_id(b"read-123/abc", 20)
    }

    #[test]
    #[should_panic(expected = "read id is blank")]
    fn assert_valid_read_id_rejects_blank_read_id() {
        assert_valid_read_id(b"", 20).unwrap();
    }

    #[test]
    #[should_panic(expected = "error in setting read id, length > 4")]
    fn assert_valid_read_id_rejects_overlength_read_id() {
        assert_valid_read_id(b"read-123", 4).unwrap();
    }

    #[test]
    #[should_panic(expected = "read_id contains strange characters and/or quotes!")]
    fn assert_valid_read_id_rejects_control_characters() {
        assert_valid_read_id(b"read\n123", 20).unwrap();
    }

    #[test]
    #[should_panic(expected = "read_id contains strange characters and/or quotes!")]
    fn assert_valid_read_id_rejects_non_ascii_bytes() {
        assert_valid_read_id(b"read\x7f123", 20).unwrap();
    }

    #[test]
    #[should_panic(expected = "read_id contains strange characters and/or quotes!")]
    fn assert_valid_read_id_rejects_quote_like_characters() {
        assert_valid_read_id(b"read`123", 20).unwrap();
    }

    #[test]
    #[should_panic(expected = "we do not accept read ids starting with a # symbol")]
    fn assert_valid_read_id_rejects_leading_hash() {
        assert_valid_read_id(b"#read123", 20).unwrap();
    }

    #[test]
    #[should_panic(
        expected = "we do not accept read ids starting with reserved leading characters"
    )]
    fn assert_valid_read_id_rejects_leading_equals() {
        assert_valid_read_id(b"=read123", 20).unwrap();
    }

    #[test]
    #[should_panic(
        expected = "we do not accept read ids starting with reserved leading characters"
    )]
    fn assert_valid_read_id_rejects_leading_plus() {
        assert_valid_read_id(b"+read123", 20).unwrap();
    }

    #[test]
    #[should_panic(
        expected = "we do not accept read ids starting with reserved leading characters"
    )]
    fn assert_valid_read_id_rejects_leading_hyphen() {
        assert_valid_read_id("-read123".as_bytes(), 20).unwrap();
    }

    #[test]
    #[should_panic(
        expected = "we do not accept read ids starting with reserved leading characters"
    )]
    fn assert_valid_read_id_rejects_leading_at_sign() {
        assert_valid_read_id(b"@read123", 20).unwrap();
    }

    #[test]
    #[should_panic(expected = "read_id contains strange characters and/or quotes!")]
    fn assert_valid_read_id_rejects_leading_tab() {
        assert_valid_read_id(b"\tread123", 20).unwrap();
    }

    #[test]
    #[should_panic(expected = "read_id contains strange characters and/or quotes!")]
    fn assert_valid_read_id_rejects_leading_nul() {
        assert_valid_read_id(b"\0read123", 20).unwrap();
    }

    #[test]
    fn assert_valid_read_id_accepts_hash_in_middle() -> Result<(), Error> {
        assert_valid_read_id(b"read#123", 20)
    }
}
