//! Helpers for (1) protecting loops against runaway iteration.
//! and for (2) protecting against memory exhaustion.

use crate::Error;

/// Increment a counter and ensure that it stays within a configured bound.
///
/// # Errors
/// Returns an error if the counter exceeds `max`.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "no overflow as bounded by max"
)]
pub fn ensure_bounded_counter(idx: &mut u32, max: u32, what: &str) -> Result<(), Error> {
    if *idx < max {
        *idx += 1;
    } else {
        return Err(Error::InvalidState(format!("{what} limit exceeded: {max}")));
    }
    Ok(())
}

/// Ensure that a counter is non-zero.
///
/// # Errors
/// Returns an error if `idx` is zero.
pub fn ensure_nonzero_counter(idx: u32, msg: &str) -> Result<(), Error> {
    if idx == 0 {
        return Err(Error::InvalidState(msg.to_owned()));
    }

    Ok(())
}

/// Ensure that a BAM record's internal data capacity stays within a configured bound.
///
/// # Errors
/// Returns an error if `m_data` exceeds `max`.
pub fn ensure_record_data_capacity(m_data: u32, max: u32, what: &str) -> Result<(), Error> {
    if m_data > max {
        return Err(Error::InvalidState(format!(
            "{what} record capacity limit exceeded: {max}"
        )));
    }

    Ok(())
}

/// Ensure that a flag is true.
///
/// # Errors
/// Returns an error if `flag` is false.
pub fn ensure_flag(flag: bool, msg: &str) -> Result<(), Error> {
    if !flag {
        return Err(Error::InvalidState(msg.to_owned()));
    }

    Ok(())
}

/// Verify that a read-id-like or contig-like identifier is safe for downstream use.
///
/// This helper deliberately applies a stricter shared policy than the underlying
/// BAM or reference-name specifications so that callers can avoid downstream
/// issues such as comment parsing, spreadsheet formula injection, and awkward
/// punctuation in text-based exports.
#[expect(clippy::else_if_without_else, reason = "simple enough structure")]
fn ensure_valid_identifier(value: &[u8], max_len: u8, what: &str) -> Result<(), Error> {
    #[expect(
        clippy::indexing_slicing,
        reason = "the first branch returns on empty input, so later `value[0]` is guarded"
    )]
    if value.is_empty() {
        return Err(Error::InvalidState(format!("{what} is blank")));
    } else if value.len() > usize::from(max_len) {
        return Err(Error::InvalidState(format!(
            "error in setting {what}, length > {max_len}"
        )));
    } else if matches!(value[0], b'#' | b'*' | b'=' | b'+' | b'-' | b'@') {
        // These are reserved either by our downstream safety rules or by the
        // reference-name specification for the first character.
        return Err(Error::InvalidState(format!(
            "we do not accept {what} values starting with reserved leading characters"
        )));
    }
    for byte in value {
        if (0..33).contains(byte)
            || (127..).contains(byte)
            || matches!(
                *byte,
                b'#' | b'`'
                    | b'"'
                    | b'\''
                    | b'\\'
                    | b','
                    | b'('
                    | b')'
                    | b'['
                    | b']'
                    | b'{'
                    | b'}'
                    | b'<'
                    | b'>'
            )
        {
            return Err(Error::InvalidState(format!(
                "{what} contains forbidden characters"
            )));
        }
    }
    Ok(())
}

/// Verify that a read id is valid under the crate's shared safe-identifier policy.
///
/// This policy is intentionally stricter than the BAM standard so that read ids
/// remain safe to reuse in downstream text formats and user-facing outputs.
///
/// # Errors
/// - if empty
/// - if above a max read id length
/// - if it contains non-ASCII, control, or forbidden punctuation characters
/// - if it starts with a reserved leading character
pub fn ensure_valid_read_id(qname: &[u8], max_len: u8) -> Result<(), Error> {
    ensure_valid_identifier(qname, max_len, "read_id").map_err(|err| {
        let Error::InvalidState(msg) = err else {
            unreachable!("shared identifier validator must return InvalidState")
        };
        Error::InvalidReadID(msg)
    })
}

/// Verify that a contig name is valid under the crate's shared safe-identifier policy.
///
/// This policy is intentionally at least as strict as the reference-name
/// specification and may reject additional punctuation for downstream safety.
///
/// # Errors
/// - if empty
/// - if above a max contig length
/// - if it contains non-ASCII, control, or forbidden punctuation characters
/// - if it starts with a reserved leading character
pub fn ensure_valid_contig(contig: &[u8], max_len: u8) -> Result<(), Error> {
    ensure_valid_identifier(contig, max_len, "contig").map_err(|err| {
        let Error::InvalidState(msg) = err else {
            unreachable!("shared identifier validator must return InvalidState")
        };
        Error::InvalidContig(msg)
    })
}

#[cfg(test)]
mod tests {
    use super::{
        ensure_bounded_counter, ensure_flag, ensure_nonzero_counter, ensure_record_data_capacity,
        ensure_valid_contig, ensure_valid_read_id,
    };
    use crate::Error;
    use crate::constants::shared::MAX_RECORD_CAPACITY_BYTES;

    #[test]
    fn ensure_bounded_counter_increments_within_bound() -> Result<(), Error> {
        let mut idx = 0;

        ensure_bounded_counter(&mut idx, 2, "peek")?;
        assert_eq!(idx, 1);

        ensure_bounded_counter(&mut idx, 2, "peek")?;
        assert_eq!(idx, 2);

        Ok(())
    }

    #[test]
    fn ensure_bounded_counter_accepts_exact_limit() -> Result<(), Error> {
        let mut idx = 0;

        ensure_bounded_counter(&mut idx, 1, "peek")?;
        assert_eq!(idx, 1);

        Ok(())
    }

    #[test]
    #[should_panic(expected = "peek record limit exceeded: 1")]
    fn ensure_bounded_counter_errors_when_limit_exceeded() {
        let mut idx = 1;

        ensure_bounded_counter(&mut idx, 1, "peek record").unwrap();
    }

    #[test]
    #[should_panic(expected = "peek record limit exceeded: 4294967295")]
    fn ensure_bounded_counter_errors_on_overflow() {
        let mut idx = u32::MAX;
        ensure_bounded_counter(&mut idx, u32::MAX, "peek record").unwrap();
    }

    #[test]
    fn ensure_nonzero_counter_accepts_nonzero() -> Result<(), Error> {
        ensure_nonzero_counter(1, "records")
    }

    #[test]
    #[should_panic(expected = "no records found!")]
    fn ensure_nonzero_counter_errors_on_zero() {
        ensure_nonzero_counter(0, "no records found!").unwrap();
    }

    #[test]
    fn ensure_record_data_capacity_accepts_exact_limit() -> Result<(), Error> {
        ensure_record_data_capacity(MAX_RECORD_CAPACITY_BYTES, MAX_RECORD_CAPACITY_BYTES, "peek")
    }

    #[test]
    #[should_panic(expected = "peek record capacity limit exceeded")]
    fn ensure_record_data_capacity_errors_when_limit_exceeded() {
        ensure_record_data_capacity(
            MAX_RECORD_CAPACITY_BYTES + 1,
            MAX_RECORD_CAPACITY_BYTES,
            "peek",
        )
        .unwrap();
    }

    #[test]
    fn ensure_flag_accepts_true() -> Result<(), Error> {
        ensure_flag(true, "should not fail")
    }

    #[test]
    #[should_panic(expected = "flag failed")]
    fn ensure_flag_errors_on_false() {
        ensure_flag(false, "flag failed").unwrap();
    }

    #[test]
    fn ensure_valid_read_id_accepts_valid_read_id() -> Result<(), Error> {
        ensure_valid_read_id(b"read-123/abc", 20)
    }

    #[test]
    #[should_panic(expected = "read_id is blank")]
    fn ensure_valid_read_id_rejects_blank_read_id() {
        ensure_valid_read_id(b"", 20).unwrap();
    }

    #[test]
    #[should_panic(expected = "error in setting read_id, length > 4")]
    fn ensure_valid_read_id_rejects_overlength_read_id() {
        ensure_valid_read_id(b"read-123", 4).unwrap();
    }

    #[test]
    fn ensure_valid_read_id_rejects_forbidden_characters() {
        for rejected in [
            b'"', b'\'', b'`', b'\t', b'\n', b'\r', b'\0', b'\\', b',', b'(', b')', b'[', b']',
            b'{', b'}', b'<', b'>',
        ] {
            let read_id = [b'r', b'e', rejected, b'i', b'd', b'_', b'1'];
            assert!(ensure_valid_read_id(&read_id, 20).is_err());
        }
    }

    #[test]
    #[should_panic(expected = "read_id contains forbidden characters")]
    fn ensure_valid_read_id_rejects_non_ascii_bytes() {
        ensure_valid_read_id(b"read\x7f123", 20).unwrap();
    }

    #[test]
    #[should_panic(
        expected = "we do not accept read_id values starting with reserved leading characters"
    )]
    fn ensure_valid_read_id_rejects_leading_hash() {
        ensure_valid_read_id(b"#read123", 20).unwrap();
    }

    #[test]
    #[should_panic(
        expected = "we do not accept read_id values starting with reserved leading characters"
    )]
    fn ensure_valid_read_id_rejects_leading_equals() {
        ensure_valid_read_id(b"=read123", 20).unwrap();
    }

    #[test]
    #[should_panic(
        expected = "we do not accept read_id values starting with reserved leading characters"
    )]
    fn ensure_valid_read_id_rejects_leading_plus() {
        ensure_valid_read_id(b"+read123", 20).unwrap();
    }

    #[test]
    #[should_panic(
        expected = "we do not accept read_id values starting with reserved leading characters"
    )]
    fn ensure_valid_read_id_rejects_leading_asterisk() {
        ensure_valid_read_id(b"*read123", 20).unwrap();
    }

    #[test]
    #[should_panic(
        expected = "we do not accept read_id values starting with reserved leading characters"
    )]
    fn ensure_valid_read_id_rejects_leading_hyphen() {
        ensure_valid_read_id("-read123".as_bytes(), 20).unwrap();
    }

    #[test]
    #[should_panic(
        expected = "we do not accept read_id values starting with reserved leading characters"
    )]
    fn ensure_valid_read_id_rejects_leading_at_sign() {
        ensure_valid_read_id(b"@read123", 20).unwrap();
    }

    #[test]
    #[should_panic(expected = "read_id contains forbidden characters")]
    fn ensure_valid_read_id_rejects_leading_tab() {
        ensure_valid_read_id(b"\tread123", 20).unwrap();
    }

    #[test]
    #[should_panic(expected = "read_id contains forbidden characters")]
    fn ensure_valid_read_id_rejects_leading_nul() {
        ensure_valid_read_id(b"\0read123", 20).unwrap();
    }

    #[test]
    #[should_panic(expected = "read_id contains forbidden characters")]
    fn ensure_valid_read_id_rejects_hash_in_middle() {
        ensure_valid_read_id(b"read#123", 20).unwrap();
    }

    #[test]
    fn ensure_valid_contig_accepts_valid_contig() -> Result<(), Error> {
        ensure_valid_contig(b"chr1_alt", 20)
    }

    #[test]
    #[should_panic(expected = "contig is blank")]
    fn ensure_valid_contig_rejects_blank_contig() {
        ensure_valid_contig(b"", 20).unwrap();
    }

    #[test]
    #[should_panic(expected = "error in setting contig, length > 4")]
    fn ensure_valid_contig_rejects_overlength_contig() {
        ensure_valid_contig(b"chr123", 4).unwrap();
    }

    #[test]
    fn ensure_valid_contig_rejects_forbidden_characters() {
        for rejected in [
            b'"', b'\'', b'`', b'\t', b'\n', b'\r', b'\0', b'\\', b',', b'(', b')', b'[', b']',
            b'{', b'}', b'<', b'>',
        ] {
            let contig = [b'r', b'e', rejected, b'i', b'd', b'_', b'1'];
            assert!(ensure_valid_contig(&contig, 20).is_err());
        }
    }

    #[test]
    #[should_panic(
        expected = "we do not accept contig values starting with reserved leading characters"
    )]
    fn ensure_valid_contig_rejects_leading_hash() {
        ensure_valid_contig(b"#chr1", 20).unwrap();
    }

    #[test]
    #[should_panic(
        expected = "we do not accept contig values starting with reserved leading characters"
    )]
    fn ensure_valid_contig_rejects_reserved_leading_character() {
        ensure_valid_contig(b"=chr1", 20).unwrap();
    }

    #[test]
    #[should_panic(
        expected = "we do not accept contig values starting with reserved leading characters"
    )]
    fn ensure_valid_contig_rejects_leading_asterisk() {
        ensure_valid_contig(b"*chr1", 20).unwrap();
    }
}
