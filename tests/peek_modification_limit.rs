//! Cross-record modification limits for the public peek command.

#[cfg(test)]
mod tests {
    use std::rc::Rc;

    use nanalogue_core::{Error, peek};
    use rust_htslib::bam;

    /// Supplies a real target array for rust-htslib's header accessors.
    fn header() -> bam::HeaderView {
        bam::HeaderView::from_header(
            bam::Header::new().push_record(
                bam::header::HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", "contig")
                    .push_tag(b"LN", 1),
            ),
        )
    }

    /// Creates a one-base unmapped record with one independently named modification.
    fn modified_record(code: u32) -> Rc<bam::Record> {
        let mut record = bam::Record::new();
        record.set(b"limit_read", None, b"C", &[30]);
        record.set_flags(4);
        record.set_tid(-1);
        record.set_pos(-1);
        record
            .push_aux(b"MM", bam::record::Aux::String(&format!("C+{code}?,0;")))
            .expect("valid MM tag");
        record
            .push_aux(b"ML", bam::record::Aux::ArrayU8((&[200u8][..]).into()))
            .expect("one probability for one call");
        Rc::new(record)
    }

    /// Exactly 100 distinct identities fit, and repeated calls do not consume slots.
    #[test]
    fn distinct_limit_accepts_duplicates_and_sorts_output() {
        let header = header();
        // Codes above the ASCII range remain numeric when displayed. Each record
        // has only one type, so this targets peek's aggregate limit rather than
        // the separate per-record parser limit. Reverse order checks sorting.
        let records = (1000..1100)
            .rev()
            .chain([1099, 1000, 1057])
            .map(|code| Ok(modified_record(code)));
        let mut output = Vec::new();
        peek::run(&mut output, &header, records).expect("100 types plus duplicates fit");
        let text = String::from_utf8(output).expect("ASCII output");
        let (preamble, modifications) = text
            .split_once("modifications:\n")
            .expect("modifications heading");
        assert_eq!(preamble, "contigs_and_lengths:\ncontig\t1\n\n");
        let expected: Vec<String> = (1000..1100).map(|code| format!("C+{code}")).collect();
        assert_eq!(
            modifications.lines().collect::<Vec<_>>(),
            expected,
            "every distinct type must appear exactly once in lexical order"
        );
    }

    /// The first new identity beyond the aggregate limit fails before printing mods.
    #[test]
    fn first_excess_type_stops_before_consuming_later_records() {
        let header = header();
        let mut consumed = 0;
        let records = (1000..1100).chain([1000, 1100, 1101]).map(|code| {
            consumed += 1;
            Ok(modified_record(code))
        });
        let mut output = Vec::new();
        let error = peek::run(&mut output, &header, records)
            .expect_err("the 101st distinct type must exceed the aggregate limit");
        assert!(
            matches!(error, Error::InvalidState(message)
                if message == "peek modification limit exceeded: > 100"),
            "peek, not per-record parsing, must reject the excess type"
        );
        assert_eq!(
            consumed, 102,
            "the duplicate is allowed and the record after the excess is not consumed"
        );
        assert_eq!(output, b"contigs_and_lengths:\ncontig\t1\n\n");
        // Once insertion rejects the 101st type, the later sorted_mods length
        // assertion cannot fail through this API. Do not manufacture an invalid
        // HashSet or change production limits to reach that impossible state.
    }

    /// The same numeric code on the other strand consumes a distinct slot.
    #[test]
    fn opposite_strand_counts_as_a_distinct_type() {
        let header = header();
        let mut extra = modified_record(1000);
        let record = Rc::get_mut(&mut extra).expect("new record has a single owner");
        record.remove_aux(b"MM").expect("fixture has an MM tag");
        record
            .push_aux(b"MM", bam::record::Aux::String("C-1000?,0;"))
            .expect("opposite strand MM tag");
        let records = (1000..1100)
            .map(|code| Ok(modified_record(code)))
            .chain([Ok(extra)]);
        let mut output = Vec::new();
        let error = peek::run(&mut output, &header, records)
            .expect_err("opposite strand must not be deduplicated with the same code");
        assert!(
            matches!(error, Error::InvalidState(message)
                if message == "peek modification limit exceeded: > 100"),
            "strand is part of the modification identity"
        );
        assert_eq!(output, b"contigs_and_lengths:\ncontig\t1\n\n");
    }
}
