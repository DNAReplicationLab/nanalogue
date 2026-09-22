#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Public `read_info::run` boundary, formatting, and error-output contracts.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::{Error, InputMods, OptionalTag, OrdPair, ThresholdState, read_info};
    use rust_htslib::bam::Record;
    use rust_htslib::bam::record::Aux;
    use rust_htslib::errors::Error as HtslibError;
    use std::io;
    use std::rc::Rc;

    /// Builds an unmapped record without requiring a header or alignment data.
    fn unmapped_record(read_id: &str, sequence: &[u8], mapq: u8) -> Record {
        let mut record = Record::new();
        record.set(
            read_id.as_bytes(),
            None,
            sequence,
            &vec![30; sequence.len()],
        );
        record.set_flags(4);
        record.set_tid(-1);
        record.set_pos(-1);
        record.set_mapq(mapq);
        record
    }

    /// Adds four asymmetric modification probabilities to an all-T record.
    fn tagged_record() -> Record {
        let mut record = unmapped_record("threshold_read", b"TTTT", 17);
        record
            .push_aux(b"MM", Aux::String("T+T,0,0,0,0;"))
            .expect("the MM tag is valid");
        record
            .push_aux(b"ML", Aux::ArrayU8((&[100u8, 150, 210, 230][..]).into()))
            .expect("one probability is supplied for each call");
        record
    }

    /// Wraps records in the public iterator item type.
    fn records(input: Vec<Record>) -> Vec<Result<Rc<Record>, HtslibError>> {
        input
            .into_iter()
            .map(|record| Ok(Rc::new(record)))
            .collect()
    }

    /// Runs the condensed mode with a selected incoming probability filter.
    fn condensed_output(filter: ThresholdState) -> String {
        let mut mods = InputMods::<OptionalTag>::default();
        mods.mod_prob_filter = filter;
        let mut output = Vec::new();
        read_info::run(&mut output, records(vec![tagged_record()]), mods, None)
            .expect("the tagged record is valid");
        String::from_utf8(output).expect("read-info output is UTF-8")
    }

    /// Condensed mode always imposes the 0.5 counting floor while preserving
    /// a stricter floor and an independently supplied exclusion interval.
    #[test]
    fn condensed_mode_combines_probability_boundaries() {
        let uncertain = OrdPair::new(170, 220).expect("ordered uncertainty interval");
        let cases = [
            (
                ThresholdState::GtEq(0),
                "T+T:3;(probabilities >= 0.5020, PHRED base qual >= 0)",
            ),
            (
                ThresholdState::GtEq(200),
                "T+T:2;(probabilities >= 0.7843, PHRED base qual >= 0)",
            ),
            (
                ThresholdState::InvertGtEqLtEq(uncertain),
                "T+T:2;(probabilities >= 0.5020 and (probabilities < 0.6667 or > 0.8627), PHRED base qual >= 0)",
            ),
            (
                ThresholdState::Both((200, uncertain)),
                "T+T:1;(probabilities >= 0.7843 and (probabilities < 0.6667 or > 0.8627), PHRED base qual >= 0)",
            ),
        ];

        for (filter, mod_count) in cases {
            let expected = format!(
                "[\n{{\n\t\"read_id\": \"threshold_read\",\n\t\"sequence_length\": 4,\n\t\"mapq\": 17,\n\t\"alignment_type\": \"unmapped\",\n\t\"mod_count\": \"{mod_count}\"\n}}\n]\n"
            );
            assert_eq!(condensed_output(filter), expected, "filter: {filter}");
        }
    }

    /// Compact and pretty detailed modes preserve record order and every
    /// asymmetric read field while differing only in JSON whitespace.
    #[test]
    fn detailed_modes_have_exact_json_layouts() {
        let compact = "[\n\
{\"alignment_type\":\"unmapped\",\"mod_table\":[],\"read_id\":\"short\",\"mapq\":7,\"seq_len\":3},\n\
{\"alignment_type\":\"unmapped\",\"mod_table\":[],\"read_id\":\"longer\",\"mapq\":29,\"seq_len\":5}\n\
]\n";
        let pretty = concat!(
            "[\n",
            "{\n",
            "  \"alignment_type\": \"unmapped\",\n",
            "  \"mod_table\": [],\n",
            "  \"read_id\": \"short\",\n",
            "  \"mapq\": 7,\n",
            "  \"seq_len\": 3\n",
            "},\n",
            "{\n",
            "  \"alignment_type\": \"unmapped\",\n",
            "  \"mod_table\": [],\n",
            "  \"read_id\": \"longer\",\n",
            "  \"mapq\": 29,\n",
            "  \"seq_len\": 5\n",
            "}\n",
            "]\n",
        );

        for (detailed, expected) in [(Some(false), compact), (Some(true), pretty)] {
            let input = records(vec![
                unmapped_record("short", b"ACG", 7),
                unmapped_record("longer", b"TGCAT", 29),
            ]);
            let mut output = Vec::new();
            read_info::run(
                &mut output,
                input,
                InputMods::<OptionalTag>::default(),
                detailed,
            )
            .expect("both unmodified records are serializable");

            assert_eq!(
                String::from_utf8(output).expect("read-info output is UTF-8"),
                expected,
                "detailed mode {detailed:?}"
            );
        }
    }

    /// Iterator failures before the first record leave output untouched.
    #[test]
    fn iterator_error_before_first_record_writes_nothing() {
        let input = vec![Err(HtslibError::BamInvalidIndex {
            target: "absent.bai".to_owned(),
        })];
        let mut output = Vec::new();
        let error = read_info::run(
            &mut output,
            input,
            InputMods::<OptionalTag>::default(),
            None,
        )
        .expect_err("an htslib iterator error must abort read-info");

        let Error::RustHtslibError(source) = error else {
            unreachable!("expected RustHtslibError, got {error:?}")
        };
        assert!(
            matches!(*source, HtslibError::BamInvalidIndex { target }
                if target == "absent.bai"),
            "the source error must be preserved"
        );
        assert!(
            output.is_empty(),
            "the JSON array must not begin before the first record is available"
        );
    }

    /// Record validation finishes before the opener or next delimiter is
    /// emitted, so a failing record contributes no bytes to streamed output.
    #[test]
    fn validation_error_does_not_emit_a_record_delimiter() {
        let mut first_output = Vec::new();
        let first_error = read_info::run(
            &mut first_output,
            records(vec![unmapped_record("invalid_first", b"", 0)]),
            InputMods::<OptionalTag>::default(),
            Some(false),
        )
        .expect_err("a zero-length first record must fail validation");
        assert!(matches!(first_error, Error::ZeroSeqLen(_)));
        assert!(
            first_output.is_empty(),
            "first-record validation must precede the array opener"
        );

        let mut later_output = Vec::new();
        let later_error = read_info::run(
            &mut later_output,
            records(vec![
                unmapped_record("valid", b"ACG", 7),
                unmapped_record("invalid_second", b"", 0),
            ]),
            InputMods::<OptionalTag>::default(),
            Some(false),
        )
        .expect_err("a zero-length second record must fail validation");
        assert!(matches!(later_error, Error::ZeroSeqLen(_)));
        assert_eq!(
            later_output,
            b"[\n{\"alignment_type\":\"unmapped\",\"mod_table\":[],\"read_id\":\"valid\",\"mapq\":7,\"seq_len\":3}",
            "the failing second record must not add a comma or newline"
        );
    }

    /// A writer that accepts a fixed prefix and then fails once.
    struct FailingWriter {
        output: Vec<u8>,
        remaining: usize,
    }

    impl io::Write for FailingWriter {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            if self.remaining == 0 {
                return Err(io::Error::other("body write failed"));
            }
            let accepted = buf.len().min(self.remaining);
            self.output.extend_from_slice(
                buf.get(..accepted)
                    .expect("accepted bytes are within the input"),
            );
            self.remaining = self
                .remaining
                .checked_sub(accepted)
                .expect("accepted bytes never exceed the remaining budget");
            Ok(accepted)
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    /// Failure while writing the first JSON body is returned as the original
    /// I/O error, after exactly the array opener and first-record newline.
    #[test]
    fn body_write_error_preserves_source_and_prefix() {
        let mut writer = FailingWriter {
            output: Vec::new(),
            remaining: 2,
        };
        let error = read_info::run(
            &mut writer,
            records(vec![unmapped_record("write_fail", b"ACG", 11)]),
            InputMods::<OptionalTag>::default(),
            None,
        )
        .expect_err("the body write must fail after the two-byte prefix");

        let Error::InputOutputError(source) = error else {
            unreachable!("expected InputOutputError, got {error:?}")
        };
        assert_eq!(source.to_string(), "body write failed");
        assert_eq!(writer.output, b"[\n");
    }
}
