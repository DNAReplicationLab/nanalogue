#![cfg_attr(coverage_nightly, feature(coverage_attribute))]

//! Focused `run`-level tests for the public subcommands, built from in-memory
//! `rust_htslib::bam::Record` values instead of fixture files.
//!
//! Values are deliberately asymmetric (different read lengths, positions,
//! qualities, and window counts) so a swapped field or a misplaced fallback
//! cannot satisfy the assertions.

#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use nanalogue_core::analysis::threshold_and_mean;
    use nanalogue_core::{
        Error, F32Bw0and1, InputMods, InputWindowing, InputWindowingBuilder, OptionalTag,
        RequiredTag, SeqDisplayOptions, find_modified_reads, peek, read_stats, reads_table,
        window_reads,
    };
    use rust_htslib::bam::ext::BamRecordExtensions as _;
    use rust_htslib::bam::record::{Aux, Cigar, CigarString};
    use rust_htslib::bam::{Header, HeaderView, Record, header::HeaderRecord};
    use rust_htslib::errors::Error as HtslibError;
    use std::io;
    use std::rc::Rc;
    use std::sync::Arc;

    /// Name of the single contig used by every in-memory header.
    const CONTIG: &str = "ctg1";
    /// Length stored in the header for [`CONTIG`].
    const CONTIG_LEN: u32 = 100;
    /// Position of the first aligned base of every mapped record.
    const START: i64 = 10;

    /// Builds the one-contig header that records need for `Record::contig()`.
    fn header() -> HeaderView {
        HeaderView::from_header(
            Header::new().push_record(
                HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", CONTIG)
                    .push_tag(b"LN", CONTIG_LEN),
            ),
        )
    }

    /// Builds a mapped record on [`CONTIG`] starting at [`START`].
    fn mapped_record(read_id: &str, cigar: Cigar, seq: &[u8], qual: &[u8]) -> Record {
        let mut record = Record::new();
        record.set(
            read_id.as_bytes(),
            Some(&CigarString(vec![cigar])),
            seq,
            qual,
        );
        record.set_header(Arc::new(header()));
        record.set_flags(0); // `Record::new` starts out unmapped
        record.set_tid(0);
        record.set_pos(START);
        record
    }

    /// Builds an unmapped record whose alignment fields stay unset.
    fn unmapped_record(read_id: &str, seq: &[u8]) -> Record {
        let qual = vec![255u8; seq.len()];
        let mut record = Record::new();
        record.set(read_id.as_bytes(), None, seq, &qual);
        record.set_header(Arc::new(header()));
        record.set_flags(4);
        record.set_tid(-1);
        record.set_pos(-1);
        record
    }

    /// Builds an all-`T` mapped record carrying one `T+T` group and its `ML` calls.
    fn modified_record(read_id: &str, mm: &str, ml: &[u8], seq: &[u8], qual: &[u8]) -> Record {
        let cigar = Cigar::Match(u32::try_from(seq.len()).expect("test sequences fit in u32"));
        let mut record = mapped_record(read_id, cigar, seq, qual);
        record
            .push_aux(b"MM", Aux::String(mm))
            .expect("MM tag must fit in the record");
        record
            .push_aux(b"ML", Aux::ArrayU8(ml.into()))
            .expect("one probability per T+T call must fit in the record");
        record
    }

    /// Wraps owned records the way the subcommand `run` iterators expect them.
    fn records(input: Vec<Record>) -> Vec<Result<Rc<Record>, HtslibError>> {
        input
            .into_iter()
            .map(|record| Ok(Rc::new(record)))
            .collect()
    }

    /// Four-position windows slid one position at a time.
    fn window_options() -> InputWindowing {
        InputWindowingBuilder::default()
            .win(4)
            .step(1)
            .build()
            .expect("a four-position window slid one position at a time is valid")
    }

    /// Output written before any record is consumed, shared by the peek tests.
    fn contig_preamble() -> String {
        format!("contigs_and_lengths:\n{CONTIG}\t{CONTIG_LEN}\n\n")
    }

    /// An unsupported flag aborts `peek` with `NotImplemented` naming the read.
    #[test]
    fn peek_propagates_unsupported_flag() {
        let mut record = mapped_record("flag_read", Cigar::Match(5), b"ACGTA", &[30u8; 5]);
        record.set_flags(0x1); // paired reads are explicitly unsupported

        let mut output = Vec::new();
        let error = peek::run(&mut output, &header(), records(vec![record]).into_iter())
            .expect_err("a paired record must abort peek");

        let rendered = format!("{error:?}");
        assert!(
            matches!(error, Error::NotImplemented(message) if message.contains("flag_read")),
            "a paired record must be NotImplemented naming flag_read, got {rendered}"
        );
        assert_eq!(
            String::from_utf8(output).expect("peek writes UTF-8"),
            contig_preamble(),
            "contig output is written before records are consumed"
        );
    }

    /// An htslib iterator error becomes the project `RustHtslibError`.
    #[test]
    fn peek_converts_htslib_iterator_error() {
        let mut output = Vec::new();
        let error = peek::run(
            &mut output,
            &header(),
            vec![Err(HtslibError::BamInvalidIndex {
                target: "missing_index".to_owned(),
            })]
            .into_iter(),
        )
        .expect_err("an iterator error must abort peek");

        let rendered = format!("{error:?}");
        let Error::RustHtslibError(inner) = error else {
            unreachable!("expected RustHtslibError, got {rendered}")
        };
        assert!(
            matches!(*inner, HtslibError::BamInvalidIndex { target }
                if target.as_str() == "missing_index"),
            "the original htslib error must be preserved unchanged"
        );
        assert_eq!(
            String::from_utf8(output).expect("peek writes UTF-8"),
            contig_preamble(),
            "contig output survives the record-level error"
        );
    }

    /// Zero-length sequences are skipped while contig output still appears.
    #[test]
    fn peek_skips_zero_length_records() {
        let record = unmapped_record("empty_read", b"");

        let mut output = Vec::new();
        peek::run(&mut output, &header(), records(vec![record]).into_iter())
            .expect("zero-length records are skipped, not fatal");

        assert_eq!(
            String::from_utf8(output).expect("peek writes UTF-8"),
            format!("{}modifications:\nNone\n", contig_preamble()),
            "a skipped record still counts as input, so peek reports no modifications"
        );
    }

    /// `read_stats` counts primary, reversed, and unmapped records exactly.
    #[test]
    fn read_stats_counts_mapped_and_unmapped_records() {
        let mut mapped = mapped_record("mapped_read", Cigar::Match(5), b"ACGTA", &[30u8; 5]);
        mapped.set_flags(16); // reverse strand
        let unmapped = unmapped_record("unmapped_read", b"ACGTACG"); // seven bases

        let mut output = Vec::new();
        read_stats::run(&mut output, records(vec![mapped, unmapped]))
            .expect("a mapped and an unmapped record are both countable");

        let expected = "key\tvalue\n\
n_primary_alignments\t1\n\
n_secondary_alignments\t0\n\
n_supplementary_alignments\t0\n\
n_unmapped_reads\t1\n\
n_reversed_reads\t1\n\
align_len_mean\t5\n\
align_len_max\t5\n\
align_len_min\t5\n\
align_len_median\t5\n\
align_len_n50\t5\n\
seq_len_mean\t6\n\
seq_len_max\t7\n\
seq_len_min\t5\n\
seq_len_median\t7\n\
seq_len_n50\t7\n";
        assert_eq!(
            String::from_utf8(output).expect("read_stats writes UTF-8"),
            expected,
            "counts, alignment lengths (5) and sequence lengths (5 and 7) must not be swapped"
        );
    }

    /// Characterises htslib's one-base floor for zero-reference-span CIGARs.
    ///
    /// This is a characterization of current behavior, not an assertion that
    /// `5S` and `5I` semantically align one reference base. Htslib's
    /// `bam_endpos` returns `pos + 1` for these records so that they can be
    /// indexed, and `read_stats` currently observes that synthetic span through
    /// `reference_end`. Rejecting such CIGARs would require inspecting them
    /// before using `reference_end`; `InvalidAlignLength` remains covered by
    /// [`read_stats_rejects_negative_alignment_start`].
    #[test]
    fn read_stats_exposes_htslib_zero_span_floor() {
        let expected = "key\tvalue\n\
n_primary_alignments\t1\n\
n_secondary_alignments\t0\n\
n_supplementary_alignments\t0\n\
n_unmapped_reads\t0\n\
n_reversed_reads\t0\n\
align_len_mean\t1\n\
align_len_max\t1\n\
align_len_min\t1\n\
align_len_median\t1\n\
align_len_n50\t1\n\
seq_len_mean\t5\n\
seq_len_max\t5\n\
seq_len_min\t5\n\
seq_len_median\t5\n\
seq_len_n50\t5\n";

        for (label, cigar) in [("5S", Cigar::SoftClip(5)), ("5I", Cigar::Ins(5))] {
            let record = mapped_record("zero_ref_span", cigar, b"ACGTA", &[30u8; 5]);
            assert_eq!(record.pos(), START, "{label}: test precondition");
            assert_eq!(
                record.reference_end(),
                START + 1,
                "{label}: htslib floors the zero reference span for indexing"
            );

            let mut output = Vec::new();
            read_stats::run(&mut output, records(vec![record])).unwrap_or_else(|error| {
                unreachable!("{label}: expected the current one-base-floor behavior, got {error:?}")
            });

            assert_eq!(
                String::from_utf8(output).expect("read_stats writes UTF-8"),
                expected,
                "{label}: current statistics expose htslib's synthetic one-base span"
            );
        }
    }

    /// A mapped record with a negative start cannot yield a span and is rejected.
    #[test]
    fn read_stats_rejects_negative_alignment_start() {
        let mut record = mapped_record("negative_start", Cigar::Match(5), b"ACGTA", &[30u8; 5]);
        record.set_pos(-1);

        let mut output = Vec::new();
        let error = read_stats::run(&mut output, records(vec![record]))
            .expect_err("a negative alignment start must be rejected");

        let rendered = format!("{error:?}");
        assert!(
            matches!(error, Error::InvalidAlignLength(message)
                if message.contains("negative_start")),
            "a negative start must be InvalidAlignLength naming negative_start, got {rendered}"
        );
        assert!(
            output.is_empty(),
            "statistics are written only after every record is processed"
        );
    }

    /// Too few calls window to nothing and produce no line; a passing record
    /// produces exactly one.
    #[test]
    fn find_modified_reads_reports_only_windowed_records() {
        // Three `T+T` calls cannot fill a four-call window, five calls can.
        // MM distances are one per call, so three calls need three entries.
        let short = modified_record(
            "short_read",
            "T+T,0,0,0;",
            &[200, 200, 200],
            b"TTT",
            &[30u8; 3],
        );
        let long = modified_record(
            "long_read",
            "T+T,0,0,0,0,0;",
            &[200, 200, 200, 200, 200],
            b"TTTTT",
            &[30u8; 5],
        );

        let mut mod_options = InputMods::<RequiredTag>::default();
        mod_options.tag = "T".parse::<RequiredTag>().expect("T is a valid tag");

        let mut output = Vec::new();
        find_modified_reads::run(
            &mut output,
            records(vec![short, long]),
            window_options(),
            &mod_options,
            threshold_and_mean,
            |_: &Vec<F32Bw0and1>| true,
        )
        .expect("a record with no windows is skipped, not fatal");

        assert_eq!(
            String::from_utf8(output).expect("find_modified_reads writes UTF-8"),
            "long_read\n",
            "only the record with full four-call windows may be reported, exactly once"
        );
    }

    /// All-255 qualities keep the sentinel; asymmetric qualities use the
    /// probability-space mean of each distinct window.
    #[test]
    fn window_reads_preserves_sentinel_and_averages_qualities() {
        let sentinel = modified_record(
            "qual_255",
            "T+T,0,0,0,0,0;",
            &[200u8; 5],
            b"TTTTT",
            &[255u8; 5],
        );
        let control = modified_record(
            "mixed_qual",
            "T+T,0,0,0,0,0;",
            &[200u8; 5],
            b"TTTTT",
            &[10, 20, 30, 40, 50],
        );

        let mods = InputMods::<OptionalTag>::default();
        let mut output = Vec::new();
        window_reads::run(
            &mut output,
            records(vec![sentinel, control]),
            window_options(),
            &mods,
            |values: &[u8]| threshold_and_mean(values).map(Into::into),
        )
        .expect("both records produce windowed output");

        let text = String::from_utf8(output).expect("window_reads writes UTF-8");
        let mut lines = text.lines();
        assert_eq!(
            lines.next().expect("header line"),
            "#contig\tref_win_start\tref_win_end\tread_id\twin_val\tstrand\tbase\tmod_strand\t\
mod_type\twin_start\twin_end\tbasecall_qual",
            "the documented header must come first"
        );
        let rows: Vec<Vec<&str>> = lines.map(|line| line.split('\t').collect()).collect();
        assert_eq!(
            rows.len(),
            4,
            "two records with five calls each give two four-call windows"
        );

        let expected = [
            ("qual_255", "0", "4", "255"),
            ("qual_255", "1", "5", "255"),
            ("mixed_qual", "0", "4", "16"),
            ("mixed_qual", "1", "5", "26"),
        ];
        for (row, (read_id, win_start, win_end, basecall_qual)) in rows.iter().zip(expected) {
            assert_eq!(row.len(), 12, "every window row has twelve columns");
            assert_eq!(row.first(), Some(&CONTIG), "contig name");
            assert_eq!(row.get(3), Some(&read_id), "read id");
            assert_eq!(row.get(9), Some(&win_start), "window start");
            assert_eq!(row.get(10), Some(&win_end), "window end");
            assert_eq!(row.get(11), Some(&basecall_qual), "basecall quality");
        }
    }

    /// Full display shows sequence and qualities, and empty `SEQ` uses fallbacks.
    #[test]
    fn reads_table_displays_sequence_and_empty_seq_fallbacks() {
        let full = mapped_record("full_seq", Cigar::Match(5), b"ACGTA", &[10, 20, 30, 40, 50]);
        let empty = unmapped_record("empty_seq", b"");

        let mut output = Vec::new();
        reads_table::run(
            &mut output,
            records(vec![full, empty]),
            None,
            SeqDisplayOptions::Full {
                show_base_qual: true,
            },
            "",
        )
        .expect("both records are displayable");

        let text = String::from_utf8(output).expect("reads_table writes UTF-8");
        let mut lines = text.lines();
        assert_eq!(
            lines.next().expect("header line"),
            "read_id\talign_length\tsequence_length_template\talignment_type\tsequence\tqualities",
            "the documented header must come first"
        );
        let mut rows: Vec<&str> = lines.collect();
        rows.sort_unstable();
        assert_eq!(
            rows,
            vec![
                "empty_seq\t0\t0\tunmapped\t*\t255",
                "full_seq\t5\t5\tprimary_forward\tACGTA\t10.20.30.40.50",
            ],
            "empty SEQ and QUAL must fall back to `*` and 255 while full data is shown verbatim"
        );
    }

    /// An output handle that either captures bytes or fails after a budget.
    struct ProbeWriter {
        /// Bytes captured in both modes.
        buffer: Vec<u8>,
        /// Bytes still accepted; `None` captures output instead of failing.
        bytes_left: Option<usize>,
        /// Captured-byte threshold at which the next `flush` fails.
        fail_flush_after: Option<usize>,
    }

    impl ProbeWriter {
        /// Captures all output.
        fn capture() -> Self {
            Self {
                buffer: Vec::new(),
                bytes_left: None,
                fail_flush_after: None,
            }
        }

        /// Accepts `bytes_left` bytes, fails one write, then accepts later writes.
        fn failing(bytes_left: usize) -> Self {
            Self {
                buffer: Vec::new(),
                bytes_left: Some(bytes_left),
                fail_flush_after: None,
            }
        }

        /// Fails one flush after at least `captured_bytes` have been written.
        fn with_failing_flush_after(mut self, captured_bytes: usize) -> Self {
            self.fail_flush_after = Some(captured_bytes);
            self
        }

        /// Returns the captured bytes as UTF-8 text.
        fn captured(&self) -> String {
            String::from_utf8(self.buffer.clone()).expect("subcommands write UTF-8")
        }
    }

    impl io::Write for ProbeWriter {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            match self.bytes_left {
                None => {
                    self.buffer.extend_from_slice(buf);
                    Ok(buf.len())
                }
                Some(0) => {
                    self.bytes_left = None;
                    Err(io::Error::other("synthetic write failure"))
                }
                Some(bytes_left) => {
                    let accepted = buf.len().min(bytes_left);
                    self.buffer.extend_from_slice(
                        buf.get(..accepted)
                            .expect("accepted bytes are within the buffer"),
                    );
                    self.bytes_left = Some(
                        bytes_left
                            .checked_sub(accepted)
                            .expect("accepted bytes never exceed the budget"),
                    );
                    Ok(accepted)
                }
            }
        }

        fn flush(&mut self) -> io::Result<()> {
            if self
                .fail_flush_after
                .is_some_and(|threshold| self.buffer.len() >= threshold)
            {
                self.fail_flush_after = None;
                Err(io::Error::other("synthetic flush failure"))
            } else {
                Ok(())
            }
        }
    }

    /// Asserts that a synthetic writer failure is propagated unchanged.
    fn assert_io_error(error: Error, label: &str, marker: &str) {
        let Error::InputOutputError(source) = error else {
            unreachable!("{label}: expected InputOutputError, got {error:?}")
        };
        assert_eq!(
            source.to_string(),
            marker,
            "{label}: the originating writer failure must be preserved"
        );
    }

    /// Fails `run_once` at every output position and expects `InputOutputError`.
    ///
    /// A successful run first captures the exact output, so each line start and
    /// the final flush pin one `?` error path inside the subcommand.
    fn assert_write_failures<F>(label: &str, mut run_once: F)
    where
        F: FnMut(&mut ProbeWriter) -> Result<(), Error>,
    {
        let mut capture = ProbeWriter::capture();
        run_once(&mut capture).unwrap_or_else(|error| {
            unreachable!("{label}: expected a successful run, got {error:?}")
        });
        let successful_output = capture.captured();

        let mut budget = 0usize;
        for line in successful_output.lines() {
            let mut writer = ProbeWriter::failing(budget);
            let error = run_once(&mut writer)
                .expect_err("a failing writer must surface an error instead of panicking");
            assert_io_error(
                error,
                &format!("{label}: write failure at byte {budget}"),
                "synthetic write failure",
            );
            budget = budget
                .checked_add(line.len().checked_add(1).expect("lines are short"))
                .expect("outputs are small");
        }

        let mut immediate_writer = ProbeWriter::capture().with_failing_flush_after(0);
        let immediate_error = run_once(&mut immediate_writer)
            .expect_err("an immediate flush failure must surface instead of panicking");
        assert_io_error(
            immediate_error,
            &format!("{label}: immediate flush failure"),
            "synthetic flush failure",
        );

        let mut final_writer =
            ProbeWriter::capture().with_failing_flush_after(successful_output.len());
        let final_error = run_once(&mut final_writer)
            .expect_err("the final flush failure must surface instead of panicking");
        assert_eq!(
            final_writer.buffer,
            successful_output.as_bytes(),
            "{label}: final flush must fail only after the complete output was written"
        );
        assert_io_error(
            final_error,
            &format!("{label}: final flush failure"),
            "synthetic flush failure",
        );
    }

    /// Every subcommand reports write and flush failures instead of panicking.
    #[test]
    fn subcommands_propagate_output_write_failures() {
        let mut required_mods = InputMods::<RequiredTag>::default();
        required_mods.tag = "T".parse::<RequiredTag>().expect("T is a valid tag");
        let optional_mods = InputMods::<OptionalTag>::default();

        assert_write_failures("peek", |writer| {
            peek::run(
                writer,
                &header(),
                records(vec![mapped_record(
                    "write_fail",
                    Cigar::Match(5),
                    b"ACGTA",
                    &[30u8; 5],
                )])
                .into_iter(),
            )
        });

        assert_write_failures("read_stats", |writer| {
            read_stats::run(
                writer,
                records(vec![mapped_record(
                    "write_fail",
                    Cigar::Match(5),
                    b"ACGTA",
                    &[30u8; 5],
                )]),
            )
        });

        assert_write_failures("find_modified_reads", |writer| {
            find_modified_reads::run(
                writer,
                records(vec![modified_record(
                    "write_fail",
                    "T+T,0,0,0,0,0;",
                    &[200u8; 5],
                    b"TTTTT",
                    &[30u8; 5],
                )]),
                window_options(),
                &required_mods,
                threshold_and_mean,
                |windows: &Vec<F32Bw0and1>| !windows.is_empty(),
            )
        });

        assert_write_failures("window_reads", |writer| {
            window_reads::run(
                writer,
                records(vec![modified_record(
                    "write_fail",
                    "T+T,0,0,0,0,0;",
                    &[200u8; 5],
                    b"TTTTT",
                    &[30u8; 5],
                )]),
                window_options(),
                &optional_mods,
                |values: &[u8]| threshold_and_mean(values).map(Into::into),
            )
        });

        assert_write_failures("reads_table", |writer| {
            reads_table::run(
                writer,
                records(vec![mapped_record(
                    "write_fail",
                    Cigar::Match(5),
                    b"ACGTA",
                    &[30u8; 5],
                )]),
                None,
                SeqDisplayOptions::Full {
                    show_base_qual: true,
                },
                "",
            )
        });
    }
}
