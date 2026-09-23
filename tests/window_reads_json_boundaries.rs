#![cfg_attr(coverage_nightly, feature(coverage_attribute))]
//! Public command routing and window-output error contracts.
#[cfg(test)]
#[cfg_attr(coverage_nightly, coverage(off))]
mod tests {
    use clap::Parser as _;
    use nanalogue_core::{
        Error, InputMods, InputWindowingBuilder, OptionalTag, SeqDisplayOptions, commands,
        read_info, reads_table, window_reads,
    };
    use rust_htslib::bam::record::{Cigar, CigarString};
    use rust_htslib::bam::{Header, HeaderView, Record, header::HeaderRecord};
    use rust_htslib::errors::Error as HtslibError;
    use std::io::{BufRead as _, BufReader, Write as _};
    use std::net::TcpListener;
    use std::rc::Rc;
    use std::sync::{
        Arc, Mutex, MutexGuard,
        atomic::{AtomicBool, Ordering},
    };
    use std::thread;
    use std::time::Duration;

    const WINDOW_HEADER: &str = concat!(
        "#contig\tref_win_start\tref_win_end\tread_id\twin_val\tstrand\t",
        "base\tmod_strand\tmod_type\twin_start\twin_end\tbasecall_qual\n"
    );
    static INDEX_LOCK: Mutex<()> = Mutex::new(());
    struct BamServer {
        stop: Arc<AtomicBool>,
        thread: Option<thread::JoinHandle<()>>,
        url: String,
        _guard: MutexGuard<'static, ()>,
    }
    impl BamServer {
        fn start(serve_index: bool) -> Self {
            let guard = INDEX_LOCK.lock().expect("index fixture lock is available");
            let listener = TcpListener::bind("127.0.0.1:0").expect("test server should bind");
            listener
                .set_nonblocking(true)
                .expect("listener should become nonblocking");
            let address = listener.local_addr().expect("listener has an address");
            let bam = std::fs::read("./examples/example_1.bam").expect("example BAM is readable");
            let bai = std::fs::read("./examples/example_1.bam.bai")
                .expect("example BAM index is readable");
            let stop = Arc::new(AtomicBool::new(false));
            let thread_stop = Arc::clone(&stop);
            let thread = thread::spawn(move || {
                let index: &[u8] = if serve_index { &bai } else { b"invalid index" };
                while !thread_stop.load(Ordering::Relaxed) {
                    match listener.accept() {
                        Ok((mut stream, _)) => {
                            stream
                                .set_read_timeout(Some(Duration::from_secs(2)))
                                .expect("read timeout can be set");
                            stream
                                .set_write_timeout(Some(Duration::from_secs(2)))
                                .expect("write timeout can be set");
                            let mut request_text = String::new();
                            let mut reader = BufReader::new(&mut stream);
                            loop {
                                let mut line = String::new();
                                let size = reader.read_line(&mut line).expect("request line reads");
                                if size == 0 || line == "\r\n" {
                                    break;
                                }
                                request_text.push_str(&line);
                            }
                            drop(reader);
                            let first_line = request_text.lines().next().unwrap_or_default();
                            let requested_content = if first_line.contains(" /input.bam ") {
                                Some(bam.as_slice())
                            } else if first_line.contains(" /input.bam.bai ")
                                || first_line.contains(" /input.bai ")
                            {
                                Some(index)
                            } else {
                                None
                            };
                            if let Some(content) = requested_content {
                                let range_start = request_text.lines().find_map(|line| {
                                    line.strip_prefix("Range: bytes=")?
                                        .split('-')
                                        .next()?
                                        .parse::<usize>()
                                        .ok()
                                });
                                let (status, body, content_range) = range_start.map_or_else(
                                    || ("200 OK", content, String::new()),
                                    |start| {
                                        let body = content.get(start..).unwrap_or_default();
                                        let end = content.len().saturating_sub(1);
                                        (
                                            "206 Partial Content",
                                            body,
                                            format!(
                                                "Content-Range: bytes {start}-{end}/{}\r\n",
                                                content.len()
                                            ),
                                        )
                                    },
                                );
                                let header = format!(
                                    "HTTP/1.1 {status}\r\nContent-Length: {}\r\n{content_range}Accept-Ranges: bytes\r\nConnection: close\r\n\r\n",
                                    body.len()
                                );
                                stream.write_all(header.as_bytes()).expect("header writes");
                                if first_line.starts_with("GET ") {
                                    stream.write_all(body).expect("BAM body writes");
                                }
                            } else {
                                stream
                                    .write_all(
                                        b"HTTP/1.1 404 Not Found\r\nContent-Length: 0\r\nConnection: close\r\n\r\n",
                                    )
                                    .expect("404 writes");
                            }
                        }
                        Err(error) if error.kind() == std::io::ErrorKind::WouldBlock => {
                            thread::sleep(Duration::from_millis(1));
                        }
                        Err(error) => unreachable!("test server accept failed: {error}"),
                    }
                }
            });
            Self {
                stop,
                thread: Some(thread),
                url: format!("http://{address}/input.bam"),
                _guard: guard,
            }
        }
    }

    impl Drop for BamServer {
        fn drop(&mut self) {
            self.stop.store(true, Ordering::Relaxed);
            self.thread
                .take()
                .expect("server thread exists")
                .join()
                .expect("server thread should stop");
            std::fs::remove_file("input.bam.bai").unwrap_or_else(|error| {
                assert_eq!(error.kind(), std::io::ErrorKind::NotFound);
            });
        }
    }

    #[test]
    fn subcommands_preserve_unsupported_record_errors() {
        let header = HeaderView::from_header(
            Header::new().push_record(
                HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", "ctg")
                    .push_tag(b"LN", 100),
            ),
        );
        let mut record = Record::new();
        record.set(
            b"paired_route",
            Some(&CigarString(vec![Cigar::Match(4)])),
            b"ACGT",
            &[30; 4],
        );
        record.set_header(Arc::new(header));
        record.set_flags(1);
        record.set_tid(0);
        record.set_pos(7);
        let record_rc = Rc::new(record);
        let one = || vec![Ok(Rc::clone(&record_rc))];
        let assert_error = |route: &str, error: Error| {
            assert!(
                matches!(error, Error::NotImplemented(message)
                    if message.contains("paired_route")),
                "{route} must preserve the unsupported flag and read ID"
            );
        };
        let mut info_output = Vec::new();
        assert_error(
            "read-info",
            read_info::run(
                &mut info_output,
                one(),
                InputMods::<OptionalTag>::default(),
                None,
            )
            .expect_err("paired input must fail read-info"),
        );
        assert!(info_output.is_empty());
        let mut table_output = Vec::new();
        assert_error(
            "read-table",
            reads_table::run(&mut table_output, one(), None, SeqDisplayOptions::No, "")
                .expect_err("paired input must fail read-table"),
        );
        assert!(table_output.is_empty());

        let windowing = InputWindowingBuilder::default()
            .win(2)
            .step(1)
            .build()
            .expect("valid test windowing");
        for route in ["window text", "window JSON"] {
            let mut output = Vec::new();
            let error = if route == "window text" {
                window_reads::run(
                    &mut output,
                    one(),
                    windowing,
                    &InputMods::<OptionalTag>::default(),
                    |_| unreachable!("record parsing fails before windowing"),
                )
            } else {
                window_reads::run_json(
                    &mut output,
                    one(),
                    windowing,
                    &InputMods::<OptionalTag>::default(),
                    |_| unreachable!("record parsing fails before windowing"),
                )
            }
            .expect_err("paired input must fail windowing");
            assert_error(route, error);
            assert_eq!(
                output,
                if route == "window text" {
                    WINDOW_HEADER.as_bytes()
                } else {
                    b"["
                },
                "{route} must stop at its format-specific prefix"
            );
        }
    }

    #[test]
    #[expect(clippy::too_many_lines, reason = "routing matrix belongs together")]
    fn command_routes_propagate_read_id_list_open_errors() {
        const BAM: &str = "./examples/example_1.bam";
        const MISSING: &str = "./missing_window_routing_read_ids.txt";
        let cases = [
            vec!["read-table-show-mods", "--read-id-list", MISSING, BAM],
            vec!["read-table-hide-mods", "--read-id-list", MISSING, BAM],
            vec!["read-stats", "--read-id-list", MISSING, BAM],
            vec!["read-info", "--read-id-list", MISSING, BAM],
            vec![
                "find-modified-reads",
                "all-dens-between",
                "--win",
                "2",
                "--step",
                "1",
                "--tag",
                "T",
                "--dens-limits",
                "0.2,0.8",
                "--read-id-list",
                MISSING,
                BAM,
            ],
            vec![
                "find-modified-reads",
                "any-dens-above",
                "--win",
                "2",
                "--step",
                "1",
                "--tag",
                "T",
                "--high",
                "0.8",
                "--read-id-list",
                MISSING,
                BAM,
            ],
            vec![
                "find-modified-reads",
                "any-dens-below",
                "--win",
                "2",
                "--step",
                "1",
                "--tag",
                "T",
                "--low",
                "0.2",
                "--read-id-list",
                MISSING,
                BAM,
            ],
            vec![
                "find-modified-reads",
                "any-dens-below-and-any-dens-above",
                "--win",
                "2",
                "--step",
                "1",
                "--tag",
                "T",
                "--low",
                "0.2",
                "--high",
                "0.8",
                "--read-id-list",
                MISSING,
                BAM,
            ],
            vec![
                "find-modified-reads",
                "dens-range-above",
                "--win",
                "2",
                "--step",
                "1",
                "--tag",
                "T",
                "--min-range",
                "0.6",
                "--read-id-list",
                MISSING,
                BAM,
            ],
            vec![
                "find-modified-reads",
                "any-abs-grad-above",
                "--win",
                "2",
                "--step",
                "1",
                "--tag",
                "T",
                "--min-grad",
                "0.3",
                "--read-id-list",
                MISSING,
                BAM,
            ],
            vec![
                "window-dens",
                "--win",
                "2",
                "--step",
                "1",
                "--read-id-list",
                MISSING,
                BAM,
            ],
            vec![
                "window-grad",
                "--win",
                "2",
                "--step",
                "1",
                "--read-id-list",
                MISSING,
                BAM,
            ],
        ];

        for args in cases {
            let route = args.first().copied().unwrap_or("unknown");
            let cli =
                commands::Cli::parse_from(std::iter::once("nanalogue").chain(args.iter().copied()));
            let mut output = Vec::new();
            let error = commands::run(cli, &mut output)
                .expect_err("a missing read-ID list must abort routing");

            assert!(
                matches!(&error, Error::InputOutputError(source)
                    if source.kind() == std::io::ErrorKind::NotFound),
                "{route} must preserve the list-open NotFound error, got {error:?}"
            );
            assert!(output.is_empty(), "{route} must fail before writing output");
        }
    }

    #[test]
    fn window_commands_propagate_record_decode_errors() {
        for subcommand in ["window-dens", "window-grad"] {
            let cli = commands::Cli::parse_from([
                "nanalogue",
                subcommand,
                "--win",
                "2",
                "--step",
                "1",
                "./examples/example_4_invalid_basequal_len.sam",
            ]);
            let mut output = Vec::new();
            let error = commands::run(cli, &mut output)
                .expect_err("the malformed SAM record must fail decoding");

            assert!(
                matches!(&error, Error::RustHtslibError(source)
                    if matches!(source.as_ref(), HtslibError::BamTruncatedRecord)),
                "{subcommand} must preserve the htslib truncated-record error, got {error:?}"
            );
            assert_eq!(
                String::from_utf8(output).expect("window header is UTF-8"),
                WINDOW_HEADER,
                "{subcommand} must fail before writing a data row"
            );
        }
    }

    #[test]
    fn window_commands_distinguish_zero_length_records() {
        for subcommand in ["window-dens", "window-grad"] {
            let cli = commands::Cli::parse_from([
                "nanalogue",
                subcommand,
                "--win",
                "2",
                "--step",
                "1",
                "--read-id",
                "read1",
                "./examples/example_2_zero_len.sam",
            ]);
            let mut output = Vec::new();
            let error = commands::run(cli, &mut output)
                .expect_err("zero-length input cannot produce windowed data");

            assert!(
                matches!(error, Error::InvalidState(message)
                if message == concat!(
                    "No records found as input for analysis. This could mean ",
                    "the input genomics file (SAM/BAM/CRAM) has no records or that filtering\n",
                    "removed records or some other possibility."
                )),
                "{subcommand} must distinguish filtered input from a decode failure"
            );
            assert_eq!(
                String::from_utf8(output).expect("window header is UTF-8"),
                WINDOW_HEADER,
                "{subcommand} must emit only its header before validation fails"
            );
        }
    }

    #[test]
    fn regional_url_falls_back_to_unindexed_windowing() {
        let server = BamServer::start(false);
        let cli = commands::Cli::parse_from([
            "nanalogue",
            "window-dens",
            "--region",
            "dummyI:1-22",
            "--win",
            "2",
            "--step",
            "1",
            &server.url,
        ]);
        let mut output = Vec::new();

        commands::run(cli, &mut output).expect("missing remote index should fall back");

        let text = String::from_utf8(output).expect("window output is UTF-8");
        assert!(text.starts_with("#contig\tref_win_start"));
        assert!(
            text.lines()
                .skip(1)
                .all(|line| line.starts_with("dummyI\t")),
            "fallback must retain the requested contig filter: {text}"
        );
        assert_eq!(
            text.lines().count(),
            4,
            "the selected read has exactly three two-call windows"
        );
    }

    #[test]
    fn regional_url_uses_remote_index_when_available() {
        let server = BamServer::start(true);
        let cli = commands::Cli::parse_from([
            "nanalogue",
            "read-stats",
            "--region",
            "dummyI:1-22",
            &server.url,
        ]);
        let mut output = Vec::new();

        commands::run(cli, &mut output).expect("the remote index should support a region fetch");

        let text = String::from_utf8(output).expect("statistics are UTF-8");
        let values: std::collections::HashMap<_, _> = text
            .lines()
            .skip(1)
            .map(|line| line.split_once('\t').expect("statistics have two columns"))
            .collect();
        assert_eq!(values.len(), 15, "every statistics field is present");
        assert_eq!(values.get("n_primary_alignments"), Some(&"1"));
        assert_eq!(values.get("n_secondary_alignments"), Some(&"0"));
        assert_eq!(values.get("align_len_mean"), Some(&"8"));
        assert_eq!(values.get("seq_len_n50"), Some(&"8"));
    }
}
