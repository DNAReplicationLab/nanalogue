//! Utility functions for file I/O operations with BAM and FASTA files.

use crate::{Error, GetDNARestrictive};
use rust_htslib::{bam, htslib};
use std::ffi::CString;
use std::fs::File;
use std::io::Write as _;
use std::num::NonZeroU32;
use std::path::{Path, PathBuf};
use url::Url;

#[cfg(unix)]
use std::os::unix::fs::MetadataExt as _;

/// Opens BAM file.
///
/// # Errors
///
/// Returns an error if the BAM file cannot be opened or read.
///
/// ```
/// use nanalogue_core::{Error, file_utils::nanalogue_bam_reader};
/// use rust_htslib::bam::Read;
/// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
/// // the above file should contain four reads, so we are checking
/// // if we load four records.
/// let mut count = 0;
/// for r in reader.records() {
///     count = count + 1;
/// }
/// assert_eq!(count, 4);
/// # Ok::<(), Error>(())
/// ```
pub fn nanalogue_bam_reader<J>(bam_path: &J) -> Result<bam::Reader, Error>
where
    J: AsRef<Path> + ?Sized,
{
    Ok(bam::Reader::from_path(bam_path)?)
}

/// Receives BAM data from stdin.
///
/// We are not writing tests for this function.
///
/// # Errors
///
/// Returns an error if the BAM data cannot be read.
pub fn nanalogue_bam_reader_from_stdin() -> Result<bam::Reader, Error> {
    Ok(bam::Reader::from_stdin()?)
}

/// Receives BAM data from a url.
///
/// We are not writing tests for this function as we have to rely on an external URL,
/// which may be unreliable.
///
/// # SSL setup
///
/// For HTTPS URLs, libcurl needs to know where the system CA bundle lives.
/// Library callers have two options:
///
/// - Set `SSL_CERT_FILE` (and optionally `CURL_CA_BUNDLE` / `SSL_CERT_DIR`) in
///   the process environment before running the program. No call into this
///   crate is needed in that case.
/// - Otherwise, invoke [`crate::init_ssl_certificates`] once at program
///   startup, before any threads are spawned, to auto-detect and populate
///   those variables. The binaries shipped with this crate already do this.
///
/// # Errors
///
/// Returns an error if the BAM data cannot be read.
pub fn nanalogue_bam_reader_from_url(url: &Url) -> Result<bam::Reader, Error> {
    Ok(bam::Reader::from_url(url)?)
}

/// Opens indexed BAM file, fetching according to input instructions.
///
/// `FetchDefinition` is a struct used by `rust_htslib` to retrieve
/// a region, a contig, all reads, unmapped reads etc. Have a look at
/// their documentation for the different variants.
///
/// # Errors
///
/// Returns an error if the BAM file cannot be opened or read or fetching does not work
///
/// # Examples
///
/// Retrieve all reads
///
/// ```
/// use nanalogue_core::{Error, file_utils::nanalogue_indexed_bam_reader};
/// use rust_htslib::bam::{Read, FetchDefinition};
/// let mut reader = nanalogue_indexed_bam_reader(&"examples/example_1.bam", FetchDefinition::All)?;
/// // the above file should contain four reads, so we are checking
/// // if we load four records.
/// assert_eq!(reader.records().count(), 4);
/// # Ok::<(), Error>(())
/// ```
///
/// Retrieve all reads from a contig.
///
/// ```
/// # use nanalogue_core::{Error, file_utils::nanalogue_indexed_bam_reader};
/// # use rust_htslib::bam::{Read, FetchDefinition};
/// let mut reader = nanalogue_indexed_bam_reader(&"examples/example_1.bam",
///     FetchDefinition::String(b"dummyI"))?;
/// // the above file should contain only one read passing through this contig.
/// assert_eq!(reader.records().count(), 1);
/// # Ok::<(), Error>(())
/// ```
///
/// Retrieve all reads overlapping with a given region.
///
/// ```
/// # use nanalogue_core::{Error, file_utils::nanalogue_indexed_bam_reader};
/// # use rust_htslib::bam::{Read, FetchDefinition};
/// // this file has a read passing through this contig so we're testing we
/// // get 1 read when we use coordinates that the read passes through and
/// // 0 when we don't.
/// let mut reader = nanalogue_indexed_bam_reader(&"examples/example_1.bam",
///     FetchDefinition::RegionString(b"dummyIII", 10, 20))?;
/// assert_eq!(reader.records().count(), 0);
/// let mut reader = nanalogue_indexed_bam_reader(&"examples/example_1.bam",
///     FetchDefinition::RegionString(b"dummyIII", 20, 30))?;
/// assert_eq!(reader.records().count(), 1);
/// // this contig doesn't have 30000bp, so we are testing if using
/// // a coordinate beyond the bounds of the contig is o.k.
/// let mut reader = nanalogue_indexed_bam_reader(&"examples/example_1.bam",
///     FetchDefinition::RegionString(b"dummyIII", 20, 30000))?;
/// assert_eq!(reader.records().count(), 1);
/// # Ok::<(), Error>(())
/// ```
///
/// Error when accessing a contig with numeric index larger than number of contigs.
///
/// ```should_panic
/// # use nanalogue_core::{Error, file_utils::nanalogue_indexed_bam_reader};
/// # use rust_htslib::bam::{Read, FetchDefinition};
/// let mut reader = nanalogue_indexed_bam_reader(&"examples/example_1.bam",
///     FetchDefinition::CompleteTid(10))?;
/// // the above panics as this file has much fewer than 10 contigs.
/// # Ok::<(), Error>(())
/// ```
///
/// Check we get an error when the index is missing.
///
/// ```
/// # use nanalogue_core::{Error, file_utils::nanalogue_indexed_bam_reader};
/// # use rust_htslib::bam::{Read, FetchDefinition};
/// use rust_htslib::errors::Error as OtherError;
/// use std::fs;
/// assert!(!(fs::exists("examples/example_1_copy_no_index.bam.bai").unwrap()));
/// // above tells us index does not exist.
/// let err = nanalogue_indexed_bam_reader(&"examples/example_1_copy_no_index.bam",
///     FetchDefinition::All).unwrap_err();
/// // we should get a missing index error, which is same as an invalid index in `rust_htslib`.
/// assert!(matches!(err, Error::RustHtslibError(err)
///     if matches!(err.as_ref(), OtherError::BamInvalidIndex{..})));
/// ```
///
/// Check we get an error when the index is malformed.
///
/// ```
/// # use nanalogue_core::{Error, file_utils::nanalogue_indexed_bam_reader};
/// # use rust_htslib::bam::{Read, FetchDefinition};
/// use rust_htslib::errors::Error as OtherError;
/// use std::fs;
/// assert!(fs::exists("examples/example_1_copy_invalid_index.bam.bai").unwrap());
/// // above tells us index exists.
/// let err = nanalogue_indexed_bam_reader(&"examples/example_1_copy_invalid_index.bam",
///     FetchDefinition::All).unwrap_err();
/// // we should get an invalid index error, as the index exists here but is just a blank file.
/// assert!(matches!(err, Error::RustHtslibError(err)
///     if matches!(err.as_ref(), OtherError::BamInvalidIndex{..})));
/// ```
pub fn nanalogue_indexed_bam_reader<J>(
    bam_path: &J,
    fetch_definition: bam::FetchDefinition,
) -> Result<bam::IndexedReader, Error>
where
    J: AsRef<Path> + ?Sized,
{
    let mut bam_reader = bam::IndexedReader::from_path(bam_path)?;
    bam_reader.fetch(fetch_definition)?;
    Ok(bam_reader)
}

/// Receives indexed BAM data from a url.
///
/// We are not writing tests for this function as we have to rely on an external URL,
/// which may be unreliable.
///
/// # SSL setup
///
/// For HTTPS URLs, libcurl needs to know where the system CA bundle lives.
/// Library callers have two options:
///
/// - Set `SSL_CERT_FILE` (and optionally `CURL_CA_BUNDLE` / `SSL_CERT_DIR`) in
///   the process environment before running the program. No call into this
///   crate is needed in that case.
/// - Otherwise, invoke [`crate::init_ssl_certificates`] once at program
///   startup, before any threads are spawned, to auto-detect and populate
///   those variables. The binaries shipped with this crate already do this.
///
/// # Errors
///
/// Returns an error if the BAM data cannot be read.
pub fn nanalogue_indexed_bam_reader_from_url(
    url: &Url,
    fetch_definition: bam::FetchDefinition,
) -> Result<bam::IndexedReader, Error> {
    let mut bam_reader = bam::IndexedReader::from_url(url)?;
    bam_reader.fetch(fetch_definition)?;
    Ok(bam_reader)
}
/// Writes contigs to a FASTA file
///
/// # Errors
///
/// Returns an error if the file cannot be created or written to.
///
/// # Examples
///
/// ```
/// use nanalogue_core::{DNARestrictive, Error, file_utils::write_fasta, uuid};
/// use std::fs;
/// use std::str::FromStr;
///
/// let contigs = vec![
///     ("seq1".to_string(), DNARestrictive::from_str("ACGTACGT").expect("no error")),
///     ("seq2".to_string(), DNARestrictive::from_str("TGCATGCA").expect("no error")),
/// ];
///
/// let temp_path = std::env::temp_dir().join(format!("{}.fa", uuid::v4_random()));
/// write_fasta(contigs, &temp_path)?;
///
/// let content = fs::read_to_string(&temp_path)?;
/// assert!(content.contains(">seq1"));
/// assert!(content.contains("ACGTACGT"));
/// assert!(content.contains(">seq2"));
/// assert!(content.contains("TGCATGCA"));
///
/// fs::remove_file(&temp_path)?;
/// # Ok::<(), Error>(())
/// ```
pub fn write_fasta<I, J, K>(sequences: I, output_path: &J) -> Result<(), Error>
where
    I: IntoIterator<Item = (String, K)>,
    J: AsRef<Path> + ?Sized,
    K: GetDNARestrictive,
{
    let mut file = File::create(output_path)?;
    for seq in sequences {
        writeln!(file, ">{}", seq.0)?;
        file.write_all(seq.1.get_dna_restrictive().get())?;
        writeln!(file)?;
    }
    Ok(())
}

/// Build the alignment header shared by de novo BAM and CRAM output.
fn denovo_alignment_header<J, K, L>(contigs: J, read_groups: K, comments: L) -> bam::Header
where
    J: IntoIterator<Item = (String, usize)>,
    K: IntoIterator<Item = String>,
    L: IntoIterator<Item = String>,
{
    let mut header = bam::Header::new();
    for (name, length) in contigs {
        let _: &mut _ = header.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", name)
                .push_tag(b"LN", length),
        );
    }
    for read_group in read_groups {
        let _: &mut _ = header.push_record(
            bam::header::HeaderRecord::new(b"RG")
                .push_tag(b"ID", read_group)
                .push_tag(b"PL", "ONT")
                .push_tag(b"LB", "blank")
                .push_tag(b"SM", "blank")
                .push_tag(b"PU", "blank"),
        );
    }
    for comment in comments {
        let _: &mut _ = header.push_comment(comment.as_bytes());
    }
    header
}

/// Validate nanalogue's deterministic coordinate-and-strand order while writing records.
fn write_sorted_reads<I, F>(reads: I, mut write: F) -> Result<(), Error>
where
    I: IntoIterator<Item = bam::Record>,
    F: FnMut(&bam::Record) -> Result<(), Error>,
{
    let mut previous_key = (false, -1, -1, false);
    for read in reads {
        let current_key = (
            read.is_unmapped(),
            read.tid(),
            read.pos(),
            read.is_reverse(),
        );
        if previous_key > current_key {
            return Err(Error::InvalidSorting(
                "reads input to the de novo alignment writer have not been sorted properly"
                    .to_string(),
            ));
        }
        previous_key = current_key;
        write(&read)?;
    }
    Ok(())
}

/// Convert a path to the null-terminated representation expected by `HTSlib`.
fn path_to_c_string(path: &Path) -> Result<CString, Error> {
    CString::new(path.as_os_str().as_encoded_bytes())
        .map_err(|error| Error::InvalidState(error.to_string()))
}

/// Return the conventional sidecar path beside an output file.
pub(crate) fn alignment_sidecar_path(output_path: &Path, suffix: &str) -> PathBuf {
    let mut path = output_path.as_os_str().to_os_string();
    path.push(suffix);
    PathBuf::from(path)
}

/// Resolve a path to the filesystem location used for output-collision checks.
///
/// Canonicalizing the parent separately supports outputs that do not exist yet while still
/// resolving relative paths, `..` components, and symlinked directories.
fn output_path_identity(path: &Path) -> Result<PathBuf, Error> {
    assert!(
        !path.as_os_str().as_encoded_bytes().is_empty(),
        "output path must not be empty"
    );
    assert!(
        path.as_os_str().as_encoded_bytes().len() <= 10_000,
        "output path must not exceed 10,000 bytes"
    );
    match std::fs::canonicalize(path) {
        Ok(identity) => Ok(identity),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            match std::fs::symlink_metadata(path) {
                Ok(metadata) if metadata.file_type().is_symlink() => {
                    return Err(Error::InvalidState(
                        "output paths must not be symbolic links to files that do not exist".into(),
                    ));
                }
                Ok(_) => return Err(error.into()),
                Err(metadata_error) if metadata_error.kind() == std::io::ErrorKind::NotFound => {}
                Err(metadata_error) => return Err(metadata_error.into()),
            }
            let file_name = path.file_name().ok_or(error)?;
            let parent = path
                .parent()
                .filter(|parent| !parent.as_os_str().is_empty());
            Ok(std::fs::canonicalize(parent.unwrap_or_else(|| Path::new(".")))?.join(file_name))
        }
        Err(error) => Err(error.into()),
    }
}

/// Return the filesystem identity of an existing path without opening it.
#[cfg(unix)]
fn existing_path_identity(path: &Path) -> Result<Option<(u64, u64)>, Error> {
    assert!(
        !path.as_os_str().as_encoded_bytes().is_empty(),
        "path must not be empty"
    );
    assert!(
        path.as_os_str().as_encoded_bytes().len() <= 10_000,
        "path must not exceed 10,000 bytes"
    );

    let metadata = match std::fs::metadata(path) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(error) => return Err(error.into()),
    };

    Ok(Some((metadata.dev(), metadata.ino())))
}

/// Return whether paths identify pairwise-distinct filesystem locations.
pub(crate) fn output_paths_are_distinct(paths: &[&Path]) -> Result<bool, Error> {
    assert!(!paths.is_empty(), "paths must not be empty");
    assert!(
        paths.len() <= 20,
        "paths must not contain more than 20 items"
    );
    for path in paths {
        assert!(
            !path.as_os_str().as_encoded_bytes().is_empty(),
            "output paths must not contain an empty path"
        );
        assert!(
            path.as_os_str().as_encoded_bytes().len() <= 10_000,
            "output paths must not exceed 10,000 bytes"
        );
    }

    let mut lexical_paths = Vec::with_capacity(paths.len());
    for path in paths {
        if lexical_paths.contains(path) {
            return Ok(false);
        }
        lexical_paths.push(*path);
    }

    #[cfg(unix)]
    {
        let mut existing_identities = Vec::with_capacity(paths.len());
        for path in paths {
            if let Some(identity) = existing_path_identity(path)? {
                if existing_identities.contains(&identity) {
                    return Ok(false);
                }
                existing_identities.push(identity);
            }
        }
    }

    let mut identities = Vec::with_capacity(paths.len());
    for path in paths {
        let identity = output_path_identity(path)?;
        if identities.contains(&identity) {
            return Ok(false);
        }
        identities.push(identity);
    }
    Ok(true)
}

/// Return the conventional CRAI path beside a CRAM output.
fn crai_path(output_path: &Path) -> PathBuf {
    alignment_sidecar_path(output_path, ".crai")
}

/// A CRAM stream whose index and final compression are explicitly completed.
#[derive(Debug)]
struct CramWriter {
    /// Raw `HTSlib` output handle.
    file: *mut htslib::htsFile,
    /// Index path retained for the lifetime required by `HTSlib`.
    index_path: CString,
}

impl CramWriter {
    /// Open a CRAM 3.1 stream using an external FASTA reference.
    fn open(output_path: &Path, reference_path: &Path, threads: NonZeroU32) -> Result<Self, Error> {
        let output = path_to_c_string(output_path)?;
        let index = path_to_c_string(&crai_path(output_path))?;
        let reference = path_to_c_string(reference_path)?;
        let mode = CString::new("wc").map_err(|error| Error::InvalidState(error.to_string()))?;
        let version =
            CString::new("3.1").map_err(|error| Error::InvalidState(error.to_string()))?;

        // SAFETY: `output` and `mode` are valid NUL-terminated strings that outlive the call.
        let file = unsafe { htslib::hts_open(output.as_ptr(), mode.as_ptr()) };
        if file.is_null() {
            return Err(Error::WriteOutput("failed to open CRAM output".into()));
        }
        let writer = Self {
            file,
            index_path: index,
        };

        // SAFETY: `writer.file` is a live handle from `hts_open`, and `version` lives across
        // the FFI call.
        if unsafe {
            htslib::hts_set_opt(
                writer.file,
                htslib::hts_fmt_option_CRAM_OPT_VERSION,
                version.as_ptr(),
            )
        } != 0
        {
            return Err(Error::WriteOutput(
                "failed to select CRAM version 3.1".into(),
            ));
        }
        #[expect(
            clippy::undocumented_unsafe_blocks,
            reason = "`writer.file` remains valid for both of these option-setting calls"
        )]
        if unsafe { htslib::hts_set_opt(writer.file, htslib::hts_fmt_option_CRAM_OPT_EMBED_REF, 0) }
            != 0
            || unsafe {
                htslib::hts_set_opt(writer.file, htslib::hts_fmt_option_CRAM_OPT_NO_REF, 0)
            } != 0
        {
            return Err(Error::WriteOutput(
                "failed to require external-reference CRAM compression".into(),
            ));
        }
        // SAFETY: `writer.file` is valid and `reference` is a valid NUL-terminated string that
        // outlives the call.
        if unsafe { htslib::hts_set_fai_filename(writer.file, reference.as_ptr()) } != 0 {
            return Err(Error::WriteOutput("failed to load FASTA reference".into()));
        }
        let thread_count = i32::try_from(threads.get())?;
        // SAFETY: `writer.file` is valid and `thread_count` is a plain integer argument.
        if unsafe { htslib::hts_set_threads(writer.file, thread_count) } != 0 {
            return Err(Error::WriteOutput(
                "failed to configure CRAM compression threads".into(),
            ));
        }
        Ok(writer)
    }

    /// Write the header and begin constructing the CRAI.
    fn write_header_and_start_index(&mut self, header: &mut bam::HeaderView) -> Result<(), Error> {
        // SAFETY: `self.file` is a live HTSlib handle opened for CRAM output and `header`
        // points to a valid header owned by Rust for the duration of the call.
        if unsafe { htslib::sam_hdr_write(self.file, header.inner_ptr()) } != 0 {
            return Err(Error::WriteOutput("failed to write CRAM header".into()));
        }
        // SAFETY: `self.file` is valid, `header` remains alive across the call, and
        // `self.index_path` is retained in `self` precisely so its pointer stays valid.
        if unsafe {
            htslib::sam_idx_init(
                self.file,
                header.inner_ptr_mut(),
                0,
                self.index_path.as_ptr(),
            )
        } != 0
        {
            return Err(Error::WriteOutput("failed to initialize CRAI index".into()));
        }
        Ok(())
    }

    /// Write one record to the CRAM stream.
    fn write(&mut self, header: &bam::HeaderView, record: &bam::Record) -> Result<(), Error> {
        // SAFETY: `self.file` is valid, `header` matches the output stream, and `record`
        // points to a live BAM record for the duration of the call.
        if unsafe { htslib::sam_write1(self.file, header.inner_ptr(), record.inner()) } < 0 {
            return Err(Error::WriteOutput("failed to write CRAM record".into()));
        }
        Ok(())
    }

    /// Save the CRAI and close the CRAM, reporting delayed write failures.
    fn finish(mut self) -> Result<(), Error> {
        // SAFETY: `self.file` is still a live handle whose index was initialized earlier.
        if unsafe { htslib::sam_idx_save(self.file) } != 0 {
            return Err(Error::WriteOutput("failed to save CRAI index".into()));
        }
        // SAFETY: `self.file` is a live HTSlib handle and is closed exactly once here.
        let result = unsafe { htslib::hts_close(self.file) };
        self.file = std::ptr::null_mut();
        if result != 0 {
            return Err(Error::WriteOutput("failed to finish CRAM output".into()));
        }
        Ok(())
    }
}

impl Drop for CramWriter {
    fn drop(&mut self) {
        if !self.file.is_null() {
            // SAFETY: `self.file` is a live HTSlib handle and Drop only reaches this branch if
            // ownership of the handle was not already consumed by `finish`.
            let _: i32 = unsafe { htslib::hts_close(self.file) };
        }
    }
}

/// Writes a new BAM file with reads. Input reads have to be sorted.
///
/// Although this function can be used for other tasks like subsetting
/// BAM files, this is not advised as the history of the BAM file
/// stored in its header would be lost as we are generating a new header here.
/// We don't do many checks here like checking if the number of contigs matches
/// the number of tids in the BAM records, whether BAM records contain valid
/// sequences etc.
///
/// # Errors
///
/// Returns an error if the BAM file cannot be created or written to,
/// if index creation fails, or if input reads are not sorted.
///
/// # Examples
///
/// ```
/// use nanalogue_core::{
///     Error,
///     file_utils::{nanalogue_bam_reader, write_bam_denovo},
///     uuid,
/// };
/// use rust_htslib::bam;
/// use rust_htslib::bam::Read;
///
/// let contigs = vec![("chr1".to_string(), 1000)];
/// let read_groups = vec!["rg1".to_string()];
/// let comments = vec!["test comment".to_string()];
/// let reads: Vec<bam::Record> = vec![];
///
/// let temp_path = std::env::temp_dir().join(format!("{}.bam", uuid::v4_random()));
/// write_bam_denovo(reads, contigs, read_groups, comments, &temp_path)?;
///
/// // Verify the file was created and can be read
/// let reader = nanalogue_bam_reader(temp_path.to_str().unwrap())?;
/// assert_eq!(reader.header().target_count(), 1);
///
/// std::fs::remove_file(&temp_path)?;
/// std::fs::remove_file(format!("{}.bai", temp_path.display()))?;
/// # Ok::<(), Error>(())
/// ```
pub fn write_bam_denovo<I, J, K, L, M>(
    reads: I,
    contigs: J,
    read_groups: K,
    comments: L,
    output_path: &M,
) -> Result<(), Error>
where
    I: IntoIterator<Item = bam::Record>,
    J: IntoIterator<Item = (String, usize)>,
    K: IntoIterator<Item = String>,
    L: IntoIterator<Item = String>,
    M: AsRef<Path> + ?Sized,
{
    let header = denovo_alignment_header(contigs, read_groups, comments);

    // Write BAM file ensuring reads are already sorted
    let mut writer = bam::Writer::from_path(output_path, &header, bam::Format::Bam)?;
    write_sorted_reads(reads, |read| writer.write(read).map_err(Error::from))?;
    drop(writer); // Close BAM file before creating index

    bam::index::build(output_path, None, bam::index::Type::Bai, 2)?;

    Ok(())
}

/// Write coordinate-and-strand sorted records directly as CRAM 3.1 and create a CRAI index.
///
/// The FASTA reference must describe the same contigs used to construct the alignment header.
/// MD and NM tags are not generated or explicitly stored; CRAM readers may reconstruct them.
///
/// # Errors
///
/// Returns an error if the CRAM or CRAI cannot be written, the reference cannot be loaded, or
/// the input records do not follow nanalogue's deterministic coordinate-and-strand order.
pub fn write_cram_denovo<I, J, K, L, M, N>(
    reads: I,
    contigs: J,
    read_groups: K,
    comments: L,
    output_path: &M,
    reference_path: &N,
    threads: NonZeroU32,
) -> Result<(), Error>
where
    I: IntoIterator<Item = bam::Record>,
    J: IntoIterator<Item = (String, usize)>,
    K: IntoIterator<Item = String>,
    L: IntoIterator<Item = String>,
    M: AsRef<Path> + ?Sized,
    N: AsRef<Path> + ?Sized,
{
    let output = output_path.as_ref();
    let reference = reference_path.as_ref();
    let crai_path = crai_path(output);
    let fai_path = alignment_sidecar_path(reference, ".fai");
    if !output_paths_are_distinct(&[output, &crai_path, reference, &fai_path])? {
        return Err(Error::InvalidState(
            "CRAM, CRAI, and FASTA outputs must use different paths".into(),
        ));
    }
    let header_template = denovo_alignment_header(contigs, read_groups, comments);
    let mut header = bam::HeaderView::from_header(&header_template);
    let mut writer = CramWriter::open(output, reference, threads)?;
    writer.write_header_and_start_index(&mut header)?;
    write_sorted_reads(reads, |read| writer.write(&header, read))?;
    writer.finish()
}

/// Read lines with a hard raw-byte cap.
///
/// The returned `u16` is the number of raw bytes consumed from `reader`,
/// including any trailing line terminator bytes.
///
/// Contract for successful returns:
/// - `bytes_read == 0` means EOF was reached before any bytes were consumed,
///   and `line` is empty.
/// - `bytes_read > 0 && line.is_empty()` can only occur for a blank physical
///   line consisting solely of `\n` or `\r\n`.
/// - Otherwise, `line` contains the line content with any trailing `\n` or
///   `\r\n` removed.
///
/// # Errors
///
/// Returns an error if reading from `reader` fails, if the internal reader
/// buffer is unexpectedly larger than `u16::MAX`, if a `\r` is not immediately
/// followed by `\n`, if the line contains non-ASCII control characters other
/// than `\t`, `\r`, and `\n`, or if consuming the next chunk would exceed
/// `line_cap`.
///
/// # Panics
///
/// Panics only if internal invariants about buffer sizing, trim bookkeeping,
/// or output-length accounting are violated.
#[expect(
    clippy::arithmetic_side_effects,
    clippy::indexing_slicing,
    reason = "(1) line & buffer lengths are checked for smallness, \
(2) `bytes_to_consume` never exceeds buffer length and indexing only happens when it is > 0"
)]
pub(crate) fn read_line_capped<R: std::io::BufRead>(
    reader: &mut R,
    line: &mut String,
    line_cap: u16,
) -> Result<u16, Error> {
    let mut total_bytes_read: u16 = 0;
    let mut is_dangling_slash_r = false;

    line.clear();
    let is_newline_found = loop {
        // `std::io` fills this buffer with data without stopping at new lines,
        // so we get many lines into this.
        // With `std::io::BufReader`, this buffered slice is typically modest in
        // size rather than a huge chunk of memory, so using `fill_buf` here is
        // acceptable for our defensive line-length checks i.e. we are not in
        // danger of loading a huge amount of data like 1GB into this buffer.
        let buffered = reader.fill_buf()?;
        if u16::try_from(buffered.len()).is_err() {
            return Err(Error::InvalidState(
                "internal buffer unexpectedly large".to_owned(),
            ));
        }
        let (bytes_to_consume, bytes_to_trim) = {
            let mut counter: u16 = 0;
            let mut bytes_to_trim: u16 = 0;
            assert!(
                u16::try_from(buffered.len()).is_ok(),
                "repeating this check to ensure `buffered` is bounded and `counter` won't overflow"
            );
            for k in buffered {
                counter += 1;
                if *k == b'\n' {
                    (is_dangling_slash_r, bytes_to_trim) =
                        match (is_dangling_slash_r, bytes_to_trim) {
                            (false | true, 0) => (false, 1),
                            (true, 1) => (false, 2),
                            _ => unreachable!(),
                        };
                    break;
                } else if is_dangling_slash_r {
                    return Err(Error::InvalidState(
                        "\\r must be followed by \\n".to_owned(),
                    ));
                } else if *k == b'\r' {
                    is_dangling_slash_r = true;
                    bytes_to_trim = 1;
                } else if !(*k == b'\t' || (32..127).contains(k)) {
                    return Err(Error::InvalidState(
                        "line has unusual characters!".to_owned(),
                    ));
                } else {
                    // pass
                }
            }
            assert!(
                counter >= bytes_to_trim,
                "`counter` cannot be smaller than `bytes_to_trim`!"
            );
            (counter, bytes_to_trim)
        };
        match bytes_to_consume {
            0 => {
                if is_dangling_slash_r {
                    return Err(Error::InvalidState(
                        "\\r must be followed by \\n".to_owned(),
                    ));
                }
                break false;
            }
            v => {
                // Reject the chunk before addition if it would push the running total
                // beyond `line_cap`; this avoids silently saturating on `u16` overflow.
                if v > line_cap.saturating_sub(total_bytes_read) {
                    return Err(Error::InvalidState(format!(
                        "line is too long (>{line_cap} bytes)"
                    )));
                }
                total_bytes_read += v;
                assert!(v >= bytes_to_trim, "`v` is less than `bytes_to_trim`");
                assert!(
                    usize::from(v) <= buffered.len(),
                    "`v` cannot be greater than `buffered.len()`"
                );
                // SAFETY: the buffered bytes were already checked to be a
                // subset of valid ASCII before this unchecked UTF-8 conversion.
                unsafe {
                    line.push_str(str::from_utf8_unchecked(
                        &buffered[..usize::from(v - bytes_to_trim)],
                    ));
                }
                reader.consume(usize::from(v));
                if bytes_to_trim > 0 && !is_dangling_slash_r {
                    break true;
                }
            }
        }
    };
    assert!(!line.ends_with('\n'), "line cannot end with a new line!");
    assert!(!line.ends_with('\r'), "line cannot end with a new line!");
    assert!(
        line.len() <= usize::from(line_cap),
        "`line` longer than `line_cap` bytes read",
    );
    assert!(
        (is_newline_found && line.len() < usize::from(total_bytes_read))
            || (!is_newline_found && line.len() == usize::from(total_bytes_read)),
        "`line` must contain fewer or equal bytes than `total_bytes_read`",
    );
    Ok(total_bytes_read)
}

#[expect(clippy::panic, reason = "panic on error is standard practice in tests")]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{DNARestrictive, uuid};
    use rust_htslib::{bam::Read as _, faidx};
    use std::{
        io::{BufReader, Cursor},
        str::FromStr as _,
    };

    fn temp_output_dir(label: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!("nanalogue_{label}_{}", uuid::v4_random()));
        std::fs::create_dir_all(&path).expect("temp dir should be creatable");
        path
    }

    fn write_test_reference(dir: &Path) -> PathBuf {
        let fasta_path = dir.join("reference.fa");
        write_fasta(
            [(
                "chr1".to_string(),
                DNARestrictive::from_str("ACGTACGTACGT").expect("valid DNA"),
            )],
            &fasta_path,
        )
        .expect("reference FASTA should be written");
        faidx::build(&fasta_path).expect("FASTA index should be written");
        fasta_path
    }

    /// Tests writing to a fasta file and check its contents
    #[test]
    fn write_fasta_works() {
        let contigs = vec![
            (
                "test_contig_0".to_string(),
                DNARestrictive::from_str("ACGT").expect("no error"),
            ),
            (
                "test_contig_1".to_string(),
                DNARestrictive::from_str("TGCA").expect("no error"),
            ),
        ];

        let temp_path = std::env::temp_dir().join(format!("{}.fa", uuid::v4_random()));
        write_fasta(contigs, &temp_path).expect("no error");

        let content = std::fs::read_to_string(&temp_path).expect("no error");
        assert_eq!(content, ">test_contig_0\nACGT\n>test_contig_1\nTGCA\n");

        std::fs::remove_file(&temp_path).expect("no error");
    }

    /// Tests reading BAM file with valid path
    #[test]
    fn nanalogue_bam_reader_valid_path() {
        let mut reader = match nanalogue_bam_reader("examples/example_1.bam") {
            Ok(r) => r,
            Err(e) => panic!("Failed to read BAM file: {e:?}"),
        };
        assert_eq!(reader.records().count(), 4);
    }

    /// Tests reading BAM file with invalid path
    #[test]
    fn nanalogue_bam_reader_invalid_path() {
        let result = nanalogue_bam_reader("nonexistent_file.bam");
        let _err = result.unwrap_err();
    }

    /// Tests indexed BAM reader with `FetchDefinition::All`
    #[test]
    fn nanalogue_indexed_bam_reader_fetch_all() {
        let mut reader =
            match nanalogue_indexed_bam_reader("examples/example_1.bam", bam::FetchDefinition::All)
            {
                Ok(r) => r,
                Err(e) => panic!("Failed to read indexed BAM file: {e:?}"),
            };
        assert_eq!(reader.records().count(), 4);
    }

    /// Tests indexed BAM reader with specific contig
    #[test]
    fn nanalogue_indexed_bam_reader_fetch_contig() {
        let mut reader = match nanalogue_indexed_bam_reader(
            "examples/example_1.bam",
            bam::FetchDefinition::String(b"dummyI"),
        ) {
            Ok(r) => r,
            Err(e) => panic!("Failed to read indexed BAM file: {e:?}"),
        };
        assert_eq!(reader.records().count(), 1);
    }

    /// Tests indexed BAM reader with nonexistent contig
    #[test]
    fn nanalogue_indexed_bam_reader_fetch_nonexistent_contig() {
        let result = nanalogue_indexed_bam_reader(
            "examples/example_1.bam",
            bam::FetchDefinition::String(b"nonexistent_contig"),
        );
        let _err = result.unwrap_err();
    }

    /// Tests indexed BAM reader with region that has reads
    #[test]
    fn nanalogue_indexed_bam_reader_fetch_region_with_reads() {
        let mut reader = match nanalogue_indexed_bam_reader(
            "examples/example_1.bam",
            bam::FetchDefinition::RegionString(b"dummyIII", 20, 30),
        ) {
            Ok(r) => r,
            Err(e) => panic!("Failed to read indexed BAM file: {e:?}"),
        };
        assert_eq!(reader.records().count(), 1);
    }

    /// Tests indexed BAM reader with region that has no reads
    #[test]
    fn nanalogue_indexed_bam_reader_fetch_region_without_reads() {
        let mut reader = match nanalogue_indexed_bam_reader(
            "examples/example_1.bam",
            bam::FetchDefinition::RegionString(b"dummyIII", 10, 20),
        ) {
            Ok(r) => r,
            Err(e) => panic!("Failed to read indexed BAM file: {e:?}"),
        };
        assert_eq!(reader.records().count(), 0);
    }

    /// Tests `write_bam_denovo` with empty reads
    #[test]
    fn write_bam_denovo_empty_reads() {
        let contigs = vec![("chr1".to_string(), 1000)];
        let read_groups = vec!["rg1".to_string()];
        let comments = vec!["test comment".to_string()];
        let reads: Vec<bam::Record> = vec![];

        let temp_path = std::env::temp_dir().join(format!("{}.bam", uuid::v4_random()));
        match write_bam_denovo(reads, contigs, read_groups, comments, &temp_path) {
            Ok(()) => (),
            Err(e) => panic!("Failed to write BAM file: {e:?}"),
        }

        // Verify the file was created
        let reader = nanalogue_bam_reader(temp_path.to_str().unwrap()).expect("no error");
        assert_eq!(reader.header().target_count(), 1);

        std::fs::remove_file(&temp_path).expect("no error");
        std::fs::remove_file(format!("{}.bai", temp_path.display())).expect("no error");
    }

    /// Tests `write_bam_denovo` with unsorted reads returns error
    #[test]
    #[expect(
        clippy::similar_names,
        reason = "read1, read2, and reads are clear in this test context"
    )]
    fn write_bam_denovo_unsorted_reads_error() {
        let contigs = vec![("chr1".to_string(), 1000), ("chr2".to_string(), 1000)];
        let read_groups = vec!["rg1".to_string()];
        let comments = vec![];

        // Create two unsorted reads
        let mut read1 = bam::Record::new();
        read1.set_tid(1);
        read1.set_pos(100);

        let mut read2 = bam::Record::new();
        read2.set_tid(0);
        read2.set_pos(50);

        let reads = vec![read1, read2]; // Unsorted: tid 1 before tid 0

        let temp_path = std::env::temp_dir().join(format!("{}.bam", uuid::v4_random()));
        let result = write_bam_denovo(reads, contigs, read_groups, comments, &temp_path);

        let err = result.unwrap_err();
        assert!(matches!(err, Error::InvalidSorting(_)));

        // Clean up temporary files (ignore errors if files don't exist)
        drop(std::fs::remove_file(&temp_path));
        drop(std::fs::remove_file(format!("{}.bai", temp_path.display())));
    }

    /// Tests the forward-before-reverse tie-break in de novo alignment ordering.
    #[test]
    fn sorted_reads_reject_reverse_before_forward_at_same_coordinate() {
        let mut reverse = bam::Record::new();
        reverse.set_tid(0);
        reverse.set_pos(50);
        reverse.set_reverse();

        let mut forward = bam::Record::new();
        forward.set_tid(0);
        forward.set_pos(50);

        let result = write_sorted_reads(vec![reverse, forward], |_| Ok(()));
        assert!(matches!(result, Err(Error::InvalidSorting(_))));
    }

    /// Tests `write_bam_denovo` with sorted reads on same contig
    #[test]
    #[expect(
        clippy::similar_names,
        reason = "read1, read2, and reads are clear in this test context"
    )]
    fn write_bam_denovo_sorted_reads_same_contig() {
        let contigs = vec![("chr1".to_string(), 1000)];
        let read_groups = vec!["rg1".to_string()];
        let comments = vec![];

        // Create sorted reads on same contig
        let mut read1 = bam::Record::new();
        read1.set_tid(0);
        read1.set_pos(50);
        read1.set(b"read1", None, b"ACGT", &[30, 30, 30, 30]);

        let mut read2 = bam::Record::new();
        read2.set_tid(0);
        read2.set_pos(100);
        read2.set(b"read2", None, b"TGCA", &[30, 30, 30, 30]);

        let reads = vec![read1, read2];

        let temp_path = std::env::temp_dir().join(format!("{}.bam", uuid::v4_random()));
        match write_bam_denovo(reads, contigs, read_groups, comments, &temp_path) {
            Ok(()) => (),
            Err(e) => panic!("Failed to write BAM file: {e:?}"),
        }

        // Verify reads were written
        let mut reader = nanalogue_bam_reader(temp_path.to_str().unwrap()).expect("no error");
        assert_eq!(reader.records().count(), 2);

        std::fs::remove_file(&temp_path).expect("no error");
        std::fs::remove_file(format!("{}.bai", temp_path.display())).expect("no error");
    }

    /// Tests `write_bam_denovo` with multiple contigs and read groups
    #[test]
    fn write_bam_denovo_multiple_contigs_and_read_groups() {
        let contigs = vec![
            ("chr1".to_string(), 1000),
            ("chr2".to_string(), 2000),
            ("chr3".to_string(), 1500),
        ];
        let read_groups = vec!["rg1".to_string(), "rg2".to_string()];
        let comments = vec![
            "comment1".to_string(),
            "comment2".to_string(),
            "comment3".to_string(),
        ];
        let reads: Vec<bam::Record> = vec![];

        let temp_path = std::env::temp_dir().join(format!("{}.bam", uuid::v4_random()));
        match write_bam_denovo(reads, contigs, read_groups, comments, &temp_path) {
            Ok(()) => (),
            Err(e) => panic!("Failed to write BAM file: {e:?}"),
        }

        let reader = nanalogue_bam_reader(temp_path.to_str().unwrap()).expect("no error");
        assert_eq!(reader.header().target_count(), 3);

        // Verify read groups exist in header
        let header = reader.header();
        let header_text = std::str::from_utf8(header.as_bytes()).expect("no error");
        assert!(header_text.contains("@RG\tID:rg1"));
        assert!(header_text.contains("@RG\tID:rg2"));

        std::fs::remove_file(&temp_path).expect("no error");
        std::fs::remove_file(format!("{}.bai", temp_path.display())).expect("no error");
    }

    /// Tests `write_bam_denovo` with unmapped reads
    #[test]
    fn write_bam_denovo_with_unmapped_reads() {
        let contigs = vec![("chr1".to_string(), 1000)];
        let read_groups = vec!["rg1".to_string()];
        let comments = vec![];

        let mut read = bam::Record::new();
        read.set(b"unmapped_read", None, b"ACGT", &[30, 30, 30, 30]);
        read.set_unmapped();

        let reads = vec![read];

        let temp_path = std::env::temp_dir().join(format!("{}.bam", uuid::v4_random()));
        match write_bam_denovo(reads, contigs, read_groups, comments, &temp_path) {
            Ok(()) => (),
            Err(e) => panic!("Failed to write BAM file: {e:?}"),
        }

        let mut reader = nanalogue_bam_reader(temp_path.to_str().unwrap()).expect("no error");
        let records: Vec<_> = reader
            .records()
            .collect::<Result<Vec<_>, _>>()
            .expect("no error");
        assert_eq!(records.len(), 1);
        let record = records.first().expect("record exists");
        assert!(record.is_unmapped());
        assert_eq!(record.mapq(), 0);

        std::fs::remove_file(&temp_path).expect("no error");
        std::fs::remove_file(format!("{}.bai", temp_path.display())).expect("no error");
    }

    #[test]
    fn write_cram_denovo_rejects_colliding_reference_sidecar_paths() {
        let temp_dir = temp_output_dir("cram_collision");
        let expected_error = "CRAM, CRAI, and FASTA outputs must use different paths";

        // `write_cram_denovo` derives the CRAI as `<output>.crai`, so using
        // `output.cram.crai` as the reference path makes the reference collide with
        // the writer's own CRAI sidecar.
        let crai_collision = write_cram_denovo(
            Vec::<bam::Record>::new(),
            [("chr1".to_string(), 12)],
            ["rg1".to_string()],
            Vec::<String>::new(),
            &temp_dir.join("output.cram"),
            &temp_dir.join("output.cram.crai"),
            NonZeroU32::MIN,
        );

        assert!(matches!(
            crai_collision,
            Err(Error::InvalidState(msg)) if msg == expected_error
        ));

        // `write_cram_denovo` also rejects using the reference `.fai` path as the CRAM output,
        // because that would make the CRAM overwrite the FASTA index sidecar.
        let fai_collision = write_cram_denovo(
            Vec::<bam::Record>::new(),
            [("chr1".to_string(), 12)],
            ["rg1".to_string()],
            Vec::<String>::new(),
            &temp_dir.join("reference.fa.fai"),
            &temp_dir.join("reference.fa"),
            NonZeroU32::MIN,
        );

        assert!(matches!(
            fai_collision,
            Err(Error::InvalidState(msg)) if msg == expected_error
        ));
        std::fs::remove_dir_all(temp_dir).expect("temp dir should be removable");
    }

    #[test]
    fn write_cram_denovo_rejects_dot_dot_alias_before_creating_output() {
        let temp_dir = temp_output_dir("cram_dot_dot_collision");
        let child_dir = temp_dir.join("child");
        std::fs::create_dir_all(&child_dir).expect("child dir should be creatable");
        let reference_path = temp_dir.join("output.cram.crai");
        std::fs::write(&reference_path, b"reference sentinel")
            .expect("reference sentinel should be writable");
        let output_path = child_dir.join("..").join("output.cram");

        let result = write_cram_denovo(
            Vec::<bam::Record>::new(),
            [("chr1".to_string(), 12)],
            ["rg1".to_string()],
            Vec::<String>::new(),
            &output_path,
            &reference_path,
            NonZeroU32::MIN,
        );

        assert!(matches!(
            result,
            Err(Error::InvalidState(msg))
                if msg == "CRAM, CRAI, and FASTA outputs must use different paths"
        ));
        assert!(!temp_dir.join("output.cram").exists());
        assert_eq!(
            std::fs::read(&reference_path).expect("reference sentinel should remain readable"),
            b"reference sentinel"
        );
        std::fs::remove_dir_all(temp_dir).expect("temp dir should be removable");
    }

    #[cfg(unix)]
    #[test]
    fn output_paths_check_lexical_duplicates_before_filesystem_identity() {
        let temp_dir = temp_output_dir("lexical_collision_precedence");
        let path = temp_dir.join("self-referential");
        std::os::unix::fs::symlink(&path, &path)
            .expect("self-referential symlink should be creatable");

        let result = output_paths_are_distinct(&[&path, &path]);

        assert!(matches!(result, Ok(false)));
        std::fs::remove_dir_all(temp_dir).expect("temp dir should be removable");
    }

    #[cfg(unix)]
    #[test]
    fn write_cram_denovo_rejects_hard_link_to_reference_before_writing() {
        let temp_dir = temp_output_dir("cram_hard_link_collision");
        let reference_path = write_test_reference(&temp_dir);
        let original_reference = std::fs::read(&reference_path)
            .expect("reference contents should be readable before writing");
        let output_path = temp_dir.join("output.cram");
        std::fs::hard_link(&reference_path, &output_path)
            .expect("hard link to reference should be creatable");

        let result = write_cram_denovo(
            Vec::<bam::Record>::new(),
            [("chr1".to_string(), 12)],
            ["rg1".to_string()],
            Vec::<String>::new(),
            &output_path,
            &reference_path,
            NonZeroU32::MIN,
        );

        assert!(matches!(
            result,
            Err(Error::InvalidState(msg))
                if msg == "CRAM, CRAI, and FASTA outputs must use different paths"
        ));
        assert_eq!(
            std::fs::read(&reference_path).expect("reference contents should remain readable"),
            original_reference,
            "reference contents must remain unchanged"
        );
        std::fs::remove_dir_all(temp_dir).expect("temp dir should be removable");
    }

    #[cfg(unix)]
    #[test]
    fn write_cram_denovo_rejects_symlinked_crai_alias_before_creating_output() {
        let temp_dir = temp_output_dir("cram_symlink_collision");
        let reference_path = temp_dir.join("reference.fa");
        std::fs::write(&reference_path, b"reference sentinel")
            .expect("reference sentinel should be writable");
        let output_path = temp_dir.join("output.cram");
        std::os::unix::fs::symlink(&reference_path, crai_path(&output_path))
            .expect("CRAI symlink should be creatable");

        let result = write_cram_denovo(
            Vec::<bam::Record>::new(),
            [("chr1".to_string(), 12)],
            ["rg1".to_string()],
            Vec::<String>::new(),
            &output_path,
            &reference_path,
            NonZeroU32::MIN,
        );

        assert!(matches!(
            result,
            Err(Error::InvalidState(msg))
                if msg == "CRAM, CRAI, and FASTA outputs must use different paths"
        ));
        assert!(!output_path.exists());
        assert_eq!(
            std::fs::read(&reference_path).expect("reference sentinel should remain readable"),
            b"reference sentinel"
        );
        std::fs::remove_dir_all(temp_dir).expect("temp dir should be removable");
    }

    #[cfg(unix)]
    #[test]
    fn write_cram_denovo_rejects_dangling_crai_symlink_to_output() {
        let temp_dir = temp_output_dir("cram_dangling_symlink_collision");
        let reference_path = temp_dir.join("reference.fa");
        std::fs::write(&reference_path, b"reference sentinel")
            .expect("reference sentinel should be writable");
        let output_path = temp_dir.join("output.cram");
        let crai_path = crai_path(&output_path);
        std::os::unix::fs::symlink(&output_path, &crai_path)
            .expect("dangling CRAI symlink should be creatable");

        let result = write_cram_denovo(
            Vec::<bam::Record>::new(),
            [("chr1".to_string(), 12)],
            ["rg1".to_string()],
            Vec::<String>::new(),
            &output_path,
            &reference_path,
            NonZeroU32::MIN,
        );

        assert!(matches!(
            result,
            Err(Error::InvalidState(msg))
                if msg == "output paths must not be symbolic links to files that do not exist"
        ));
        assert!(!output_path.exists());
        assert!(
            std::fs::symlink_metadata(crai_path).is_ok(),
            "CRAI symlink should not be replaced"
        );
        assert_eq!(
            std::fs::read(&reference_path).expect("reference sentinel should remain readable"),
            b"reference sentinel"
        );
        std::fs::remove_dir_all(temp_dir).expect("temp dir should be removable");
    }

    #[test]
    fn cram_writer_open_reports_output_open_failure() {
        let temp_dir = temp_output_dir("cram_open_failure");
        let fasta_path = write_test_reference(&temp_dir);
        let missing_dir_output = temp_dir.join("missing").join("output.cram");

        let result = CramWriter::open(&missing_dir_output, &fasta_path, NonZeroU32::MIN);

        assert!(
            matches!(result, Err(Error::WriteOutput(msg)) if msg == "failed to open CRAM output")
        );
        std::fs::remove_dir_all(temp_dir).expect("temp dir should be removable");
    }

    #[test]
    fn cram_writer_open_reports_missing_reference_failure() {
        let temp_dir = temp_output_dir("cram_reference_failure");
        let missing_reference = temp_dir.join("missing.fa");
        let cram_path = temp_dir.join("output.cram");

        let result = CramWriter::open(&cram_path, &missing_reference, NonZeroU32::MIN);

        assert!(
            matches!(result, Err(Error::WriteOutput(msg)) if msg == "failed to load FASTA reference")
        );
        drop(std::fs::remove_file(cram_path));
        std::fs::remove_dir_all(temp_dir).expect("temp dir should be removable");
    }

    #[test]
    fn cram_writer_reports_index_init_failure() {
        let temp_dir = temp_output_dir("cram_header_failure");
        let fasta_path = write_test_reference(&temp_dir);
        let cram_path = temp_dir.join("output.cram");
        let mut writer =
            CramWriter::open(&cram_path, &fasta_path, NonZeroU32::MIN).expect("writer should open");
        // HTSlib does not require the CRAI name to match the CRAM name, only that it be a
        // writable file path. We deliberately replace it with a directory path here so
        // `sam_idx_init` fails deterministically and exercises the error branch.
        writer.index_path = path_to_c_string(&temp_dir).expect("temp dir path should convert");
        let header_template = denovo_alignment_header(
            [("chr1".to_string(), 12)],
            ["rg1".to_string()],
            Vec::<String>::new(),
        );
        let mut header = bam::HeaderView::from_header(&header_template);

        let result = writer.write_header_and_start_index(&mut header);

        assert!(
            matches!(result, Err(Error::WriteOutput(msg)) if msg == "failed to initialize CRAI index")
        );
        std::fs::remove_dir_all(temp_dir).expect("temp dir should be removable");
    }

    #[test]
    fn read_line_capped_rejects_slash_r_split_across_buffer_boundary() {
        let payload = format!("{}\rX", "a".repeat(7));
        let cursor = Cursor::new(payload.into_bytes());
        let mut reader = BufReader::with_capacity(8, cursor);
        let mut line = String::new();

        let result = read_line_capped(&mut reader, &mut line, 64);

        assert!(
            matches!(result, Err(Error::InvalidState(msg)) if msg == "\\r must be followed by \\n"),
            "\\r split across buffer boundary should fail"
        );
    }

    #[test]
    fn read_line_capped_accepts_slash_r_slash_n_split_across_buffer_boundary() {
        let payload = format!("{}\r\n", "a".repeat(7));
        let cursor = Cursor::new(payload.into_bytes());
        let mut reader = BufReader::with_capacity(8, cursor);
        let mut line = String::new();

        let bytes_read = read_line_capped(&mut reader, &mut line, 64)
            .expect("\\r\\n across buffer boundary should be accepted");

        assert_eq!(line, "aaaaaaa");
        assert_eq!(bytes_read, 9);
    }
}
