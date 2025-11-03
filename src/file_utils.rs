//! Utility functions for file I/O operations with BAM and FASTA files.

use crate::Error;
use rust_htslib::bam;
use std::fs::File;
use std::io::Write as _;
use std::path::Path;

/// Opens BAM file, copied and edited from fibertools-rs repo.
///
/// The fibertools repo (<https://github.com/fiberseq/fibertools-rs>)
/// is under the MIT license (please see their Cargo.toml).
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
pub fn nanalogue_bam_reader(bam_path: &str) -> Result<bam::Reader, Error> {
    if bam_path == "-" {
        Ok(bam::Reader::from_stdin()?)
    } else {
        Ok(bam::Reader::from_path(bam_path)?)
    }
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
/// assert!(matches!(err, Error::RustHtslibError(OtherError::BamInvalidIndex{..})));
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
/// assert!(matches!(err, Error::RustHtslibError(OtherError::BamInvalidIndex{..})));
/// ```
pub fn nanalogue_indexed_bam_reader(
    bam_path: &str,
    fetch_definition: bam::FetchDefinition,
) -> Result<bam::IndexedReader, Error> {
    let mut bam_reader = bam::IndexedReader::from_path(bam_path)?;
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
/// use nanalogue_core::{Error, file_utils::write_fasta};
/// use std::fs;
/// use uuid::Uuid;
///
/// let contigs = vec![
///     ("seq1".to_string(), b"ACGTACGT".to_vec()),
///     ("seq2".to_string(), b"TGCATGCA".to_vec()),
/// ];
///
/// let temp_path = std::env::temp_dir().join(format!("{}.fa", Uuid::new_v4()));
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
pub fn write_fasta<I, J>(contigs: I, output_path: &J) -> Result<(), Error>
where
    I: IntoIterator<Item = (String, Vec<u8>)>,
    J: AsRef<Path> + ?Sized,
{
    let mut file = File::create(output_path)?;
    for contig in contigs {
        writeln!(file, ">{}", contig.0)?;
        file.write_all(&contig.1)?;
        writeln!(file)?;
    }
    Ok(())
}

/// Writes a new BAM file with reads. Input reads have to be sorted.
/// Although this function can be used for other tasks like subsetting
/// BAM files, this is not advised as the history of the BAM file
/// stored in its header would be lost as we are generating a new header here.
///
/// # Errors
///
/// Returns an error if the BAM file cannot be created or written to,
/// if index creation fails, or if input reads are not sorted.
///
/// # Examples
///
/// ```
/// use nanalogue_core::{Error, file_utils::{write_bam_denovo, nanalogue_bam_reader}};
/// use rust_htslib::bam;
/// use rust_htslib::bam::Read;
/// use uuid::Uuid;
///
/// let contigs = vec![("chr1".to_string(), 1000)];
/// let read_groups = vec!["rg1".to_string()];
/// let comments = vec!["test comment".to_string()];
/// let reads: Vec<bam::Record> = vec![];
///
/// let temp_path = std::env::temp_dir().join(format!("{}.bam", Uuid::new_v4()));
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
    let header = {
        let mut header = bam::Header::new();
        for k in contigs {
            let _: &mut _ = header.push_record(
                bam::header::HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", k.0)
                    .push_tag(b"LN", k.1),
            );
        }
        for k in read_groups {
            let _: &mut _ = header.push_record(
                bam::header::HeaderRecord::new(b"RG")
                    .push_tag(b"ID", k)
                    .push_tag(b"PL", "ONT")
                    .push_tag(b"LB", "blank")
                    .push_tag(b"SM", "blank")
                    .push_tag(b"PU", "blank"),
            );
        }
        for k in comments {
            let _: &mut _ = header.push_comment(k.as_bytes());
        }
        header
    };

    // Write BAM file ensuring reads are already sorted
    let mut writer = bam::Writer::from_path(output_path, &header, bam::Format::Bam)?;
    let mut curr_read_key: (bool, i32, i64, bool);
    let mut prev_read_key: (bool, i32, i64, bool) = (false, -1, -1, false);
    for read in reads {
        curr_read_key = (
            read.is_unmapped(),
            read.tid(),
            read.pos(),
            read.is_reverse(),
        );
        if prev_read_key > curr_read_key {
            return Err(Error::InvalidSorting(
                "reads input to write bam denovo have not been sorted properly".to_string(),
            ));
        }
        prev_read_key = curr_read_key;
        writer.write(&read)?;
    }
    drop(writer); // Close BAM file before creating index

    bam::index::build(output_path, None, bam::index::Type::Bai, 2)?;

    Ok(())
}

#[expect(clippy::panic, reason = "panic on error is standard practice in tests")]
#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::Read as _;
    use uuid::Uuid;

    /// Tests writing to a fasta file and check its contents
    #[test]
    fn write_fasta_works() {
        let contigs = vec![
            ("test_contig_0".to_string(), b"ACGT".to_vec()),
            ("test_contig_1".to_string(), b"TGCA".to_vec()),
        ];

        let temp_path = std::env::temp_dir().join(format!("{}.fa", Uuid::new_v4()));
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

        let temp_path = std::env::temp_dir().join(format!("{}.bam", Uuid::new_v4()));
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

        let temp_path = std::env::temp_dir().join(format!("{}.bam", Uuid::new_v4()));
        let result = write_bam_denovo(reads, contigs, read_groups, comments, &temp_path);

        let err = result.unwrap_err();
        assert!(matches!(err, Error::InvalidSorting(_)));
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

        let temp_path = std::env::temp_dir().join(format!("{}.bam", Uuid::new_v4()));
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

        let temp_path = std::env::temp_dir().join(format!("{}.bam", Uuid::new_v4()));
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

        let temp_path = std::env::temp_dir().join(format!("{}.bam", Uuid::new_v4()));
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
        assert!(records.first().expect("record exists").is_unmapped());

        std::fs::remove_file(&temp_path).expect("no error");
        std::fs::remove_file(format!("{}.bai", temp_path.display())).expect("no error");
    }
}
