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
/// ```
/// use nanalogue_core::{Error, file_utils::nanalogue_indexed_bam_reader};
/// use rust_htslib::bam::{Read, FetchDefinition};
/// let mut reader = nanalogue_indexed_bam_reader(&"examples/example_1.bam", FetchDefinition::All)?;
/// // the above file should contain four reads, so we are checking
/// // if we load four records.
/// assert_eq!(reader.records().count(), 4);
/// # Ok::<(), Error>(())
/// ```
/// ```
/// # use nanalogue_core::{Error, file_utils::nanalogue_indexed_bam_reader};
/// # use rust_htslib::bam::{Read, FetchDefinition};
/// let mut reader = nanalogue_indexed_bam_reader(&"examples/example_1.bam",
///     FetchDefinition::String(b"dummyI"))?;
/// // the above file should contain only one read passing through this contig.
/// assert_eq!(reader.records().count(), 1);
/// # Ok::<(), Error>(())
/// ```
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
/// # Ok::<(), Error>(())
/// ```
/// ```should_panic
/// # use nanalogue_core::{Error, file_utils::nanalogue_indexed_bam_reader};
/// # use rust_htslib::bam::{Read, FetchDefinition};
/// let mut reader = nanalogue_indexed_bam_reader(&"examples/example_1.bam",
///     FetchDefinition::CompleteTid(10))?;
/// // the above panics as this file has much fewer than 10 contigs.
/// # Ok::<(), Error>(())
/// ```
/// ```
/// # use nanalogue_core::{Error, file_utils::nanalogue_indexed_bam_reader};
/// # use rust_htslib::bam::{Read, FetchDefinition};
/// use rust_htslib::errors::Error as OtherError;
/// let err = nanalogue_indexed_bam_reader(&"examples/example_1_copy_no_index.bam",
///     FetchDefinition::All).unwrap_err();
/// // we should get a missing index error.
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

#[cfg(test)]
mod tests {
    use super::*;
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
}
