//! # FileUtils
//!
//! Utility functions for file I/O operations with BAM and FASTA files.

use crate::Error;
use rust_htslib::bam;
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// Opens BAM file, copied and edited from fiberseq repo.
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

/// Writes contigs to a FASTA file
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
            let _: &mut bam::Header = header.push_record(
                bam::header::HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", k.0)
                    .push_tag(b"LN", k.1),
            );
        }
        for k in read_groups {
            let _: &mut bam::Header = header.push_record(
                bam::header::HeaderRecord::new(b"RG")
                    .push_tag(b"ID", k)
                    .push_tag(b"PL", "ONT")
                    .push_tag(b"LB", "blank")
                    .push_tag(b"SM", "blank")
                    .push_tag(b"PU", "blank"),
            );
        }
        for k in comments {
            let _: &mut bam::Header = header.push_comment(k.as_bytes());
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
        } else {
            prev_read_key = curr_read_key;
            writer.write(&read)?;
        }
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
    fn test_write_fasta() {
        let contigs = vec![
            ("test_contig_0".to_string(), b"ACGT".to_vec()),
            ("test_contig_1".to_string(), b"TGCA".to_vec()),
        ];

        let temp_path = std::env::temp_dir().join(format!("{}.fa", Uuid::new_v4()));
        write_fasta(contigs.into_iter(), &temp_path).expect("no error");

        let content = std::fs::read_to_string(&temp_path).expect("no error");
        assert_eq!(content, ">test_contig_0\nACGT\n>test_contig_1\nTGCA\n");

        std::fs::remove_file(&temp_path).expect("no error");
    }
}
