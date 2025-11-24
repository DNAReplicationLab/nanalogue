//! # Gets statistics on reads
//!
//! This module calculates statistics such as mean, median, N50
//! read and alignment lengths etc. from a BAM file.

use crate::{CurrRead, Error, ReadState};
use rust_htslib::bam;
use std::collections::BinaryHeap;
use std::rc::Rc;

#[expect(
    clippy::arithmetic_side_effects,
    reason = "sums can overflow but only if _huge_ amounts of data i.e. total sequence bp = (2^64-1)"
)]
/// Receives a heap of numbers and calculates statistics
/// (count, mean, longest, shortest, median, N50).
/// We assume the heap is long enough that we use approximate
/// formulae for e.g. the median. We assume there is enough read-to-read
/// variation that integer precision is ok for these quantities.
///
/// # Panics
/// Should not panic as heap size is not expected to exceed `u64::MAX` (~2^64),
/// and running total will not disagree with total length.
fn get_stats_from_heap(
    mut input: BinaryHeap<u64>,
    total_length: u64,
) -> (u64, u64, u64, u64, u64, u64) {
    // process heaps to get statistics
    let mut counter: u64 = 0;
    let mut longest: u64 = 0;
    let mut shortest: u64 = 0;
    let mut median: u64 = 0;
    let mut n50: u64 = 0;
    let mut running_total_length: u64 = 0;
    let heap_size: u64 = input.len().try_into().expect("heap cannot be this large");

    while let Some(v) = input.pop() {
        if counter == 0 {
            longest = v;
        }

        running_total_length += v;

        // In both median and N50 calculations, we assume lots of reads
        // and a continuous distribution. If not, the answers can be slightly
        // different... for example, for the median, technically speaking,
        // if the number of reads is even, we are supposed to take the mean
        // of the two lengths in the middle of the pack. We don't do this as
        // we assume lots of reads so there's no point in making such an accurate
        // calculation.
        if median == 0 && counter > heap_size.div_ceil(2).saturating_sub(1) {
            median = v;
        }

        if n50 == 0 && running_total_length > total_length.div_ceil(2) {
            n50 = v;
        }

        shortest = v;

        counter += 1;
    }

    assert_eq!(
        running_total_length, total_length,
        "running_total and pre-calculated total are expected to be equal"
    );

    // This should be a floating point number, but in any reasonable dataset
    // we will not need a precision < 1 bp.
    let mean = total_length.checked_div(counter).unwrap_or(0u64);

    (counter, mean, longest, shortest, median, n50)
}

/// Reads the input BAM file and prints statistics
/// such as mean and median read lengths, N50s etc.
///
/// # Errors
/// Returns an error if BAM record reading, parsing, or calculating or writing statistics fails.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "sums can overflow but only if _huge_ amounts of data i.e. total sequence bp = (2^64-1)"
)]
pub fn run<W, D>(handle: &mut W, bam_records: D) -> Result<(), Error>
where
    W: std::io::Write,
    D: IntoIterator<Item = Result<Rc<bam::Record>, rust_htslib::errors::Error>>,
{
    // declare counts for different types of reads
    let mut primary_count: u64 = 0;
    let mut secondary_count: u64 = 0;
    let mut supplementary_count: u64 = 0;
    let mut unmapped_count: u64 = 0;
    let mut reversed_count: u64 = 0;

    // create two read length heaps
    let mut seq_len_heap = BinaryHeap::<u64>::new();
    let mut align_len_heap = BinaryHeap::<u64>::new();

    // total read length variables
    let mut seq_len_total = 0;
    let mut align_len_total = 0;

    // Go record by record in the BAM file,
    for r in bam_records {
        // read records
        let record = r?;

        let curr_read = CurrRead::default()
            .set_read_state(&record)?
            .set_seq_len(&record)?;

        if curr_read.read_state().strand() == '-' {
            reversed_count += 1;
        }

        match curr_read.read_state() {
            ReadState::PrimaryFwd | ReadState::PrimaryRev => {
                primary_count += 1;
            }
            ReadState::SecondaryFwd | ReadState::SecondaryRev => {
                secondary_count += 1;
            }
            ReadState::SupplementaryFwd | ReadState::SupplementaryRev => {
                supplementary_count += 1;
            }
            ReadState::Unmapped => unmapped_count += 1,
        }

        // get length of sequence.
        if let Ok(v) = curr_read.seq_len() {
            seq_len_total += v;
            seq_len_heap.push(v);
        }

        // get length of alignment
        match curr_read.set_align_len(&record) {
            Ok(cr) => {
                let v = cr.align_len()?;
                align_len_total += v;
                align_len_heap.push(v);
            }
            Err(Error::Unmapped) => {}
            Err(e) => return Err(e),
        }
    }

    // process heaps to get statistics
    let (_, seq_len_mean, seq_len_max, seq_len_min, seq_len_median, seq_len_n50) =
        get_stats_from_heap(seq_len_heap, seq_len_total);
    let (_, align_len_mean, align_len_max, align_len_min, align_len_median, align_len_n50) =
        get_stats_from_heap(align_len_heap, align_len_total);

    writeln!(handle, "key\tvalue")?;
    writeln!(handle, "n_primary_alignments\t{primary_count}")?;
    writeln!(handle, "n_secondary_alignments\t{secondary_count}")?;
    writeln!(handle, "n_supplementary_alignments\t{supplementary_count}")?;
    writeln!(handle, "n_unmapped_reads\t{unmapped_count}")?;
    writeln!(handle, "n_reversed_reads\t{reversed_count}")?;
    writeln!(handle, "align_len_mean\t{align_len_mean}")?;
    writeln!(handle, "align_len_max\t{align_len_max}")?;
    writeln!(handle, "align_len_min\t{align_len_min}")?;
    writeln!(handle, "align_len_median\t{align_len_median}")?;
    writeln!(handle, "align_len_n50\t{align_len_n50}")?;
    writeln!(handle, "seq_len_mean\t{seq_len_mean}")?;
    writeln!(handle, "seq_len_max\t{seq_len_max}")?;
    writeln!(handle, "seq_len_min\t{seq_len_min}")?;
    writeln!(handle, "seq_len_median\t{seq_len_median}")?;
    writeln!(handle, "seq_len_n50\t{seq_len_n50}")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use indoc::indoc;
    use rust_htslib::bam::Read as _;

    #[test]
    fn read_stats_example_1() -> Result<(), Error> {
        let mut output = Vec::new();

        let mut bam_reader = bam::Reader::from_path("./examples/example_1.bam")
            .expect("Failed to open example_1.bam");

        let bam_records = bam_reader.rc_records();

        run(&mut output, bam_records)?;

        let output_str = String::from_utf8(output).expect("Invalid UTF-8 output");

        let expected_output = indoc! {"key\tvalue
        n_primary_alignments\t3
        n_secondary_alignments\t0
        n_supplementary_alignments\t0
        n_unmapped_reads\t1
        n_reversed_reads\t1
        align_len_mean\t29
        align_len_max\t48
        align_len_min\t8
        align_len_median\t8
        align_len_n50\t48
        seq_len_mean\t34
        seq_len_max\t48
        seq_len_min\t8
        seq_len_median\t33
        seq_len_n50\t48\n"};

        assert_eq!(output_str, expected_output);
        Ok(())
    }

    #[test]
    fn read_stats_example_3() -> Result<(), Error> {
        let mut output = Vec::new();

        let mut bam_reader = bam::Reader::from_path("./examples/example_3.bam")
            .expect("Failed to open example_3.bam");

        let bam_records = bam_reader.rc_records();

        run(&mut output, bam_records)?;

        let output_str = String::from_utf8(output).expect("Invalid UTF-8 output");

        let expected_output = indoc! {"key\tvalue
        n_primary_alignments\t6
        n_secondary_alignments\t2
        n_supplementary_alignments\t2
        n_unmapped_reads\t0
        n_reversed_reads\t3
        align_len_mean\t12
        align_len_max\t20
        align_len_min\t6
        align_len_median\t12
        align_len_n50\t15
        seq_len_mean\t12
        seq_len_max\t20
        seq_len_min\t6
        seq_len_median\t12
        seq_len_n50\t15\n"};

        assert_eq!(output_str, expected_output);
        Ok(())
    }
}
