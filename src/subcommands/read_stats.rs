//! # Gets statistics on reads
//!
//! This module calculates statistics such as mean, median, N50
//! read and alignment lengths etc. from a BAM file.

use crate::{CurrRead, Error, ReadState};
use rust_htslib::bam;
use std::collections::BinaryHeap;
use std::rc::Rc;

fn get_stats_from_heap(
    mut input: BinaryHeap<u64>,
    total_length: u64,
) -> Result<(u64, f32, u64, u64, u64, u64), Error> {
    // process heaps to get statistics
    let mut counter: u64 = 0;
    let mut longest: u64 = 0;
    let mut shortest: u64 = 0;
    let mut median: u64 = 0;
    let mut n50: u64 = 0;
    let mut running_total_length: u64 = 0;
    let heap_size: u64 = match input.len().try_into() {
        Ok(v) => v,
        Err(_) => Err(Error::RareHeapTooLarge)?,
    };

    while let Some(v) = input.pop() {
        if counter == 0 {
            longest = v
        };

        running_total_length += v;

        if median == 0 && counter > heap_size / 2 {
            median = v;
        }

        if n50 == 0 && running_total_length > total_length / 2 {
            n50 = v;
        }

        shortest = v;

        counter += 1;
    }

    assert_eq!(running_total_length, total_length);

    let mean = match counter {
        0 => 0.0,
        v => (total_length as f32) / (v as f32 * 1.0),
    };

    Ok((counter, mean, longest, shortest, median, n50))
}

/// Reads the input BAM file and prints statistics
/// such as mean and median read lengths, N50s etc.
pub fn run<W, D>(handle: &mut W, bam_records: D) -> Result<bool, Error>
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

        // set the read state using the type of alignment
        // and increment read counter
        let mut curr_read_state = CurrRead::default();
        curr_read_state.set_read_state(&record)?;
        match curr_read_state.read_state() {
            ReadState::Unknown => {}
            ReadState::PrimaryFwd => primary_count += 1,
            ReadState::SecondaryFwd => secondary_count += 1,
            ReadState::SupplementaryFwd => supplementary_count += 1,
            ReadState::Unmapped => unmapped_count += 1,
            ReadState::PrimaryRev => {
                primary_count += 1;
                reversed_count += 1
            }
            ReadState::SecondaryRev => {
                secondary_count += 1;
                reversed_count += 1
            }
            ReadState::SupplementaryRev => {
                supplementary_count += 1;
                reversed_count += 1
            }
        };

        // get length of alignment
        if let Some(v) = curr_read_state.set_align_len(&record)? {
            align_len_total += v;
            align_len_heap.push(v);
        };

        // get length of sequence
        if let Some(v) = curr_read_state.set_seq_len(&record)? {
            seq_len_total += v;
            seq_len_heap.push(v);
        };
    }

    // process heaps to get statistics
    let (_, seq_len_mean, seq_len_max, seq_len_min, seq_len_median, seq_len_n50) =
        get_stats_from_heap(seq_len_heap, seq_len_total)?;
    let (_, align_len_mean, align_len_max, align_len_min, align_len_median, align_len_n50) =
        get_stats_from_heap(align_len_heap, align_len_total)?;

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

    Ok(true)
}
