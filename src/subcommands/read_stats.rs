use rust_htslib::{bam, bam::Read, bam::ext::BamRecordExtensions};
use std::error::Error;
use std::collections::BinaryHeap;
use crate::{ReadState, ReadTransition, CurrRead};

fn get_stats_from_heap(mut input: BinaryHeap::<u64>, total_length: u64) -> (u64, f32, u64, u64, u64, u64){
    // process heaps to get statistics
    let mut counter :u64 = 0;
    let mut longest :u64 = 0;
    let mut shortest :u64 = 0;
    let mut median :u64 = 0;
    let mut n50 :u64 = 0;
    let mut running_total_length :u64 = 0;
    let heap_size :u64 = match input.len().try_into(){
        Ok(v) => v,
        Err(_) => {
            eprintln!("Error while getting heap length. Is it too large?");
            std::process::exit(1);
        },
    };

    while let Some(v) = input.pop() {

        if counter == 0 { longest = v };

        running_total_length += v;

        if median == 0 && counter > heap_size / 2 {
            median = v;
        }

        if n50 == 0 && running_total_length > total_length / 2 {
            n50 = v;
        }

        shortest = v;

        counter += 1;
    };

    assert_eq!(running_total_length, total_length);

    let mean = match counter {
        0 => 0.0,
        v => (total_length as f32) / (v as f32 * 1.0),
    };

    (counter, mean, longest, shortest, median, n50)
}

pub fn run(bam_path: &str) -> Result<(), Box<dyn Error>> {

    // open BAM file
    let mut bam = match bam::Reader::from_path(bam_path) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("Problem opening file, error: {e}");
            std::process::exit(1)
        },
    };

    // declare counts for different types of reads
    let mut primary_count :u64 = 0;
    let mut secondary_count :u64 = 0;
    let mut supplementary_count :u64 = 0;
    let mut unmapped_count :u64 = 0;

    // create two read length heaps
    let mut seq_len_heap = BinaryHeap::<u64>::new();
    let mut align_len_heap = BinaryHeap::<u64>::new();

    // total read length variables
    let mut seq_len_total = 0;
    let mut align_len_total = 0;

    // Go record by record in the BAM file,
    for r in bam.records() {

        // read records
        let record = match r {
            Ok(v) => v,
            Err(e) => {
                eprintln!("Some error while reading records {e}");
                std::process::exit(1)
            },
        };

        // set the read state using the type of alignment
        let mut curr_read_state = CurrRead::new();
        if record.is_unmapped() { curr_read_state.transition(ReadTransition::PrimaryToUnmapped); }
        if record.is_secondary() { curr_read_state.transition(ReadTransition::PrimaryToSecondary); }
        if record.is_supplementary() { curr_read_state.transition(ReadTransition::PrimaryToSupplementary); }

        // increment read counter suitably
        match curr_read_state.state {
            ReadState::Primary => primary_count += 1,
            ReadState::Secondary => secondary_count += 1,
            ReadState::Supplementary => supplementary_count += 1,
            ReadState::Unmapped => unmapped_count += 1,
        };

        // get length of alignment
        match curr_read_state.state {
            ReadState::Unmapped => {},
            _ => {
                let align_len :u64 = match (record.reference_end() - record.pos()).try_into() {
                    Ok(v) => v,
                    Err(_) => {
                        eprintln!("Problem getting alignment length");
                        std::process::exit(1);
                    },
                };
                align_len_heap.push(align_len);
                align_len_total += align_len;
            },
        };

        // get length of sequence
        let seq_len :u64 = match record.seq_len().try_into() {
            Ok(v) => v,
            Err(_) => {
                eprintln!("Error while getting sequencing length!");
                std::process::exit(1);
            },
        };
        seq_len_heap.push(seq_len);
        seq_len_total += seq_len;
	}

    // process heaps to get statistics 
    let (_, seq_len_mean, seq_len_max, seq_len_min, seq_len_median, seq_len_n50) = 
        get_stats_from_heap(seq_len_heap, seq_len_total);
    let (_, align_len_mean, align_len_max, align_len_min, align_len_median, align_len_n50) = 
        get_stats_from_heap(align_len_heap, align_len_total);

    println!("# input bam {bam_path}");
    println!("key\tvalue");
    println!("n_primary_alignments\t{primary_count}");
    println!("n_secondary_alignments\t{secondary_count}");
    println!("n_supplementary_alignments\t{supplementary_count}");
    println!("n_unmapped_reads\t{unmapped_count}");
    println!("align_len_mean\t{align_len_mean}");
    println!("align_len_max\t{align_len_max}");
    println!("align_len_min\t{align_len_min}");
    println!("align_len_median\t{align_len_median}");
    println!("align_len_n50\t{align_len_n50}");
    println!("seq_len_mean\t{seq_len_mean}");
    println!("seq_len_max\t{seq_len_max}");
    println!("seq_len_min\t{seq_len_min}");
    println!("seq_len_median\t{seq_len_median}");
    println!("seq_len_n50\t{seq_len_n50}");

    Ok(())
}

