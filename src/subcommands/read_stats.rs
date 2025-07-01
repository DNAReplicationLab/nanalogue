use rust_htslib::{bam, bam::Read, bam::ext::BamRecordExtensions};
use std::error::Error;
use std::collections::BinaryHeap;

// A read can exist in four states
enum ReadState {
	Primary,
	Secondary,
	Supplementary,
	Unmapped,
}

// Assume a read is primary unless otherwise stated,
// so we have three possible transitions
enum ReadTransition {
	PrimaryToSecondary,
	PrimaryToSupplementary,
	PrimaryToUnmapped,
}

struct CurrRead {
	state: ReadState,
}

impl CurrRead {
	fn new() -> Self {
		Self { state: ReadState::Primary }
	}
	fn transition(&mut self, transition: ReadTransition){
		match (&self.state, transition) {
            (ReadState::Primary, ReadTransition::PrimaryToSecondary) => self.state = ReadState::Secondary,
            (ReadState::Primary, ReadTransition::PrimaryToSupplementary) => self.state = ReadState::Supplementary,
            (ReadState::Primary, ReadTransition::PrimaryToUnmapped) => self.state = ReadState::Unmapped,
            _ => {
                eprintln!("Invalid state reached!");
                std::process::exit(1);
            }
		}
	}
}

pub fn run(bam_path: &str) -> Result<(), Box<dyn Error>> {

	// open BAM file
    let mut bam = match bam::Reader::from_path(bam_path) {
        Ok(v) => {
            v
        },
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
	}

    println!("Number of primary reads: {primary_count}");
    println!("Number of secondary reads: {secondary_count}");
    println!("Number of supplementary reads: {supplementary_count}");
    println!("Number of unmapped reads: {unmapped_count}");

    if let Some(v) = align_len_heap.pop() { println!("Longest alignment length: {v}") }

    if let Some(v) = seq_len_heap.pop() { println!("Longest sequence length: {v}") }

    Ok(())
}

