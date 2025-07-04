use std::fmt;

// Declare the modules. 
pub mod subcommands;

// A read can exist in four states
pub enum ReadState {
	Primary,
	Secondary,
	Supplementary,
	Unmapped,
}

// Assume a read is primary unless otherwise stated,
// so we have three possible transitions
pub enum ReadTransition {
	PrimaryToSecondary,
	PrimaryToSupplementary,
	PrimaryToUnmapped,
}

pub struct CurrRead {
	state: ReadState,
    read_id: Option<String>,
    seq_len: Option<u64>,
    align_len: Option<u64>,
}

impl CurrRead {
	fn new() -> Self {
		Self { 
            state: ReadState::Primary,
            read_id: None,
            seq_len: None,
            align_len: None,
        }
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
    fn header_string() -> String {
        "read_id\tsequence_length_template\talign_length\talignment_type".to_string()
    }
}

impl fmt::Display for CurrRead {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {

        let mut output_string = String::from("");

        if let Some(v) = &self.read_id {
            output_string = output_string + v;
        }

        if let Some(v) = self.seq_len {
            output_string = output_string + "\t" + &v.to_string();
        } else {
            output_string += "\tNA";
        }

        if let Some(v) = self.align_len {
            output_string = output_string + "\t" + &v.to_string();
        } else {
            output_string += "\tNA";
        }

        match self.state {
            ReadState::Primary => output_string = output_string + "\t" + "primary",
            ReadState::Secondary => output_string = output_string + "\t" + "secondary",
            ReadState::Supplementary => output_string = output_string + "\t" + "supplementary",
            ReadState::Unmapped => output_string = output_string + "\t" + "unmapped",
        }

        write!(f, "{output_string}")

    }
}
