//! Utils module providing shared datatypes for nanalogue
//! Includes genomic coordinates, constrained numerics, and BAM-related types

pub mod allowed_agctn;
pub mod contains;
pub mod dna_restrictive;
pub mod f32_abs_val_at_most1;
pub mod f32_bw0and1;
pub mod filter_by_ref_coords;
pub mod genomic_region;
pub mod intersects;
pub mod mod_char;
pub mod ord_pair;
pub mod path_or_url_or_stdin;
pub mod read_state;
pub mod read_states;
pub mod restrict_mod_called_strand;
pub mod threshold_state;

// Re-export public utility types and functions to expose the utils API
pub use allowed_agctn::*;
pub use contains::*;
pub use dna_restrictive::*;
pub use f32_abs_val_at_most1::*;
pub use f32_bw0and1::*;
pub use filter_by_ref_coords::*;
pub use genomic_region::*;
pub use intersects::*;
pub use mod_char::*;
pub use ord_pair::*;
pub use path_or_url_or_stdin::*;
pub use read_state::*;
pub use read_states::*;
pub use restrict_mod_called_strand::*;
pub use threshold_state::*;
