//! Central location for crate-level constants.
//!
//! These constants are grouped by feature/module to avoid polluting the crate
//! root namespace.

/// Constants shared across multiple commands/modules.
pub mod shared {
    /// Hard cap on the maximum capacity (bytes) of a BAM record's internal
    /// `data` buffer (`bam1_t::m_data`) that commands will accept.
    pub const MAX_RECORD_CAPACITY_BYTES: u32 = 32 * 1024 * 1024 - 1;

    /// Hard cap on the maximum number of records that commands will accept
    pub const MAX_RECORDS: u32 = 500_000_000;
}

/// Constants used by the `peek` subcommand/module.
pub mod peek {
    /// Hard cap on the number of records processed by `peek`.
    /// `peek` is intended as a fast metadata inspection command,
    /// so it is reasonable that this is lower than [`crate::constants::shared::MAX_RECORDS`].
    pub const MAX_RECORDS: u32 = 100_000;

    /// Hard cap on the number of distinct modification strings shown by `peek`.
    pub const MAX_MODIFICATIONS: u32 = 100;

    /// Hard cap on the number of contigs displayed by `peek`.
    pub const MAX_CONTIGS: u32 = 1_000_000;
}
