//! Central location for crate-level constants.
//!
//! For various constants here, our goal is to protect ourselves against processing
//! pathologically large resources. So, we are stricter than the BAM format.
//! For example, the BAM format allows records whose memory footprint is several GB
//! (as you can have any number of tags that are arbitrarily long for example),
//! and allows read ids that are 255 characters long I believe.
//! We go stricter than this to protect ourselves to avoid out of memory errors (OOM) etc.
//! This is a best-faith effort as we still rely on htslib and as far as I know,
//! htslib will process pathologically long records and cause an OOM.
//!
//! These constants are grouped by feature/module to avoid polluting the crate
//! root namespace.

/// Constants shared across multiple commands/modules.
pub mod shared {
    /// Hard cap on the maximum capacity (bytes) of a BAM record's internal
    /// `data` buffer (`bam1_t::m_data`) that commands will accept.
    /// As far as I know, we cannot prevent htslib from loading a large record;
    /// we can only check after the fact. The hope here is if we run into a
    /// large record, it is likely that larger records are in the future, so
    /// we can exit safely if we want. It is also possible that we will run
    /// into an extremely large record without warning, in which case
    /// we will crash as we depend on htslib and there's no way to prevent this.
    pub const MAX_RECORD_CAPACITY_BYTES: u32 = 32 * 1024 * 1024 - 1;

    /// Hard cap on the maximum number of records that commands will accept
    pub const MAX_RECORDS: u32 = 500_000_000;

    /// Hard cap on the number of types of mods per BAM record
    pub const MAX_MOD_TYPES: u8 = 100;

    /// Hard cap on the total number of serialized modification annotations per read.
    pub const MAX_TOTAL_MOD_ANNOTATIONS_PER_READ: u32 = u32::MAX;

    /// Hard cap on the number of lines in an external text file of read ids for e.g. filtering
    pub const MAX_READ_IDS_FOR_FILTERING: u32 = 1_000_000;

    /// Hard cap on length of a read id
    pub const MAX_READ_ID_LEN: u8 = 200;

    /// Hard cap on length of a contig name.
    pub const MAX_CONTIG_NAME_LENGTH: u8 = 200;

    /// Hard cap on the number of contigs.
    pub const MAX_CONTIGS: u32 = 1_000_000;

    /// Hard cap on the length of a path e.g. to a file or a URL
    pub const MAX_PATH_LENGTH: u16 = 900;

    /// Hard cap on the length of an array from an ML tag
    pub const MAX_ML_ARRAY_LENGTH: u32 = 100_000_000;

    /// Hard cap on the length of an MM tag
    pub const MAX_MM_TAG_LENGTH: u32 = 100_000_000;

    /// Hard cap on the length of a genomic region string (e.g. chr1:1000-2000)
    pub const MAX_GENOMIC_REGION_STRING_LENGTH: u8 = 255;

    /// Shared error message used when analysis commands receive zero records.
    pub const NO_RECORDS_FOUND_FOR_ANALYSIS: &str = concat!(
        "No records found as input for analysis. This could mean ",
        "the input genomics file (SAM/BAM/CRAM) has no records or that filtering\n",
        "removed records or some other possibility.",
    );

    // a region is just contig:num-num or a shorter string, so the length of the
    // ":num-num" part is what we are asserting here. Even if num is in single-digit billions,
    // this string is just 22 characters (10 + 10 + 1 + 1). So 26 is comfortably above
    // any limit we can hit. When the program supports extremely large contigs, we would
    // need to revise this.
    const _: () = assert!(
        MAX_GENOMIC_REGION_STRING_LENGTH - MAX_CONTIG_NAME_LENGTH > 26,
        "genomic region must be large enough to accommodate coordinates"
    );
}

/// Constants used by the `peek` subcommand/module.
pub mod peek {
    /// Hard cap on the number of records processed by `peek`.
    /// `peek` is intended as a fast metadata inspection command,
    /// so it is reasonable that this is lower than [`crate::constants::shared::MAX_RECORDS`].
    pub const MAX_RECORDS: u32 = 100_000;
}

/// Constants used by the `reads_table` subcommand/module.
pub mod reads_table {
    /// Hard cap on the number of characters per line in the sequencing summary file.
    pub const MAX_SEQ_SUMM_SIZE_PER_LINE: u16 = 1_000;

    /// Hard cap on sequencing summary file size.
    pub const MAX_SEQ_SUMM_BYTES: u64 = 10u64 * 1024u64 * 1024u64 * 1024u64;
}
