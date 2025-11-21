//! Tests for `curr_reads_to_dataframe` function

use nanalogue_core::{
    AlignmentInfoBuilder, CurrReadBuilder, Error, ModTableEntryBuilder, ReadState,
    curr_reads_to_dataframe,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn curr_reads_to_dataframe_basic_example() -> Result<(), Error> {
        // Build modification table entries with two data points
        let mod_table_entry = ModTableEntryBuilder::default()
            .base('C')
            .is_strand_plus(true)
            .mod_code("m".into())
            .data([(0, 15, 200), (2, 25, 100)].into())
            .build()?;

        // Build a CurrRead<AlignAndModData>
        let read = CurrReadBuilder::default()
            .read_id("test_read".into())
            .seq_len(40)
            .alignment_type(ReadState::PrimaryFwd)
            .alignment(
                AlignmentInfoBuilder::default()
                    .start(10)
                    .end(60)
                    .contig("chr1".into())
                    .contig_id(1)
                    .build()?,
            )
            .mod_table([mod_table_entry].into())
            .build()?;

        // Convert to DataFrame
        let df = curr_reads_to_dataframe(&[read])?;

        // Verify DataFrame structure
        assert_eq!(df.height(), 2); // Two modification data points
        assert_eq!(df.width(), 13); // 13 columns

        // Check read-level fields (repeated across both rows)
        let read_id_col = df.column("read_id")?.str()?;
        assert_eq!(read_id_col.get(0), Some("test_read"));
        assert_eq!(read_id_col.get(1), Some("test_read"));

        let seq_len_col = df.column("seq_len")?.u64()?;
        assert_eq!(seq_len_col.get(0), Some(40));
        assert_eq!(seq_len_col.get(1), Some(40));

        let alignment_type_col = df.column("alignment_type")?.str()?;
        assert_eq!(alignment_type_col.get(0), Some("primary_forward"));
        assert_eq!(alignment_type_col.get(1), Some("primary_forward"));

        // Check alignment fields (repeated across both rows)
        let align_start_col = df.column("align_start")?.u64()?;
        assert_eq!(align_start_col.get(0), Some(10));
        assert_eq!(align_start_col.get(1), Some(10));

        let align_end_col = df.column("align_end")?.u64()?;
        assert_eq!(align_end_col.get(0), Some(60));
        assert_eq!(align_end_col.get(1), Some(60));

        let contig_col = df.column("contig")?.str()?;
        assert_eq!(contig_col.get(0), Some("chr1"));
        assert_eq!(contig_col.get(1), Some("chr1"));

        let contig_id_col = df.column("contig_id")?.i32()?;
        assert_eq!(contig_id_col.get(0), Some(1));
        assert_eq!(contig_id_col.get(1), Some(1));

        // Check modification entry fields (repeated across both rows)
        let base_col = df.column("base")?.str()?;
        assert_eq!(base_col.get(0), Some("C"));
        assert_eq!(base_col.get(1), Some("C"));

        let is_strand_plus_col = df.column("is_strand_plus")?.bool()?;
        assert_eq!(is_strand_plus_col.get(0), Some(true));
        assert_eq!(is_strand_plus_col.get(1), Some(true));

        let mod_code_col = df.column("mod_code")?.str()?;
        assert_eq!(mod_code_col.get(0), Some("m"));
        assert_eq!(mod_code_col.get(1), Some("m"));

        // Check modification data points (unique per row)
        let position_col = df.column("position")?.u64()?;
        assert_eq!(position_col.get(0), Some(0)); // First mod position
        assert_eq!(position_col.get(1), Some(2)); // Second mod position

        let ref_position_col = df.column("ref_position")?.i64()?;
        assert_eq!(ref_position_col.get(0), Some(15)); // First ref position
        assert_eq!(ref_position_col.get(1), Some(25)); // Second ref position

        let mod_quality_col = df.column("mod_quality")?.u32()?;
        assert_eq!(mod_quality_col.get(0), Some(200)); // First mod_quality score
        assert_eq!(mod_quality_col.get(1), Some(100)); // Second mod_quality score

        Ok(())
    }
}
