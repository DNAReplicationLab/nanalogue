//! Tests for `curr_reads_to_dataframe` function

use nanalogue_core::{
    AlignmentInfoBuilder, CurrReadBuilder, Error, ModTableEntryBuilder, ReadState,
    curr_reads_to_dataframe,
};

#[cfg(test)]
mod tests {
    use super::*;

    /// Test unmapped read with minimal information (no mods)
    #[test]
    fn unmapped_minimal() -> Result<(), Error> {
        let read = CurrReadBuilder::default().build()?;
        let df = curr_reads_to_dataframe(&[read])?;

        // No modifications means DataFrame should be empty
        assert_eq!(df.height(), 0);
        // But should still have all 13 columns
        assert_eq!(df.width(), 13);

        Ok(())
    }

    /// Test unmapped read with basic info (no mods)
    #[test]
    fn unmapped_with_info() -> Result<(), Error> {
        let read = CurrReadBuilder::default()
            .read_id("some_read".into())
            .seq_len(40)
            .build()?;
        let df = curr_reads_to_dataframe(&[read])?;

        // No modifications means DataFrame should be empty
        assert_eq!(df.height(), 0);
        assert_eq!(df.width(), 13);

        Ok(())
    }

    /// Test mapped read without modification information
    #[test]
    fn mapped_basic() -> Result<(), Error> {
        let read = CurrReadBuilder::default()
            .read_id("some_read".into())
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
            .build()?;
        let df = curr_reads_to_dataframe(&[read])?;

        // No modifications means DataFrame should be empty
        assert_eq!(df.height(), 0);
        assert_eq!(df.width(), 13);

        Ok(())
    }

    /// Test mapped read with single modification type
    #[test]
    fn mapped_with_single_mod() -> Result<(), Error> {
        let mod_table_entry = ModTableEntryBuilder::default()
            .base('C')
            .is_strand_plus(true)
            .mod_code("m".into())
            .data([(0, 15, 200), (2, 25, 100)].into())
            .build()?;

        let read = CurrReadBuilder::default()
            .read_id("some_read".into())
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

        let df = curr_reads_to_dataframe(&[read])?;

        // Two modification data points
        assert_eq!(df.height(), 2);
        assert_eq!(df.width(), 13);

        // Verify read-level fields
        let read_id_col = df.column("read_id")?.str()?;
        assert_eq!(read_id_col.get(0), Some("some_read"));
        assert_eq!(read_id_col.get(1), Some("some_read"));

        let seq_len_col = df.column("seq_len")?.u64()?;
        assert_eq!(seq_len_col.get(0), Some(40));
        assert_eq!(seq_len_col.get(1), Some(40));

        let alignment_type_col = df.column("alignment_type")?.str()?;
        assert_eq!(alignment_type_col.get(0), Some("primary_forward"));
        assert_eq!(alignment_type_col.get(1), Some("primary_forward"));

        // Verify alignment fields
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

        // Verify modification entry fields
        let base_col = df.column("base")?.str()?;
        assert_eq!(base_col.get(0), Some("C"));
        assert_eq!(base_col.get(1), Some("C"));

        let is_strand_plus_col = df.column("is_strand_plus")?.bool()?;
        assert_eq!(is_strand_plus_col.get(0), Some(true));
        assert_eq!(is_strand_plus_col.get(1), Some(true));

        let mod_code_col = df.column("mod_code")?.str()?;
        assert_eq!(mod_code_col.get(0), Some("m"));
        assert_eq!(mod_code_col.get(1), Some("m"));

        // Verify modification data points
        let position_col = df.column("position")?.u64()?;
        assert_eq!(position_col.get(0), Some(0));
        assert_eq!(position_col.get(1), Some(2));

        let ref_position_col = df.column("ref_position")?.i64()?;
        assert_eq!(ref_position_col.get(0), Some(15));
        assert_eq!(ref_position_col.get(1), Some(25));

        let mod_quality_col = df.column("mod_quality")?.u32()?;
        assert_eq!(mod_quality_col.get(0), Some(200));
        assert_eq!(mod_quality_col.get(1), Some(100));

        Ok(())
    }

    /// Test mapped read with multiple modification types
    #[test]
    fn mapped_with_multiple_mods() -> Result<(), Error> {
        let mod_table_entry_1 = ModTableEntryBuilder::default()
            .base('C')
            .is_strand_plus(true)
            .mod_code("m".into())
            .data([(0, 15, 200), (2, 25, 100)].into())
            .build()?;

        let mod_table_entry_2 = ModTableEntryBuilder::default()
            .base('A')
            .is_strand_plus(true)
            .mod_code("a".into())
            .data([(1, 20, 50), (3, 30, 225)].into())
            .build()?;

        let read = CurrReadBuilder::default()
            .read_id("some_read".into())
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
            .mod_table([mod_table_entry_1, mod_table_entry_2].into())
            .build()?;

        let df = curr_reads_to_dataframe(&[read])?;

        // Four modification data points total (2 from each mod type)
        assert_eq!(df.height(), 4);
        assert_eq!(df.width(), 13);

        // Verify read-level fields (should be same for all rows)
        let read_id_col = df.column("read_id")?.str()?;
        for i in 0..4 {
            assert_eq!(read_id_col.get(i), Some("some_read"));
        }

        let seq_len_col = df.column("seq_len")?.u64()?;
        for i in 0..4 {
            assert_eq!(seq_len_col.get(i), Some(40));
        }

        let alignment_type_col = df.column("alignment_type")?.str()?;
        for i in 0..4 {
            assert_eq!(alignment_type_col.get(i), Some("primary_forward"));
        }

        // Verify alignment fields (should be same for all rows)
        let align_start_col = df.column("align_start")?.u64()?;
        for i in 0..4 {
            assert_eq!(align_start_col.get(i), Some(10));
        }

        let align_end_col = df.column("align_end")?.u64()?;
        for i in 0..4 {
            assert_eq!(align_end_col.get(i), Some(60));
        }

        let contig_col = df.column("contig")?.str()?;
        for i in 0..4 {
            assert_eq!(contig_col.get(i), Some("chr1"));
        }

        let contig_id_col = df.column("contig_id")?.i32()?;
        for i in 0..4 {
            assert_eq!(contig_id_col.get(i), Some(1));
        }

        // Verify modification data points
        let base_col = df.column("base")?.str()?;
        let mod_code_col = df.column("mod_code")?.str()?;
        let position_col = df.column("position")?.u64()?;
        let ref_position_col = df.column("ref_position")?.i64()?;
        let mod_quality_col = df.column("mod_quality")?.u32()?;

        let mut data_points = Vec::new();
        for i in 0..4 {
            data_points.push((
                base_col.get(i).unwrap(),
                mod_code_col.get(i).unwrap(),
                position_col.get(i).unwrap(),
                ref_position_col.get(i).unwrap(),
                mod_quality_col.get(i).unwrap(),
            ));
        }

        assert_eq!(data_points.first(), Some(&("A", "a", 1, 20, 50)));
        assert_eq!(data_points.get(1), Some(&("A", "a", 3, 30, 225)));
        assert_eq!(data_points.get(2), Some(&("C", "m", 0, 15, 200)));
        assert_eq!(data_points.get(3), Some(&("C", "m", 2, 25, 100)));

        Ok(())
    }

    /// Test unmapped read with multiple modification types
    #[test]
    fn unmapped_with_multiple_mods() -> Result<(), Error> {
        let mod_table_entry_1 = ModTableEntryBuilder::default()
            .base('C')
            .is_strand_plus(true)
            .mod_code("m".into())
            .data([(0, -1, 200), (2, -1, 100)].into())
            .build()?;

        let mod_table_entry_2 = ModTableEntryBuilder::default()
            .base('A')
            .is_strand_plus(true)
            .mod_code("a".into())
            .data([(1, -1, 50), (3, -1, 225)].into())
            .build()?;

        let read = CurrReadBuilder::default()
            .read_id("some_read".into())
            .seq_len(40)
            .mod_table([mod_table_entry_1, mod_table_entry_2].into())
            .build()?;

        let df = curr_reads_to_dataframe(&[read])?;

        // Four modification data points total
        assert_eq!(df.height(), 4);
        assert_eq!(df.width(), 13);

        // Verify read-level fields
        let read_id_col = df.column("read_id")?.str()?;
        for i in 0..4 {
            assert_eq!(read_id_col.get(i), Some("some_read"));
        }

        let seq_len_col = df.column("seq_len")?.u64()?;
        for i in 0..4 {
            assert_eq!(seq_len_col.get(i), Some(40));
        }

        let alignment_type_col = df.column("alignment_type")?.str()?;
        for i in 0..4 {
            assert_eq!(alignment_type_col.get(i), Some("unmapped"));
        }

        // Verify alignment fields are None for unmapped reads
        let align_start_col = df.column("align_start")?.u64()?;
        for i in 0..4 {
            assert_eq!(align_start_col.get(i), None);
        }

        let align_end_col = df.column("align_end")?.u64()?;
        for i in 0..4 {
            assert_eq!(align_end_col.get(i), None);
        }

        let contig_col = df.column("contig")?.str()?;
        for i in 0..4 {
            assert_eq!(contig_col.get(i), None);
        }

        let contig_id_col = df.column("contig_id")?.i32()?;
        for i in 0..4 {
            assert_eq!(contig_id_col.get(i), None);
        }

        // Verify modification data points
        let base_col = df.column("base")?.str()?;
        let mod_code_col = df.column("mod_code")?.str()?;
        let position_col = df.column("position")?.u64()?;
        let ref_position_col = df.column("ref_position")?.i64()?;
        let mod_quality_col = df.column("mod_quality")?.u32()?;

        // All ref_positions should be -1 for unmapped reads
        for i in 0..4 {
            assert_eq!(ref_position_col.get(i), Some(-1));
        }

        let mut data_points = Vec::new();
        for i in 0..4 {
            data_points.push((
                base_col.get(i).unwrap(),
                mod_code_col.get(i).unwrap(),
                position_col.get(i).unwrap(),
                mod_quality_col.get(i).unwrap(),
            ));
        }

        assert_eq!(data_points.first(), Some(&("A", "a", 1, 50)));
        assert_eq!(data_points.get(1), Some(&("A", "a", 3, 225)));
        assert_eq!(data_points.get(2), Some(&("C", "m", 0, 200)));
        assert_eq!(data_points.get(3), Some(&("C", "m", 2, 100)));

        Ok(())
    }
}
