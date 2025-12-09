//! Tests for `CurrReadBuilder` extracted from doctests
//! These tests verify the builder pattern functionality for creating `CurrRead` instances

use nanalogue_core::{
    AlignmentInfoBuilder, CurrRead, CurrReadBuilder, Error, ModTableEntryBuilder, ReadState,
    read_utils::AlignAndModData,
};

#[cfg(test)]
mod tests {
    use super::*;

    /// First example, unmapped read with very little information
    #[test]
    fn unmapped_minimal() -> Result<(), Error> {
        let _read: CurrRead<AlignAndModData> = CurrReadBuilder::default().build()?;
        Ok(())
    }

    /// Add some simple information, still unmapped
    #[test]
    fn unmapped_with_info() -> Result<(), Error> {
        let _read: CurrRead<AlignAndModData> = CurrReadBuilder::default()
            .read_id("some_read".into())
            .seq_len(40)
            .build()?;
        Ok(())
    }

    /// Mapped read without alignment information should panic
    #[test]
    #[should_panic(expected = "UnknownAlignState")]
    fn mapped_without_alignment() {
        let _: CurrRead<AlignAndModData> = CurrReadBuilder::default()
            .read_id("some_read".into())
            .seq_len(40)
            .alignment_type(ReadState::PrimaryFwd)
            .build()
            .unwrap();
    }

    /// Mapped read building
    #[test]
    fn mapped_basic() -> Result<(), Error> {
        let _read: CurrRead<AlignAndModData> = CurrReadBuilder::default()
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
        Ok(())
    }

    /// Mapped read building with modification information
    #[test]
    fn mapped_with_single_mod() -> Result<(), Error> {
        let mod_table_entry = ModTableEntryBuilder::default()
            .base('C')
            .is_strand_plus(true)
            .mod_code("m".into())
            .data([(0, 15, 200), (2, 25, 100)])
            .build()?;

        let _read: CurrRead<AlignAndModData> = CurrReadBuilder::default()
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
        Ok(())
    }

    /// Mapped read building with multiple mod types
    #[test]
    fn mapped_with_multiple_mods() -> Result<(), Error> {
        let mod_table_entry_1 = ModTableEntryBuilder::default()
            .base('C')
            .is_strand_plus(true)
            .mod_code("m".into())
            .data([(0, 15, 200), (2, 25, 100)])
            .build()?;

        let mod_table_entry_2 = ModTableEntryBuilder::default()
            .base('A')
            .is_strand_plus(true)
            .mod_code("a".into())
            .data([(1, 20, 50), (3, 30, 225)])
            .build()?;

        let _read: CurrRead<AlignAndModData> = CurrReadBuilder::default()
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
        Ok(())
    }

    /// Unmapped read building with multiple mod types
    #[test]
    fn unmapped_with_multiple_mods() -> Result<(), Error> {
        let mod_table_entry_1 = ModTableEntryBuilder::default()
            .base('C')
            .is_strand_plus(true)
            .mod_code("m".into())
            .data([(0, -1, 200), (2, -1, 100)])
            .build()?;

        let mod_table_entry_2 = ModTableEntryBuilder::default()
            .base('A')
            .is_strand_plus(true)
            .mod_code("a".into())
            .data([(1, -1, 50), (3, -1, 225)])
            .build()?;

        let _read: CurrRead<AlignAndModData> = CurrReadBuilder::default()
            .read_id("some_read".into())
            .seq_len(40)
            .mod_table([mod_table_entry_1, mod_table_entry_2].into())
            .build()?;
        Ok(())
    }

    /// Coordinates are not sorted - should panic
    #[test]
    #[should_panic(expected = "InvalidModCoords")]
    fn unsorted_coords_forward() {
        let mod_table_entry = ModTableEntryBuilder::default()
            .base('C')
            .is_strand_plus(true)
            .mod_code("m".into())
            .data([(2, 25, 100), (0, 15, 200)])
            .build()
            .unwrap();

        let read_before_build = CurrReadBuilder::default()
            .read_id("some_read".into())
            .seq_len(40)
            .alignment(
                AlignmentInfoBuilder::default()
                    .start(10)
                    .end(60)
                    .contig("chr1".into())
                    .contig_id(1)
                    .build()
                    .unwrap(),
            )
            .mod_table([mod_table_entry].into());

        let _: CurrRead<AlignAndModData> = read_before_build
            .alignment_type(ReadState::PrimaryFwd)
            .build()
            .unwrap();
    }

    /// Coordinates are not sorted (reverse strand) - should panic
    #[test]
    #[should_panic(expected = "InvalidModCoords")]
    fn unsorted_coords_reverse() {
        let mod_table_entry = ModTableEntryBuilder::default()
            .base('C')
            .is_strand_plus(true)
            .mod_code("m".into())
            .data([(2, 25, 100), (0, 15, 200)])
            .build()
            .unwrap();

        let read_before_build = CurrReadBuilder::default()
            .read_id("some_read".into())
            .seq_len(40)
            .alignment(
                AlignmentInfoBuilder::default()
                    .start(10)
                    .end(60)
                    .contig("chr1".into())
                    .contig_id(1)
                    .build()
                    .unwrap(),
            )
            .mod_table([mod_table_entry].into());

        let _: CurrRead<AlignAndModData> = read_before_build
            .alignment_type(ReadState::PrimaryRev)
            .build()
            .unwrap();
    }

    /// Unmapped read but mod reference coordinates are not set to -1 - should panic
    #[test]
    #[should_panic(expected = "InvalidAlignCoords")]
    fn unmapped_with_invalid_ref_coords() {
        let mod_table_entry = ModTableEntryBuilder::default()
            .base('C')
            .is_strand_plus(true)
            .mod_code("m".into())
            .data([(0, -1, 200), (2, 11, 100)])
            .build()
            .unwrap();

        let _: CurrRead<AlignAndModData> = CurrReadBuilder::default()
            .read_id("some_read".into())
            .seq_len(40)
            .mod_table([mod_table_entry].into())
            .build()
            .unwrap();
    }

    /// Read with mod coordinates larger than sequence length (unmapped) - should panic
    #[test]
    #[should_panic(expected = "InvalidModCoords")]
    fn mod_coords_exceed_seq_len_unmapped() {
        let mod_table_entry = ModTableEntryBuilder::default()
            .base('C')
            .is_strand_plus(true)
            .mod_code("m".into())
            .data([(0, -1, 200), (42, -1, 100)])
            .build()
            .unwrap();

        let _: CurrRead<AlignAndModData> = CurrReadBuilder::default()
            .read_id("some_read".into())
            .seq_len(40)
            .mod_table([mod_table_entry].into())
            .build()
            .unwrap();
    }

    /// Read with mod coordinates larger than sequence length (mapped) - should panic
    #[test]
    #[should_panic(expected = "InvalidModCoords")]
    fn mod_coords_exceed_seq_len_mapped() {
        let mod_table_entry = ModTableEntryBuilder::default()
            .base('C')
            .is_strand_plus(true)
            .mod_code("m".into())
            .data([(0, 20, 200), (42, 30, 100)])
            .build()
            .unwrap();

        let _: CurrRead<AlignAndModData> = CurrReadBuilder::default()
            .read_id("some_read".into())
            .seq_len(40)
            .alignment_type(ReadState::PrimaryFwd)
            .alignment(
                AlignmentInfoBuilder::default()
                    .start(10)
                    .end(60)
                    .contig("chr1".into())
                    .contig_id(1)
                    .build()
                    .unwrap(),
            )
            .mod_table([mod_table_entry].into())
            .build()
            .unwrap();
    }

    /// Read with mod coordinates beyond alignment coordinates - should panic
    #[test]
    #[should_panic(expected = "InvalidAlignCoords")]
    fn mod_ref_coords_beyond_alignment() {
        let mod_table_entry = ModTableEntryBuilder::default()
            .base('C')
            .is_strand_plus(true)
            .mod_code("m".into())
            .data([(0, 1, 200), (22, 30, 100)])
            .build()
            .unwrap();

        let _: CurrRead<AlignAndModData> = CurrReadBuilder::default()
            .read_id("some_read".into())
            .seq_len(40)
            .alignment_type(ReadState::PrimaryFwd)
            .alignment(
                AlignmentInfoBuilder::default()
                    .start(10)
                    .end(60)
                    .contig("chr1".into())
                    .contig_id(1)
                    .build()
                    .unwrap(),
            )
            .mod_table([mod_table_entry].into())
            .build()
            .unwrap();
    }

    /// Reads with multiple tracks of the same kind of modification - should panic
    #[test]
    #[should_panic(expected = "InvalidDuplicates")]
    fn duplicate_modification_tracks() {
        let mod_table_entry_1 = ModTableEntryBuilder::default()
            .base('C')
            .is_strand_plus(true)
            .mod_code("m".into())
            .data([(0, -1, 200), (2, -1, 100)])
            .build()
            .unwrap();

        let mod_table_entry_2 = ModTableEntryBuilder::default()
            .base('C')
            .is_strand_plus(true)
            .mod_code("m".into())
            .data([(1, -1, 50), (3, -1, 225)])
            .build()
            .unwrap();

        let _: CurrRead<AlignAndModData> = CurrReadBuilder::default()
            .read_id("some_read".into())
            .seq_len(40)
            .mod_table([mod_table_entry_1, mod_table_entry_2].into())
            .build()
            .unwrap();
    }

    /// Alignment with end < start - should panic with `InvalidAlignCoords`
    #[test]
    #[should_panic(expected = "InvalidAlignCoords")]
    fn alignment_end_before_start() {
        let _: CurrRead<AlignAndModData> = CurrReadBuilder::default()
            .read_id("invalid_align_read".into())
            .seq_len(40)
            .alignment_type(ReadState::PrimaryFwd)
            .alignment(
                AlignmentInfoBuilder::default()
                    .start(60)
                    .end(10)
                    .contig("chr1".into())
                    .contig_id(1)
                    .build()
                    .unwrap(),
            )
            .build()
            .unwrap();
    }
}
