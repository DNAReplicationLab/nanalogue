//! Tests for `set_mod_data_restricted_options` method in `read_utils.rs`
//! Covers filtering by region, tag, strand, probability, base quality, and read end trimming

use bedrs::Bed3;
use nanalogue_core::{
    CurrRead, Error, InputModOptions, InputRegionOptions, ModChar, RestrictModCalledStrand,
    ThresholdState, curr_reads_to_dataframe, nanalogue_bam_reader,
};
use rust_htslib::bam::Read as _;
use std::str::FromStr as _;

/// Mock struct implementing both `InputModOptions` and `InputRegionOptions` for testing
#[derive(Debug, Clone)]
struct MockModOptions {
    tag: Option<ModChar>,
    mod_strand: Option<RestrictModCalledStrand>,
    mod_prob_filter: ThresholdState,
    trim_read_ends_mod: usize,
    base_qual_filter_mod: u8,
    region_filter: Option<Bed3<i32, u64>>,
}

impl MockModOptions {
    fn new() -> Self {
        Self {
            tag: None,
            mod_strand: None,
            mod_prob_filter: ThresholdState::GtEq(0),
            trim_read_ends_mod: 0,
            base_qual_filter_mod: 0,
            region_filter: None,
        }
    }

    fn with_tag(mut self, tag: ModChar) -> Self {
        self.tag = Some(tag);
        self
    }

    fn with_mod_strand(mut self, strand: RestrictModCalledStrand) -> Self {
        self.mod_strand = Some(strand);
        self
    }

    fn with_mod_prob_filter(mut self, threshold: ThresholdState) -> Self {
        self.mod_prob_filter = threshold;
        self
    }

    fn with_trim_read_ends(mut self, trim: usize) -> Self {
        self.trim_read_ends_mod = trim;
        self
    }

    fn with_base_qual_filter(mut self, qual: u8) -> Self {
        self.base_qual_filter_mod = qual;
        self
    }

    fn with_region_filter(mut self, region: Bed3<i32, u64>) -> Self {
        self.region_filter = Some(region);
        self
    }
}

impl InputModOptions for MockModOptions {
    fn tag(&self) -> Option<ModChar> {
        self.tag
    }

    fn mod_strand(&self) -> Option<RestrictModCalledStrand> {
        self.mod_strand
    }

    fn mod_prob_filter(&self) -> ThresholdState {
        self.mod_prob_filter
    }

    fn trim_read_ends_mod(&self) -> usize {
        self.trim_read_ends_mod
    }

    fn base_qual_filter_mod(&self) -> u8 {
        self.base_qual_filter_mod
    }
}

impl InputRegionOptions for MockModOptions {
    fn region_filter(&self) -> &Option<Bed3<i32, u64>> {
        &self.region_filter
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_filters_applied() -> Result<(), Error> {
        // Test: Basic case with no filters - all modifications should be retained
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let options = MockModOptions::new();

        // Use the second record which has 3 T modifications
        let record = reader.records().nth(1).unwrap()?;
        let curr_read = CurrRead::default().try_from_only_alignment(&record)?;
        let result = curr_read.set_mod_data_restricted_options(&record, &options)?;

        // Should have modifications (example_1.bam second record has 3 T mods)
        let mod_count = result.base_count_per_mod();
        assert!(
            mod_count.get(&ModChar::new('T')).copied().unwrap_or(0) > 0,
            "Expected modifications to be present"
        );

        Ok(())
    }

    #[test]
    fn high_probability_filter() -> Result<(), Error> {
        // Test: Probability threshold filtering with explicit position checking
        // Second record has T mods at positions with qualities: [3,26,221], [8,31,242], [27,50,3], [39,62,47], [47,70,239]
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let record = reader.records().nth(1).unwrap()?;

        // Test with threshold 128: should keep mods with qual >= 128 (positions 3, 8, 47)
        let curr_read_low = CurrRead::default().try_from_only_alignment(&record)?;
        let options_low = MockModOptions::new().with_mod_prob_filter(ThresholdState::GtEq(128));
        let result_low = curr_read_low.set_mod_data_restricted_options(&record, &options_low)?;

        // Convert to DataFrame to check exact positions
        let df_low = curr_reads_to_dataframe(&[result_low])?;

        // Should have 3 mods (qual 221, 242, 239 are all >= 128)
        assert_eq!(df_low.height(), 3, "Should have 3 mods with qual >= 128");

        // Check that all remaining mods have quality >= 128
        let qual_col = df_low.column("mod_quality")?.u32()?;
        for i in 0..df_low.height() {
            let qual = qual_col.get(i).unwrap();
            assert!(
                qual >= 128,
                "All mods should have quality >= 128, found {qual}"
            );
        }

        // Check specific positions that should be present (read positions 3, 8, 47)
        let pos_col = df_low.column("position")?.u64()?;
        let positions: Vec<u64> = (0..df_low.height())
            .filter_map(|i| pos_col.get(i))
            .collect();
        assert!(
            positions.contains(&3),
            "Position 3 (qual 221) should be present"
        );
        assert!(
            positions.contains(&8),
            "Position 8 (qual 242) should be present"
        );
        assert!(
            positions.contains(&47),
            "Position 47 (qual 239) should be present"
        );

        // Test with threshold 255: should filter out all mods (none have qual >= 255)
        let curr_read_high = CurrRead::default().try_from_only_alignment(&record)?;
        let options_high = MockModOptions::new().with_mod_prob_filter(ThresholdState::GtEq(255));
        let result_high = curr_read_high.set_mod_data_restricted_options(&record, &options_high)?;

        let df_high = curr_reads_to_dataframe(&[result_high])?;

        // Should have 0 mods
        assert_eq!(df_high.height(), 0, "No mods should have quality >= 255");

        Ok(())
    }

    #[test]
    fn region_filter_full_overlap() -> Result<(), Error> {
        // Test: Region filter fully contains the read - no filtering should occur
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let record = reader.records().nth(1).unwrap()?;

        // Second record in example_1.bam is on contig 2, positions 23-71
        // Get baseline count without region filter
        let curr_read_baseline = CurrRead::default().try_from_only_alignment(&record)?;
        let options_baseline = MockModOptions::new();
        let result_baseline =
            curr_read_baseline.set_mod_data_restricted_options(&record, &options_baseline)?;
        let baseline_count = result_baseline
            .base_count_per_mod()
            .get(&ModChar::new('T'))
            .copied()
            .unwrap_or(0);

        // Create a region that fully contains the read
        let curr_read = CurrRead::default().try_from_only_alignment(&record)?;
        let region = Bed3::new(2, 0, 100); // Fully contains the read
        let options = MockModOptions::new().with_region_filter(region);

        let result = curr_read.set_mod_data_restricted_options(&record, &options)?;

        // Should retain all modifications since region fully contains the read
        let mod_count = result
            .base_count_per_mod()
            .get(&ModChar::new('T'))
            .copied()
            .unwrap_or(0);
        assert_eq!(
            mod_count, baseline_count,
            "All modifications should be retained with full overlap"
        );

        Ok(())
    }

    #[test]
    fn region_filter_partial_overlap() -> Result<(), Error> {
        // Test: Partial overlap - only mods in overlapping region should be retained
        // Second record: contig 2, align 23-71, mods at ref positions [26, 31, 50, 62, 70]
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let record = reader.records().nth(1).unwrap()?;

        // Create a region that partially overlaps (region 50-100 overlaps with read 23-71)
        // Only mods at ref positions >= 50 should be kept (positions 50, 62, 70)
        let curr_read = CurrRead::default().try_from_only_alignment(&record)?;
        let region = Bed3::new(2, 50, 100); // Partial overlap on right side
        let options = MockModOptions::new().with_region_filter(region);

        let result = curr_read.set_mod_data_restricted_options(&record, &options)?;
        let df = curr_reads_to_dataframe(&[result])?;

        // Should have exactly 3 mods (at ref positions 50, 62, 70)
        assert_eq!(df.height(), 3, "Should have 3 mods in region 50-100");

        // Check that all ref_positions are within the region [50, 100)
        let ref_pos_col = df.column("ref_position")?.i64()?;
        for i in 0..df.height() {
            let ref_pos = ref_pos_col.get(i).unwrap();
            assert!(
                (50..100).contains(&ref_pos),
                "Ref position {ref_pos} should be in range [50, 100)"
            );
        }

        // Check specific ref positions that should be present
        let ref_positions: Vec<i64> = (0..df.height())
            .filter_map(|i| ref_pos_col.get(i))
            .collect();
        assert!(
            ref_positions.contains(&50),
            "Ref position 50 should be present"
        );
        assert!(
            ref_positions.contains(&62),
            "Ref position 62 should be present"
        );
        assert!(
            ref_positions.contains(&70),
            "Ref position 70 should be present"
        );

        // Verify positions NOT in region are excluded
        assert!(
            !ref_positions.contains(&26),
            "Ref position 26 should be excluded"
        );
        assert!(
            !ref_positions.contains(&31),
            "Ref position 31 should be excluded"
        );

        Ok(())
    }

    #[test]
    fn region_filter_no_overlap() -> Result<(), Error> {
        // Test: No overlap - no modifications should be retained
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let record = reader.records().nth(1).unwrap()?;

        // Second record in example_1.bam is on contig 2, positions 23-71
        // Create a region on a different contig (no overlap)
        let curr_read = CurrRead::default().try_from_only_alignment(&record)?;
        let region = Bed3::new(0, 0, 100); // Different contig, no overlap
        let options = MockModOptions::new().with_region_filter(region);

        let result = curr_read.set_mod_data_restricted_options(&record, &options)?;

        // Should have no modifications
        let mod_count = result.base_count_per_mod();
        let total_mods: u32 = mod_count.values().sum();
        assert_eq!(
            total_mods, 0,
            "No overlap should result in zero modifications"
        );

        Ok(())
    }

    #[test]
    fn region_filter_exact_boundaries() -> Result<(), Error> {
        // Test: Region exactly matches read boundaries
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let record = reader.records().nth(1).unwrap()?;

        // Second record in example_1.bam is on contig 2, positions 23-71
        // Get baseline count without region filter
        let curr_read_baseline = CurrRead::default().try_from_only_alignment(&record)?;
        let options_baseline = MockModOptions::new();
        let result_baseline =
            curr_read_baseline.set_mod_data_restricted_options(&record, &options_baseline)?;
        let baseline_count = result_baseline
            .base_count_per_mod()
            .get(&ModChar::new('T'))
            .copied()
            .unwrap_or(0);

        let curr_read = CurrRead::default().try_from_only_alignment(&record)?;
        let region = Bed3::new(2, 23, 71); // Exact match
        let options = MockModOptions::new().with_region_filter(region);

        let result = curr_read.set_mod_data_restricted_options(&record, &options)?;

        // Should retain all modifications
        let mod_count = result
            .base_count_per_mod()
            .get(&ModChar::new('T'))
            .copied()
            .unwrap_or(0);
        assert_eq!(
            mod_count, baseline_count,
            "Exact boundary match should retain all modifications"
        );

        Ok(())
    }

    #[test]
    fn trim_read_ends() -> Result<(), Error> {
        // Test: Trimming read ends should exclude modifications near ends
        // Second record: seq_len=48, mods at read positions [3, 8, 27, 39, 47]
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let record = reader.records().nth(1).unwrap()?;

        // Test with trimming 10bp from each end (keeps positions 10-37 in a 48bp read)
        let curr_read_trim = CurrRead::default().try_from_only_alignment(&record)?;
        let options_trim = MockModOptions::new().with_trim_read_ends(10);
        let result_trim = curr_read_trim.set_mod_data_restricted_options(&record, &options_trim)?;
        let df_trim = curr_reads_to_dataframe(&[result_trim])?;

        // Should only keep position 27 (positions 3, 8 are < 10; positions 39, 47 are >= 38)
        assert_eq!(
            df_trim.height(),
            1,
            "Should have 1 mod after trimming 10bp from each end"
        );

        // Check that remaining position is 27
        let pos_col = df_trim.column("position")?.u64()?;
        let position = pos_col.get(0).unwrap();
        assert_eq!(
            position, 27,
            "Only position 27 should remain after trimming"
        );

        Ok(())
    }

    #[test]
    fn base_quality_filter() -> Result<(), Error> {
        // Test: Base quality filtering
        // NOTE that we use example_5 here that has valid basecalling qualities.
        let mut reader = nanalogue_bam_reader("examples/example_5_valid_basequal.sam")?;
        let record = reader.records().next().unwrap()?;

        // Get count with no base quality filter
        let curr_read_no_filter = CurrRead::default().try_from_only_alignment(&record)?;
        let options_no_filter = MockModOptions::new().with_base_qual_filter(0);
        let result_no_filter =
            curr_read_no_filter.set_mod_data_restricted_options(&record, &options_no_filter)?;
        let count_no_filter = result_no_filter.base_count_per_mod().values().sum::<u32>();

        // Get count with high base quality filter
        let curr_read_high_filter = CurrRead::default().try_from_only_alignment(&record)?;
        let options_high_filter = MockModOptions::new().with_base_qual_filter(40);
        let result_high_filter =
            curr_read_high_filter.set_mod_data_restricted_options(&record, &options_high_filter)?;
        let count_high_filter = result_high_filter
            .base_count_per_mod()
            .values()
            .sum::<u32>();

        // Higher base quality threshold should filter out more modifications
        assert!(
            count_high_filter < count_no_filter,
            "Higher base quality threshold should result in fewer or equal mods"
        );

        Ok(())
    }

    #[test]
    fn base_quality_filter_2() -> Result<(), Error> {
        // Test: Base quality filtering
        // NOTE that we use example_8, example_9 here that has valid basecalling qualities on one
        // read and on an identical read on the next line has no basecalling qualities
        for file_name in ["examples/example_8.sam", "examples/example_9.sam"] {
            let mut reader = nanalogue_bam_reader(file_name)?;
            let record_1 = reader.records().next().unwrap()?;
            let record_2 = reader.records().next().unwrap()?;

            // Get count with no base quality filter
            let curr_read_no_filter_1 = CurrRead::default().try_from_only_alignment(&record_1)?;
            let options_no_filter = MockModOptions::new().with_base_qual_filter(0);
            let result_no_filter_1 = curr_read_no_filter_1
                .set_mod_data_restricted_options(&record_1, &options_no_filter)?;
            let count_no_filter_1 = result_no_filter_1
                .base_count_per_mod()
                .values()
                .sum::<u32>();

            // Get count with high base quality filter
            let curr_read_high_filter_1 = CurrRead::default().try_from_only_alignment(&record_1)?;
            let options_high_filter = MockModOptions::new().with_base_qual_filter(40);
            let result_high_filter_1 = curr_read_high_filter_1
                .set_mod_data_restricted_options(&record_1, &options_high_filter)?;
            let count_high_filter_1 = result_high_filter_1
                .base_count_per_mod()
                .values()
                .sum::<u32>();

            // Higher base quality threshold should filter out more modifications
            assert!(
                count_no_filter_1 == 6 && count_high_filter_1 == 3,
                "Higher base quality threshold should result in fewer or equal mods"
            );

            // Get count with no base quality filter
            let curr_read_no_filter_2 = CurrRead::default().try_from_only_alignment(&record_2)?;
            let result_no_filter_2 = curr_read_no_filter_2
                .set_mod_data_restricted_options(&record_2, &options_no_filter)?;
            let count_no_filter_2 = result_no_filter_2
                .base_count_per_mod()
                .values()
                .sum::<u32>();

            // Get count with high base quality filter
            let curr_read_high_filter_2 = CurrRead::default().try_from_only_alignment(&record_2)?;
            let result_high_filter_2 = curr_read_high_filter_2
                .set_mod_data_restricted_options(&record_2, &options_high_filter)?;
            let count_high_filter_2 = result_high_filter_2
                .base_count_per_mod()
                .values()
                .sum::<u32>();

            // Higher base quality threshold should filter out more modifications
            assert!(
                count_no_filter_2 == 6 && count_high_filter_2 == 0,
                "Another test of higher base quality threshold should result in fewer or equal mods"
            );
        }

        Ok(())
    }

    #[test]
    fn base_implicit_reading_test() -> Result<(), Error> {
        // Test: Implicit base quality reading
        // NOTE that we use example_8, example_9 here that has implicit mod quals on one
        // read and on an identical read on the third line but mods are marked explicit here
        for file_name in ["examples/example_8.sam", "examples/example_9.sam"] {
            let mut reader = nanalogue_bam_reader(file_name)?;
            let record_1 = reader.records().next().unwrap()?;
            let _record_2 = reader.records().next().unwrap()?;
            let record_3 = reader.records().next().unwrap()?;

            // Get count of two records
            let curr_read_no_filter_1 = CurrRead::default().try_from_only_alignment(&record_1)?;
            let options_no_filter = MockModOptions::new();
            let result_no_filter_1 = curr_read_no_filter_1
                .set_mod_data_restricted_options(&record_1, &options_no_filter)?;
            let count_no_filter_1 = result_no_filter_1
                .base_count_per_mod()
                .values()
                .sum::<u32>();

            let curr_read_no_filter_3 = CurrRead::default().try_from_only_alignment(&record_3)?;
            let result_no_filter_3 = curr_read_no_filter_3
                .set_mod_data_restricted_options(&record_3, &options_no_filter)?;
            let count_no_filter_3 = result_no_filter_3
                .base_count_per_mod()
                .values()
                .sum::<u32>();

            // Check that implicit bases are read correctly
            assert!(
                count_no_filter_1 == 6 && count_no_filter_3 == 4,
                "Test of implicit base counts"
            );
        }

        Ok(())
    }

    #[test]
    fn base_quality_filter_in_file_with_no_base_qual() -> Result<(), Error> {
        // Test: Base quality filtering
        // NOTE that we use example_1 here that has no basecalling qualities.
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let record = reader.records().nth(1).unwrap()?;

        // Get count with no base quality filter
        let curr_read_no_filter = CurrRead::default().try_from_only_alignment(&record)?;
        let options_no_filter = MockModOptions::new().with_base_qual_filter(0);
        let result_no_filter =
            curr_read_no_filter.set_mod_data_restricted_options(&record, &options_no_filter)?;
        let count_no_filter = result_no_filter.base_count_per_mod().values().sum::<u32>();

        // Get count with high base quality filter, should throw out all positions
        // as this read does not have base quality information
        let curr_read_high_filter = CurrRead::default().try_from_only_alignment(&record)?;
        let options_high_filter = MockModOptions::new().with_base_qual_filter(94);
        let result_high_filter =
            curr_read_high_filter.set_mod_data_restricted_options(&record, &options_high_filter)?;
        let count_high_filter = result_high_filter
            .base_count_per_mod()
            .values()
            .sum::<u32>();

        // These numbers should be equal, and the high filter count equal to zero.
        assert!(
            count_no_filter > 0 && count_high_filter == 0,
            "Base quality threshold should reject all calls as there are no base qualities in this file"
        );

        Ok(())
    }

    #[test]
    fn tag_filter_specific_tag() -> Result<(), Error> {
        // Test: Filter to specific modification tag
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;

        // Fourth record has both T and another mod type
        let record = reader.records().nth(3).unwrap()?;

        let curr_read = CurrRead::default().try_from_only_alignment(&record)?;
        let options = MockModOptions::new().with_tag(ModChar::new('T'));

        let result = curr_read.set_mod_data_restricted_options(&record, &options)?;

        // Should only have T modifications
        let mod_count = result.base_count_per_mod();

        // Check all tags are 'T' - collect to Vec first to avoid iteration over hash type
        let tags: Vec<_> = mod_count.keys().copied().collect();
        for tag in tags {
            assert!(
                tag == ModChar::new('T'),
                "Should only contain T modifications when filtered by tag"
            );
        }

        Ok(())
    }

    #[test]
    fn tag_filter_none_accepts_all() -> Result<(), Error> {
        // Test: tag=None should accept all tags
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let record = reader.records().nth(3).unwrap()?;

        let curr_read = CurrRead::default().try_from_only_alignment(&record)?;
        let options = MockModOptions::new(); // tag is None by default

        let result = curr_read.set_mod_data_restricted_options(&record, &options)?;

        // Should have modifications of potentially multiple types
        let mod_count = result.base_count_per_mod();
        assert!(
            !mod_count.is_empty(),
            "Should have modifications when no tag filter is applied"
        );

        Ok(())
    }

    #[test]
    fn combined_filters() -> Result<(), Error> {
        // Test: Multiple filters working together
        // Second record: mods at [3,26,221], [8,31,242], [27,50,3], [39,62,47], [47,70,239]
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let record = reader.records().nth(1).unwrap()?;

        // Apply multiple filters: trimming + probability + region
        let curr_read = CurrRead::default().try_from_only_alignment(&record)?;
        let region = Bed3::new(2, 50, 71); // Region [50, 71)
        let options = MockModOptions::new()
            .with_trim_read_ends(5) // Keeps read positions [5, 43)
            .with_mod_prob_filter(ThresholdState::GtEq(200)) // Keeps qual >= 200
            .with_region_filter(region); // Keeps ref positions [50, 71)

        let result = curr_read.set_mod_data_restricted_options(&record, &options)?;
        let df = curr_reads_to_dataframe(&[result])?;

        // Expected: NO mods pass all three filters
        // - Position 8 (ref 31, qual 242): passes trim & qual, but fails region (31 < 50)
        // - Position 27 (ref 50, qual 3): passes trim & region, but fails qual (3 < 200)
        // - Position 39 (ref 62, qual 47): passes trim & region, but fails qual (47 < 200)
        // - Positions 3 and 47 are trimmed
        assert_eq!(df.height(), 0, "No mods should pass all three filters");

        Ok(())
    }

    #[test]
    fn empty_interval_removes_all_mods() -> Result<(), Error> {
        // Test: Empty interval (0..0) should remove all modifications
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let record = reader.records().nth(1).unwrap()?;

        // Create a scenario with no overlap (different contig) which creates 0..0 interval
        let curr_read = CurrRead::default().try_from_only_alignment(&record)?;
        let region = Bed3::new(1, 0, 10); // Different contig than record's contig 2
        let options = MockModOptions::new().with_region_filter(region);

        let result = curr_read.set_mod_data_restricted_options(&record, &options)?;

        // Should have zero modifications
        let total_mods: u32 = result.base_count_per_mod().values().sum();
        assert_eq!(
            total_mods, 0,
            "Empty interval should result in zero modifications"
        );

        Ok(())
    }

    #[expect(
        clippy::shadow_unrelated,
        reason = "code blocks are separate enough, so we won't get confused"
    )]
    #[test]
    fn mod_strand_filter() -> Result<(), Error> {
        // Test: Mod strand filtering
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let record = reader.records().nth(3).unwrap()?;

        // Test with basecalled strand filter
        let curr_read = CurrRead::default().try_from_only_alignment(&record)?;
        let options =
            MockModOptions::new().with_mod_strand(RestrictModCalledStrand::from_str("bc")?);

        let result = curr_read.set_mod_data_restricted_options(&record, &options)?;

        // Check all tags are 'T' - collect to Vec first to avoid iteration over hash type
        let mod_count = result.base_count_per_mod();
        let tags: Vec<_> = mod_count.keys().copied().collect();
        for tag in tags {
            assert!(
                tag == ModChar::new('T'),
                "Should only contain T modifications when filtered by bc"
            );
        }

        // Test with basecalled complementary strand filter
        let curr_read = CurrRead::default().try_from_only_alignment(&record)?;
        let options =
            MockModOptions::new().with_mod_strand(RestrictModCalledStrand::from_str("bc_comp")?);

        let result = curr_read.set_mod_data_restricted_options(&record, &options)?;

        // Check all tags are 7200 - collect to Vec first to avoid iteration over hash type
        let mod_count = result.base_count_per_mod();
        let tags: Vec<_> = mod_count.keys().copied().collect();
        for tag in tags {
            assert!(
                tag == ModChar::from_str("7200")?,
                "Should only contain modification of code 7200 when filtered by bc_comp"
            );
        }

        Ok(())
    }

    #[test]
    fn aggressive_trimming_removes_all() -> Result<(), Error> {
        // Additional test: Very aggressive trimming should remove most/all mods
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let record = reader.records().nth(1).unwrap()?;

        // Record has seq_len of 48, trim 24 from each end leaves nothing
        let curr_read = CurrRead::default().try_from_only_alignment(&record)?;
        let options = MockModOptions::new().with_trim_read_ends(24);

        let result = curr_read.set_mod_data_restricted_options(&record, &options)?;

        // Should have zero or very few modifications
        let total_mods: u32 = result.base_count_per_mod().values().sum();
        assert_eq!(total_mods, 0, "Aggressive trimming should remove all mods");

        Ok(())
    }

    #[test]
    fn unmapped_read_handling() -> Result<(), Error> {
        // Additional test: Verify behavior with unmapped reads
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;

        // Fourth record (index 3) is unmapped
        let record = reader.records().nth(3).unwrap()?;
        let curr_read_result = CurrRead::default().try_from_only_alignment(&record);

        // Should successfully create CurrRead for unmapped read
        assert!(curr_read_result.is_ok(), "Should handle unmapped reads");

        Ok(())
    }
}
