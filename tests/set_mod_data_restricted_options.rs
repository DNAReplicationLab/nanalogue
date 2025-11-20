//! Tests for `set_mod_data_restricted_options` method in `read_utils.rs`
//! Covers filtering by region, tag, strand, probability, base quality, and read end trimming

use bedrs::Bed3;
use nanalogue_core::{
    CurrRead, Error, InputModOptions, InputRegionOptions, ModChar, RestrictModCalledStrand,
    ThresholdState, nanalogue_bam_reader,
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
        // Test: Probability threshold filtering
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;

        // Use second record which has 3 T mods with threshold 128
        let record = reader.records().nth(1).unwrap()?;

        // First, get count with low threshold
        let curr_read_low = CurrRead::default().try_from_only_alignment(&record)?;
        let options_low = MockModOptions::new().with_mod_prob_filter(ThresholdState::GtEq(128));
        let result_low = curr_read_low.set_mod_data_restricted_options(&record, &options_low)?;
        let count_low = result_low
            .base_count_per_mod()
            .get(&ModChar::new('T'))
            .copied()
            .unwrap_or(0);

        // Then get count with high threshold
        let curr_read_high = CurrRead::default().try_from_only_alignment(&record)?;
        let options_high = MockModOptions::new().with_mod_prob_filter(ThresholdState::GtEq(255));
        let result_high = curr_read_high.set_mod_data_restricted_options(&record, &options_high)?;
        let count_high = result_high
            .base_count_per_mod()
            .get(&ModChar::new('T'))
            .copied()
            .unwrap_or(0);

        // Higher threshold should filter out more modifications
        assert!(
            count_high < count_low && count_high == 0,
            "Higher probability threshold should result in fewer or equal mods"
        );

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
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let record = reader.records().nth(1).unwrap()?;

        // Second record in example_1.bam is on contig 2, positions 23-71
        // Get baseline count without region filter
        let curr_read_baseline = CurrRead::default().try_from_only_alignment(&record)?;
        let options_baseline = MockModOptions::new();
        let result_baseline =
            curr_read_baseline.set_mod_data_restricted_options(&record, &options_baseline)?;
        let baseline_count = result_baseline.base_count_per_mod().values().sum::<u32>();

        // Second record in example_1.bam is on contig 2, positions 23-71
        // Create a region that partially overlaps (only covers part of the read)
        let curr_read = CurrRead::default().try_from_only_alignment(&record)?;
        let region = Bed3::new(2, 50, 100); // Partial overlap on right side
        let options = MockModOptions::new().with_region_filter(region);

        let result = curr_read.set_mod_data_restricted_options(&record, &options)?;

        // With partial overlap, we should have fewer modifications
        // (exact count depends on where mods are located in the read)
        let region_count = result.base_count_per_mod().values().sum::<u32>();
        assert!(
            region_count < baseline_count,
            "Partial overlap should retain at most all original mods"
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
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let record = reader.records().nth(1).unwrap()?;

        // Test with no trimming first
        let curr_read_no_trim = CurrRead::default().try_from_only_alignment(&record)?;
        let options_no_trim = MockModOptions::new();
        let result_no_trim =
            curr_read_no_trim.set_mod_data_restricted_options(&record, &options_no_trim)?;
        let count_no_trim = result_no_trim.base_count_per_mod().values().sum::<u32>();

        // Test with trimming (trim 10bp from each end of a 48bp read)
        let curr_read_trim = CurrRead::default().try_from_only_alignment(&record)?;
        let options_trim = MockModOptions::new().with_trim_read_ends(10);
        let result_trim = curr_read_trim.set_mod_data_restricted_options(&record, &options_trim)?;
        let count_trim = result_trim.base_count_per_mod().values().sum::<u32>();

        // Trimming should result in fewer or equal modifications
        assert!(
            count_trim <= count_no_trim,
            "Trimming should result in fewer or equal modifications"
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

        // Get count with high base quality filter
        let curr_read_high_filter = CurrRead::default().try_from_only_alignment(&record)?;
        let options_high_filter = MockModOptions::new().with_base_qual_filter(94);
        let result_high_filter =
            curr_read_high_filter.set_mod_data_restricted_options(&record, &options_high_filter)?;
        let count_high_filter = result_high_filter
            .base_count_per_mod()
            .values()
            .sum::<u32>();

        // Higher base quality threshold should filter out more modifications
        assert!(
            count_high_filter == count_no_filter,
            "Base quality threshold makes no difference; there are no base qualities in this file"
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
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let record = reader.records().nth(1).unwrap()?;

        // Apply multiple filters: trimming + probability + region
        let curr_read = CurrRead::default().try_from_only_alignment(&record)?;
        let region = Bed3::new(2, 50, 71); // Partial region
        let options = MockModOptions::new()
            .with_trim_read_ends(5)
            .with_mod_prob_filter(ThresholdState::GtEq(200))
            .with_region_filter(region);

        let result = curr_read.set_mod_data_restricted_options(&record, &options)?;

        // Should have some filtering effect
        // The exact count will depend on the data, but it should not error
        let mod_count = result.base_count_per_mod();
        assert!(
            mod_count.values().sum::<u32>() < 3,
            "original count is 5 I believe, combined filters should retain < 3"
        );

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
