//! Tests for `set_mod_data_restricted_options` method in `read_utils.rs`
//! Covers filtering by region, tag, strand, probability, base quality, and read end trimming

#[cfg(feature = "polars")]
use nanalogue_core::curr_reads_to_dataframe;
use nanalogue_core::{
    CurrRead, Error, GenomicBed3, InputModOptions, InputRegionOptions, ModChar,
    RestrictModCalledStrand, ThresholdState, nanalogue_bam_reader,
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
    region_filter: Option<GenomicBed3>,
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

    fn with_region_filter(mut self, region: GenomicBed3) -> Self {
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
    fn region_filter(&self) -> &Option<GenomicBed3> {
        &self.region_filter
    }

    fn region_filter_genomic_string(&self) -> Option<nanalogue_core::GenomicRegion> {
        None
    }

    fn set_region_filter(&mut self, value: Option<GenomicBed3>) {
        self.region_filter = value;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const HIGH_PROBABILITY_THRESHOLD: u8 = 255;
    const HIGH_PROBABILITY_POSITIONS: [u32; 3] = [3, 8, 47];
    const LOW_PROBABILITY_THRESHOLD: u8 = 128;
    const PARTIAL_OVERLAP_END: u32 = 100;
    const PARTIAL_OVERLAP_REF_POSITIONS: [i64; 3] = [50, 62, 70];
    const PARTIAL_OVERLAP_START: u32 = 50;
    const TRIM_READ_ENDS: usize = 10;
    const TRIMMED_POSITION: u32 = 27;

    fn filtered_second_record(
        options: &MockModOptions,
    ) -> Result<CurrRead<nanalogue_core::read_utils::AlignAndModData>, Error> {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let record = reader.records().nth(1).unwrap()?;
        CurrRead::default()
            .try_from_only_alignment(&record)?
            .set_mod_data_restricted_options(&record, options)
    }

    fn annotations(
        read: &CurrRead<nanalogue_core::read_utils::AlignAndModData>,
    ) -> Vec<(u32, i64, u8)> {
        read.mod_data()
            .0
            .base_mods
            .iter()
            .flat_map(|base_mod| {
                base_mod.ranges.annotations.iter().map(|annotation| {
                    (
                        annotation.pos,
                        annotation.ref_pos.map_or(-1, i64::from),
                        annotation.qual,
                    )
                })
            })
            .collect()
    }

    fn combined_filter_options() -> MockModOptions {
        MockModOptions::new()
            .with_trim_read_ends(5)
            .with_mod_prob_filter(ThresholdState::GtEq(200))
            .with_region_filter(GenomicBed3::new(2, 50, 71))
    }

    fn partial_overlap_options() -> MockModOptions {
        MockModOptions::new().with_region_filter(GenomicBed3::new(
            2,
            PARTIAL_OVERLAP_START,
            PARTIAL_OVERLAP_END,
        ))
    }

    fn probability_filter_options(threshold: u8) -> MockModOptions {
        MockModOptions::new().with_mod_prob_filter(ThresholdState::GtEq(threshold))
    }

    fn trim_read_ends_options() -> MockModOptions {
        MockModOptions::new().with_trim_read_ends(TRIM_READ_ENDS)
    }

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

    #[cfg(feature = "polars")]
    #[test]
    fn high_probability_filter_with_polars() -> Result<(), Error> {
        // Test: Probability threshold filtering with explicit position checking
        // Second record has T mods at positions with qualities: [3,26,221], [8,31,242], [27,50,3], [39,62,47], [47,70,239]
        // Test with threshold 128: should keep mods with qual >= 128 (positions 3, 8, 47)
        let options_low = probability_filter_options(LOW_PROBABILITY_THRESHOLD);
        let df_low = curr_reads_to_dataframe(&[filtered_second_record(&options_low)?])?;

        // Should have 3 mods (qual 221, 242, 239 are all >= 128)
        assert_eq!(
            df_low.height(),
            HIGH_PROBABILITY_POSITIONS.len(),
            "Should have 3 mods with qual >= 128"
        );

        // Check that all remaining mods have quality >= 128
        let qual_col = df_low.column("mod_quality")?.u32()?;
        for i in 0..df_low.height() {
            let qual = qual_col.get(i).unwrap();
            assert!(
                qual >= u32::from(LOW_PROBABILITY_THRESHOLD),
                "All mods should have quality >= 128, found {qual}"
            );
        }

        // Check specific positions that should be present (read positions 3, 8, 47)
        let pos_col = df_low.column("position")?.u32()?;
        let mut positions: Vec<u32> = (0..df_low.height())
            .filter_map(|i| pos_col.get(i))
            .collect();
        positions.sort_unstable();
        assert_eq!(positions, HIGH_PROBABILITY_POSITIONS);

        // Test with threshold 255: should filter out all mods (none have qual >= 255)
        let options_high = probability_filter_options(HIGH_PROBABILITY_THRESHOLD);
        let df_high = curr_reads_to_dataframe(&[filtered_second_record(&options_high)?])?;

        // Should have 0 mods
        assert_eq!(df_high.height(), 0, "No mods should have quality >= 255");

        Ok(())
    }

    #[test]
    fn high_probability_filter_without_polars() -> Result<(), Error> {
        let low = filtered_second_record(&probability_filter_options(LOW_PROBABILITY_THRESHOLD))?;
        let low_annotations = annotations(&low);
        assert_eq!(low_annotations.len(), 3);
        assert!(
            low_annotations
                .iter()
                .all(|entry| entry.2 >= LOW_PROBABILITY_THRESHOLD)
        );
        let mut positions = low_annotations
            .iter()
            .map(|entry| entry.0)
            .collect::<Vec<_>>();
        positions.sort_unstable();
        assert_eq!(positions, HIGH_PROBABILITY_POSITIONS);

        let high = filtered_second_record(&probability_filter_options(HIGH_PROBABILITY_THRESHOLD))?;
        assert!(annotations(&high).is_empty());
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
        let region = GenomicBed3::new(2, 0, 100); // Fully contains the read
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

    #[cfg(feature = "polars")]
    #[test]
    fn region_filter_partial_overlap_with_polars() -> Result<(), Error> {
        // Test: Partial overlap - only mods in overlapping region should be retained
        // Second record: contig 2, align 23-71, mods at ref positions [26, 31, 50, 62, 70]
        // Create a region that partially overlaps (region 50-100 overlaps with read 23-71)
        // Only mods at ref positions >= 50 should be kept (positions 50, 62, 70)
        let options = partial_overlap_options();
        let df = curr_reads_to_dataframe(&[filtered_second_record(&options)?])?;

        // Should have exactly 3 mods (at ref positions 50, 62, 70)
        assert_eq!(
            df.height(),
            PARTIAL_OVERLAP_REF_POSITIONS.len(),
            "Should have 3 mods in region 50-100"
        );

        // Check that all ref_positions are within the region [50, 100)
        let ref_pos_col = df.column("ref_position")?.i64()?;
        for i in 0..df.height() {
            let ref_pos = ref_pos_col.get(i).unwrap();
            assert!(
                (i64::from(PARTIAL_OVERLAP_START)..i64::from(PARTIAL_OVERLAP_END))
                    .contains(&ref_pos),
                "Ref position {ref_pos} should be in range [50, 100)"
            );
        }

        // Check specific ref positions that should be present
        let mut ref_positions: Vec<i64> = (0..df.height())
            .filter_map(|i| ref_pos_col.get(i))
            .collect();
        ref_positions.sort_unstable();
        assert_eq!(ref_positions, PARTIAL_OVERLAP_REF_POSITIONS);

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
    fn region_filter_partial_overlap_without_polars() -> Result<(), Error> {
        let read = filtered_second_record(&partial_overlap_options())?;
        let mut ref_positions = annotations(&read)
            .into_iter()
            .map(|entry| entry.1)
            .collect::<Vec<_>>();
        ref_positions.sort_unstable();
        assert_eq!(ref_positions, PARTIAL_OVERLAP_REF_POSITIONS);
        assert!(ref_positions.iter().all(|position| {
            (i64::from(PARTIAL_OVERLAP_START)..i64::from(PARTIAL_OVERLAP_END)).contains(position)
        }));
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
        let region = GenomicBed3::new(0, 0, 100); // Different contig, no overlap
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
        let region = GenomicBed3::new(2, 23, 71); // Exact match
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

    #[cfg(feature = "polars")]
    #[test]
    fn trim_read_ends_with_polars() -> Result<(), Error> {
        // Test: Trimming read ends should exclude modifications near ends
        // Second record: seq_len=48, mods at read positions [3, 8, 27, 39, 47]
        // Test with trimming 10bp from each end (keeps positions 10-37 in a 48bp read)
        let options_trim = trim_read_ends_options();
        let df_trim = curr_reads_to_dataframe(&[filtered_second_record(&options_trim)?])?;

        // Should only keep position 27 (positions 3, 8 are < 10; positions 39, 47 are >= 38)
        assert_eq!(
            df_trim.height(),
            1,
            "Should have 1 mod after trimming 10bp from each end"
        );

        // Check that remaining position is 27
        let pos_col = df_trim.column("position")?.u32()?;
        let position = pos_col.get(0).unwrap();
        assert_eq!(
            position, TRIMMED_POSITION,
            "Only position 27 should remain after trimming"
        );

        Ok(())
    }

    #[test]
    fn trim_read_ends_without_polars() -> Result<(), Error> {
        let read = filtered_second_record(&trim_read_ends_options())?;
        assert_eq!(
            annotations(&read)
                .iter()
                .map(|entry| entry.0)
                .collect::<Vec<_>>(),
            [TRIMMED_POSITION]
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

    #[cfg(feature = "polars")]
    #[test]
    fn combined_filters_with_polars() -> Result<(), Error> {
        // Test: Multiple filters working together
        // Second record: mods at [3,26,221], [8,31,242], [27,50,3], [39,62,47], [47,70,239]
        // Apply multiple filters: trimming + probability + region
        let options = combined_filter_options();
        let df = curr_reads_to_dataframe(&[filtered_second_record(&options)?])?;

        // Expected: NO mods pass all three filters
        // - Position 8 (ref 31, qual 242): passes trim & qual, but fails region (31 < 50)
        // - Position 27 (ref 50, qual 3): passes trim & region, but fails qual (3 < 200)
        // - Position 39 (ref 62, qual 47): passes trim & region, but fails qual (47 < 200)
        // - Positions 3 and 47 are trimmed
        assert_eq!(df.height(), 0, "No mods should pass all three filters");

        Ok(())
    }

    #[test]
    fn combined_filters_without_polars() -> Result<(), Error> {
        let read = filtered_second_record(&combined_filter_options())?;
        assert!(
            annotations(&read).is_empty(),
            "No mods should pass all three filters"
        );
        Ok(())
    }

    #[test]
    fn empty_interval_due_to_no_overlap_removes_all_mods() -> Result<(), Error> {
        // Test: Empty interval (0..0) should remove all modifications
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        let record = reader.records().nth(1).unwrap()?;

        // Create a scenario with no overlap (different contig) which creates 0..0 interval
        let curr_read = CurrRead::default().try_from_only_alignment(&record)?;
        let region = GenomicBed3::new(1, 0, 10); // Different contig than record's contig 2
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
