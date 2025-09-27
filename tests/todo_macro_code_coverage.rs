#![deny(
    missing_copy_implementations,
    missing_debug_implementations,
    missing_docs,
    trivial_casts,
    trivial_numeric_casts,
    unused_extern_crates,
    unused_import_braces,
    unused_qualifications
)]
//! Test file to provide code coverage for all todo!() macros in trait default implementations
//! This ensures code coverage tools recognize that todo!() panics are intentionally tested

use bedrs::prelude::Bed3;
use nanalogue_core::{
    BamPreFilt, F32Bw0and1, ReadStates,
    cli::{InputModOptions, InputRegionOptions, TagState},
    utils::FilterByRefCoords,
};
use std::collections::HashSet;
use std::num::NonZeroU32;

/// Test struct that uses TagState default implementation
#[derive(Debug, Default)]
struct TestTagState;

impl TagState for TestTagState {}

/// Test struct that uses InputModOptions default implementation
#[derive(Debug, Default)]
struct TestInputModOptions;

impl InputModOptions for TestInputModOptions {}

/// Test struct that uses InputRegionOptions default implementation
#[derive(Debug, Default)]
struct TestInputRegionOptions;

impl InputRegionOptions for TestInputRegionOptions {}

/// Test struct that uses FilterByRefCoords default implementation
#[derive(Debug, Default)]
struct TestFilterByRefCoords;

impl FilterByRefCoords for TestFilterByRefCoords {}

/// Test struct that uses BamPreFilt default implementation
#[derive(Debug, Default)]
struct TestBamPreFilt;

impl BamPreFilt for TestBamPreFilt {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_tag_state_tag_panics() {
        let test_obj = TestTagState::default();
        test_obj.tag();
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_input_mod_options_tag_panics() {
        let test_obj = TestInputModOptions::default();
        test_obj.tag();
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_input_mod_options_mod_strand_panics() {
        let test_obj = TestInputModOptions::default();
        test_obj.mod_strand();
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_input_mod_options_mod_prob_filter_panics() {
        let test_obj = TestInputModOptions::default();
        test_obj.mod_prob_filter();
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_input_mod_options_trim_read_ends_mod_panics() {
        let test_obj = TestInputModOptions::default();
        test_obj.trim_read_ends_mod();
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_input_mod_options_base_qual_filter_mod_panics() {
        let test_obj = TestInputModOptions::default();
        test_obj.base_qual_filter_mod();
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_input_region_options_region_filter_panics() {
        let test_obj = TestInputRegionOptions::default();
        test_obj.region_filter();
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_input_region_options_region_filter_genomic_string_panics() {
        let test_obj = TestInputRegionOptions::default();
        test_obj.region_filter_genomic_string();
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_input_region_options_set_region_filter_panics() {
        let mut test_obj = TestInputRegionOptions::default();
        test_obj.set_region_filter(None);
    }

    #[test]
    fn test_input_region_options_is_full_overlap_false() {
        let test_obj = TestInputRegionOptions::default();
        assert!(!test_obj.is_full_overlap());
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_filter_by_ref_coords_filter_by_ref_pos_panics() {
        let mut test_obj = TestFilterByRefCoords::default();
        test_obj.filter_by_ref_pos(0, 100);
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_bam_pre_filt_pre_filt_panics() {
        let test_obj = TestBamPreFilt::default();
        let input_bam = nanalogue_core::cli::InputBam {
            bam_path: "test.bam".to_string(),
            min_seq_len: 0,
            min_align_len: None,
            read_id: None,
            read_id_list: None,
            read_id_set: None,
            threads: NonZeroU32::new(2).unwrap(),
            include_zero_len: false,
            read_filter: None,
            sample_fraction: F32Bw0and1::new(1.0).unwrap(),
            mapq_filter: 0,
            exclude_mapq_unavail: false,
            region: None,
            region_bed3: None,
            full_region: false,
        };
        test_obj.pre_filt(&input_bam);
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_bam_pre_filt_filt_by_len_panics() {
        let test_obj = TestBamPreFilt::default();
        test_obj.filt_by_len(100, false);
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_bam_pre_filt_filt_by_align_len_panics() {
        let test_obj = TestBamPreFilt::default();
        test_obj.filt_by_align_len(100);
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_bam_pre_filt_filt_by_read_id_panics() {
        let test_obj = TestBamPreFilt::default();
        test_obj.filt_by_read_id("test_read");
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_bam_pre_filt_filt_by_read_id_set_panics() {
        let test_obj = TestBamPreFilt::default();
        let read_ids = HashSet::new();
        test_obj.filt_by_read_id_set(&read_ids);
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_bam_pre_filt_filt_by_bitwise_or_flags_panics() {
        let test_obj = TestBamPreFilt::default();
        let states = ReadStates::default();
        test_obj.filt_by_bitwise_or_flags(&states);
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_bam_pre_filt_filt_random_subset_panics() {
        let test_obj = TestBamPreFilt::default();
        let fraction = F32Bw0and1::new(0.5).unwrap();
        test_obj.filt_random_subset(fraction);
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_bam_pre_filt_filt_by_mapq_panics() {
        let test_obj = TestBamPreFilt::default();
        test_obj.filt_by_mapq(30, false);
    }

    #[test]
    #[should_panic(expected = "not yet implemented")]
    fn test_bam_pre_filt_filt_by_region_panics() {
        let test_obj = TestBamPreFilt::default();
        let region = Bed3::new(1, 1000, 2000);
        test_obj.filt_by_region(&region, false);
    }
}
