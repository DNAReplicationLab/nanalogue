//! `GenomicRegion` struct for representing genomic coordinates
//! Handles parsing of genomic regions from standard string formats

use super::ord_pair::OrdPair;
use crate::Error;
use bedrs::prelude::Bed3;
use rust_htslib::bam;
use serde::{Deserialize, Serialize};
use std::str::FromStr;

/// Datatype holding a genomic region
#[derive(Debug, Default, Clone, PartialOrd, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub struct GenomicRegion((String, Option<OrdPair<u64>>));

/// Obtains genomic region from a string with the standard region format of name[:begin[-end]].
///
/// ```
/// use nanalogue_core::GenomicRegion;
/// use std::str::FromStr;
///
/// // Simple contig name only
/// let region = GenomicRegion::from_str("chr1")?;
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
///
/// ```
/// # use nanalogue_core::GenomicRegion;
/// # use std::str::FromStr;
/// #
/// // Contig with coordinates
/// let region = GenomicRegion::from_str("chr1:1000-2000")?;
/// let region = GenomicRegion::from_str("chr1:1000-")?;
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
///
/// ```
/// # use nanalogue_core::GenomicRegion;
/// # use std::str::FromStr;
/// #
/// // Contig name with colons (e.g., from some assemblies)
/// let region = GenomicRegion::from_str("chr1:alternate:1000-2000")?;
/// # Ok::<(), nanalogue_core::Error>(())
/// ```
impl FromStr for GenomicRegion {
    type Err = Error;

    fn from_str(val_str: &str) -> Result<Self, Self::Err> {
        let mut colon_split: Vec<&str> = val_str.split(':').collect();
        match colon_split.len() {
            0 => unreachable!(),
            1 => Ok(GenomicRegion((val_str.to_string(), None))),
            _ => {
                let interval_str = colon_split.pop().expect("no error");
                Ok(GenomicRegion((
                    colon_split.join(":").to_string(),
                    Some(OrdPair::<u64>::from_interval(interval_str)?),
                )))
            }
        }
    }
}

impl GenomicRegion {
    /// converts genomic region from genomic string representation to bed3 representation
    ///
    /// # Errors
    /// Returns an error if the contig name is not found in the header,
    /// or if the specified region is invalid (start/end position exceeds contig length).
    pub fn try_to_bed3(self, header: &bam::HeaderView) -> Result<Bed3<i32, u64>, Error> {
        let region_bed = {
            let GenomicRegion((contig_name, coords)) = self;
            let numeric_contig: i32 = header
                .tid(contig_name.as_bytes())
                .ok_or(Error::InvalidAlignCoords)?
                .try_into()?;

            let (start, end) = if let Some(c) = coords {
                let start = c.get_low();
                let end = c.get_high();

                // Check if start position exceeds contig length
                if let Some(contig_length) = header.target_len(u32::try_from(numeric_contig)?) {
                    if start >= contig_length {
                        let region_str = if end == u64::MAX {
                            format!("{contig_name}:{start}-")
                        } else {
                            format!("{contig_name}:{start}-{end}")
                        };

                        return Err(Error::InvalidRegionError {
                            region: region_str,
                            pos: start,
                            contig_length,
                        });
                    }

                    // Also check end coordinate for closed intervals
                    if end != u64::MAX && end > contig_length {
                        let region_str = format!("{contig_name}:{start}-{end}");
                        return Err(Error::InvalidRegionError {
                            region: region_str,
                            pos: end,
                            contig_length,
                        });
                    }
                }

                (start, end)
            } else {
                (u64::MIN, u64::MAX)
            };

            Bed3::<i32, u64>::new(numeric_contig, start, end)
        };
        Ok(region_bed)
    }

    /// Gets the contig
    /// ```
    /// use nanalogue_core::GenomicRegion;
    /// use std::str::FromStr;
    ///
    /// assert_eq!(GenomicRegion::from_str("chr1").unwrap().contig(), "chr1");
    /// assert_eq!(GenomicRegion::from_str("chr10_4:1000-").unwrap().contig(), "chr10_4");
    /// assert_eq!(GenomicRegion::from_str("xz_4:a:1000-2000").unwrap().contig(), "xz_4:a");
    /// ```
    #[must_use]
    pub fn contig(&self) -> &str {
        &self.0.0
    }

    /// Gets a (start, end) pair if available
    /// ```
    /// use nanalogue_core::GenomicRegion;
    /// use std::str::FromStr;
    ///
    /// assert_eq!(GenomicRegion::from_str("chr1").unwrap().start_end(), None);
    /// assert_eq!(GenomicRegion::from_str("chr10_4:1000-").unwrap().start_end(), Some((1000, u64::MAX)));
    /// assert_eq!(GenomicRegion::from_str("chr10_4:a:1000-2000").unwrap().start_end(), Some((1000, 2000)));
    /// ```
    #[must_use]
    pub fn start_end(&self) -> Option<(u64, u64)> {
        self.0.1.as_ref().map(|v| (v.get_low(), v.get_high()))
    }
}

impl<'a> TryFrom<&'a GenomicRegion> for bam::FetchDefinition<'a> {
    type Error = Error;

    /// Converts into `FetchDefinition`, a struct used by `rust_htslib` to fetch by region from
    /// indexed BAM files.
    /// ```
    /// use nanalogue_core::GenomicRegion;
    /// use rust_htslib::bam::FetchDefinition;
    /// use std::str::FromStr;
    ///
    /// // To do examples below, we need to convert every character to ASCII
    ///
    /// let region1 = GenomicRegion::from_str("chr1")?;
    /// let region1_fd = FetchDefinition::try_from(&region1)?;
    /// assert_eq!(format!("{:?}", region1_fd), "String([99, 104, 114, 49])");
    ///
    /// let region2 = GenomicRegion::from_str("chr10_4:1000-")?;
    /// let region2_fd = FetchDefinition::try_from(&region2)?;
    /// assert_eq!(format!("{:?}", region2_fd), "RegionString([99, 104, 114, 49, 48, 95, 52], \
    /// 1000, 9223372036854775807)");
    ///
    /// let region3 = GenomicRegion::from_str("chr10_4:a:1000-2000")?;
    /// // change call from `try_from` to `try_into` for sake of variety.
    /// let region3_fd :FetchDefinition = (&region3).try_into()?;
    /// assert_eq!(format!("{:?}", region3_fd), "RegionString([99, 104, 114, 49, 48, 95, 52, 58, 97], \
    /// 1000, 2000)");
    ///
    /// # Ok::<(), nanalogue_core::Error>(())
    /// ```
    fn try_from(value: &'a GenomicRegion) -> Result<bam::FetchDefinition<'a>, Error> {
        match (value.contig(), value.start_end()) {
            (c, None) => Ok(bam::FetchDefinition::from(c)),
            (c, Some((st, en))) => {
                let st_i64 = i64::try_from(st)?;
                let en_i64 = if en == u64::MAX {
                    i64::MAX
                } else {
                    i64::try_from(en)?
                };
                Ok(bam::FetchDefinition::from((c, st_i64, en_i64)))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bedrs::Coordinates as _;
    use indoc::indoc;

    /// Tests comprehensive `GenomicRegion` parsing
    #[expect(
        clippy::shadow_unrelated,
        reason = "repetition is fine; each block is clearly separated"
    )]
    #[test]
    fn genomic_region_parsing() {
        // Simple contig name only
        let region = GenomicRegion::from_str("chr1").expect("should parse");
        assert_eq!(region.0.0, "chr1");
        assert!(region.0.1.is_none());

        // Contig with coordinates
        let region = GenomicRegion::from_str("chr1:1000-2000").expect("should parse");
        assert_eq!(region.0.0, "chr1");
        let coords = region.0.1.unwrap();
        assert_eq!(coords.get_low(), 1000);
        assert_eq!(coords.get_high(), 2000);

        // Contig name with colons (e.g., from some assemblies)
        let region = GenomicRegion::from_str("chr1:alternate:1000-2000").expect("should parse");
        assert_eq!(region.0.0, "chr1:alternate");
        let coords = region.0.1.unwrap();
        assert_eq!(coords.get_low(), 1000);
        assert_eq!(coords.get_high(), 2000);

        // Complex contig names
        let region = GenomicRegion::from_str("scaffold_123:456-789").expect("should parse");
        assert_eq!(region.0.0, "scaffold_123");
        let coords = region.0.1.unwrap();
        assert_eq!(coords.get_low(), 456);
        assert_eq!(coords.get_high(), 789);
    }

    /// Tests `GenomicRegion` parsing with open-ended support
    #[test]
    fn genomic_region_open_ended() {
        // Open-ended interval support
        let region = GenomicRegion::from_str("chr1:1000-").expect("should parse");
        assert_eq!(region.0.0, "chr1");
        let coords = region.0.1.unwrap();
        assert_eq!(coords.get_low(), 1000);
        assert_eq!(coords.get_high(), u64::MAX);
    }

    /// Tests `GenomicRegion` parsing with wrong order coordinates
    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn genomic_region_parsing_wrong_order() {
        let _: GenomicRegion = GenomicRegion::from_str("chr1:2000-1000").unwrap();
    }

    /// Tests `GenomicRegion` parsing with equal start and end (strict inequality)
    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn genomic_region_parsing_equal_coordinates() {
        let _: GenomicRegion = GenomicRegion::from_str("chr1:1000-1000").unwrap();
    }

    /// Tests `GenomicRegion` parsing with invalid coordinate format
    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn genomic_region_parsing_invalid_coordinates() {
        let _: GenomicRegion = GenomicRegion::from_str("chr1:abc-def").unwrap();
    }

    /// Tests `GenomicRegion` parsing with invalid interval format (too many dashes)
    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn genomic_region_parsing_too_many_dashes() {
        let _: GenomicRegion = GenomicRegion::from_str("chr1:1000-2000-3000").unwrap();
    }

    /// Tests `GenomicRegion` parsing with missing start coordinate
    #[test]
    #[should_panic(expected = "OrdPairConversionError")]
    fn genomic_region_parsing_missing_start() {
        let _: GenomicRegion = GenomicRegion::from_str("chr1:-2000").unwrap();
    }

    #[test]
    fn genomic_region_basic_access() {
        // Test accessing the contig name and coordinates
        let region = GenomicRegion::from_str("chr22:5000-10000").expect("should parse");
        assert_eq!(region.0.0, "chr22");

        let coords = &region.0.1.unwrap();
        assert_eq!(coords.get_low(), 5000);
        assert_eq!(coords.get_high(), 10000);
    }

    #[test]
    fn genomic_region_no_coordinates() {
        let region = GenomicRegion::from_str("chrX").expect("should parse");
        assert_eq!(region.0.0, "chrX");
        assert!(region.0.1.is_none());
    }

    /// Creates a sample BAM header for testing
    fn create_test_header() -> bam::HeaderView {
        bam::HeaderView::from_bytes(
            indoc! {b"@HD\tVN:1.6\tSO:coordinate
            @SQ\tSN:chr1\tLN:248956422
            @SQ\tSN:chr2\tLN:242193529
            @SQ\tSN:chrX\tLN:156040895
            @PG\tID:test\tPN:test"
            }
            .as_ref(),
        )
    }

    /// Tests basic region conversion to BED3 format with chromosome and coordinates
    #[test]
    fn try_to_bed3_basic_region() {
        // Create test header
        let header = create_test_header();

        // Test region with coordinates
        let region = GenomicRegion::from_str("chr1:1000-2000").expect("should parse");
        let bed3 = region.try_to_bed3(&header).expect("should convert to bed3");

        // Verify the BED3 record
        assert_eq!(*bed3.chr(), 0); // chr1 should be the first chromosome (index 0)
        assert_eq!(bed3.start(), 1000);
        assert_eq!(bed3.end(), 2000);
    }

    /// Tests conversion of a region with no coordinates to BED3 format
    #[test]
    fn try_to_bed3_no_coordinates() {
        // Create test header
        let header = create_test_header();

        // Test region without coordinates (should use full chromosome)
        let region = GenomicRegion::from_str("chr2").expect("should parse");
        let bed3 = region.try_to_bed3(&header).expect("should convert to bed3");

        // Verify the BED3 record
        assert_eq!(*bed3.chr(), 1); // chr2 should be the second chromosome (index 1)
        assert_eq!(bed3.start(), u64::MIN);
        assert_eq!(bed3.end(), u64::MAX);
    }

    /// Tests error case when start position exceeds contig length
    #[test]
    #[should_panic(expected = "InvalidRegionError")]
    fn try_to_bed3_start_exceeds_contig_length() {
        let header = create_test_header();

        // Create a region with start position > contig length for chr1
        let region = GenomicRegion::from_str("chr1:300000000-400000000").expect("should parse");

        // This should panic with InvalidRegionError
        let _: Bed3<i32, u64> = region.try_to_bed3(&header).unwrap();
    }

    /// Tests error case when contig doesn't exist in the header
    #[test]
    #[should_panic(expected = "InvalidAlignCoords")]
    fn try_to_bed3_nonexistent_contig() {
        let header = create_test_header();

        // Create a region with a contig that doesn't exist in the header
        let region = GenomicRegion::from_str("nonexistent_contig:1000-2000").expect("should parse");

        // This should panic with InvalidAlignCoords error
        let _: Bed3<i32, u64> = region.try_to_bed3(&header).unwrap();
    }

    /// Tests conversion of an open-ended region to BED3 format
    #[test]
    fn try_to_bed3_open_ended() {
        let header = create_test_header();

        // Create an open-ended region
        let region = GenomicRegion::from_str("chr1:1000-").expect("should parse");

        // Convert to BED3
        let bed3 = region.try_to_bed3(&header).expect("should convert to bed3");

        // Verify the BED3 record
        assert_eq!(*bed3.chr(), 0);
        assert_eq!(bed3.start(), 1000);
        assert_eq!(bed3.end(), u64::MAX);
    }
}
