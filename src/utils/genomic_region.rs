//! GenomicRegion struct for representing genomic coordinates
//! Handles parsing of genomic regions from standard string formats

use super::ord_pair::OrdPair;
use crate::Error;
use bedrs::prelude::Bed3;
use rust_htslib::bam;
use serde::{Deserialize, Serialize};
use std::str::FromStr;

/// Datatype holding a genomic region
#[derive(Debug, Default, Clone, PartialOrd, PartialEq, Serialize, Deserialize)]
pub struct GenomicRegion(pub (String, Option<OrdPair<u64>>));

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
        let mut colon_split: Vec<&str> = val_str.split(":").collect();
        match colon_split.len() {
            0 => Err(Error::InvalidAlignCoords),
            1 => Ok(GenomicRegion((val_str.to_string(), None))),
            _ => {
                let interval_str = colon_split.pop().ok_or(Error::UnknownError)?;
                Ok(GenomicRegion((
                    colon_split.join(":").to_string(),
                    Some(OrdPair::<u64>::from_interval(interval_str)?),
                )))
            }
        }
    }
}

/// Converts genomic region from genomic string representation to bed3 representation
impl GenomicRegion {
    /// converts genomic region from genomic string representation to bed3 representation
    pub fn try_to_bed3(self, header: bam::HeaderView) -> Result<Bed3<i32, u64>, Error> {
        let region_bed = {
            let GenomicRegion((contig_name, coords)) = &self;
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
                            format!("{}:{}-", contig_name, start)
                        } else {
                            format!("{}:{}-{}", contig_name, start, end)
                        };

                        return Err(Error::InvalidRegionError {
                            region: region_str,
                            start,
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
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tests comprehensive GenomicRegion parsing
    #[test]
    fn test_genomic_region_parsing() {
        // Simple contig name only
        let region = GenomicRegion::from_str("chr1").expect("should parse");
        assert_eq!(region.0.0, "chr1");
        assert_eq!(region.0.1, None);

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

    /// Tests GenomicRegion parsing with open-ended support
    #[test]
    fn test_genomic_region_open_ended() {
        // Open-ended interval support
        let region = GenomicRegion::from_str("chr1:1000-").expect("should parse");
        assert_eq!(region.0.0, "chr1");
        let coords = region.0.1.unwrap();
        assert_eq!(coords.get_low(), 1000);
        assert_eq!(coords.get_high(), u64::MAX);
    }

    /// Tests GenomicRegion parsing error cases
    #[test]
    fn test_genomic_region_parsing_errors() {
        // Wrong order coordinates should fail
        assert!(matches!(
            GenomicRegion::from_str("chr1:2000-1000"),
            Err(Error::WrongOrder)
        ));

        // Equal start and end should now fail (strict inequality)
        assert!(matches!(
            GenomicRegion::from_str("chr1:1000-1000"),
            Err(Error::WrongOrder)
        ));

        // Invalid coordinate format should fail
        assert!(GenomicRegion::from_str("chr1:abc-def").is_err());

        // Invalid format should fail
        assert!(GenomicRegion::from_str("chr1:1000-2000-3000").is_err());
    }

    #[test]
    fn test_genomic_region_basic_access() {
        // Test accessing the contig name and coordinates
        let region = GenomicRegion::from_str("chr22:5000-10000").expect("should parse");
        assert_eq!(region.0.0, "chr22");

        if let Some(coords) = &region.0.1 {
            assert_eq!(coords.get_low(), 5000);
            assert_eq!(coords.get_high(), 10000);
        } else {
            panic!("Expected coordinates to be present");
        }
    }

    #[test]
    fn test_genomic_region_no_coordinates() {
        let region = GenomicRegion::from_str("chrX").expect("should parse");
        assert_eq!(region.0.0, "chrX");
        assert!(region.0.1.is_none());
    }
}
