//! BAM/CRAM matrix tests for simulated alignment parity.

#[cfg(test)]
mod tests {
    use nanalogue_core::simulate_mod_bam::{AlignmentFormat, SimulationConfig, TempBamSimulation};
    use nanalogue_core::{CurrRead, Error, ThresholdState};
    use rust_htslib::bam;
    use rust_htslib::bam::Read as _;
    use rust_htslib::bam::ext::BamRecordExtensions as _;
    use rust_htslib::bam::record::Aux;
    use std::fs;

    /// Stable per-record summary used to compare BAM and CRAM decoding.
    #[derive(Debug, Clone, PartialEq, Eq)]
    struct BasicRecordSummary {
        /// Query name.
        read_id: String,
        /// BAM flag bits.
        flags: u16,
        /// Numeric contig id, or -1 for unmapped records.
        tid: i32,
        /// Leftmost reference position, or -1 for unmapped records.
        pos: i64,
        /// Rightmost reference position.
        reference_end: i64,
        /// Mapping quality.
        mapq: u8,
        /// Sequence length.
        seq_len: usize,
        /// CIGAR string.
        cigar: String,
        /// Raw MM tag string, if present.
        mm: Option<String>,
        /// Raw ML probabilities, if present.
        ml: Option<Vec<u8>>,
    }

    /// Semantic summary for one mapped, mod-carrying record.
    #[derive(Debug, Clone, PartialEq, Eq)]
    struct DecodedModSummary {
        /// Query name.
        read_id: String,
        /// Internal read-state label.
        read_state: String,
        /// Sequence length.
        seq_len: u32,
        /// Alignment span on the reference.
        align_len: u32,
        /// Numeric contig id.
        contig_id: i32,
        /// Alignment start.
        start: u32,
        /// Counts per modification type.
        mod_counts: Vec<(char, u32)>,
    }

    /// Shared BAM/CRAM simulation pair reused within one test.
    struct MatrixFixture {
        /// Seeded BAM simulation.
        bam_sim: TempBamSimulation,
        /// Seeded CRAM simulation.
        cram_sim: TempBamSimulation,
    }

    /// Deterministic simulation config used for BAM/CRAM parity tests.
    fn matrix_config() -> SimulationConfig {
        let config_json = r#"{
        "contigs": {
            "number": 1,
            "len_range": [2000, 2000],
            "repeated_seq": "ACGTACGT"
        },
        "reads": [{
            "number": 200,
            "mapq_range": [30, 30],
            "base_qual_range": [35, 35],
            "len_range": [0.4, 0.4],
            "mods": [{
                "base": "C",
                "is_strand_plus": true,
                "mod_code": "m",
                "win": [3],
                "mod_range": [[0.9, 0.9]]
            }]
        }],
        "seed": 12345
    }"#;
        serde_json::from_str(config_json).expect("matrix simulation config should deserialize")
    }

    /// Build one seeded BAM/CRAM pair for a single test.
    fn matrix_fixture() -> Result<MatrixFixture, Error> {
        let config = matrix_config();
        Ok(MatrixFixture {
            bam_sim: TempBamSimulation::new(config.clone(), AlignmentFormat::Bam)?,
            cram_sim: TempBamSimulation::new(config, AlignmentFormat::Cram)?,
        })
    }

    /// Open a plain BAM/CRAM reader.
    fn open_reader(sim: &TempBamSimulation) -> Result<bam::Reader, Error> {
        Ok(bam::Reader::from_path(sim.bam_path())?)
    }

    /// Open an indexed BAM/CRAM reader.
    fn open_indexed_reader(sim: &TempBamSimulation) -> Result<bam::IndexedReader, Error> {
        Ok(bam::IndexedReader::from_path(sim.bam_path())?)
    }

    /// Summarize one BAM/CRAM record using stable comparison fields.
    fn summarize_record(record: &bam::Record) -> Result<BasicRecordSummary, Error> {
        let mm = match record.aux(b"MM") {
            Ok(Aux::String(value)) => Some(value.to_owned()),
            Ok(_) => return Err(Error::InvalidState("MM tag had unexpected type".into())),
            Err(_) => None,
        };
        let ml = match record.aux(b"ML") {
            Ok(Aux::ArrayU8(values)) => Some(values.iter().collect()),
            Ok(_) => return Err(Error::InvalidState("ML tag had unexpected type".into())),
            Err(_) => None,
        };

        Ok(BasicRecordSummary {
            read_id: String::from_utf8_lossy(record.qname()).into_owned(),
            flags: record.flags(),
            tid: record.tid(),
            pos: record.pos(),
            reference_end: record.reference_end(),
            mapq: record.mapq(),
            seq_len: record.seq_len(),
            cigar: record.cigar().to_string(),
            mm,
            ml,
        })
    }

    /// Summarize every record yielded by a BAM/CRAM reader.
    fn summarize_reader_records<R: bam::Read>(
        reader: &mut R,
    ) -> Result<Vec<BasicRecordSummary>, Error> {
        reader
            .records()
            .map(|result| summarize_record(&result?))
            .collect()
    }

    /// Collect stable per-record summaries from a plain reader.
    fn summarize_records(sim: &TempBamSimulation) -> Result<Vec<BasicRecordSummary>, Error> {
        let mut reader = open_reader(sim)?;
        summarize_reader_records(&mut reader)
    }

    /// Collect stable per-record summaries from an indexed regional fetch.
    fn summarize_region_records(
        sim: &TempBamSimulation,
        tid: u32,
        start: i64,
        end: i64,
    ) -> Result<Vec<BasicRecordSummary>, Error> {
        let mut reader = open_indexed_reader(sim)?;
        reader.fetch((tid, start, end))?;
        summarize_reader_records(&mut reader)
    }

    /// Decode one mapped, mod-bearing record into a semantic summary.
    fn summarize_mod_record(record: &bam::Record) -> Result<Option<DecodedModSummary>, Error> {
        if record.tid() < 0 || record.aux(b"MM").is_err() {
            return Ok(None);
        }

        let curr_read = CurrRead::default()
            .set_read_state_and_id(record)?
            .set_seq_len(record)?
            .set_align_len(record)?
            .set_contig_id_and_start(record)?
            .set_mod_data(record, ThresholdState::default(), 0)?;
        let (contig_id, start) = curr_read.contig_id_and_start()?;
        let mut mod_counts: Vec<(char, u32)> = curr_read
            .base_count_per_mod()
            .into_iter()
            .map(|(mod_char, count)| (mod_char.val(), count))
            .collect();
        mod_counts.sort_unstable();

        Ok(Some(DecodedModSummary {
            read_id: curr_read.read_id().to_owned(),
            read_state: curr_read.read_state().to_string(),
            seq_len: curr_read.seq_len()?,
            align_len: curr_read.align_len()?,
            contig_id,
            start,
            mod_counts,
        }))
    }

    /// Decode all mapped, mod-bearing records into semantic summaries.
    fn summarize_mod_records(sim: &TempBamSimulation) -> Result<Vec<DecodedModSummary>, Error> {
        let mut reader = open_reader(sim)?;
        let mut summaries = Vec::new();
        for result in reader.records() {
            if let Some(summary) = summarize_mod_record(&result?)? {
                summaries.push(summary);
            }
        }
        Ok(summaries)
    }

    /// Extract a fetchable region from the first mapped record summary.
    fn first_mapped_region(summaries: &[BasicRecordSummary]) -> (u32, i64, i64) {
        let mapped = summaries
            .iter()
            .find(|summary| summary.tid >= 0 && summary.pos >= 0)
            .expect("matrix simulation should yield at least one mapped record");
        let tid = u32::try_from(mapped.tid).expect("non-negative tid should fit into u32");
        (tid, mapped.pos, mapped.reference_end)
    }

    /// Simulated BAM and CRAM should decode to the same stable per-record summaries.
    #[test]
    fn simulated_bam_and_cram_have_matching_record_summaries() -> Result<(), Error> {
        let fixture = matrix_fixture()?;

        let bam_summaries = summarize_records(&fixture.bam_sim)?;
        let cram_summaries = summarize_records(&fixture.cram_sim)?;

        assert!(
            !bam_summaries.is_empty(),
            "matrix BAM should contain records"
        );
        assert_eq!(bam_summaries, cram_summaries);
        Ok(())
    }

    /// Simulated BAM and CRAM should support identical indexed regional fetches.
    #[test]
    fn simulated_bam_and_cram_have_matching_indexed_region_summaries() -> Result<(), Error> {
        let fixture = matrix_fixture()?;

        let bam_summaries = summarize_records(&fixture.bam_sim)?;
        let (tid, start, end) = first_mapped_region(&bam_summaries);
        let bam_region = summarize_region_records(&fixture.bam_sim, tid, start, end)?;
        let cram_region = summarize_region_records(&fixture.cram_sim, tid, start, end)?;

        assert!(
            !bam_region.is_empty(),
            "indexed BAM fetch should return records"
        );
        assert_eq!(bam_region, cram_region);
        Ok(())
    }

    /// Mapped mod-bearing records should decode to the same semantic summaries in BAM and CRAM.
    #[test]
    fn simulated_bam_and_cram_decode_matching_mod_record_summaries() -> Result<(), Error> {
        let fixture = matrix_fixture()?;

        let bam_mod_summaries = summarize_mod_records(&fixture.bam_sim)?;
        let cram_mod_summaries = summarize_mod_records(&fixture.cram_sim)?;

        assert!(
            !bam_mod_summaries.is_empty(),
            "matrix BAM should contain mapped mod-bearing records"
        );
        assert_eq!(bam_mod_summaries, cram_mod_summaries);
        Ok(())
    }

    /// Simulated BAM and CRAM should emit measurable, non-empty file sizes.
    #[test]
    fn simulated_bam_and_cram_have_non_empty_file_sizes() -> Result<(), Error> {
        let fixture = matrix_fixture()?;

        let bam_len = fs::metadata(fixture.bam_sim.bam_path())?.len();
        let cram_len = fs::metadata(fixture.cram_sim.bam_path())?.len();

        assert!(bam_len > 0, "simulated BAM should be non-empty");
        assert!(cram_len > 0, "simulated CRAM should be non-empty");
        assert_ne!(
            bam_len, cram_len,
            "BAM and CRAM files should differ in size"
        );
        Ok(())
    }
}
