//! # Nanalogue Core (Nanalogue = Nucleic Acid Analogue)
//!
//! We process and calculate data associated with DNA molecules, their alignments to
//! reference genomes, modification information on them, and other miscellaneous
//! information.

use bedrs::{Bed3, Coordinates as _};
use bio::alphabets::dna::revcomp;
use bio_types::sequence::SequenceRead as _;
use fibertools_rs::utils::{
    bamannotations::{FiberAnnotation, Ranges},
    basemods::{BaseMod, BaseMods},
    bio_io::{convert_seq_uppercase, get_u8_tag},
};
use rand::random;
use regex::Regex;
use rust_htslib::{bam, bam::ext::BamRecordExtensions as _, bam::record::Aux, tpool};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead as _, BufReader};
use std::sync::LazyLock;

// Declare the modules.
pub mod analysis;
pub mod cli;
pub mod commands;
pub mod error;
pub mod file_utils;
pub mod read_utils;
pub mod simulate_mod_bam;
pub mod subcommands;
pub mod utils;

// Re-exports
pub use cli::{
    InputBam, InputModOptions, InputMods, InputRegionOptions, InputWindowing, OptionalTag,
    RequiredTag, SeqDisplayOptions,
};
pub use error::Error;
pub use file_utils::{
    nanalogue_bam_reader, nanalogue_bam_reader_from_stdin, nanalogue_bam_reader_from_url,
    nanalogue_indexed_bam_reader, nanalogue_indexed_bam_reader_from_url, write_bam_denovo,
    write_fasta,
};
pub use read_utils::{
    AlignmentInfoBuilder, CurrRead, CurrReadBuilder, ModTableEntryBuilder, curr_reads_to_dataframe,
};
pub use simulate_mod_bam::SimulationConfig;
pub use subcommands::{find_modified_reads, read_info, read_stats, reads_table, window_reads};
pub use utils::{
    AllowedAGCTN, Contains, DNARestrictive, F32AbsValAtMost1, F32Bw0and1, FilterByRefCoords,
    GenomicRegion, GetDNARestrictive, Intersects, ModChar, OrdPair, PathOrURLOrStdin, ReadState,
    ReadStates, RestrictModCalledStrand, SeqCoordCalls, ThresholdState,
};

/// Extracts mod information from BAM record to the `fibertools-rs` `BaseMods` Struct.
///
/// We are copying and modifying code from the fibertools-rs repository
/// (<https://github.com/fiberseq/fibertools-rs>) which is under the MIT license
/// (please see their Cargo.toml) to create this function.
///
/// Function should cover almost all mod bam cases, but will fail in the following scenarios:
/// - If multiple mods are present on the same base e.g. methylation and hydroxymethylation,
///   most BAM files the author has come across use the notation MM:Z:C+m,...;C+h,...;,
///   which this function can parse. But the notation MM:Z:C+mh,...; is also allowed.
///   We do not parse this for now, please contribute code if you want to add this functionality!
/// - We do not know any technology/technique which, for instance, replaces both C and T
///   with `BrdU` i.e. two different bases have the same substitution. This also leads to failure.
///   If you deal with this, please notify us! On the other hand, it is possible that both C
///   and A are methylated e.g. 5-Methylcytosine and 6-Methyladenine. Our program can deal
///   with this as the "tags" corresponding to these are 'm' and 'a' respectively i.e. the
///   tags are different. For this failure mode, the modifications have to be identical,
///   not just conceptually identical i.e. although 5-Methylcytosine and 6-Methyladenine
///   fall under "methylation", they are not chemically identical and thus have different
///   tags associated with them in the mod BAM format, and thus will not lead to failure.
///
/// Following is an example of reading and parsing modification data, with
/// no filters applied.
/// We are using an example mod BAM file which has very short reads and very few
/// modified positions, and just examining two reads in it below.
/// ```
/// use nanalogue_core::{Error, nanalogue_bam_reader, nanalogue_mm_ml_parser};
/// use rust_htslib::bam::Read;
/// use fibertools_rs::utils::basemods::{BaseMods, BaseMod};
/// use fibertools_rs::utils::bamannotations::{FiberAnnotation, Ranges};
/// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
/// let mut count = 0;
/// for record in reader.records(){
///     let r = record?;
///     let Ok(BaseMods{base_mods: v}) = nanalogue_mm_ml_parser(&r, |&_| true, |&_| true,
///         |&_, &_, &_| true, 0) else { unreachable!() };
///     match count {
///     0 => assert_eq!(v, vec![BaseMod{
///             modified_base: b'T',
///             strand: '+',
///             modification_type: 'T',
///             ranges: Ranges {
///                 annotations: vec![
///                     FiberAnnotation {
///                         start: 0, end: 1, length: 1, qual: 4,
///                         reference_start: Some(9), reference_end: Some(9),
///                         reference_length: Some(0), extra_columns: None,
///                     },
///                     FiberAnnotation {
///                         start: 3, end: 4, length: 1, qual: 7,
///                         reference_start: Some(12), reference_end: Some(12),
///                         reference_length: Some(0), extra_columns: None,
///                     },
///                     FiberAnnotation {
///                         start: 4, end: 5, length: 1, qual: 9,
///                         reference_start: Some(13), reference_end: Some(13),
///                         reference_length: Some(0), extra_columns: None,
///                     },
///                     FiberAnnotation {
///                         start: 7, end: 8, length: 1, qual: 6,
///                         reference_start: Some(16), reference_end: Some(16),
///                         reference_length: Some(0), extra_columns: None,
///                     },
///                 ],
///                 seq_len: 8,
///                 reverse: false,
///             },
///             record_is_reverse: false,
///     }]),
///     2 => assert_eq!(v, vec![BaseMod{
///             modified_base: b'T',
///             strand: '+',
///             modification_type: 'T',
///             ranges: Ranges {
///                 annotations: vec![
///                     FiberAnnotation {
///                         start: 12, end: 13, length: 1, qual: 3,
///                         reference_start: Some(15), reference_end: Some(15),
///                         reference_length: Some(0), extra_columns: None,
///                     },
///                     FiberAnnotation {
///                         start: 13, end: 14, length: 1, qual: 3,
///                         reference_start: Some(16), reference_end: Some(16),
///                         reference_length: Some(0), extra_columns: None,
///                     },
///                     FiberAnnotation {
///                         start: 16, end: 17, length: 1, qual: 4,
///                         reference_start: Some(19), reference_end: Some(19),
///                         reference_length: Some(0), extra_columns: None,
///                     },
///                     FiberAnnotation {
///                         start: 19, end: 20, length: 1, qual: 3,
///                         reference_start: Some(22), reference_end: Some(22),
///                         reference_length: Some(0), extra_columns: None,
///                     },
///                     FiberAnnotation {
///                         start: 20, end: 21, length: 1, qual: 182,
///                         reference_start: Some(23), reference_end: Some(23),
///                         reference_length: Some(0), extra_columns: None,
///                     },
///                 ],
///                 seq_len: 33,
///                 reverse: true,
///             },
///             record_is_reverse: true,
///     }]),
///     _ => {},
///     }
///     count = count + 1;
/// }
/// # Ok::<(), Error>(())
/// ```
///
/// # Errors
/// If MM/ML BAM tags are malformed, you will get `InvalidModProbs` or `InvalidModCoords`.
/// Most integer overflows are dealt with using `except`, except one which gives `ArithmeticError`.
/// `InvalidDuplicates` occurs if the same tag, strand combination occurs many times.
/// Please read the function documentation above as well, which explains some scenarios where
/// even valid tags can be marked as malformed.
///
#[expect(clippy::too_many_lines, reason = "Complex BAM MM/ML tag parsing logic")]
#[expect(
    clippy::missing_panics_doc,
    reason = "either impossible scenarios or integer overflows \
which are also unlikely as genomic coordinates are much less than ~2^63"
)]
pub fn nanalogue_mm_ml_parser<F, G, H>(
    record: &bam::Record,
    filter_mod_prob: F,
    filter_mod_pos: G,
    filter_mod_base_strand_tag: H,
    min_qual: u8,
) -> Result<BaseMods, Error>
where
    F: Fn(&u8) -> bool,
    G: Fn(&usize) -> bool,
    H: Fn(&u8, &char, &ModChar) -> bool,
{
    // Regular expression for matching modification data in the MM tag
    static MM_RE: LazyLock<Regex> = LazyLock::new(|| {
        Regex::new(r"((([ACGTUN])([-+])([A-Za-z]+|[0-9]+)([.?]?))((,[0-9]+)*;)*)")
            .expect("no error")
    });
    // Array to store all the different modifications within the MM tag
    let mut rtn: Vec<BaseMod> = Vec::new();

    let ml_tag: Vec<u8> = get_u8_tag(record, b"ML");

    let mut num_mods_seen: usize = 0;

    let is_reverse = record.is_reverse();

    // if there is an MM tag iterate over all the regex matches
    if let Ok(Aux::String(mm_text)) = record.aux(b"MM") {
        // base qualities; must reverse if rev comp.
        // NOTE this is always equal to number of bases in the sequence, otherwise
        // rust_htslib will throw an error, so we don't have to check this.
        let base_qual: Vec<u8> = match (min_qual, is_reverse) {
            (0, _) => Vec::new(),
            (_, true) => record.qual().iter().rev().copied().collect(),
            (_, false) => record.qual().to_vec(),
        };

        // get forward sequence bases from the bam record
        let forward_bases = {
            let seq = convert_seq_uppercase(record.seq().as_bytes());
            if is_reverse { revcomp(seq) } else { seq }
        };

        let seq_len = forward_bases.len();

        let pos_map = {
            let temp: Vec<Option<i64>> = {
                if record.is_unmapped() {
                    std::iter::repeat_n(None, seq_len).collect()
                } else {
                    record
                        .aligned_pairs_full()
                        .filter(|x| x[0].is_some())
                        .map(|x| x[1])
                        .collect()
                }
            };
            if temp.len() == seq_len {
                temp
            } else {
                return Err(Error::InvalidState(format!(
                    "rust_htslib failure! seq coordinates malformed {} != {}",
                    temp.len(),
                    seq_len
                )));
            }
        };

        for cap in MM_RE.captures_iter(mm_text) {
            let mod_base = cap
                .get(3)
                .map(|m| {
                    *m.as_str()
                        .as_bytes()
                        .first()
                        .expect("regex match guaranteed to have at least 1 byte")
                })
                .expect("no error");
            let mod_strand = cap
                .get(4)
                .map_or("", |m| m.as_str())
                .chars()
                .next()
                .expect("no error");

            // get modification type
            let modification_type: ModChar = cap
                .get(5)
                .map_or("", |m| m.as_str())
                .parse()
                .expect("no error");

            let is_implicit = match cap.get(6).map_or("", |m| m.as_str()).as_bytes() {
                b"" | b"." => true,
                b"?" => false,
                _ => unreachable!("our regex expression must have prevented other possibilities"),
            };
            let mod_dists_str = cap.get(7).map_or("", |m| m.as_str());
            // parse the string containing distances between modifications into a vector of i64

            let mod_dists: Vec<usize> = mod_dists_str
                .trim_end_matches(';')
                .split(',')
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(|s| s.parse().unwrap())
                .collect();

            // do we include bases with zero probabilities?
            let is_include_zero_prob = filter_mod_prob(&0);

            // find real positions in the forward sequence
            let mut cur_mod_idx: usize = 0;
            let mut dist_from_last_mod_base: usize = 0;

            // declare vectors with an approximate with_capacity
            let (mut modified_positions, mut modified_probabilities) = {
                #[expect(
                    clippy::arithmetic_side_effects,
                    reason = "in rare chance of overflow, vectors are initially missized but will be resized anyway \
as they are populated. This will result in a small performance hit but we are ok as this will probably never happen \
when usize is 64-bit as genomic sequences are not that long"
                )]
                let mod_data_len_approx = if is_implicit {
                    // In implicit mode, there may be any number of bases
                    // after the MM data is over, which must be assumed as unmodified.
                    // So we cannot know the exact length of the data before actually
                    // parsing it, and this is just a lower bound of the length.
                    mod_dists.len() + mod_dists.iter().sum::<usize>()
                } else {
                    mod_dists.len()
                };
                (
                    Vec::<usize>::with_capacity(mod_data_len_approx),
                    Vec::<u8>::with_capacity(mod_data_len_approx),
                )
            };

            #[expect(
                clippy::arithmetic_side_effects,
                reason = "one counter is checked for overflow and the other is incremented only when below a ceiling"
            )]
            for (cur_seq_idx, &_) in forward_bases
                .iter()
                .enumerate()
                .filter(|&(_, &k)| mod_base == b'N' || k == mod_base)
            {
                let is_seq_pos_pass: bool = filter_mod_pos(&cur_seq_idx)
                    && (min_qual > 0).then(|| {
                        base_qual
                            .get(cur_seq_idx)
                            .filter(|&x| *x >= min_qual && *x != 255u8)
                            .is_some()
                    }) != Some(false);
                if cur_mod_idx < mod_dists.len()
                    && dist_from_last_mod_base
                        == *mod_dists
                            .get(cur_mod_idx)
                            .expect("cur_mod_idx < mod_dists.len()")
                {
                    let prob = ml_tag
                        .get(cur_mod_idx..)
                        .expect("cur_mod_idx < mod_dists.len() and ml_tag has same length")
                        .get(num_mods_seen)
                        .ok_or(Error::InvalidModProbs(
                            "ML tag appears to be insufficiently long!".into(),
                        ))?;
                    if filter_mod_prob(prob) && is_seq_pos_pass {
                        modified_positions.push(cur_seq_idx);
                        modified_probabilities.push(*prob);
                    }
                    dist_from_last_mod_base = 0;
                    cur_mod_idx += 1;
                } else if cur_mod_idx < mod_dists.len()
                    && dist_from_last_mod_base
                        > *mod_dists
                            .get(cur_mod_idx)
                            .expect("cur_mod_idx < mod_dists.len()")
                {
                    return Err(Error::InvalidModCoords(String::from(
                        "Problem with parsing distances in MM/ML data",
                    )));
                } else {
                    if is_include_zero_prob && is_implicit && is_seq_pos_pass {
                        modified_positions.push(cur_seq_idx);
                        modified_probabilities.push(0);
                    }
                    dist_from_last_mod_base += 1;
                }
            }

            if cur_mod_idx == mod_dists.len() {
                num_mods_seen = num_mods_seen
                    .checked_add(cur_mod_idx)
                    .ok_or(Error::Arithmetic("in MM ML parsing".to_owned()))?;
            } else {
                return Err(Error::InvalidModCoords(String::from(
                    "Problem with parsing MM/ML data, counts do not match",
                )));
            }

            // if data matches filters, add to struct.
            modified_positions.shrink_to(0);
            modified_probabilities.shrink_to(0);

            #[expect(
                clippy::arithmetic_side_effects,
                reason = "`seq_len - 1 - k` (protected as mod pos cannot exceed seq_len), \
`k.0 + 1` (overflow unlikely as genomic coords << 2^63)"
            )]
            #[expect(
                clippy::indexing_slicing,
                reason = "`pos_map[k.0]`; neither pos_map's len nor mod pos entry can exceed seq_len"
            )]
            if filter_mod_base_strand_tag(&mod_base, &mod_strand, &modification_type) {
                let annotations: Vec<FiberAnnotation> = modified_positions
                    .iter()
                    .map(|k| if is_reverse { seq_len - 1 - k } else { *k })
                    .zip(modified_probabilities.iter())
                    .map(|k| {
                        Ok(FiberAnnotation {
                            start: i64::try_from(k.0)?,
                            end: i64::try_from(k.0 + 1)?,
                            length: 1,
                            qual: *k.1,
                            reference_start: pos_map[k.0],
                            reference_end: pos_map[k.0],
                            reference_length: pos_map[k.0].is_some().then_some(0),
                            extra_columns: None,
                        })
                    })
                    .collect::<Result<Vec<FiberAnnotation>, Error>>()?;
                let mods = BaseMod {
                    modified_base: mod_base,
                    strand: mod_strand,
                    modification_type: modification_type.val(),
                    ranges: Ranges {
                        annotations: if is_reverse {
                            annotations.into_iter().rev().collect()
                        } else {
                            annotations
                        },
                        seq_len: i64::try_from(seq_len)?,
                        reverse: is_reverse,
                    },
                    record_is_reverse: is_reverse,
                };
                rtn.push(mods);
            }
        }
    }

    if num_mods_seen == ml_tag.len() {
        // needed so I can compare methods
        rtn.sort();

        // Check for duplicate strand, modification_type combinations
        let mut seen_combinations = HashSet::new();
        for base_mod in &rtn {
            let combination = (base_mod.strand, base_mod.modification_type);
            if seen_combinations.contains(&combination) {
                return Err(Error::InvalidDuplicates(format!(
                    "Duplicate strand '{}' and modification_type '{}' combination found",
                    base_mod.strand, base_mod.modification_type
                )));
            }
            let _: bool = seen_combinations.insert(combination);
        }

        Ok(BaseMods { base_mods: rtn })
    } else {
        Err(Error::InvalidModProbs(
            "MM and ML tag lengths do not match!".to_owned(),
        ))
    }
}

/// A global struct which contains BAM records for further usage.
/// NOTE: we don't derive many traits here as the `RcRecords` object
/// does not have many traits.
///
/// ```
/// use nanalogue_core::{Error, nanalogue_bam_reader, BamRcRecords, InputBam, InputMods};
/// use rust_htslib::bam::Read;
/// let mut reader = nanalogue_bam_reader(&"examples/example_1.bam")?;
/// let BamRcRecords = BamRcRecords::new(&mut reader, &mut InputBam::default(),
///     &mut InputMods::default())?;
/// assert_eq!(BamRcRecords.header.tid(b"dummyI"), Some(0));
/// assert_eq!(BamRcRecords.header.tid(b"dummyII"), Some(1));
/// assert_eq!(BamRcRecords.header.tid(b"dummyIII"), Some(2));
/// # Ok::<(), Error>(())
/// ```
#[derive(Debug)]
#[non_exhaustive]
pub struct BamRcRecords<'a, R>
where
    R: bam::Read,
{
    /// `RcRecords` object output by rust htslib which we can iterate over
    pub rc_records: bam::RcRecords<'a, R>,
    /// Header of the bam file
    pub header: bam::HeaderView,
}

impl<'a, R: bam::Read> BamRcRecords<'a, R> {
    /// Extracts `RcRecords` from a BAM Reader
    ///
    /// # Errors
    /// Returns an error if thread pool creation, BAM region fetching/processing,
    /// or read ID file processing fails.
    pub fn new<T: InputRegionOptions>(
        bam_reader: &'a mut R,
        bam_opts: &mut InputBam,
        mod_opts: &mut T,
    ) -> Result<Self, Error> {
        let header = bam_reader.header().clone();
        let tp = tpool::ThreadPool::new(bam_opts.threads.get())?;
        bam_reader.set_thread_pool(&tp)?;

        // Load read ID list if specified
        bam_opts.read_id_set = if let Some(file_path) = bam_opts.read_id_list.as_ref() {
            let file = File::open(file_path).map_err(Error::InputOutputError)?;
            let reader = BufReader::new(file);
            let mut read_ids = HashSet::new();

            for raw_line in reader.lines() {
                let temp_line = raw_line.map_err(Error::InputOutputError)?;
                let line = temp_line.trim();
                if !line.is_empty() && !line.starts_with('#') {
                    let _: bool = read_ids.insert(line.to_string());
                }
            }
            Some(read_ids)
        } else {
            None
        };

        bam_opts.convert_region_to_bed3(header.clone())?;
        mod_opts.convert_region_to_bed3(header.clone())?;
        let rc_records = bam_reader.rc_records();
        Ok(BamRcRecords::<R> { rc_records, header })
    }
}

/// Trait that performs filtration
pub trait BamPreFilt {
    /// apply default filtration
    fn pre_filt(&self, _bam_opts: &InputBam) -> bool {
        unimplemented!()
    }
    /// filtration by length
    fn filt_by_len(&self, _min_seq_len: u64, _include_zero_len: bool) -> bool {
        unimplemented!()
    }
    /// filtration by alignment length
    fn filt_by_align_len(&self, _min_align_len: i64) -> bool {
        unimplemented!()
    }
    /// filtration by read id
    fn filt_by_read_id(&self, _read_id: &str) -> bool {
        unimplemented!()
    }
    /// filtration by read id set
    fn filt_by_read_id_set(&self, _read_id_set: &HashSet<String>) -> bool {
        unimplemented!()
    }
    /// filtration using flags
    fn filt_by_bitwise_or_flags(&self, _states: &ReadStates) -> bool {
        unimplemented!()
    }
    /// random filtration
    fn filt_random_subset(&self, _fraction: F32Bw0and1) -> bool {
        unimplemented!()
    }
    /// filtration by mapq
    fn filt_by_mapq(&self, _min_mapq: u8, _exclude_mapq_unavail: bool) -> bool {
        unimplemented!()
    }
    /// filtration by region
    fn filt_by_region(&self, _region: &Bed3<i32, u64>, _full_region: bool) -> bool {
        unimplemented!()
    }
}

/// Trait that performs filtration on `rust_htslib` Record
///
/// Filter in action below, 100 bp read passed with a 20 bp filter
/// but fails with a 120 bp filter
/// ```
/// use nanalogue_core::{BamPreFilt, InputBam, Error};
/// use rust_htslib::bam::record;
/// use rust_htslib::bam::record::{Cigar, CigarString};
/// let mut bam_record = record::Record::new();
/// bam_record.set(&vec![b'r',b'e',b'a',b'd'], Some(&CigarString(vec![Cigar::Match(100)])),
///       &vec![ b'A' as u8; 100], &vec![255 as u8; 100]);
/// let mut bam_opts = InputBam::default();
/// bam_opts.min_seq_len = 20;
/// assert_eq!(bam_record.pre_filt(&bam_opts), true);
/// bam_opts.min_seq_len = 120;
/// assert_eq!(bam_record.pre_filt(&bam_opts), false);
/// # Ok::<(), Error>(())
/// ```
/// By default, zero-length records are excluded to avoid processing errors.
/// ```
/// # use nanalogue_core::{BamPreFilt, InputBam, Error};
/// # use rust_htslib::bam::record;
/// # use rust_htslib::bam::record::{Cigar, CigarString};
/// let mut bam_record = record::Record::new();
/// let mut bam_opts = InputBam::default();
/// bam_opts.min_seq_len = 20;
/// bam_opts.include_zero_len = false;  // default behavior
/// assert_eq!(bam_record.pre_filt(&bam_opts), false);  // excluded
/// # Ok::<(), Error>(())
/// ```
/// Zero-length records can be included if explicitly requested.
/// ```
/// # use nanalogue_core::{BamPreFilt, InputBam, Error};
/// # use rust_htslib::bam::record;
/// # use rust_htslib::bam::record::{Cigar, CigarString};
/// let mut bam_record = record::Record::new();
/// let mut bam_opts = InputBam::default();
/// bam_opts.min_seq_len = 20;
/// bam_opts.include_zero_len = true;
/// assert_eq!(bam_record.pre_filt(&bam_opts), true);  // included
/// # Ok::<(), Error>(())
/// ```
impl BamPreFilt for bam::Record {
    /// apply default filtration by read length
    fn pre_filt(&self, bam_opts: &InputBam) -> bool {
        self.filt_by_len(bam_opts.min_seq_len, bam_opts.include_zero_len)
            & self.filt_random_subset(bam_opts.sample_fraction)
            & self.filt_by_mapq(bam_opts.mapq_filter, bam_opts.exclude_mapq_unavail)
            & {
                if let Some(v) = bam_opts.read_id.as_ref() {
                    self.filt_by_read_id(v)
                } else if let Some(read_id_set) = bam_opts.read_id_set.as_ref() {
                    self.filt_by_read_id_set(read_id_set)
                } else {
                    true
                }
            }
            & {
                if let Some(v) = bam_opts.min_align_len {
                    self.filt_by_align_len(v)
                } else {
                    true
                }
            }
            & {
                if let Some(v) = bam_opts.read_filter.as_ref() {
                    self.filt_by_bitwise_or_flags(v)
                } else {
                    true
                }
            }
            & {
                if let Some(v) = bam_opts.region_filter().as_ref() {
                    self.filt_by_region(v, bam_opts.is_full_overlap())
                } else {
                    true
                }
            }
    }
    /// filtration by read length
    fn filt_by_len(&self, min_seq_len: u64, include_zero_len: bool) -> bool {
        // self.len() returns a usize which we convert to u64
        match (min_seq_len, self.len() as u64, include_zero_len) {
            (_, 0, false) => false, // Exclude zero-length by default
            (_, 0, true) => true,   // Include zero-length when explicitly requested
            (l_min, v, _) => v >= l_min,
        }
    }
    /// filtration by alignment length
    fn filt_by_align_len(&self, min_align_len: i64) -> bool {
        !self.is_unmapped() && {
            let ref_end = self.reference_end();
            let ref_start = self.pos();
            ref_end >= ref_start
                && ref_start >= 0
                && ref_end
                    .checked_sub(ref_start)
                    .expect("ref_end >= ref_start so overflow is impossible")
                    >= min_align_len
        }
    }
    /// filtration by read id
    fn filt_by_read_id(&self, read_id: &str) -> bool {
        read_id.as_bytes() == self.name()
    }
    /// filtration by read id set
    fn filt_by_read_id_set(&self, read_id_set: &HashSet<String>) -> bool {
        if let Ok(name_str) = std::str::from_utf8(self.name()) {
            read_id_set.contains(name_str)
        } else {
            false
        }
    }
    /// filtration by flag list
    fn filt_by_bitwise_or_flags(&self, states: &ReadStates) -> bool {
        states.bam_flags().contains(&self.flags())
    }
    /// random filtration
    fn filt_random_subset(&self, fraction: F32Bw0and1) -> bool {
        match fraction.val() {
            1.0 => true,
            0.0 => false,
            v => {
                let random_number: f32 = random();
                random_number < v
            }
        }
    }
    /// filtration by mapq
    fn filt_by_mapq(&self, min_mapq: u8, exclude_mapq_unavail: bool) -> bool {
        !(exclude_mapq_unavail && self.mapq() == 255) && self.mapq() >= min_mapq
    }
    /// filtration by region
    fn filt_by_region(&self, region: &Bed3<i32, u64>, full_region: bool) -> bool {
        !self.is_unmapped() && (self.tid() == *region.chr()) && {
            let region_start = region.start();
            let region_end = region.end();
            (region_start == u64::MIN && region_end == u64::MAX) || {
                let start: u64 = self.pos().try_into().expect("no error");
                let end: u64 = self.reference_end().try_into().expect("no error");
                if full_region {
                    (start..end).contains(&region_start) && (start..=end).contains(&region_end)
                } else {
                    (start..end).intersects(&(region_start..region_end))
                }
            }
        }
    }
}

#[cfg(test)]
mod mod_parse_tests {
    use super::*;
    use rust_htslib::bam::Read as _;

    /// Tests if Mod BAM modification parsing is alright,
    /// some of the test cases here may be a repeat of the doctest above.
    #[test]
    #[expect(clippy::too_many_lines, reason = "Comprehensive integration test")]
    fn mod_bam_parsing_from_example_1_bam() -> Result<(), Error> {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam")?;
        for (count, record) in reader.records().enumerate() {
            let r = record?;
            let BaseMods { base_mods: v } =
                nanalogue_mm_ml_parser(&r, |&_| true, |&_| true, |&_, &_, &_| true, 0).unwrap();
            match count {
                0 => assert_eq!(
                    v,
                    vec![BaseMod {
                        modified_base: b'T',
                        strand: '+',
                        modification_type: 'T',
                        ranges: Ranges {
                            annotations: vec![
                                FiberAnnotation {
                                    start: 0,
                                    end: 1,
                                    length: 1,
                                    qual: 4,
                                    reference_start: Some(9),
                                    reference_end: Some(9),
                                    reference_length: Some(0),
                                    extra_columns: None,
                                },
                                FiberAnnotation {
                                    start: 3,
                                    end: 4,
                                    length: 1,
                                    qual: 7,
                                    reference_start: Some(12),
                                    reference_end: Some(12),
                                    reference_length: Some(0),
                                    extra_columns: None,
                                },
                                FiberAnnotation {
                                    start: 4,
                                    end: 5,
                                    length: 1,
                                    qual: 9,
                                    reference_start: Some(13),
                                    reference_end: Some(13),
                                    reference_length: Some(0),
                                    extra_columns: None,
                                },
                                FiberAnnotation {
                                    start: 7,
                                    end: 8,
                                    length: 1,
                                    qual: 6,
                                    reference_start: Some(16),
                                    reference_end: Some(16),
                                    reference_length: Some(0),
                                    extra_columns: None,
                                },
                            ],
                            seq_len: 8,
                            reverse: false,
                        },
                        record_is_reverse: false,
                    }]
                ),
                1 => assert_eq!(
                    v,
                    vec![BaseMod {
                        modified_base: b'T',
                        strand: '+',
                        modification_type: 'T',
                        ranges: Ranges {
                            annotations: vec![
                                FiberAnnotation {
                                    start: 3,
                                    end: 4,
                                    length: 1,
                                    qual: 221,
                                    reference_start: Some(26),
                                    reference_end: Some(26),
                                    reference_length: Some(0),
                                    extra_columns: None,
                                },
                                FiberAnnotation {
                                    start: 8,
                                    end: 9,
                                    length: 1,
                                    qual: 242,
                                    reference_start: Some(31),
                                    reference_end: Some(31),
                                    reference_length: Some(0),
                                    extra_columns: None,
                                },
                                FiberAnnotation {
                                    start: 27,
                                    end: 28,
                                    length: 1,
                                    qual: 3,
                                    reference_start: Some(50),
                                    reference_end: Some(50),
                                    reference_length: Some(0),
                                    extra_columns: None,
                                },
                                FiberAnnotation {
                                    start: 39,
                                    end: 40,
                                    length: 1,
                                    qual: 47,
                                    reference_start: Some(62),
                                    reference_end: Some(62),
                                    reference_length: Some(0),
                                    extra_columns: None,
                                },
                                FiberAnnotation {
                                    start: 47,
                                    end: 48,
                                    length: 1,
                                    qual: 239,
                                    reference_start: Some(70),
                                    reference_end: Some(70),
                                    reference_length: Some(0),
                                    extra_columns: None,
                                },
                            ],
                            seq_len: 48,
                            reverse: false,
                        },
                        record_is_reverse: false,
                    }]
                ),
                2 => assert_eq!(
                    v,
                    vec![BaseMod {
                        modified_base: b'T',
                        strand: '+',
                        modification_type: 'T',
                        ranges: Ranges {
                            annotations: vec![
                                FiberAnnotation {
                                    start: 12,
                                    end: 13,
                                    length: 1,
                                    qual: 3,
                                    reference_start: Some(15),
                                    reference_end: Some(15),
                                    reference_length: Some(0),
                                    extra_columns: None,
                                },
                                FiberAnnotation {
                                    start: 13,
                                    end: 14,
                                    length: 1,
                                    qual: 3,
                                    reference_start: Some(16),
                                    reference_end: Some(16),
                                    reference_length: Some(0),
                                    extra_columns: None,
                                },
                                FiberAnnotation {
                                    start: 16,
                                    end: 17,
                                    length: 1,
                                    qual: 4,
                                    reference_start: Some(19),
                                    reference_end: Some(19),
                                    reference_length: Some(0),
                                    extra_columns: None,
                                },
                                FiberAnnotation {
                                    start: 19,
                                    end: 20,
                                    length: 1,
                                    qual: 3,
                                    reference_start: Some(22),
                                    reference_end: Some(22),
                                    reference_length: Some(0),
                                    extra_columns: None,
                                },
                                FiberAnnotation {
                                    start: 20,
                                    end: 21,
                                    length: 1,
                                    qual: 182,
                                    reference_start: Some(23),
                                    reference_end: Some(23),
                                    reference_length: Some(0),
                                    extra_columns: None,
                                },
                            ],
                            seq_len: 33,
                            reverse: true,
                        },
                        record_is_reverse: true,
                    }]
                ),
                3 => assert_eq!(
                    v,
                    vec![
                        BaseMod {
                            modified_base: b'G',
                            strand: '-',
                            modification_type: '\u{1C20}',
                            ranges: Ranges {
                                annotations: vec![
                                    FiberAnnotation {
                                        start: 28,
                                        end: 29,
                                        length: 1,
                                        qual: 0,
                                        reference_start: None,
                                        reference_end: None,
                                        reference_length: None,
                                        extra_columns: None,
                                    },
                                    FiberAnnotation {
                                        start: 29,
                                        end: 30,
                                        length: 1,
                                        qual: 0,
                                        reference_start: None,
                                        reference_end: None,
                                        reference_length: None,
                                        extra_columns: None,
                                    },
                                    FiberAnnotation {
                                        start: 30,
                                        end: 31,
                                        length: 1,
                                        qual: 0,
                                        reference_start: None,
                                        reference_end: None,
                                        reference_length: None,
                                        extra_columns: None,
                                    },
                                    FiberAnnotation {
                                        start: 32,
                                        end: 33,
                                        length: 1,
                                        qual: 0,
                                        reference_start: None,
                                        reference_end: None,
                                        reference_length: None,
                                        extra_columns: None,
                                    },
                                    FiberAnnotation {
                                        start: 43,
                                        end: 44,
                                        length: 1,
                                        qual: 77,
                                        reference_start: None,
                                        reference_end: None,
                                        reference_length: None,
                                        extra_columns: None,
                                    },
                                    FiberAnnotation {
                                        start: 44,
                                        end: 45,
                                        length: 1,
                                        qual: 0,
                                        reference_start: None,
                                        reference_end: None,
                                        reference_length: None,
                                        extra_columns: None,
                                    },
                                ],
                                seq_len: 48,
                                reverse: false,
                            },
                            record_is_reverse: false,
                        },
                        BaseMod {
                            modified_base: b'T',
                            strand: '+',
                            modification_type: 'T',
                            ranges: Ranges {
                                annotations: vec![
                                    FiberAnnotation {
                                        start: 3,
                                        end: 4,
                                        length: 1,
                                        qual: 221,
                                        reference_start: None,
                                        reference_end: None,
                                        reference_length: None,
                                        extra_columns: None,
                                    },
                                    FiberAnnotation {
                                        start: 8,
                                        end: 9,
                                        length: 1,
                                        qual: 242,
                                        reference_start: None,
                                        reference_end: None,
                                        reference_length: None,
                                        extra_columns: None,
                                    },
                                    FiberAnnotation {
                                        start: 27,
                                        end: 28,
                                        length: 1,
                                        qual: 0,
                                        reference_start: None,
                                        reference_end: None,
                                        reference_length: None,
                                        extra_columns: None,
                                    },
                                    FiberAnnotation {
                                        start: 39,
                                        end: 40,
                                        length: 1,
                                        qual: 47,
                                        reference_start: None,
                                        reference_end: None,
                                        reference_length: None,
                                        extra_columns: None,
                                    },
                                    FiberAnnotation {
                                        start: 47,
                                        end: 48,
                                        length: 1,
                                        qual: 239,
                                        reference_start: None,
                                        reference_end: None,
                                        reference_length: None,
                                        extra_columns: None,
                                    },
                                ],
                                seq_len: 48,
                                reverse: false,
                            },
                            record_is_reverse: false,
                        }
                    ]
                ),
                _ => {}
            }
        }
        Ok(())
    }

    fn create_test_record_and_parse(mm_value: &str, ml_values: &[u8]) -> Result<BaseMods, Error> {
        let mut record = bam::Record::new();

        // Create a sequence - ATCG repeated
        let seq_bytes = b"ATCGATCGATCGATCGATCG";
        let qname = b"test_read";
        let cigar = bam::record::CigarString::from(vec![bam::record::Cigar::Match(20)]);
        let quals = vec![30u8; 20]; // Quality scores of 30 for all bases
        record.set(qname, Some(&cigar), seq_bytes, &quals);

        // Set MM tag
        record.push_aux(b"MM", Aux::String(mm_value))?;

        // Set ML tag (modification probabilities)
        record.push_aux(b"ML", Aux::ArrayU8((&ml_values).into()))?;

        // Call the parser and return result
        nanalogue_mm_ml_parser(
            &record,
            |&_| true,         // Accept all probabilities
            |&_| true,         // Accept all positions
            |&_, &_, &_| true, // Accept all base/strand/tag combinations
            0,                 // No quality threshold
        )
    }

    #[test]
    #[should_panic(expected = "InvalidDuplicates")]
    fn nanalogue_mm_ml_parser_detects_duplicates() {
        // Two T+ modifications with same strand and modification_type
        let mm_value = "T+T,0,3;T+T,1,2;";
        let ml_values = Vec::from([100u8, 200u8, 150u8, 180u8]);
        let _: BaseMods = create_test_record_and_parse(mm_value, &ml_values).unwrap();
    }

    #[test]
    fn nanalogue_mm_ml_parser_accepts_unique_combinations() {
        // T+, C+, T- (all different combinations)
        let mm_value = "T+T,0,3;C+m,0,1;T-T,1,2;";
        let ml_values = Vec::from([100u8, 200u8, 150u8, 180u8, 120u8, 140u8]);
        let _: BaseMods = create_test_record_and_parse(mm_value, &ml_values).unwrap();
    }

    #[test]
    #[should_panic(expected = "InvalidModCoords")]
    fn nanalogue_mm_ml_parser_detects_invalid_mod_coords_1() {
        // Invalid coordinates: seq does not have 11 As (4 + 1 + 5 + 1)
        let mm_value = "T+T,0,3;A-T,4,5;";
        let ml_values = Vec::from([100u8, 200u8, 150u8, 180u8]);
        let _: BaseMods = create_test_record_and_parse(mm_value, &ml_values).unwrap();
    }

    #[test]
    #[should_panic(expected = "InvalidModCoords")]
    fn nanalogue_mm_ml_parser_detects_invalid_mod_coords_2() {
        // Invalid coords: there aren't 11 C in the sequence
        let mm_value = "T+T,0,3;C+m,4,5;T-T,1,2;";
        let ml_values = Vec::from([100u8, 200u8, 150u8, 180u8, 120u8, 140u8]);
        let _: BaseMods = create_test_record_and_parse(mm_value, &ml_values).unwrap();
    }

    #[test]
    #[should_panic(expected = "InvalidModProbs")]
    fn nanalogue_mm_ml_parser_detects_ml_tag_longer_than_mm_data() {
        // ML tag has more values than modifications in MM tag
        let mm_value = "T+T,0,3;"; // 2 modifications
        let ml_values = Vec::from([100u8, 200u8, 150u8, 180u8]); // 4 values - too many!
        let _: BaseMods = create_test_record_and_parse(mm_value, &ml_values).unwrap();
    }

    #[test]
    #[should_panic(expected = "InvalidModProbs")]
    fn nanalogue_mm_ml_parser_detects_ml_tag_shorter_than_mm_data() {
        // ML tag has fewer values than modifications in MM tag
        let mm_value = "T+T,0,3;"; // 2 modifications
        let ml_values = Vec::from([100u8]); // 1 value - too few!
        let _: BaseMods = create_test_record_and_parse(mm_value, &ml_values).unwrap();
    }
}

#[cfg(test)]
mod zero_length_filtering_tests {
    use super::*;
    use rust_htslib::bam::Read as _;

    #[test]
    fn zero_length_filtering_with_example_2_zero_len_sam() -> Result<(), Error> {
        // Test with include_zero_len = false and min_seq_len = 1 (should get 1 read)
        let mut reader = nanalogue_bam_reader("examples/example_2_zero_len.sam")?;
        let bam_opts_exclude_zero: InputBam =
            serde_json::from_str(r#"{"min_seq_len": 1}"#).unwrap();

        let mut count_exclude_zero = 0;
        for record_result in reader.records() {
            let record = record_result?;
            if record.pre_filt(&bam_opts_exclude_zero) {
                count_exclude_zero += 1;
            }
        }

        // Test with include_zero_len = true and min_seq_len = 1 (should get 2 reads)
        let mut reader_same_file = nanalogue_bam_reader("examples/example_2_zero_len.sam")?;
        let bam_opts_include_zero: InputBam =
            serde_json::from_str(r#"{"min_seq_len": 1, "include_zero_len": true}"#).unwrap();

        let mut count_include_zero = 0;
        for record_result in reader_same_file.records() {
            let record = record_result?;
            if record.pre_filt(&bam_opts_include_zero) {
                count_include_zero += 1;
            }
        }

        assert_eq!(count_exclude_zero, 1);
        assert_eq!(count_include_zero, 2);

        Ok(())
    }
}

#[cfg(test)]
mod invalid_seq_length_tests {
    use super::*;
    use rust_htslib::bam::Read as _;

    #[test]
    fn invalid_seq_length_error_with_example_4() {
        let mut reader =
            nanalogue_bam_reader("examples/example_4_invalid_basequal_len.sam").unwrap();

        let mut cnt = 0;
        for record_result in reader.records() {
            cnt += 1;
            let _: rust_htslib::errors::Error = record_result.unwrap_err();
        }
        assert_eq!(cnt, 1);
    }
}

#[cfg(test)]
mod base_qual_filtering_tests {
    use super::*;
    use rust_htslib::bam::Read as _;

    #[test]
    fn base_qual_filtering_with_example_5_valid_basequal_sam() -> Result<(), Error> {
        let mut reader = nanalogue_bam_reader("examples/example_5_valid_basequal.sam")?;

        for record_result in reader.records() {
            let record = record_result?;
            let BaseMods { base_mods: v } =
                nanalogue_mm_ml_parser(&record, |&_| true, |&_| true, |&_, &_, &_| true, 60)
                    .unwrap();

            assert_eq!(v.len(), 1);
            let base_mod = v.first().expect("v has exactly 1 element");
            assert_eq!(base_mod.modified_base, b'T');
            assert_eq!(base_mod.strand, '+');
            assert_eq!(base_mod.modification_type, 'T');
            assert_eq!(base_mod.ranges.qual(), vec![7, 9]);
            assert_eq!(base_mod.ranges.starts(), vec![3, 4]);
            assert_eq!(base_mod.ranges.ends(), vec![4, 5]);
        }

        Ok(())
    }
}

#[cfg(test)]
mod bam_rc_record_tests {
    use super::*;
    use rand::random_range;
    use rust_htslib::bam::record;
    use rust_htslib::bam::record::{Cigar, CigarString};

    #[test]
    fn bam_rc_records() {
        let mut reader = nanalogue_bam_reader("examples/example_1.bam").unwrap();
        let bam_rc_records = BamRcRecords::new(
            &mut reader,
            &mut InputBam::default(),
            &mut InputMods::default(),
        )
        .unwrap();
        assert_eq!(bam_rc_records.header.tid(b"dummyI"), Some(0));
        assert_eq!(bam_rc_records.header.tid(b"dummyII"), Some(1));
        assert_eq!(bam_rc_records.header.tid(b"dummyIII"), Some(2));
    }

    #[test]
    fn example_3_read_list_1() {
        // Test with example_3_subset_1 - should see 2 reads
        let mut reader = nanalogue_bam_reader("examples/example_3.bam").unwrap();
        let mut bam_opts: InputBam =
            serde_json::from_str(r#"{"read_id_list": "examples/example_3_subset_1"}"#).unwrap();

        let bam_rc_records =
            BamRcRecords::new(&mut reader, &mut bam_opts, &mut InputMods::default()).unwrap();

        // Count lines in the read ID file
        let file = File::open("examples/example_3_subset_1").unwrap();
        let reader_file = BufReader::new(file);
        let mut line_count = 0;
        for raw_line in reader_file.lines() {
            let temp_line = raw_line.unwrap();
            let line = temp_line.trim();
            if !line.is_empty() && !line.starts_with('#') {
                line_count += 1;
            }
        }
        assert_eq!(line_count, 2);

        let mut count = 0;
        for record_result in bam_rc_records.rc_records {
            let record = record_result.unwrap();
            if record.pre_filt(&bam_opts) {
                count += 1;
            }
        }
        assert_eq!(count, 2);
    }

    #[test]
    fn example_3_read_list_2() {
        // Test with example_3_subset_w_invalid - should see 2 reads
        // even though file contains 3 read IDs
        let mut reader = nanalogue_bam_reader("examples/example_3.bam").unwrap();
        let mut bam_opts: InputBam =
            serde_json::from_str(r#"{"read_id_list": "examples/example_3_subset_w_invalid"}"#)
                .unwrap();
        let bam_rc_records =
            BamRcRecords::new(&mut reader, &mut bam_opts, &mut InputMods::default()).unwrap();

        // Count lines in the read ID file
        let file = File::open("examples/example_3_subset_w_invalid").unwrap();
        let reader_file = BufReader::new(file);
        let mut line_count = 0;
        for raw_line in reader_file.lines() {
            let temp_line = raw_line.unwrap();
            let line = temp_line.trim();
            if !line.is_empty() && !line.starts_with('#') {
                line_count += 1;
            }
        }
        assert_eq!(line_count, 3);

        let mut count = 0;
        for record_result in bam_rc_records.rc_records {
            let record = record_result.unwrap();
            if record.pre_filt(&bam_opts) {
                count += 1;
            }
        }
        assert_eq!(count, 2);
    }

    #[test]
    fn random_retrieval() {
        let mut count_retained = 0;
        let bam_opts: InputBam =
            serde_json::from_str(r#"{"sample_fraction": 0.5, "include_zero_len": true}"#).unwrap();

        for _ in 0..10000 {
            let record = record::Record::new();
            if record.pre_filt(&bam_opts) {
                count_retained += 1;
            }
        }

        // 50% retention rate => 5000 reads, so we test if 4500-5500 reads come through
        // (this is quite lax actually)
        assert!((4500..=5500).contains(&count_retained));
    }

    #[test]
    fn single_read_id_filtering() {
        let mut count_retained = 0;

        let bam_opts: InputBam =
            serde_json::from_str(r#"{"read_id": "read2", "include_zero_len": true}"#).unwrap();

        for i in 1..=1000 {
            let mut record = record::Record::new();
            let read_id = format!("read{i}");
            record.set_qname(read_id.as_bytes());

            if record.pre_filt(&bam_opts) {
                count_retained += 1;
            }
        }

        assert_eq!(count_retained, 1);
    }

    #[test]
    fn seq_and_align_len_filtering() {
        let bam_opts_min_len: InputBam = serde_json::from_str(r#"{"min_seq_len": 5000}"#).unwrap();
        let bam_opts_min_align_len: InputBam =
            serde_json::from_str(r#"{"min_align_len": 2500}"#).unwrap();

        let mut count_retained = (0, 0);

        for _ in 1..=10000 {
            let mut record = record::Record::new();
            let seq_len: usize = random_range(2..=10000);
            let match_len = seq_len / 2;
            let hard_clip_len = seq_len - match_len;

            record.set(
                b"read",
                Some(&CigarString(vec![
                    Cigar::Match(u32::try_from(match_len).unwrap()),
                    Cigar::HardClip(u32::try_from(hard_clip_len).unwrap()),
                ])),
                &vec![b'A'; seq_len],
                &vec![50; seq_len],
            );
            record.unset_flags();
            record.set_flags(loop {
                let random_state: ReadState = random();
                match random_state {
                    ReadState::Unmapped => {}
                    v @ (ReadState::PrimaryFwd
                    | ReadState::PrimaryRev
                    | ReadState::SecondaryFwd
                    | ReadState::SecondaryRev
                    | ReadState::SupplementaryFwd
                    | ReadState::SupplementaryRev) => break u16::from(v),
                }
            });
            if record.pre_filt(&bam_opts_min_align_len) {
                count_retained.0 += 1;
            }
            if record.pre_filt(&bam_opts_min_len) {
                count_retained.1 += 1;
            }
        }

        assert!(count_retained.0 >= 4500 && count_retained.0 <= 5500);
        assert!(count_retained.1 >= 4500 && count_retained.1 <= 5500);
    }

    #[test]
    fn filt_by_bitwise_or_flags() {
        // Create random subset of read states (ensure at least one)
        let selected_states = {
            let mut selected_states = Vec::new();
            let all_states = vec!["0", "4", "16", "256", "272", "2048", "2064"];
            for state in &all_states {
                if random::<f32>() < 0.5 {
                    selected_states.push(*state);
                }
            }
            if selected_states.is_empty() {
                let random_idx = random_range(0..all_states.len());
                selected_states.push(
                    *all_states
                        .get(random_idx)
                        .expect("random_idx is within all_states range"),
                );
            }
            selected_states
        };

        let states_string = selected_states.join(",");
        let bam_opts: InputBam = serde_json::from_str(
            format!("{{\"read_filter\": [{states_string}], \"include_zero_len\": true}}").as_str(),
        )
        .unwrap();

        let mut count_retained = 0;

        for _ in 1..=70000 {
            let mut record = record::Record::new();
            record.unset_flags();
            record.set_flags({
                let random_state: ReadState = random();
                u16::from(random_state)
            });
            if record.pre_filt(&bam_opts) {
                count_retained += 1;
            }
        }

        let expected_count = selected_states.len() * 10000;
        let tolerance = expected_count.div_ceil(5);
        let min_count = expected_count.saturating_sub(tolerance);
        let max_count = expected_count + tolerance;

        assert!(count_retained >= min_count && count_retained <= max_count);
    }

    #[test]
    fn filt_by_region() {
        let mut count_retained = (0, 0);

        for _ in 1..=10000 {
            let mut record = record::Record::new();

            // Draw four numbers from 0 to 10000
            let four_nums = [
                random_range(0..=10000),
                random_range(0..=10000),
                random_range(0..=10000),
                random_range(0..=10000u64),
            ];

            // First two for record
            let mut pos_nums = [four_nums[0], four_nums[1]];
            pos_nums.sort_unstable();
            let start_pos = pos_nums[0];
            let end_pos = pos_nums[1];

            // Use the last two to set up region, with a random contig chosen
            // from a set of two.
            let region_tid = i32::from(random::<bool>());
            let mut region_nums = [four_nums[2], four_nums[3]];
            region_nums.sort_unstable();
            let region_start = region_nums[0];
            let region_end = region_nums[1];
            let bam_opts_no_full_region: InputBam = serde_json::from_str(
                format!("{{\"region_bed3\": [{region_tid}, {region_start}, {region_end}]}}")
                    .as_str(),
            )
            .unwrap();
            let bam_opts_full_region: InputBam = serde_json::from_str(
                format!("{{\"full_region\": true, \"region_bed3\": [{region_tid}, {region_start}, {region_end}]}}")
                    .as_str(),
            )
            .unwrap();

            // Calculate sequence length (ensure at least 1)
            let seq_len = if end_pos == start_pos {
                1
            } else {
                end_pos - start_pos
            };

            // Set sequence details first
            record.set(
                b"read",
                Some(&CigarString(vec![Cigar::Match(
                    u32::try_from(seq_len).unwrap(),
                )])),
                &vec![b'A'; usize::try_from(seq_len).unwrap()],
                &vec![255; usize::try_from(seq_len).unwrap()],
            );

            // Set tid and position, contig chosen from a random set of two.
            let tid = i32::from(random::<bool>());
            record.set_tid(tid);
            record.set_pos(i64::try_from(start_pos).unwrap());

            // Verify reference_end calculation
            let expected_ref_end = i64::try_from(start_pos + seq_len).unwrap();
            assert_eq!(record.reference_end(), expected_ref_end);

            if record.pre_filt(&bam_opts_no_full_region) {
                count_retained.0 += 1;
            }
            if record.pre_filt(&bam_opts_full_region) {
                count_retained.1 += 1;
            }
        }

        // the chance that two randomly chosen intervals on the interval [0, l]
        // intersect is 2/3, and then we have the added random element of one
        // of two contigs, so the probability is 2/3 * 1/2 = 1/3.
        let expected_count = 10000usize.div_ceil(3); // ~3333
        let tolerance = expected_count.div_ceil(5);
        let min_count = expected_count.saturating_sub(tolerance);
        let max_count = expected_count + tolerance;
        assert!(count_retained.0 >= min_count && count_retained.0 <= max_count);

        // the chance that two randomly chosen intervals are such that one is
        // contained within the other is 1/3, and we are looking for the region
        // being contained completely within the read and NOT the read being contained
        // within the region, so the probability is 1/6. This is one fourth of the
        // 2/3 used above. So the same criterion will work with the second count quadrupled.
        assert!(4 * count_retained.1 >= min_count && 4 * count_retained.1 <= max_count);
    }
}
