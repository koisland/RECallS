use std::{collections::HashMap, ops::Bound, path::Path};

use eyre::{ContextCompat, bail};
use itertools::Itertools;
use noodles::{
    bam::{self},
    core::{Position, Region},
    sam::alignment::{
        Record,
        record::{
            Cigar, Flags,
            cigar::{Op, op::Kind},
            data::field::Value,
        },
    },
};
use rust_lapper::{Interval, Lapper};

use crate::{
    baseline::ReadIndelSummaryStats,
    events::{Event, MismatchSignal, SV, SVType, Signal},
    unbalanced_aln::is_unbalanced_alignment,
    utils::{Direction, overlap_length},
};

const CLIP_OVL_DST: usize = 2000;

/// Look-up table for nts
/// https://github.com/Daniel-Liu-c0deb0t/frfr/blob/e07781cd4f43877a3aba8d0bc42217b4419a8a67/src/lib.rs#L249-L260
static LUT: [u8; 128] = {
    let mut l = [0u8; 128];
    l[b'A' as usize] = 0b00;
    l[b'a' as usize] = 0b00;
    l[b'C' as usize] = 0b01;
    l[b'c' as usize] = 0b01;
    l[b'G' as usize] = 0b10;
    l[b'g' as usize] = 0b10;
    l[b'T' as usize] = 0b11;
    l[b't' as usize] = 0b11;
    l
};

/// Order of nucleotides in [`NtCounts`]
static NT_ORDER: [u8; 4] = *b"ACGT";

/// Nucleotide counts for each base position
/// * Order is [`NT_ORDER`]
pub(crate) struct NtCounts(pub(crate) [usize; 4]);

impl NtCounts {
    /// Top nucleotide and its count for a given position.
    #[inline]
    fn top_nt(&self) -> (u8, usize) {
        self.0
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.cmp(b.1))
            .map(|(idx, cnt)| (NT_ORDER[idx], *cnt))
            .expect("SAFETY: Never empty")
    }
}

/// Record elements for `SA` tag for `minimap2` `v2.28`.
/// * Corresponds to `[chrom, start, strand, cigar, mapq, num_mismatches_gaps]`
/// * ex. `chr7,2441699,-,12760S18196M69I,60,120`
type SARecord<'a> = (&'a str, &'a str, &'a str, &'a str, &'a str, &'a str);

/// Get clip direction for CIGAR string.
fn get_clip_direction(cg: impl Iterator<Item = Op>) -> Option<Direction> {
    let mut left_op = None;
    let mut right_op = None;
    for (i, op) in cg
        .enumerate()
        .filter(|(_, op)| matches!(op.kind(), Kind::SoftClip | Kind::HardClip))
    {
        if i == 0 {
            left_op = Some(op)
        } else {
            right_op = Some(op)
        }
    }
    match (left_op, right_op) {
        (None, None) => None,
        (None, Some(_)) => Some(Direction::Right),
        (Some(_), None) => Some(Direction::Left),
        (Some(_), Some(_)) => Some(Direction::Both),
    }
}

/// Check if given position overlaps a given interval and is not an ignored position.
/// 
/// # Arguments
/// * ref_pos: evaluated position in reference
/// * st: start position of checked interval
/// * end: end position of checked interval
/// * itree_ignore: interval tree of ignored intervals
/// 
/// # Returns
/// * If valid.
#[inline]
fn check_valid_pos(
    ref_pos: usize,
    st: usize,
    end: usize,
    itree_ignore: &Lapper<usize, String>,
) -> bool {
    ref_pos >= st && ref_pos <= end && itree_ignore.count(ref_pos, ref_pos) == 0
}

/// Read marker position information used in detection of unbalanced reads and indel outlier reads
pub struct ReadMarkers {
    /// Positions as float for kde calculation
    pub pos: Vec<f64>,
    /// Length of indel (I/D) operations.
    pub n_indels: usize,
    /// Length of mismatch (X) operations.
    pub n_mismatches: usize,
}

/// Collect [`ReadMarkers`] for a given BAM record.
/// * Sets a hardcoded QV filter of 30 (99.9% accurate)
///     * Avoids picking up basecalling artifacts.
/// * Only gets mismatch positions ([`ReadMarkers::pos`]) that are the top count in `secondary_snv_counts`.
///     * Secondary alignments in segdups are often other reads from paralogs.
///     * Therefore, if the read has mismatches matching the top secondary read count, indicates a potential inversion where an alignment was forced.
/// 
/// # Argument
/// * rec: BAM record
/// * st: Start position of BAM query
/// * end: End position of BAM query
/// * itree_ignore: Interval tree of ignored regions
/// * secondary_snv_counts: map of reference position to secondary read nucleotide counts
/// 
/// # Returns
/// * [`ReadMarkers`] for a read.
pub fn collect_read_markers(
    rec: &bam::Record,
    st: usize,
    end: usize,
    itree_ignore: &Lapper<usize, String>,
    secondary_snv_counts: &HashMap<usize, NtCounts>,
) -> eyre::Result<ReadMarkers> {
    let cg = rec.cigar();
    let seq = rec.sequence();
    let qscores = rec.quality_scores().as_bytes();
    let mut marker_qpos = vec![];
    let mut n_mismatches: usize = 0;
    let mut n_indels: usize = 0;

    let mut ref_pos: usize = rec
        .alignment_start()
        .with_context(|| format!("No alignment start for {rec:?}"))??
        .get();
    let mut qpos: usize = 0;
    for (op, l) in cg.iter().flatten().map(|op| (op.kind(), op.len())) {
        match op {
            Kind::Match | Kind::SequenceMatch => {
                for _ in ref_pos..(ref_pos + l) {
                    qpos += 1
                }
                ref_pos += l
            }
            Kind::SequenceMismatch => {
                for _ in ref_pos..(ref_pos + l) {
                    qpos += 1
                }
                // Must not be in ignored region
                // Base QV must be 99.9% accurate
                // Mismatched base must be top secondary nt
                if check_valid_pos(ref_pos, st, end, itree_ignore)
                    && qscores.get(qpos).cloned().unwrap_or_default() > 30
                    && seq
                        .get(qpos)
                        .and_then(|nt| {
                            secondary_snv_counts
                                .get(&ref_pos)
                                .map(|nt_cnts| (nt, nt_cnts))
                        })
                        .is_some_and(|(nt, nt_cnts)| {
                            let (top_nt, _top_nt_cnt) = nt_cnts.top_nt();
                            top_nt == nt
                        })
                {
                    n_mismatches += 1;
                    marker_qpos.push(qpos as f64);
                }
                ref_pos += l
            }
            Kind::Pad | Kind::SoftClip => {
                qpos += l;
            }
            Kind::Insertion => {
                if check_valid_pos(ref_pos, st, end, itree_ignore) {
                    for _ in ref_pos..(ref_pos + l) {
                        n_indels += 1;
                    }
                }
                qpos += l;
            }
            Kind::Deletion => {
                if check_valid_pos(ref_pos, st, end, itree_ignore) {
                    for _ in ref_pos..(ref_pos + l) {
                        n_indels += 1;
                    }
                }
                ref_pos += l
            }
            Kind::HardClip => {
                continue;
            }
            Kind::Skip => ref_pos += l,
        }
    }
    Ok(ReadMarkers {
        pos: marker_qpos,
        n_indels,
        n_mismatches,
    })
}

/// Supplementary interval signals. Represents the paired intervals of supplementary alignment and their directions. 
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SupplIntervalSignals {
    /// Chromosome name
    pub chrom: String,
    /// Primary interval and its direction
    pub itv: Interval<usize, Direction>,
    /// Supplementary interval and its direction
    pub suppl_itv: Interval<usize, Direction>,
    /// Read name
    pub read: String,
}

impl SupplIntervalSignals {
    /// Get breakpoints based on direction.
    /// * If supplementary alignments on both ends, default to the end that minimizes the length of the interval.
    /// 
    /// # Returns
    /// * tuple of start and stop position
    pub fn breakpoints(&self) -> (usize, usize) {
        match (self.itv.val, self.suppl_itv.val) {
            (Direction::Left, Direction::Left) | (Direction::Left, Direction::Both) => (
                std::cmp::min(self.itv.start, self.suppl_itv.start),
                std::cmp::max(self.itv.start, self.suppl_itv.start),
            ),
            (Direction::Left, Direction::Right) => (
                std::cmp::min(self.itv.start, self.suppl_itv.stop),
                std::cmp::max(self.itv.start, self.suppl_itv.stop),
            ),
            (Direction::Right, Direction::Left)
            | (Direction::Right, Direction::Both)
            | (Direction::Both, Direction::Left)
            | (Direction::Both, Direction::Both) => (
                std::cmp::min(self.itv.stop, self.suppl_itv.start),
                std::cmp::max(self.itv.stop, self.suppl_itv.start),
            ),
            (Direction::Right, Direction::Right) | (Direction::Both, Direction::Right) => (
                std::cmp::min(self.itv.stop, self.suppl_itv.stop),
                std::cmp::max(self.itv.stop, self.suppl_itv.stop),
            ),
        }
    }

    /// Get span of supplementary alignment breakpoints.
    /// 
    /// # Returns
    /// * [`Interval`] with [`Interval::val`] representing spanned length.
    pub fn span(&self) -> Interval<usize, usize> {
        let (start, stop) = self.breakpoints();
        let span = start.saturating_sub(stop);
        Interval {
            start,
            stop,
            val: span,
        }
    }
}

/// Assign a [`SVType`] by read clipping direction.
/// 
/// * For inversion, expect both ends to be in same clipping direction
/// ```ignore
/// |x >|              |< *|
///     |< *|      |x >|
/// ```
/// 
/// * For deletion, ends will be in opposing directions.
/// ```ignore
/// |x >|      |< *|
/// |* >|      |< x|
/// ```
/// 
/// * If not one of these, omit.
/// 
/// # Arguments
/// * itv: primary alignment and direction
/// * suppl_itv: supplementary alignment interval and direction
/// 
/// # Returns
/// * Clip direction, if any. Invalid direction pairs are omitted.
fn assign_type_by_clip_direction(
    itv: &Interval<usize, Direction>,
    suppl_itv: &Interval<usize, Direction>,
) -> Option<SVType> {
    // What is the order of the current alignment relative to the suppl alignment?
    let sa_upstream = suppl_itv.start > itv.stop;

    match (sa_upstream, &itv.val, &suppl_itv.val) {
            (_, Direction::Left, Direction::Left) |
            (_, Direction::Right, Direction::Right) => {
                Some(SVType::Inversion)
            },
            // Invalid
            (_, Direction::Both, _) |
            (_, _, Direction::Both) |
            //       *
            // |<  | |  >|
            (true, Direction::Left, Direction::Right) |
            // *     
            // |<  | |  >|
            (false, Direction::Right, Direction::Left) => {
                None
            },
            //       *
            // |  >| |<  |
            // *     
            // |  >| |<  |
            (true, Direction::Right, Direction::Left) |
            (false, Direction::Left, Direction::Right) => {
                Some(SVType::Deletion)
            }
        }
}

/// Get secondary read alignment SNV nucleotide counts.
/// * Enables filtering of primary alignment SNVs that may be due to a inversion
/// * Has hard-coded QV filter of 30 on mismatches
/// 
/// # Arguments
/// * bam: BAM file
/// * region: region to query
/// * itree_ignore: Regions to ignore. If read has majority overlap with a region, is filtered.
/// 
/// # Returns
/// * map of reference position to SNV counts for secondary alignments
fn get_secondary_snv_counts(
    bam: &Path,
    region: &Region,
    itree_ignore: &Lapper<usize, String>,
) -> eyre::Result<HashMap<usize, NtCounts>> {
    let mut fh = bam::io::indexed_reader::Builder::default().build_from_path(bam)?;
    let header = fh.read_header()?;
    // Collect SNVs across secondary alignments
    // Effectively SNVs of other repeats that can help detect if recombination event occurs
    let mut mismatch_pos_cnts = HashMap::new();
    let query = fh.query(&header, region)?;
    let (Bound::Included(st), Bound::Included(end)) = (
        region.start().map(|v| v.get()),
        region.end().map(|v| v.get()),
    ) else {
        unreachable!("Cannot get region bounds: {region:?}")
    };

    for rec in query
        .records()
        .flatten()
        .filter(|rec| rec.flags().contains(Flags::SECONDARY))
    {
        let cg = rec.cigar();
        let seq = rec.sequence();
        let qscores = rec.quality_scores().as_bytes();
        let (rst, rend) = (
            rec.alignment_start().unwrap().map(|p| p.get())?,
            rec.alignment_end().unwrap().map(|p| p.get())?,
        );
        let aln_len = rend - rst;
        let aln_len_f = aln_len as f64;

        // Ignore if read overlaps majority of any ignored region.
        if itree_ignore.find(rst, rend).any(|ovl| {
            let ovl_len = overlap_length(rst, rend, ovl.start, ovl.stop) as f64;
            (ovl_len / aln_len_f) > 0.5
        }) {
            continue;
        }

        let mut ref_pos: usize = rec
            .alignment_start()
            .with_context(|| format!("No alignment start for {rec:?}"))??
            .get();
        let mut qpos: usize = 0;

        for (op, l) in cg.iter().flatten().map(|op| (op.kind(), op.len())) {
            match op {
                Kind::Match | Kind::SequenceMatch => {
                    for _ in ref_pos..(ref_pos + l) {
                        qpos += 1
                    }
                    ref_pos += l
                }
                Kind::SequenceMismatch => {
                    for _ in ref_pos..(ref_pos + l) {
                        qpos += 1
                    }
                    // Must have mismatch
                    if check_valid_pos(ref_pos, st, end, itree_ignore)
                        && qscores.get(qpos).cloned().unwrap_or_default() > 30
                        && let Some(nt) = seq.get(qpos)
                    {
                        mismatch_pos_cnts
                            .entry(ref_pos)
                            .and_modify(|nt_cnts: &mut NtCounts| {
                                // Get LUT value for nt (A=0=00, C=1=01, G=2=10, T=3=11)
                                nt_cnts.0[LUT[nt as usize] as usize] += 1;
                            })
                            .or_insert_with(|| NtCounts([0; 4]));
                    }
                    ref_pos += l
                }
                Kind::Pad | Kind::SoftClip => {
                    qpos += l;
                }
                Kind::Insertion => {
                    qpos += l;
                }
                Kind::Deletion => ref_pos += l,
                Kind::HardClip => {
                    continue;
                }
                Kind::Skip => ref_pos += l,
            }
        }
    }
    Ok(mismatch_pos_cnts)
}

/// Main loop for initial event detection.
/// 1. Get secondary SNV counts across the query window
/// 2. Get mismatch markers in every read
/// 3. Calculate z-score for every read to see if it's indel proportion is an outlier
/// 4. Check if read mismatches are "unbalanced" with SNVs pre-dominantly on one end of the read.
/// 5. Filter mismatch markers that correspond to the largest secondary SNV counts at a given position
/// 6. If the read is part of a supplementary alignment, check to make sure that orientation of clipping make sense for an inversion/deletion.
/// 7. Find overlaps of supplementary alignment breakpoints.
/// 8. Filter supplementary event pairs that are not reciprocal. i.e. both ends need to have a breakpoint spanning each other + 2kbp 
/// 
/// # Arguments
/// * bam: BAM file
/// * itv: interval to query
/// * indel_read_stats: indel read summary stats for primary and secondary alignments
/// * itree_ignore: interval tree of ignored regions.
/// * inv_indel_zscore: z-score to detect outlier indel reads for inversions
/// * inv_min_aln_len: minimum aligned read length to consider read for inversions
/// * inv_thr_unbalanced: threshold in the absolute difference in mismatch distribution between the left/right
///   side of a read to consider a read unbalanced for inversions
/// * inv_min_num_snvs: threshold in minimum required number of SNVs to consider read for inversions
/// * del_min_mapq: threshold MAPQ of read for deletions
/// 
/// # Returns
/// * All putative [`Event`]s in this window.
#[allow(clippy::too_many_arguments)]
pub fn detect_events(
    bam: &Path,
    itv: &Interval<usize, String>,
    indel_read_stats: &ReadIndelSummaryStats,
    itree_ignore: &Lapper<usize, String>,
    inv_indel_zscore: f64,
    inv_min_aln_len: usize,
    inv_thr_unbalanced: f64,
    inv_min_num_snvs: usize,
    del_min_mapq: u8,
) -> eyre::Result<Vec<Event>> {
    let mut fh = bam::io::indexed_reader::Builder::default().build_from_path(bam)?;
    let header = fh.read_header()?;
    let chrom = &itv.val;
    let (st, end) = (itv.start, itv.stop);
    let region = Region::new(
        chrom.to_owned(),
        Position::new(itv.start.clamp(1, usize::MAX)).unwrap()..=Position::new(itv.stop).unwrap(),
    );
    let inv_min_aln_len = inv_min_aln_len as f64;

    let indel_read_stats = [&indel_read_stats.primary, &indel_read_stats.secondary];
    let mut events = vec![];
    let mut sa_itvs = vec![];

    // Get secondary read SNVs in region
    let secondary_snv_counts = get_secondary_snv_counts(bam, &region, itree_ignore)?;

    let query = fh.query(&header, &region)?;
    for rec in query.records().flatten() {
        let read = rec.name().unwrap();
        let cg = rec.cigar();
        let is_suppl = rec.data().get(b"SA").is_some();
        let is_sec = rec.flags().contains(Flags::SECONDARY);
        let typ_read_stats = &indel_read_stats[is_sec as usize];
        let (rst, rend) = (
            rec.alignment_start().unwrap().map(|p| p.get())?,
            rec.alignment_end().unwrap().map(|p| p.get())?,
        );
        let aln_len = rend - rst;
        let aln_len_f = aln_len as f64;

        // Ignore if read overlaps majority of any ignored region.
        if itree_ignore.find(rst, rend).any(|ovl| {
            let ovl_len = overlap_length(rst, rend, ovl.start, ovl.stop) as f64;
            (ovl_len / aln_len_f) > 0.5
        }) {
            continue;
        }

        // Look for:
        // * unbalanced reads bordered by large indels. check secondary alignment
        // * supplementary alignments on same chrom (for now)
        let read_markers =
            collect_read_markers(&rec, st, end, itree_ignore, &secondary_snv_counts)?;
        let indel_rate_zscore = typ_read_stats.zscore(read_markers.n_indels as f64 / aln_len_f);
        let unbalanced_summary = is_unbalanced_alignment(
            &read_markers.pos,
            (rst, rend),
            inv_min_num_snvs,
            inv_thr_unbalanced,
        )?;

        // Inversion signal based on mismatches
        if indel_rate_zscore > inv_indel_zscore
            && aln_len_f > inv_min_aln_len
            && read_markers.n_mismatches >= inv_min_num_snvs
            && let Some(unbalanced_dir) =
                unbalanced_summary.as_ref().and_then(|summ| summ.direction)
            && let Some(unbalanced_span) = unbalanced_summary.and_then(|summ| summ.span)
        {
            let event = MismatchSignal {
                itv: Interval {
                    start: rst,
                    stop: rend,
                    val: chrom.to_owned(),
                },
                read: String::from_utf8(read.to_vec())?,
                n_indels: read_markers.n_indels,
                aln_len,
                is_secondary: is_sec,
                unbalanced_dir,
                unbalanced_span,
            };
            events.push(Event::Signal(Signal::MismatchSignal(event)));
        }
        // Deletion signal based on directional supplementary aln signals
        if is_suppl {
            let Value::String(sa_tag) = rec
                .data()
                .get(b"SA")
                .with_context(|| format!("Must have SA tag for {read}."))??
            else {
                bail!("Invalid type for SA tag for {read}.")
            };

            // Store suppl inversion signal intervals
            let clip_direction = get_clip_direction(cg.iter().flatten()).with_context(|| {
                format!("Read {read} must have soft/hardclipped operation in cigar: {cg:?}")
            })?;
            for (_sa_chrom, sa_start, _sa_strand, sa_cigar, sa_mapq, _sa_nm) in
                str::from_utf8(sa_tag)?
                    .split(';')
                    .flat_map(|rec| rec.splitn(6, ',').collect_tuple::<SARecord>())
                    // Must be same chromosome
                    .filter(|sa_rec| sa_rec.0 == chrom)
            {
                let sa_cigar = noodles::sam::record::Cigar::new(sa_cigar.as_bytes());
                let sa_start: usize = sa_start.parse()?;
                let sa_aln_len = sa_cigar.alignment_span()?;
                let sa_end = sa_start + sa_aln_len;

                let sa_mapq: u8 = sa_mapq.parse()?;
                if sa_mapq < del_min_mapq {
                    continue;
                }

                let sa_clip_direction = get_clip_direction(sa_cigar.iter().flatten())
                    .with_context(|| format!("Read {read} must have soft/hardclipped operation in SA cigar: {sa_cigar:?}"))?;

                // Add a small overlap distance of 2 kbp in case of alignment forcing suppl align signal slightly upstream/downstream
                match clip_direction {
                    // |<  |
                    Direction::Left => sa_itvs.push(Interval {
                        start: rst.saturating_sub(CLIP_OVL_DST),
                        stop: rst + CLIP_OVL_DST,
                        val: SupplIntervalSignals {
                            chrom: chrom.to_owned(),
                            itv: Interval {
                                start: rst,
                                stop: rend,
                                val: clip_direction,
                            },
                            suppl_itv: Interval {
                                start: sa_start,
                                stop: sa_end,
                                val: sa_clip_direction,
                            },
                            read: String::from_utf8(read.to_vec())?,
                        },
                    }),
                    // |  >|
                    Direction::Right => sa_itvs.push(Interval {
                        start: rend.saturating_sub(CLIP_OVL_DST),
                        stop: rend + CLIP_OVL_DST,
                        val: SupplIntervalSignals {
                            chrom: chrom.to_owned(),
                            itv: Interval {
                                start: rst,
                                stop: rend,
                                val: clip_direction,
                            },
                            suppl_itv: Interval {
                                start: sa_start,
                                stop: sa_end,
                                val: sa_clip_direction,
                            },
                            read: String::from_utf8(read.to_vec())?,
                        },
                    }),
                    // |< >|
                    Direction::Both => {
                        sa_itvs.push(Interval {
                            start: rst.saturating_sub(CLIP_OVL_DST),
                            stop: rst + CLIP_OVL_DST,
                            val: SupplIntervalSignals {
                                chrom: chrom.to_owned(),
                                itv: Interval {
                                    start: rst,
                                    stop: rend,
                                    val: clip_direction,
                                },
                                suppl_itv: Interval {
                                    start: sa_start,
                                    stop: sa_end,
                                    val: sa_clip_direction,
                                },
                                read: String::from_utf8(read.to_vec())?,
                            },
                        });
                        sa_itvs.push(Interval {
                            start: rend.saturating_sub(CLIP_OVL_DST),
                            stop: rend + CLIP_OVL_DST,
                            val: SupplIntervalSignals {
                                chrom: chrom.to_owned(),
                                itv: Interval {
                                    start: rst,
                                    stop: rend,
                                    val: clip_direction,
                                },
                                suppl_itv: Interval {
                                    start: sa_start,
                                    stop: sa_end,
                                    val: sa_clip_direction,
                                },
                                read: String::from_utf8(read.to_vec())?,
                            },
                        })
                    }
                };
            }
        }
    }

    // Inversion supplementary read intersections
    let itree_suppl_inv = Lapper::new(sa_itvs);
    for sitv in itree_suppl_inv.iter() {
        let ovl_itvs = itree_suppl_inv
            .find(sitv.start, sitv.stop)
            // Omit self
            .filter(|oitv| oitv.val != sitv.val)
            // Each overlapping interval should also overlap with its supplementary itv
            // And should be different reads
            .filter(|oitv| {
                oitv.val
                    .suppl_itv
                    .overlap(sitv.val.suppl_itv.start, sitv.val.suppl_itv.stop)
                    && sitv.val.read != oitv.val.read
            })
            .collect_vec();

        if ovl_itvs.is_empty() {
            continue;
        }
        if ovl_itvs.len() > 2 {
            log::debug!("More overlaps than expected for itv: {itv:?}");
            log::debug!("{ovl_itvs:#?}");
            continue;
        }

        let all_contained = ovl_itvs.iter().all(|oitv| {
            itv.overlap(oitv.val.itv.start, oitv.val.itv.stop)
                && itv.overlap(oitv.val.suppl_itv.start, oitv.val.suppl_itv.stop)
        });
        // If all contained in this window, means we're done
        if all_contained {
            if let (Some(itv_1), Some(itv_2)) = (ovl_itvs.first(), ovl_itvs.get(1)) {
                let putative_type_1 =
                    assign_type_by_clip_direction(&itv_1.val.itv, &itv_1.val.suppl_itv);
                let putative_type_2 =
                    assign_type_by_clip_direction(&itv_2.val.itv, &itv_2.val.suppl_itv);
                log::debug!(
                    "{:?},{itv:?}\n\t{putative_type_1:?},{putative_type_2:?},{},{}",
                    (sitv.start, sitv.stop),
                    itv_1.val.read,
                    itv_2.val.read
                )
            }
            if let Some(itv_1) = ovl_itvs.first() {
                let Some(putative_type_1) =
                    assign_type_by_clip_direction(&itv_1.val.itv, &itv_1.val.suppl_itv)
                else {
                    continue;
                };

                let itv_span = itv_1.val.span();
                events.push(Event::SV(SV {
                    itv: Interval {
                        start: itv_span.start,
                        stop: itv_span.stop,
                        val: chrom.to_owned(),
                    },
                    reads: [Some(itv_1.val.read.clone()), None],
                    typ: putative_type_1,
                }));
            }
        } else if let Some(putative_type) =
            assign_type_by_clip_direction(&sitv.val.itv, &sitv.val.suppl_itv)
        {
            // Otherwise pass
            // SupplSignal
            match putative_type {
                SVType::Deletion => {
                    let span = sitv.val.span();
                    events.push(Event::SV(SV {
                        itv: Interval {
                            start: span.start,
                            stop: span.stop,
                            val: chrom.to_owned(),
                        },
                        reads: [Some(sitv.val.read.clone()), None],
                        typ: putative_type,
                    }));
                }
                SVType::Inversion => {
                    events.push(Event::Signal(Signal::SupplSignal(sitv.val.clone())))
                }
            }
        }
    }

    Ok(events)
}
