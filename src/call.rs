use std::path::Path;

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
    self_align::Paf,
    unbalanced_aln::is_unbalanced_alignment,
    utils::overlap_length,
};

const CLIP_OVL_DST: usize = 2000;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum ClipDirection {
    Left,
    Right,
    Both,
}

type SARecord<'a> = (&'a str, &'a str, &'a str, &'a str, &'a str, &'a str);

fn get_clip_direction(cg: impl Iterator<Item = Op>) -> Option<ClipDirection> {
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
        (None, Some(_)) => Some(ClipDirection::Right),
        (Some(_), None) => Some(ClipDirection::Left),
        (Some(_), Some(_)) => Some(ClipDirection::Both),
    }
}

#[inline]
fn check_valid_itv(
    ref_pos: usize,
    st: usize,
    end: usize,
    itree_ignore: &Lapper<usize, String>,
) -> bool {
    ref_pos >= st && ref_pos <= end && itree_ignore.count(ref_pos, ref_pos) == 0
}

#[allow(unused)]
pub struct ReadMarkers {
    pub pos: Vec<f64>,
    pub n_indels: usize,
    pub n_mismatches: usize,
}

pub fn collect_read_markers(
    rec: &bam::Record,
    st: usize,
    end: usize,
    itree_ignore: &Lapper<usize, String>,
) -> eyre::Result<ReadMarkers> {
    let cg = rec.cigar();
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
                // Must have mismatch
                if check_valid_itv(ref_pos, st, end, itree_ignore)
                    && qscores.get(qpos).cloned().unwrap_or_default() > 30
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
                if check_valid_itv(ref_pos, st, end, itree_ignore) {
                    for _ in ref_pos..(ref_pos + l) {
                        n_indels += 1;
                        // marker_qpos.push(qpos as f64);
                    }
                }
                qpos += l;
            }
            Kind::Deletion => {
                if check_valid_itv(ref_pos, st, end, itree_ignore) {
                    for _ in ref_pos..(ref_pos + l) {
                        n_indels += 1;
                        // marker_qpos.push(qpos as f64);
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SupplIntervalSignals {
    pub chrom: String,
    pub itv: Interval<usize, ClipDirection>,
    pub suppl_itv: Interval<usize, ClipDirection>,
    pub read: String,
}

impl SupplIntervalSignals {
    pub fn new(
        chrom: String,
        itv: Interval<usize, ClipDirection>,
        suppl_itv: Interval<usize, ClipDirection>,
        read: String,
    ) -> Self {
        SupplIntervalSignals {
            chrom,
            itv,
            suppl_itv,
            read,
        }
    }

    pub fn breakpoints(&self) -> (usize, usize) {
        match (self.itv.val, self.suppl_itv.val) {
            (ClipDirection::Left, ClipDirection::Left)
            | (ClipDirection::Left, ClipDirection::Both) => (
                std::cmp::min(self.itv.start, self.suppl_itv.start),
                std::cmp::max(self.itv.start, self.suppl_itv.start),
            ),
            (ClipDirection::Left, ClipDirection::Right) => (
                std::cmp::min(self.itv.start, self.suppl_itv.stop),
                std::cmp::max(self.itv.start, self.suppl_itv.stop),
            ),
            (ClipDirection::Right, ClipDirection::Left)
            | (ClipDirection::Right, ClipDirection::Both)
            | (ClipDirection::Both, ClipDirection::Left)
            | (ClipDirection::Both, ClipDirection::Both) => (
                std::cmp::min(self.itv.stop, self.suppl_itv.start),
                std::cmp::max(self.itv.stop, self.suppl_itv.start),
            ),
            (ClipDirection::Right, ClipDirection::Right)
            | (ClipDirection::Both, ClipDirection::Right) => (
                std::cmp::min(self.itv.stop, self.suppl_itv.stop),
                std::cmp::max(self.itv.stop, self.suppl_itv.stop),
            ),
        }
    }
    // if both, take the minimal distance
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

// For inversion, expect both ends to be in same clipping direction
// |x >|              |< *|
//     |< *|      |x >|
// For deletion, ends will be in opposing directions.
// |x >|      |< *|
// |* >|      |< x|
// If not one of these, omit.
fn assign_type_by_clip_direction(
    itv: &Interval<usize, ClipDirection>,
    suppl_itv: &Interval<usize, ClipDirection>,
) -> Option<SVType> {
    // What is the order of the current alignment relative to the suppl alignment?
    let sa_upstream = suppl_itv.start > itv.stop;

    match (sa_upstream, &itv.val, &suppl_itv.val) {
            (_, ClipDirection::Left, ClipDirection::Left) |
            (_, ClipDirection::Right, ClipDirection::Right) => {
                Some(SVType::Inversion)
            },
            // Invalid
            (_, ClipDirection::Both, _) |
            (_, _, ClipDirection::Both) |
            //       *
            // |<  | |  >|
            (true, ClipDirection::Left, ClipDirection::Right) |
            // *     
            // |<  | |  >|
            (false, ClipDirection::Right, ClipDirection::Left) => {
                None
            },
            //       *
            // |  >| |<  |
            // *     
            // |  >| |<  |
            (true, ClipDirection::Right, ClipDirection::Left) |
            (false, ClipDirection::Left, ClipDirection::Right) => {
                // let n_similar = itree_self_similar.count(sa_start, sa_start+sa_aln_len);
                // if n_similar != 0 {
                //     let event = SupplSignal {
                //         chrom: chrom.to_owned(),
                //         start: rst,
                //         stop: rend,
                //         suppl_start: sa_start,
                //         suppl_stop: sa_end,
                //         rname: String::from_utf8(rname.to_vec())?,
                //     };
                //     // TODO: Dedup due to two suppl alns
                //     events.push(Event::Deletion(event));
                // }
                Some(SVType::Deletion)
            }
        }
}

#[allow(clippy::too_many_arguments)]
pub fn detect_events(
    bam: &Path,
    itv: &Interval<usize, String>,
    read_stats: &ReadIndelSummaryStats,
    itree_ignore: &Lapper<usize, String>,
    itree_self_similar: &Lapper<usize, Paf>,
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
    // Get intervaltree of ignored regions
    let query = fh.query(&header, &region)?;

    let indel_read_stats = [&read_stats.primary, &read_stats.secondary];
    let mut events = vec![];
    let mut sa_itvs = vec![];

    for rec in query.records().flatten() {
        let read = rec.name().unwrap();
        let cg = rec.cigar();
        let is_suppl = rec.data().get(b"SA").is_some();
        let is_sec = rec.flags().contains(Flags::SECONDARY);
        let typ_read_stats = &indel_read_stats[is_sec as usize];
        let aln_len = noodles::sam::alignment::Record::alignment_span(&rec).unwrap()?;
        let aln_len_f = aln_len as f64;
        let (rst, rend) = (
            rec.alignment_start().unwrap().map(|p| p.get())?,
            rec.alignment_end().unwrap().map(|p| p.get())?,
        );

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
        let read_markers = collect_read_markers(&rec, st, end, itree_ignore)?;
        let indel_rate_zscore = typ_read_stats.zscore(read_markers.n_indels as f64 / aln_len_f);
        let unbalanced_summary = is_unbalanced_alignment(
            &read_markers.pos,
            aln_len as f64,
            inv_min_num_snvs,
            inv_thr_unbalanced,
        )?;

        // Inversion signal based on mismatches
        let is_unbalanced = unbalanced_summary
            .as_ref()
            .map(|s| s.is_unbalanced)
            .unwrap_or_default();
        if indel_rate_zscore > inv_indel_zscore
            && aln_len_f > inv_min_aln_len
            && read_markers.n_mismatches >= inv_min_num_snvs
            && is_unbalanced
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
                is_unbalanced,
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
            // minimap2 v2.28 SA tag format
            //
            // chrom,start,strand,cigar,mapq,num_mismatches_gaps
            // chr7,2441699,-,12760S18196M69I,60,120
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

                match clip_direction {
                    // |<  |
                    ClipDirection::Left => sa_itvs.push(Interval {
                        start: rst.saturating_sub(CLIP_OVL_DST),
                        stop: rst + CLIP_OVL_DST,
                        val: SupplIntervalSignals::new(
                            chrom.to_owned(),
                            Interval {
                                start: rst,
                                stop: rend,
                                val: clip_direction,
                            },
                            Interval {
                                start: sa_start,
                                stop: sa_end,
                                val: sa_clip_direction,
                            },
                            String::from_utf8(read.to_vec())?,
                        ),
                    }),
                    // |  >|
                    ClipDirection::Right => sa_itvs.push(Interval {
                        start: rend.saturating_sub(CLIP_OVL_DST),
                        stop: rend + CLIP_OVL_DST,
                        val: SupplIntervalSignals::new(
                            chrom.to_owned(),
                            Interval {
                                start: rst,
                                stop: rend,
                                val: clip_direction,
                            },
                            Interval {
                                start: sa_start,
                                stop: sa_end,
                                val: sa_clip_direction,
                            },
                            String::from_utf8(read.to_vec())?,
                        ),
                    }),
                    // |< >|
                    ClipDirection::Both => {
                        sa_itvs.push(Interval {
                            start: rst.saturating_sub(CLIP_OVL_DST),
                            stop: rst + CLIP_OVL_DST,
                            val: SupplIntervalSignals::new(
                                chrom.to_owned(),
                                Interval {
                                    start: rst,
                                    stop: rend,
                                    val: clip_direction,
                                },
                                Interval {
                                    start: sa_start,
                                    stop: sa_end,
                                    val: sa_clip_direction,
                                },
                                String::from_utf8(read.to_vec())?,
                            ),
                        });
                        sa_itvs.push(Interval {
                            start: rend.saturating_sub(CLIP_OVL_DST),
                            stop: rend + CLIP_OVL_DST,
                            val: SupplIntervalSignals::new(
                                chrom.to_owned(),
                                Interval {
                                    start: rst,
                                    stop: rend,
                                    val: clip_direction,
                                },
                                Interval {
                                    start: sa_start,
                                    stop: sa_end,
                                    val: sa_clip_direction,
                                },
                                String::from_utf8(read.to_vec())?,
                            ),
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
            eprintln!("More overlaps than expected for itv: {itv:?}");
            eprintln!("{ovl_itvs:#?}");
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
                eprintln!(
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
                    read: itv_1.val.read.clone(),
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
                        read: sitv.val.read.clone(),
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
