use std::{collections::HashMap, ops::Add, path::Path};

use itertools::Itertools;
use noodles::{
    bam,
    core::{Position, Region},
    sam::alignment::{
        Record,
        record::{Flags, cigar::op::Kind},
    },
};
use rayon::prelude::*;
use rust_lapper::{Interval, Lapper};

#[derive(Debug, Clone)]
pub struct ReadIndelSummaryStats {
    pub primary: IndelSummaryStats,
    pub secondary: IndelSummaryStats,
}

impl Add for ReadIndelSummaryStats {
    type Output = ReadIndelSummaryStats;

    fn add(mut self, rhs: Self) -> Self::Output {
        self.primary.mean = self.primary.mean.algebraic_add(rhs.primary.mean);
        self.primary.var = self.primary.var.algebraic_add(rhs.primary.var);
        self.primary.n += rhs.primary.n;
        // secondary
        self.secondary.mean = self.secondary.mean.algebraic_add(rhs.secondary.mean);
        self.secondary.var = self.secondary.var.algebraic_add(rhs.secondary.var);
        self.secondary.n += rhs.secondary.n;
        self
    }
}

#[derive(Debug, Default, Clone)]
pub struct IndelSummaryStats {
    pub mean: f64,
    pub var: f64,
    pub stdev: f64,
    pub n: usize,
}

// https://rust-lang-nursery.github.io/rust-cookbook/science/mathematics/statistics.html
impl IndelSummaryStats {
    pub fn new(data: &[f64]) -> Self {
        let n = data.len() as f64;
        if n == 0.0 {
            return IndelSummaryStats::default();
        }
        let mean = data.iter().sum::<f64>() / n;
        let var = data
            .iter()
            .map(|value| {
                let diff = mean - *value;
                diff * diff
            })
            .sum::<f64>()
            / n;
        Self {
            mean,
            var,
            stdev: var.sqrt(),
            n: data.len(),
        }
    }

    pub fn zscore(&self, x: f64) -> f64 {
        let diff = x - self.mean;
        diff / self.stdev
    }
}

/// Calculate mean indel rate for primary/supplementary and secondary alignments as prior
///
/// # Arguments
/// * bam: path to bam
/// * itv: whole-contig interval
/// * ignore_bed: regions to ignore
///
/// # Returns
/// * tuple of `IndelSummaryStats` where first is of primary alignments and second is of secondary alignments.
pub fn calculate_stats_indel_rate(
    bam: &Path,
    itv: &Interval<usize, String>,
    ignore_bed: &HashMap<String, Lapper<usize, String>>,
) -> eyre::Result<(IndelSummaryStats, IndelSummaryStats)> {
    let mut fh = bam::io::indexed_reader::Builder::default().build_from_path(bam)?;
    let header = fh.read_header()?;
    let chrom = &itv.val;
    // Cannot be 0
    let (st, stop) = (
        itv.start.clamp(1, usize::MAX),
        itv.stop.clamp(1, usize::MAX),
    );
    let region = Region::new(
        chrom.to_owned(),
        Position::new(st).unwrap()..=Position::new(stop).unwrap(),
    );
    // Get intervaltree of ignored regions
    let null_itree_ignore = Lapper::new(vec![]);
    let itree_ignore = ignore_bed.get(chrom).unwrap_or(&null_itree_ignore);
    let query = fh.query(&header, &region)?;

    // Primary, secondary
    let mut both_perc_indel: [Vec<f64>; 2] = [vec![], vec![]];

    // Iter thru all records
    for rec in query.records().flatten() {
        let cg = rec.cigar();
        // Index into arr
        let idx = rec.flags().contains(Flags::SECONDARY) as usize;

        // Can be 0?
        let aln_len = rec.alignment_span().transpose()?.unwrap_or_default() as f64;
        let mut indel_len = 0;
        let mut omit_len = 0;

        let mut pos: usize = rec.alignment_start().unwrap()?.get();
        for op in cg.iter() {
            let op = op?;
            let kind = op.kind();
            let length = op.len();
            match kind {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    if itree_ignore.count(pos, pos + length) != 0 {
                        omit_len += length;
                    }
                    pos += length;
                }
                Kind::Insertion => {
                    // Even though doesn't consume reference position, we want to track length.
                    // May result in indel rate over 100%?
                    if itree_ignore.count(pos, pos + length) == 0 {
                        indel_len += length;
                    }
                }
                Kind::Pad | Kind::SoftClip | Kind::HardClip => {}
                Kind::Deletion => {
                    if itree_ignore.count(pos, pos + length) != 0 {
                        omit_len += length;
                    } else {
                        indel_len += length;
                    }
                    pos += length;
                }
                Kind::Skip => {
                    if itree_ignore.count(pos, pos + length) != 0 {
                        omit_len += length;
                    }
                    pos += length
                }
            }
        }
        // divide by zero
        let adj_aln_len = (aln_len - omit_len as f64).clamp(0.0, f64::MAX);
        if adj_aln_len != 0.0 {
            let perc_indel = (indel_len as f64) / adj_aln_len;
            // Add to 0 (primary/suppl) or 1 (secondary)
            both_perc_indel[idx].push(perc_indel);
        }
    }

    let stats_prim_perc_indel = IndelSummaryStats::new(&both_perc_indel[0]);
    let stats_sec_perc_indel = IndelSummaryStats::new(&both_perc_indel[1]);
    Ok((stats_prim_perc_indel, stats_sec_perc_indel))
}

pub fn aggregate_stats_indel_rate(
    regions: &HashMap<String, Vec<Interval<usize, String>>>,
    bam: &Path,
    ignore_bed: &HashMap<String, Lapper<usize, String>>,
) -> HashMap<String, ReadIndelSummaryStats> {
    // https://stats.stackexchange.com/a/26647
    let regions: Vec<&Interval<usize, String>> = regions.values().flatten().collect();
    let mut chrom_read_stats: HashMap<String, ReadIndelSummaryStats> = regions
        .into_par_iter()
        .map(|region| {
            (
                region.val.clone(),
                calculate_stats_indel_rate(bam, region, ignore_bed).unwrap(),
            )
        })
        .fold(
            HashMap::new,
            |mut acc: HashMap<String, ReadIndelSummaryStats>, (chrom, (prim_stats, sec_stats))| {
                if let Some(read_stats) = acc.get_mut(&chrom) {
                    // Sum up stats (mean, stdev, and n) across windows
                    let new_read_stats = read_stats.clone()
                        + ReadIndelSummaryStats {
                            primary: prim_stats,
                            secondary: sec_stats,
                        };
                    *read_stats = new_read_stats
                } else {
                    acc.insert(
                        chrom.to_owned(),
                        ReadIndelSummaryStats {
                            primary: prim_stats,
                            secondary: sec_stats,
                        },
                    );
                }
                acc
            },
        )
        .reduce(HashMap::new, |mut a_stats, mut b_stats| {
            let mut all_stats = HashMap::new();
            let chroms = a_stats
                .keys()
                .chain(b_stats.keys())
                .unique()
                .cloned()
                .collect_vec();
            for chrom in chroms {
                let a_chrom_stats = a_stats.remove(&chrom);
                let b_chrom_stats = b_stats.remove(&chrom);
                // Merge hashmaps
                match (a_chrom_stats, b_chrom_stats) {
                    (None, Some(b_chrom_stats)) => {
                        all_stats.insert(chrom, b_chrom_stats);
                    }
                    (Some(a_chrom_stats), None) => {
                        all_stats.insert(chrom, a_chrom_stats);
                    }
                    (Some(a_chrom_stats), Some(b_chrom_stats)) => {
                        let new_read_stats = a_chrom_stats + b_chrom_stats;
                        all_stats.insert(chrom, new_read_stats);
                    }
                    _ => unreachable!("Keys derived from a_stats so not possible"),
                }
            }
            all_stats
        });

    // Then average stats
    for val in chrom_read_stats.values_mut() {
        let prim = &mut val.primary;
        let sec = &mut val.secondary;
        // Update mean
        prim.mean /= prim.n as f64;
        sec.mean /= sec.n as f64;
        // Update variance
        prim.var /= prim.n as f64;
        sec.var /= sec.n as f64;
        // Update stdev
        prim.stdev = prim.var.sqrt();
        sec.stdev = sec.var.sqrt();
    }
    chrom_read_stats
}
