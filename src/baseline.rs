use std::{collections::HashMap, path::Path};

use noodles::{
    bam,
    core::{Position, Region},
    sam::alignment::{
        Record,
        record::{Flags, cigar::op::Kind},
    },
};
use rust_lapper::{Interval, Lapper};

#[derive(Debug)]
pub struct ReadSummaryStats {
    pub primary: SummaryStats,
    pub secondary: SummaryStats,
}

#[derive(Debug, Default)]
pub struct SummaryStats {
    pub mean: f64,
    pub var: f64,
    pub stdev: f64,
    pub n: usize,
}

// https://rust-lang-nursery.github.io/rust-cookbook/science/mathematics/statistics.html
impl SummaryStats {
    pub fn new(data: &[f64]) -> Self {
        let n = data.len() as f64;
        if n == 0.0 {
            return SummaryStats::default();
        }
        let mean = data.iter().sum::<f64>() / n;
        let var = data
            .iter()
            .map(|value| {
                let diff = mean - (*value as f64);
                diff * diff
            })
            .sum::<f64>()
            / n;
        return Self {
            mean,
            var,
            stdev: var.sqrt(),
            n: data.len(),
        };
    }

    pub fn zscore(&self, x: f64) -> f64 {
        let diff = x - self.mean;
        diff / self.stdev
    }
}

/// Calculate mean indel rate for primary/supplementary and secondary alignments as prior
pub fn calculate_stats_indel_rate(
    bam: &Path,
    itv: &Interval<usize, String>,
    ignore_bed: &HashMap<String, Lapper<usize, String>>,
) -> eyre::Result<(SummaryStats, SummaryStats)> {
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
                    indel_len += length;
                }
                Kind::Pad | Kind::SoftClip | Kind::HardClip => {}
                Kind::Deletion => {
                    if itree_ignore.count(pos, pos + length) != 0 {
                        omit_len += length;
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

    let stats_prim_perc_indel = SummaryStats::new(&both_perc_indel[0]);
    let stats_sec_perc_indel = SummaryStats::new(&both_perc_indel[1]);
    Ok((stats_prim_perc_indel, stats_sec_perc_indel))
}
