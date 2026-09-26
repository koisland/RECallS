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

/// Calculate mean indel rate for primary/supplementary and secondary alignments as prior
pub fn calculate_mean_indel_rate(
    bam: &Path,
    itv: &Interval<usize, String>,
    ignore_bed: &HashMap<String, Lapper<usize, String>>,
) -> eyre::Result<(f64, f64)> {
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
    let mut n_reads: [usize; 2] = [0, 0];
    let mut sum_perc_indel: [f64; 2] = [0f64, 0f64];

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
            sum_perc_indel[idx] = sum_perc_indel[idx].algebraic_add(perc_indel);
            n_reads[idx] += 1;
        }
    }
    let prim_n_reads_indel = n_reads[0] as f64;
    let prim_sum_perc_indel = sum_perc_indel[0];
    let sec_n_reads_indel = n_reads[1] as f64;
    let sec_sum_perc_indel = sum_perc_indel[1];
    let prim_avg_indel = if prim_n_reads_indel == 0.0 {
        0.0
    } else {
        prim_sum_perc_indel / prim_n_reads_indel
    };
    let sec_avg_indel = if sec_n_reads_indel == 0.0 {
        0.0
    } else {
        sec_sum_perc_indel / sec_n_reads_indel
    };
    Ok((prim_avg_indel, sec_avg_indel))
}
