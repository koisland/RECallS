use std::{collections::HashMap, ops::Bound};

use eyre::bail;
use noodles::{
    bam::{self, Record},
    core::{Position, Region},
    sam::alignment::record::{Flags, cigar::op::Kind},
};

mod cli;
mod dotplot;

/// Convert cigar string to operations.
/// * Adapted from <https://github.com/pysam-developers/pysam/blob/3e3c8b0b5ac066d692e5c720a85d293efc825200/pysam/libcalignedsegment.pyx#L2009>
pub(crate) fn get_aligned_pairs(
    cg: impl Iterator<Item = (Kind, usize)>,
    pos: usize,
) -> eyre::Result<Vec<(usize, usize, Kind)>> {
    let mut pos: usize = pos;
    let mut qpos: usize = 0;
    let mut pairs = vec![];
    // Matches only
    for (op, l) in cg {
        match op {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in pos..(pos + l) {
                    pairs.push((qpos, i, op));
                    qpos += 1
                }
                pos += l
            }
            // Track indels and softclips.
            Kind::Pad | Kind::Insertion | Kind::SoftClip => {
                qpos += l;
                continue;
            }
            Kind::Deletion => {
                for i in pos..(pos + l) {
                    pairs.push((qpos, i, op));
                }
                pos += l
            }
            Kind::HardClip => {
                continue;
            }
            Kind::Skip => pos += l,
        }
    }
    Ok(pairs)
}

fn generate_mismatch_dag(aln: &str, region: Region) -> eyre::Result<()> {
    let mut indexed_reader = bam::io::indexed_reader::Builder::default().build_from_path(&aln)?;
    let header = indexed_reader.read_header()?;
    let (Bound::Included(st), Bound::Included(end)) = (
        region.start().map(|b| b.get()),
        region.end().map(|b| b.get()),
    ) else {
        bail!("Invalid st or end")
    };
    let length = end - st;
    let query = indexed_reader.query(&header, &region)?;

    // Global metrics
    let mut mapq: Vec<usize> = vec![0; length + 1];
    let mut cov: Vec<usize> = vec![0; length + 1];

    for rec in query
        .records()
        .flatten()
        .filter(|aln| !aln.flags().contains(Flags::SECONDARY))
    {
        let cg: bam::record::Cigar<'_> = rec.cigar();
        let aln_pairs = get_aligned_pairs(
            cg.iter().flatten().map(|op| (op.kind(), op.len())),
            rec.alignment_start().unwrap()?.get(),
        )?;
        let qscores = rec.quality_scores().as_bytes();
        let rmapq = rec.mapping_quality().unwrap_or_default().get() as usize;
        let seq = rec.sequence();

        for (qpos, refpos, kind) in aln_pairs
            .into_iter()
            .filter(|(_, refpos, _)| *refpos >= st && *refpos <= end)
        {
            let ipos = refpos - st;
            let cnt = match kind {
                Kind::SoftClip => 0,
                Kind::Insertion => 0,
                Kind::SequenceMatch | Kind::Deletion => 1,
                Kind::SequenceMismatch => {
                    // 0-93 ASCII+33 for pacbio
                    let qscore = qscores[qpos];
                    if qscore > 30 {
                        let nt = seq.get(qpos).unwrap();
                    }
                    1
                }
                _ => 1,
            };
            // Store coverage and MAPQ (total)
            cov[ipos] += cnt;
            mapq[ipos] += rmapq;
        }
    }

    Ok(())
}

fn main() -> eyre::Result<()> {
    let bam = "test/CT22_ENA_CBCUDK010000011_CBCUDK010000011.1_6335921-6341074.bam";
    // let fa = "test/CT22_ENA_CBCUDK010000011_CBCUDK010000011.1.fa.gz";

    // ENA_CBCUDK010000011_CBCUDK010000011.1:6334613-6342169
    let region = Region::new(
        "ENA_CBCUDK010000011_CBCUDK010000011.1",
        Position::new(6317808).unwrap()..=Position::new(6359969).unwrap(),
    );
    generate_mismatch_dag(bam, region)?;
    Ok(())
}
