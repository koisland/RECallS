use std::path::Path;

use clap::Parser;
use eyre::ContextCompat;
use noodles::{
    bam,
    core::Region,
    sam::alignment::record::{Flags, cigar::op::Kind},
};

use crate::{
    cli::Args,
    io::{aligned_intervals_windows, read_bed},
    utils::{get_aligned_pairs, get_coords_from_region},
};

mod cli;
mod dotplot;
mod io;
mod unbalanced_aln;
mod utils;

fn detect_events(bam: &Path, _fa: &Path, region: Region) -> eyre::Result<()> {
    let mut fh = bam::io::indexed_reader::Builder::default().build_from_path(bam)?;
    let header = fh.read_header()?;
    let (st, end) = get_coords_from_region(&region)?;
    let length = end - st;
    let query = fh.query(&header, &region)?;

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
        let seq = rec.sequence();

        // Look for:
        // * unbalanced reads bordered by large indels. check secondary alignment
        // * supplementary alignments on same chrom (for now)
        for (qpos, refpos, kind) in aln_pairs
            .into_iter()
            .filter(|(_, refpos, _)| *refpos >= st && *refpos <= end)
        {
            match kind {
                Kind::Insertion => {},
                Kind::Deletion => {},
                Kind::SequenceMismatch => {
                    // 0-93 ASCII+33 for pacbio
                    let qscore = qscores[qpos];
                    if qscore > 30 {
                        let nt = seq.get(qpos).unwrap();
                        eprintln!("{nt}")
                    }
                }
                _ => {},
            };
        }
    }

    Ok(())
}

fn main() -> eyre::Result<()> {
    let args = Args::parse();

    let bam = &args.bam;
    let fa = &args.fa;

    // ENA_CBCUDK010000011_CBCUDK010000011.1:6334613-6342169
    let regions = if let Some(bed) = &args.bed {
        read_bed(bed.as_ref()).with_context(|| format!("No valid intervals in {bed:?}"))
    } else {
        let mut fh = bam::io::indexed_reader::Builder::default().build_from_path(bam)?;
        aligned_intervals_windows(&mut fh, args.wg_window)
    }?;

    for region in regions {
        detect_events(bam, fa, region)?;
    }

    Ok(())
}
