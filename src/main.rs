use std::{collections::HashMap, fs::File, io::BufWriter, path::Path};

use clap::Parser;
use eyre::ContextCompat;
use noodles::{
    bam::{self},
    core::{Position, Region},
    sam::alignment::{
        RecordBuf,
        io::Write,
        record::{Flags, cigar::op::Kind, data::field::Tag},
        record_buf::data::field::Value,
    },
};
use rust_lapper::{Interval, Lapper};

use crate::{
    baseline::calculate_mean_indel_rate,
    cli::Args,
    io::{aligned_intervals_windows, read_bed},
    unbalanced_aln::is_unbalanced_alignment,
    utils::get_aligned_pairs,
};

mod baseline;
mod cli;
mod dotplot;
mod io;
mod unbalanced_aln;
mod utils;

fn detect_events(
    bam: &Path,
    _fa: &Path,
    itv: &Interval<usize, String>,
    indel_rate: (f64, f64),
    ignore_bed: &HashMap<String, Lapper<usize, String>>,
) -> eyre::Result<()> {
    let mut fh = bam::io::indexed_reader::Builder::default().build_from_path(bam)?;
    let header = fh.read_header()?;
    let chrom = &itv.val;
    let (st, end) = (itv.start, itv.stop);
    let region = Region::new(
        chrom.to_owned(),
        Position::new(itv.start.clamp(1, usize::MAX)).unwrap()..=Position::new(itv.stop).unwrap(),
    );
    // Get intervaltree of ignored regions
    let null_itree_ignore = Lapper::new(vec![]);
    let itree_ignore = ignore_bed.get(chrom).unwrap_or(&null_itree_ignore);
    let query = fh.query(&header, &region)?;

    // let outfile_name = format!("{chrom}:{st}-{end}.bam");
    // let outfile = BufWriter::new(File::create_new(&outfile_name)?);
    // let mut out_bam = bam::io::Writer::new(outfile);
    // out_bam.write_header(&header)?;

    for rec in query
        .records()
        .flatten()
        .filter(|aln| !aln.flags().contains(Flags::SECONDARY))
    {
        let rname = rec.name().unwrap();
        let cg: bam::record::Cigar<'_> = rec.cigar();
        let aln_pairs = get_aligned_pairs(
            cg.iter().flatten().map(|op| (op.kind(), op.len())),
            rec.alignment_start().unwrap()?.get(),
        )?;
        let qscores = rec.quality_scores().as_bytes();
        let seq = rec.sequence();
        let is_suppl = rec.flags().contains(Flags::SUPPLEMENTARY);
        // if is_suppl {
        // if verbose {
        //     eprintln!("{:?}\n{:?}", rec.flags(), rec.data())
        // }
        let aln_len = noodles::sam::alignment::Record::alignment_span(&rec).unwrap()? as f64;

        // Look for:
        // * unbalanced reads bordered by large indels. check secondary alignment
        // * supplementary alignments on same chrom (for now)
        let mut mismatch_qpos = vec![];
        for (qpos, refpos, kind) in aln_pairs.into_iter().filter(|(_, refpos, _)| {
            *refpos >= st && *refpos <= end && itree_ignore.count(*refpos, *refpos) == 0
        }) {
            match kind {
                Kind::Insertion => {}
                Kind::Deletion => {}
                Kind::SequenceMismatch => {
                    // 0-93 ASCII+33 for pacbio
                    let qscore = qscores[qpos];
                    if qscore > 30 {
                        let nt = seq.get(qpos).unwrap();
                        mismatch_qpos.push(qpos as f64);
                    }
                }
                _ => {}
            };
        }
        let is_unbalanced = is_unbalanced_alignment(&mismatch_qpos, aln_len, 5, false)?;
        // if is_unbalanced {
        //     // let (st, end) = (
        //     //     rec.alignment_start().unwrap().map(|p| p.get())?,
        //     //     rec.alignment_end().unwrap().map(|p| p.get())?,
        //     // );
        //     // let rname = str::from_utf8(rname)?;
        //     let mut record_buf = RecordBuf::try_from_alignment_record(&header, &rec)?;
        //     let data = record_buf.data_mut();
        //     data.insert(Tag::new(b'U', b'B'), Value::from(1));
        //     out_bam.write_alignment_record(&header, &record_buf)?;
        // } else {
        //     out_bam.write_record(&header, &rec)?;
        // }
    }

    // bam::fs::index(outfile_name)?;
    Ok(())
}

fn main() -> eyre::Result<()> {
    let args = Args::parse();

    let bam = &args.bam;
    let fa = &args.fa;

    let ignore_bed: HashMap<String, Lapper<usize, String>> = args
        .ignore_bed
        .as_ref()
        .map(|bed| {
            let itvs = read_bed(&bed).unwrap_or_default();
            itvs.into_iter()
                .map(|(chrom, itvs)| (chrom, Lapper::new(itvs)))
                .collect()
        })
        .unwrap_or_default();

    if args.ignore_bed.is_some() {
        eprintln!(
            "Loaded {} ignored region(s) across {} chromosome(s).",
            ignore_bed.values().map(|b| b.len()).sum::<usize>(),
            ignore_bed.len()
        );
    }

    let regions = if let Some(bed) = &args.bed {
        read_bed(bed.as_ref()).with_context(|| format!("No valid intervals in {bed:?}"))
    } else {
        let mut fh = bam::io::indexed_reader::Builder::default().build_from_path(bam)?;
        aligned_intervals_windows(&mut fh, args.wg_window)
    }?;

    eprintln!(
        "Computing indel rates across {} chromosome(s).",
        ignore_bed.len()
    );
    let mut indel_rates: HashMap<String, (Vec<f64>, Vec<f64>)> = HashMap::new();
    for region in regions.values().flatten() {
        let (prim_indel_rate, sec_indel_rate) =
            calculate_mean_indel_rate(bam, region, &ignore_bed)?;
        if let Some((prim_indel_rates, sec_indel_rates)) = indel_rates.get_mut(&region.val) {
            prim_indel_rates.push(prim_indel_rate);
            sec_indel_rates.push(sec_indel_rate);
        } else {
            indel_rates.insert(
                region.val.to_owned(),
                (vec![prim_indel_rate], vec![sec_indel_rate]),
            );
        }
    }
    let indel_rates: HashMap<String, (f64, f64)> = indel_rates
        .into_iter()
        .map(|(chrom, (prim_indel_rates, sec_indel_rates))| {
            let prim_indel_rate =
                prim_indel_rates.iter().sum::<f64>() / prim_indel_rates.len() as f64;
            let sec_indel_rate = sec_indel_rates.iter().sum::<f64>() / sec_indel_rates.len() as f64;
            (chrom, (prim_indel_rate, sec_indel_rate))
        })
        .collect();
    eprintln!("{indel_rates:?}");

    eprintln!(
        "Detecting events across {} window(s) in {} chromosome(s).",
        ignore_bed.values().map(|b| b.len()).sum::<usize>(),
        ignore_bed.len()
    );

    for region in regions.values().flatten() {
        let indel_rate = indel_rates[&region.val];
        eprintln!("On {region:?}...");
        detect_events(bam, fa, region, indel_rate, &ignore_bed)?;
    }

    Ok(())
}
