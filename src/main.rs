use std::{collections::HashMap, fs::File, io::BufWriter, path::Path};

use clap::Parser;
use eyre::ContextCompat;
use noodles::{
    bam::{self},
    core::{Position, Region},
    sam::alignment::{
        Record, RecordBuf,
        io::Write,
        record::{Flags, cigar::op::Kind, data::field::Tag},
        record_buf::data::field::Value,
    },
};
use rust_lapper::{Interval, Lapper};

use crate::{
    baseline::{ReadSummaryStats, calculate_stats_indel_rate}, cli::Args, io::{aligned_intervals_windows, read_bed}, unbalanced_aln::{UnbalancedSummary, is_unbalanced_alignment}, utils::get_aligned_pairs,
};

mod baseline;
mod cli;
mod dotplot;
mod io;
mod unbalanced_aln;
mod utils;

struct Event {
    chrom: String,
    start: usize,
    stop: usize,
    rname: String,
    n_indels: usize,
    aln_len: f64,
    is_secondary: bool,
    unbalanced_summary: Option<UnbalancedSummary>
}
impl Event {
    fn as_bed(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:?}",
            self.chrom, self.start, self.stop, self.rname, self.n_indels, self.aln_len, self.is_secondary, self.unbalanced_summary
        )
    }
}

fn detect_events(
    bam: &Path,
    _fa: &Path,
    itv: &Interval<usize, String>,
    read_stats: &ReadSummaryStats,
    ignore_bed: &HashMap<String, Lapper<usize, String>>,
) -> eyre::Result<Vec<Event>> {
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

    let indel_read_stats = [&read_stats.primary, &read_stats.secondary];
    let mut events = vec![];

    for rec in query
        .records()
        .flatten()
    {
        let rname = rec.name().unwrap();
        let cg: bam::record::Cigar<'_> = rec.cigar();
        let aln_pairs = get_aligned_pairs(
            cg.iter().flatten().map(|op| (op.kind(), op.len())),
            rec.alignment_start().unwrap()?.get(),
        )?;
        let qscores = rec.quality_scores().as_bytes();
        let is_suppl = rec.flags().contains(Flags::SUPPLEMENTARY);
        let is_sec = rec.flags().contains(Flags::SECONDARY);
        let typ_read_stats = &indel_read_stats[is_sec as usize];
        // if is_suppl {
        //     eprintln!("{:?}\n{:?}", rec.flags(), rec.data())
        // }
        let aln_len = noodles::sam::alignment::Record::alignment_span(&rec).unwrap()? as f64;

        // Look for:
        // * unbalanced reads bordered by large indels. check secondary alignment
        // * supplementary alignments on same chrom (for now)
        let mut marker_qpos = vec![];
        let mut n_indels: usize = 0;
        for (qpos, _, kind) in aln_pairs.into_iter().filter(|(_, refpos, _)| {
            *refpos >= st && *refpos <= end && itree_ignore.count(*refpos, *refpos) == 0
        }) {
            match kind {
                Kind::Insertion | Kind::Deletion => {
                    n_indels += 1;
                }
                Kind::SequenceMismatch => {
                    // 0-93 ASCII+33 for pacbio
                    let Some(qscore) = qscores.get(qpos) else {
                        continue;
                    };
                    if *qscore > 30 {
                        marker_qpos.push(qpos as f64);
                    }
                }
                _ => {}
            };
        }

        let indel_rate_zscore = typ_read_stats.zscore(n_indels as f64 / aln_len);
        let is_unbalanced = is_unbalanced_alignment(&marker_qpos, aln_len, 5)?;

        if indel_rate_zscore > 3.4 && aln_len > 10_000.0 {
            let (rst, rend) = (
                rec.alignment_start().unwrap().map(|p| p.get())?,
                rec.alignment_end().unwrap().map(|p| p.get())?,
            );
            let event = Event {
                chrom: chrom.to_owned(),
                start: rst,
                stop: rend,
                rname: String::from_utf8(rname.to_vec())?,
                n_indels,
                aln_len,
                is_secondary: is_sec,
                unbalanced_summary: is_unbalanced,
            };
            events.push(event);
        }
        // if is_unbalanced {
        //     let mut record_buf = RecordBuf::try_from_alignment_record(&header, &rec)?;
        //     let data = record_buf.data_mut();
        //     data.insert(Tag::new(b'U', b'B'), Value::from(1));
        //     out_bam.write_alignment_record(&header, &record_buf)?;
        // } else {
        //     out_bam.write_record(&header, &rec)?;
        // }
    }

    // bam::fs::index(outfile_name)?;
    Ok(events)
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

    // https://stats.stackexchange.com/a/26647
    let mut chrom_read_stats = regions
        .values()
        .flatten()
        .map(|region| {
            (
                region.val.clone(),
                calculate_stats_indel_rate(bam, region, &ignore_bed).unwrap(),
            )
        })
        .fold(
            HashMap::new(),
            |mut acc: HashMap<String, ReadSummaryStats>, (chrom, (prim_stats, sec_stats))| {
                if let Some(read_stats) = acc.get_mut(&chrom) {
                    read_stats.primary.mean =
                        read_stats.primary.mean.algebraic_add(prim_stats.mean);
                    read_stats.primary.var = read_stats.primary.var.algebraic_add(prim_stats.var);
                    read_stats.primary.n += prim_stats.n;
                    read_stats.secondary.mean =
                        read_stats.secondary.mean.algebraic_add(sec_stats.mean);
                    read_stats.secondary.var =
                        read_stats.secondary.var.algebraic_add(sec_stats.var);
                    read_stats.secondary.n += sec_stats.n;
                } else {
                    acc.insert(
                        chrom.to_owned(),
                        ReadSummaryStats {
                            primary: prim_stats,
                            secondary: sec_stats,
                        },
                    );
                }
                acc
            },
        );
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

    eprintln!(
        "Detecting events across {} window(s) in {} chromosome(s).",
        ignore_bed.values().map(|b| b.len()).sum::<usize>(),
        ignore_bed.len()
    );

    let mut read_events: HashMap<String, Vec<Event>> = HashMap::new();
    for region in regions.values().flatten() {
        let read_stats = &chrom_read_stats[&region.val];
        eprintln!("On {region:?}...");
        let events = detect_events(bam, fa, region, read_stats, &ignore_bed)?;
        for event in events {
            if let Some(read_events) = read_events.get_mut(&event.rname) {
                read_events.push(event);
            } else {
                read_events.insert(event.rname.to_owned(), vec![event]);
            }
        }
    }
    // Must have more than one event per read 
    // Must have at least one primary alignment
    read_events.retain(|_, v| v.len() > 1 && v.iter().any(|e| !e.is_secondary));

    for (_, events) in read_events {
        for event in events {
            println!("{}", event.as_bed())
        }
    }

    Ok(())
}
