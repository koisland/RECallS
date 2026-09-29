use std::{collections::HashMap, path::Path};

use clap::Parser;
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
    baseline::{ReadIndelSummaryStats, aggregate_stats_indel_rate},
    cli::Args,
    dotplot::generate_whole_contig_dotplot,
    io::{
        FastaHandle, aligned_intervals_windows, read_bed, read_indel_read_stats, read_paf,
        write_indel_read_stats, write_itvs_self_similar_paf,
    },
    unbalanced_aln::{UnbalancedSummary, is_unbalanced_alignment},
};

mod baseline;
mod cli;
mod dotplot;
mod io;
mod unbalanced_aln;

pub struct InversionEvent {
    chrom: String,
    start: usize,
    stop: usize,
    rname: String,
    n_indels: usize,
    aln_len: f64,
    is_secondary: bool,
    unbalanced_summary: Option<UnbalancedSummary>,
}
impl InversionEvent {
    pub fn as_bed(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:?}",
            self.chrom,
            self.start,
            self.stop,
            self.rname,
            self.n_indels,
            self.aln_len,
            self.is_secondary,
            self.unbalanced_summary
        )
    }
}

pub struct DeletionEvent {}

pub enum Event {
    Deletion(DeletionEvent),
    Inversion(InversionEvent),
}

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

pub struct ReadMarkers {
    pub pos: Vec<f64>,
    pub n_indels: usize,
    pub n_mismatches: usize,
}

fn collect_read_markers(
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
                        marker_qpos.push(qpos as f64);
                    }
                }
                qpos += l;
            }
            Kind::Deletion => {
                if check_valid_itv(ref_pos, st, end, itree_ignore) {
                    for _ in ref_pos..(ref_pos + l) {
                        n_indels += 1;
                        marker_qpos.push(qpos as f64);
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

fn detect_events(
    bam: &Path,
    _fa: &Path,
    itv: &Interval<usize, String>,
    read_stats: &ReadIndelSummaryStats,
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

    let indel_read_stats = [&read_stats.primary, &read_stats.secondary];
    let mut events = vec![];

    for rec in query.records().flatten() {
        let rname = rec.name().unwrap();
        let cg = rec.cigar();
        let is_suppl = rec.flags().contains(Flags::SUPPLEMENTARY);
        let is_sec = rec.flags().contains(Flags::SECONDARY);
        let typ_read_stats = &indel_read_stats[is_sec as usize];
        let aln_len = noodles::sam::alignment::Record::alignment_span(&rec).unwrap()? as f64;

        // Look for:
        // * unbalanced reads bordered by large indels. check secondary alignment
        // * supplementary alignments on same chrom (for now)
        let read_markers = collect_read_markers(&rec, st, end, itree_ignore)?;
        // TODO: use number of snv as filter
        let indel_rate_zscore = typ_read_stats.zscore(read_markers.n_indels as f64 / aln_len);
        let is_unbalanced = is_unbalanced_alignment(&read_markers.pos, aln_len, 5, 0.33)?;
        let (rst, rend) = (
            rec.alignment_start().unwrap().map(|p| p.get())?,
            rec.alignment_end().unwrap().map(|p| p.get())?,
        );
        if indel_rate_zscore > 3.4 && aln_len > 10_000.0 {
            let event = InversionEvent {
                chrom: chrom.to_owned(),
                start: rst,
                stop: rend,
                rname: String::from_utf8(rname.to_vec())?,
                n_indels: read_markers.n_indels,
                aln_len,
                is_secondary: is_sec,
                unbalanced_summary: is_unbalanced,
            };
            events.push(Event::Inversion(event));
        }
        if is_suppl {
            let Value::String(sa_tag) = rec
                .data()
                .get(b"SA")
                .with_context(|| format!("Must have SA tag for {rname}."))??
            else {
                bail!("Invalid type for SA tag for {rname}.")
            };
            let clip_direction = get_clip_direction(cg.iter().flatten()).with_context(|| {
                format!("Read {rname} must have soft/hardclipped operation in cigar: {cg:?}")
            })?;

            // minimap2 v2.28 SA tag format
            //
            // chrom,start,strand,cigar,mapq,num_mismatches_gaps
            // chr7,2441699,-,12760S18196M69I,60,120
            for (sa_chrom, sa_start, sa_strand, sa_cigar, sa_mapq, sa_nm) in str::from_utf8(sa_tag)?
                .split(';')
                .flat_map(|rec| rec.splitn(6, ',').collect_tuple::<SARecord>())
                // Must be same chromosome
                .filter(|sa_rec| sa_rec.0 == chrom)
            {
                let sa_start: usize = sa_start.parse()?;

                // What is the order of the current alignment relative to the suppl alignment?
                let sa_upstream = sa_start > rend;

                let sa_cigar = noodles::sam::record::Cigar::new(sa_cigar.as_bytes());
                let sa_aln_len = sa_cigar.alignment_span()?;
                let sa_clip_direction = get_clip_direction(sa_cigar.iter().flatten())
                    .with_context(|| format!("Read {rname} must have soft/hardclipped operation in SA cigar: {sa_cigar:?}"))?;

                match (sa_upstream, &clip_direction, sa_clip_direction) {
                    // Invalid
                    (_, ClipDirection::Left, ClipDirection::Left) |
                    (_, ClipDirection::Right, ClipDirection::Right) |
                    (_, ClipDirection::Both, _) |
                    (_, _, ClipDirection::Both) |
                    //       *
                    // |<  | |  >|
                    (true, ClipDirection::Left, ClipDirection::Right) |
                    // *     
                    // |<  | |  >|
                    (false, ClipDirection::Right, ClipDirection::Left) => {
                        continue;
                    },
                    //       *
                    // |  >| |<  |
                    (true, ClipDirection::Right, ClipDirection::Left) => {
                        // println!("{chrom}\t{rst}\t{rend}\t{sa_chrom}:{sa_start}-{}", sa_start+sa_aln_len)
                    },
                    // *     
                    // |  >| |<  |
                    (false, ClipDirection::Left, ClipDirection::Right) => {
                        // println!("{chrom}\t{rst}\t{rend}\t{sa_chrom}:{sa_start}-{}", sa_start+sa_aln_len)
                    }
                }
            }
        }
    }

    Ok(events)
}

fn main() -> eyre::Result<()> {
    let args = Args::parse();

    let bam = &args.bam;
    let fa = &args.fa;
    let output_dir = &args.output_dir;
    std::fs::create_dir_all(output_dir)?;

    let ignore_bed: HashMap<String, Lapper<usize, String>> = args
        .ignore_bed
        .as_ref()
        .map(|bed| {
            let itvs = read_bed(bed).unwrap_or_default();
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
        "Detecting homologous regions from self-alignment of {} chromosome(s).",
        ignore_bed.len()
    );
    // Generate dotplot per contig
    let mut fh = FastaHandle::new(fa)?;
    let seq_lens: HashMap<String, usize> = fh
        .fai
        .as_ref()
        .iter()
        .map(|rec| {
            let ctg_name = str::from_utf8(rec.name()).unwrap();
            let ctg_len = rec.length() as usize;
            (ctg_name.to_owned(), ctg_len)
        })
        .collect();

    let paf_self_align = output_dir.join("chrom_self_align.paf");
    let itvs_self_similar = if !paf_self_align.exists() {
        let itvs_self_similar = generate_whole_contig_dotplot(&mut fh, &seq_lens)?;
        write_itvs_self_similar_paf(&itvs_self_similar, &paf_self_align)?;
        itvs_self_similar
    } else {
        read_paf(&paf_self_align)?
    };

    eprintln!(
        "Computing indel rates across {} chromosome(s).",
        ignore_bed.len()
    );
    let tsv_indel_stats = output_dir.join("chrom_indel_stats.tsv");
    let chrom_read_stats = if !tsv_indel_stats.exists() {
        let indel_stats = aggregate_stats_indel_rate(&regions, bam, &ignore_bed);
        write_indel_read_stats(&indel_stats, &tsv_indel_stats)?;
        indel_stats
    } else {
        read_indel_read_stats(&tsv_indel_stats)?
    };

    eprintln!(
        "Detecting events across {} window(s) in {} chromosome(s).",
        ignore_bed.values().map(|b| b.len()).sum::<usize>(),
        ignore_bed.len()
    );

    let mut read_inv_events: HashMap<String, Vec<InversionEvent>> = HashMap::new();
    for region in regions.values().flatten() {
        let read_stats = &chrom_read_stats[&region.val];
        eprintln!("On {region:?}...");
        let events = detect_events(bam, fa, region, read_stats, &ignore_bed)?;
        for event in events {
            match event {
                Event::Deletion(deletion_event) => todo!(),
                Event::Inversion(inversion_event) => {
                    if let Some(read_events) = read_inv_events.get_mut(&inversion_event.rname) {
                        read_events.push(inversion_event);
                    } else {
                        read_inv_events
                            .insert(inversion_event.rname.to_owned(), vec![inversion_event]);
                    }
                }
            }
        }
    }
    // Must have more than one event per read
    // Must have at least one primary alignment
    read_inv_events.retain(|_, v| v.len() > 1 && v.iter().any(|e| !e.is_secondary));

    for (_, events) in read_inv_events {
        for event in events {
            println!("{}", event.as_bed())
        }
    }

    Ok(())
}
