use std::{
    collections::HashMap,
    fs::File,
    io::{BufWriter, Write},
};

use clap::Parser;
use eyre::ContextCompat;
use noodles::bam::{self};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rust_lapper::{Interval, Lapper};

use crate::{
    baseline::aggregate_stats_indel_rate,
    call::detect_events,
    cli::Args,
    events::{DeletionEvent, Event, InversionEvent},
    io::{
        FastaHandle, aligned_intervals_windows, read_bed, read_indel_read_stats, read_paf,
        write_indel_read_stats, write_itvs_self_similar_paf,
    },
    self_align::generate_contig_self_alignment,
};

mod baseline;
mod call;
mod cli;
mod events;
mod io;
mod self_align;
mod unbalanced_aln;

fn main() -> eyre::Result<()> {
    let args = Args::parse();

    let bam = &args.bam;
    let fa = &args.fa;
    let output_dir = &args.output_dir;
    std::fs::create_dir_all(output_dir)?;

    // set up global threadpool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;

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
    let null_itree_ignore = Lapper::new(vec![]);

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
    let fh = FastaHandle::new(fa)?;
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
    std::mem::drop(fh);

    let paf_self_align = output_dir.join("chrom_self_align.paf");
    let itvs_self_similar = if !paf_self_align.exists() {
        let itvs_self_similar = generate_contig_self_alignment(fa, &seq_lens, args.del_max_rgn_dv)?;
        write_itvs_self_similar_paf(&itvs_self_similar, &paf_self_align)?;
        itvs_self_similar
    } else {
        read_paf(&paf_self_align, args.del_max_rgn_dv)?
    };
    let null_itree_self_similar = Lapper::new(vec![]);

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

    let intervals: Vec<&Interval<usize, String>> = regions.values().flatten().collect();
    // let mut read_inv_events: HashMap<String, Vec<InversionEvent>> = HashMap::new();
    let read_events: Vec<Event> = intervals
        .into_par_iter()
        .flat_map(|region| {
            let read_stats = &chrom_read_stats[&region.val];
            let itree_ignore_chrom = ignore_bed.get(&region.val).unwrap_or(&null_itree_ignore);
            let itree_self_similar_chrom = &itvs_self_similar
                .get(&region.val)
                .unwrap_or(&null_itree_self_similar);

            eprintln!("On {}:{}-{}...", region.val, region.start, region.stop);
            let events = detect_events(
                bam,
                region,
                read_stats,
                itree_ignore_chrom,
                itree_self_similar_chrom,
                args.inv_indel_zscore as f64,
                args.inv_min_aln_len,
                args.del_min_mapq,
            );
            match events {
                Ok(events) => Some(events),
                Err(err) => {
                    eprintln!("Failed on {region:?}: {err}");
                    None
                }
            }
        })
        .flatten()
        .collect();

    let (mut read_inv_events, read_itvs_del_events) = read_events.into_iter().fold(
        (HashMap::<String, Vec<InversionEvent>>::new(), Vec::new()),
        |(mut read_inv_events, mut itvs_del_events), event| {
            match event {
                Event::Deletion(deletion_event) => {
                    // Add both the main aligned position and the supplementary position
                    itvs_del_events.push(Interval {
                        start: deletion_event.start,
                        stop: deletion_event.stop,
                        val: deletion_event.clone(),
                    });
                    itvs_del_events.push(Interval {
                        start: deletion_event.suppl_start,
                        stop: deletion_event.suppl_stop,
                        val: deletion_event.clone(),
                    })
                }
                Event::Inversion(inversion_event) => {
                    if let Some(read_events) = read_inv_events.get_mut(&inversion_event.rname) {
                        read_events.push(inversion_event);
                    } else {
                        read_inv_events
                            .insert(inversion_event.rname.to_owned(), vec![inversion_event]);
                    }
                }
            }
            (read_inv_events, itvs_del_events)
        },
    );

    eprintln!(
        "Filtering {} candidate inversion and {} candidate deletion events.",
        read_inv_events
            .values()
            .map(|events| events.len())
            .sum::<usize>(),
        read_itvs_del_events.len() / 2
    );

    // Must have more than one event per read
    // Must have at least one primary alignment
    let outfile_inv = output_dir.join("calls_inv.bed");
    let mut outfile_inv_fh = BufWriter::new(File::create(outfile_inv)?);

    writeln!(&mut outfile_inv_fh, "{}", InversionEvent::header())?;

    read_inv_events.retain(|_, v| v.len() > 1 && v.iter().any(|e| !e.is_secondary));

    for (_, events) in read_inv_events {
        for event in events {
            writeln!(&mut outfile_inv_fh, "{}", event.as_bed())?;
        }
    }

    // Then check overlaps, to be confident require that other end also produces suppl on same side so at least 2 ovl
    let outfile_del = output_dir.join("calls_del.bed");
    let mut outfile_del_fh = BufWriter::new(File::create(outfile_del)?);
    let read_del_events = Lapper::new(read_itvs_del_events);

    writeln!(&mut outfile_del_fh, "{}", DeletionEvent::header())?;

    for itv in &read_del_events {
        let ovl_cnt = read_del_events.count(itv.start, itv.stop);
        if ovl_cnt >= args.del_min_ovl_cnt {
            writeln!(&mut outfile_del_fh, "{}", itv.val.as_bed())?;
        }
    }

    eprintln!("Done!");

    Ok(())
}
