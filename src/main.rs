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
    events::{Event, Signal},
    io::{aligned_intervals_windows, read_bed, read_indel_read_stats, write_indel_read_stats},
    postprocess::postprocess_events,
    tag::create_tagged_bam,
};

mod baseline;
mod call;
mod cli;
mod events;
mod io;
mod postprocess;
mod tag;
mod unbalanced_aln;
mod utils;

fn main() -> eyre::Result<()> {
    let args = Args::parse();
    simple_logger::init_with_level(args.log_level)?;

    let bam = &args.bam;
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
        log::info!(
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

    let tsv_indel_stats = output_dir.join("chrom_indel_stats.tsv");
    let chrom_read_stats = if !tsv_indel_stats.exists() {
        log::info!(
            "Computing indel rates across {} chromosome(s).",
            regions.len()
        );
        let indel_stats = aggregate_stats_indel_rate(bam, &regions, &ignore_bed);
        write_indel_read_stats(&indel_stats, &tsv_indel_stats)?;
        indel_stats
    } else {
        log::info!("Loaded existing indel rates from {tsv_indel_stats:?}.");
        read_indel_read_stats(&tsv_indel_stats)?
    };

    log::info!(
        "Detecting events across {} window(s) in {} chromosome(s).",
        regions.values().map(|b| b.len()).sum::<usize>(),
        regions.len()
    );

    let intervals: Vec<&Interval<usize, String>> = regions.values().flatten().collect();
    // let mut read_inv_events: HashMap<String, Vec<InversionEvent>> = HashMap::new();
    let read_events: HashMap<String, Vec<Event>> = intervals
        .into_par_iter()
        .flat_map(|region| {
            let read_stats = &chrom_read_stats[&region.val];
            let itree_ignore_chrom = ignore_bed.get(&region.val).unwrap_or(&null_itree_ignore);

            log::info!("On {}:{}-{}...", region.val, region.start, region.stop);
            let events = detect_events(
                bam,
                region,
                read_stats,
                itree_ignore_chrom,
                args.inv_indel_zscore as f64,
                args.inv_min_aln_len,
                args.inv_thr_unbalanced as f64,
                args.inv_min_num_snvs,
                args.del_min_mapq,
            );
            match events {
                Ok(events) => Some(events),
                Err(err) => {
                    log::error!("Failed on {region:?}: {err}");
                    None
                }
            }
        })
        .fold(
            HashMap::new,
            |mut acc: HashMap<String, Vec<Event>>, events| {
                // empty window with no events
                if events.is_empty() {
                    return acc;
                }
                // Store by chromosome
                let chrom = match &events[0] {
                    Event::SV(sv) => &sv.itv.val,
                    Event::Signal(Signal::MismatchSignal(msignal)) => &msignal.itv.val,
                    Event::Signal(Signal::SupplSignal(ssignal)) => &ssignal.chrom,
                };

                if let Some(chrom_events) = acc.get_mut(chrom) {
                    chrom_events.extend(events);
                } else {
                    acc.insert(chrom.to_owned(), events);
                }

                acc
            },
        )
        .reduce(HashMap::new, |mut a, b| {
            // Reduce hashmap per window of chrom events into one hashmap
            for (chrom, other_events) in b.into_iter() {
                if let Some(events) = a.get_mut(&chrom) {
                    events.extend(other_events);
                } else {
                    a.insert(chrom, other_events);
                }
            }
            a
        });

    log::info!(
        "Postprocessing and filtering {} candidate events.",
        read_events.values().map(|b| b.len()).sum::<usize>(),
    );
    let sv_events = postprocess_events(read_events);
    let output_bed = output_dir.join("calls.bed");
    let mut output_bed_fh = BufWriter::new(File::create(output_bed)?);
    for events in sv_events.values() {
        for event in events {
            writeln!(&mut output_bed_fh, "{}", event.as_bed())?;
        }
    }

    if let Some(out_bam) = args.output_bam {
        log::info!("Generating tagged BAM.",);
        create_tagged_bam(&args.bam, &out_bam, &sv_events)?;
    }

    log::info!("Done!");

    Ok(())
}
