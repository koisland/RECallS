use std::collections::{HashMap, HashSet};

use clap::Parser;
use eyre::ContextCompat;
use itertools::Itertools;
use noodles::bam::{self};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rust_lapper::{Interval, Lapper};

use crate::{
    baseline::aggregate_stats_indel_rate,
    call::detect_events,
    cli::Args,
    events::{Event, MismatchSignal, SVType, Signal},
    io::{aligned_intervals_windows, read_bed, read_indel_read_stats, write_indel_read_stats},
    utils::merge_intervals,
};

mod baseline;
mod call;
mod cli;
mod events;
mod io;
mod tag;
mod unbalanced_aln;
mod utils;

fn main() -> eyre::Result<()> {
    let args = Args::parse();

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
        regions.len()
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

            eprintln!("On {}:{}-{}...", region.val, region.start, region.stop);
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
                    eprintln!("Failed on {region:?}: {err}");
                    None
                }
            }
        })
        .fold(HashMap::new, |mut acc: HashMap<String, Vec<Event>>, x| {
            // empty window with no events
            if x.is_empty() {
                return acc;
            }
            let chrom = match &x[0] {
                Event::SV(sv) => &sv.itv.val,
                Event::Signal(Signal::MismatchSignal(msignal)) => &msignal.itv.val,
                Event::Signal(Signal::SupplSignal(ssignal)) => &ssignal.chrom,
            };

            if let Some(events) = acc.get_mut(chrom) {
                events.extend(x);
            } else {
                acc.insert(chrom.to_owned(), x);
            }

            acc
        })
        .reduce(HashMap::new, |mut a, b| {
            for (chrom, other_events) in b.into_iter() {
                if let Some(events) = a.get_mut(&chrom) {
                    events.extend(other_events);
                } else {
                    a.insert(chrom, other_events);
                }
            }
            a
        });

    // Merge overlapping reads
    for (chrom, events) in read_events.iter() {
        // Complete events
        let mut inversion_events: Vec<Interval<usize, HashSet<&str>>> = vec![];
        let mut deletion_events: Vec<Interval<usize, HashSet<&str>>> = vec![];

        // Create itree of events split by type
        let mut suppl_signals: Vec<Interval<usize, &str>> = vec![];
        let mut mismatch_signals: Vec<Interval<usize, &MismatchSignal>> = vec![];
        let mut read_mismatch_signals: HashMap<&str, Vec<Interval<usize, &MismatchSignal>>> =
            HashMap::new();
        for event in events.iter() {
            match event {
                Event::SV(sv) => match sv.typ {
                    SVType::Deletion => {
                        deletion_events.push(Interval {
                            start: sv.itv.start,
                            stop: sv.itv.stop,
                            val: HashSet::from_iter([&sv.read[..]]),
                        });
                    }
                    SVType::Inversion => {
                        inversion_events.push(Interval {
                            start: sv.itv.start,
                            stop: sv.itv.stop,
                            val: HashSet::from_iter([&sv.read[..]]),
                        });
                    }
                },
                Event::Signal(Signal::MismatchSignal(msignal)) => {
                    mismatch_signals.push(Interval {
                        start: msignal.itv.start,
                        stop: msignal.itv.stop,
                        val: msignal,
                    });
                    read_mismatch_signals
                        .entry(&msignal.read)
                        .and_modify(|itv| {
                            itv.push(Interval {
                                start: msignal.itv.start,
                                stop: msignal.itv.stop,
                                val: msignal,
                            })
                        })
                        .or_insert_with(|| {
                            vec![Interval {
                                start: msignal.itv.start,
                                stop: msignal.itv.stop,
                                val: msignal,
                            }]
                        });
                }
                Event::Signal(Signal::SupplSignal(ssignal)) => {
                    suppl_signals.push(Interval {
                        start: ssignal.itv.start,
                        stop: ssignal.itv.stop,
                        val: &ssignal.read,
                    });
                }
            }
        }
        // Merge complete events
        let merged_inversion_events = merge_intervals(
            inversion_events.into_iter(),
            1,
            |_a, _b| true,
            |a, b| a.val.iter().chain(b.val.iter()).cloned().collect(),
        );
        let merged_deletion_events = merge_intervals(
            deletion_events.into_iter(),
            1,
            |_a, _b| true,
            |a, b| a.val.iter().chain(b.val.iter()).cloned().collect(),
        );
        // for itv in merged_inversion_events.iter() {
        //     println!(
        //         "{chrom}\t{}\t{}\t{}",
        //         itv.start,
        //         itv.stop,
        //         itv.val.iter().sorted().join(",")
        //     )
        // }
        // for itv in merged_deletion_events.iter() {
        //     println!(
        //         "{chrom}\t{}\t{}\t{}",
        //         itv.start,
        //         itv.stop,
        //         itv.val.iter().sorted().join(",")
        //     )
        // }

        // Partition by self-similar reads
        let itree_mismatch_events = Lapper::new(mismatch_signals);
        let mut eval_reads: HashSet<&str> = HashSet::new();
        for itv_msignal in itree_mismatch_events.iter() {
            if eval_reads.contains(&itv_msignal.val.read[..]) {
                continue;
            }

            let other_msignals_same_read = &read_mismatch_signals[&itv_msignal.val.read[..]];

            // Check overlap of current position
            // Do they have no unbalanced overlap
            let mut ovl = itree_mismatch_events
                .find(itv_msignal.start, itv_msignal.stop)
                .collect_vec();
            let ovl_reads: HashSet<&str> = ovl.iter().map(|itv| &itv.val.read[..]).collect();

            // Check overlap with other positions
            // Do they have no unbalanced overlap
            let mut ovl_sec = other_msignals_same_read
                .iter()
                .filter(|itv| *itv != itv_msignal)
                .flat_map(|itv| itree_mismatch_events.find(itv.start, itv.stop))
                .collect_vec();
            let ovl_sec_reads: HashSet<&str> =
                ovl_sec.iter().map(|itv| &itv.val.read[..]).collect();

            let shared_reads: HashSet<_> = ovl_reads.intersection(&ovl_sec_reads).collect();
            if shared_reads.len() != 2 {
                continue;
            }
            let shared_reads_str = shared_reads.iter().join("/");

            ovl.retain(|itv| shared_reads.contains(&&itv.val.read[..]));
            ovl_sec.retain(|itv| shared_reads.contains(&&itv.val.read[..]));

            eval_reads.extend(shared_reads);

            let n_dir_uniq = ovl
                .iter()
                .map(|itv| itv.val.unbalanced_dir)
                .unique()
                .count();
            let n_dir_uniq_sec = ovl
                .iter()
                .map(|itv| itv.val.unbalanced_dir)
                .unique()
                .count();

            if n_dir_uniq != 2 || n_dir_uniq_sec != 2 {
                continue;
            }

            let has_prim = ovl.iter().any(|itv| !itv.val.is_secondary);
            if !has_prim {
                continue;
            }

            println!(
                "{chrom}\t{}\t{}\t{chrom}\t{}\t{}\t{}",
                ovl[0].start, ovl[0].stop, ovl_sec[0].start, ovl_sec[0].stop, shared_reads_str,
            )

            // // println!("1.");
            // for itv in ovl.iter() {
            //     println!(
            //         "{}\t{}\t{}\t{:?}-{}-{}",
            //         chrom,
            //         itv.start,
            //         itv.stop,
            //         itv.val.unbalanced_dir,
            //         itv.val.is_secondary,
            //         itv.val.read
            //     )
            // }
            // // println!("2.");
            // for itv in ovl_sec.iter() {
            //     println!(
            //         "{}\t{}\t{}\t{:?}-{}-{}",
            //         chrom,
            //         itv.start,
            //         itv.stop,
            //         itv.val.unbalanced_dir,
            //         itv.val.is_secondary,
            //         itv.val.read
            //     )
            // }
        }

        // let itree_suppl_signals = Lapper::new(suppl_signals);
    }

    // if let Some(out_bam) = args.output_bam {
    //     eprintln!("Generating tagged BAM.",);
    //     tag_bam(
    //         &args.bam,
    //         &out_bam,
    //         &final_read_inv_events,
    //         &final_read_del_events,
    //     )?;
    //     // bam::fs::index(out_bam)?;
    // }

    eprintln!("Done!");

    Ok(())
}
