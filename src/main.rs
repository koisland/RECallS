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
    events::{Event, MismatchSignal, SupplJuncSignal, SupplSignal},
    io::{
        FastaHandle, aligned_intervals_windows, read_bed, read_indel_read_stats, read_paf,
        write_indel_read_stats, write_itvs_self_similar_paf,
    },
    self_align::generate_contig_self_alignment,
    tag::tag_bam,
    utils::overlap_length,
};

mod baseline;
mod call;
mod cli;
mod events;
mod io;
mod self_align;
mod tag;
mod unbalanced_aln;
mod utils;

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

    eprintln!(
        "Detecting homologous regions from self-alignment of {} chromosome(s).",
        seq_lens.len()
    );
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
        seq_lens.len()
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
        seq_lens.len()
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
        .flatten()
        .collect();

    let (read_inv_events, chrom_inv_events, chrom_inv_junc_events, chrom_del_events) =
        read_events.into_iter().fold(
            (
                HashMap::<String, Vec<MismatchSignal>>::new(),
                HashMap::<String, Vec<Interval<usize, MismatchSignal>>>::new(),
                HashMap::<String, Vec<Interval<usize, SupplJuncSignal>>>::new(),
                HashMap::<String, Vec<Interval<usize, SupplSignal>>>::new(),
            ),
            |(
                mut read_inv_events,
                mut chrom_inv_events,
                mut chrom_inv_junc_events,
                mut chrom_del_events,
            ),
             event| {
                match event {
                    Event::Deletion(suppl_signal) => {
                        // Add both the main aligned position and the supplementary position
                        let itv_main = Interval {
                            start: suppl_signal.start,
                            stop: suppl_signal.stop,
                            val: suppl_signal.clone(),
                        };
                        let itv_suppl = Interval {
                            start: suppl_signal.suppl_start,
                            stop: suppl_signal.suppl_stop,
                            val: suppl_signal.clone(),
                        };
                        if let Some(chrom_events) = chrom_del_events.get_mut(&suppl_signal.chrom) {
                            chrom_events.push(itv_main);
                            chrom_events.push(itv_suppl)
                        } else {
                            chrom_del_events
                                .insert(suppl_signal.chrom.to_owned(), vec![itv_main, itv_suppl]);
                        }
                    }
                    Event::Inversion(suppl_junc_signal) => {
                        let itv = Interval {
                            start: suppl_junc_signal.start_1,
                            stop: suppl_junc_signal.start_2,
                            val: suppl_junc_signal.clone(),
                        };
                        if let Some(chrom_events) =
                            chrom_inv_junc_events.get_mut(&suppl_junc_signal.chrom)
                        {
                            chrom_events.push(itv);
                        } else {
                            chrom_inv_junc_events
                                .insert(suppl_junc_signal.chrom.to_owned(), vec![itv]);
                        }
                    }
                    Event::InversionInferred(mismatch_signal) => {
                        if let Some(read_events) = read_inv_events.get_mut(&mismatch_signal.rname) {
                            read_events.push(mismatch_signal.clone());
                        } else {
                            read_inv_events.insert(
                                mismatch_signal.rname.to_owned(),
                                vec![mismatch_signal.clone()],
                            );
                        };
                        if let Some(chrom_events) = chrom_inv_events.get_mut(&mismatch_signal.chrom)
                        {
                            chrom_events.push(Interval {
                                start: mismatch_signal.start,
                                stop: mismatch_signal.stop,
                                val: mismatch_signal,
                            });
                        } else {
                            chrom_inv_events.insert(
                                mismatch_signal.chrom.to_owned(),
                                vec![Interval {
                                    start: mismatch_signal.start,
                                    stop: mismatch_signal.stop,
                                    val: mismatch_signal,
                                }],
                            );
                        }
                    }
                }
                (
                    read_inv_events,
                    chrom_inv_events,
                    chrom_inv_junc_events,
                    chrom_del_events,
                )
            },
        );

    let n_del_events = chrom_del_events
        .values()
        .map(|events| events.len())
        .sum::<usize>()
        / 2;
    let itrees_chrom_inv_events: HashMap<String, Lapper<usize, MismatchSignal>> = chrom_inv_events
        .into_iter()
        .map(|(chrom, itvs)| (chrom, Lapper::new(itvs)))
        .collect();
    let itrees_chrom_inv_junc_events: HashMap<String, Lapper<usize, SupplJuncSignal>> =
        chrom_inv_junc_events
            .into_iter()
            .map(|(chrom, itvs)| (chrom, Lapper::new(itvs)))
            .collect();
    let itrees_chrom_del_events: HashMap<String, Lapper<usize, SupplSignal>> = chrom_del_events
        .into_iter()
        .map(|(chrom, itvs)| (chrom, Lapper::new(itvs)))
        .collect();

    eprintln!(
        "Filtering {} candidate inversion and {n_del_events} candidate deletion events.",
        read_inv_events
            .values()
            .map(|events| events.len())
            .sum::<usize>(),
    );

    let mut final_read_inv_events = HashMap::new();
    let mut final_read_del_events = HashMap::new();

    // Must have more than one event per read (ex. sec and primary)
    // Must have at least one primary alignment
    // Must overlap with another self similar region
    let outfile_inv = output_dir.join("calls_inv.bed");
    let mut outfile_inv_fh = BufWriter::new(File::create(outfile_inv)?);
    writeln!(&mut outfile_inv_fh, "{}", MismatchSignal::header())?;

    for (chrom, itree_inv_events) in itrees_chrom_inv_events.iter() {
        let itree_self_similar = itvs_self_similar
            .get(chrom)
            .unwrap_or(&null_itree_self_similar);
        for itv in itree_inv_events.iter() {
            let itv_len = (itv.stop - itv.start) as f64;
            let all_read_events = &read_inv_events[&itv.val.rname];
            // secondary aln check.
            let sec_check =
                all_read_events.len() < 2 || all_read_events.iter().all(|e| e.is_secondary);
            let same_chrom = all_read_events.iter().all(|e| e.chrom == *chrom);
            if sec_check || !same_chrom {
                continue;
            }
            // self-similar regions in genome to this event
            let n_itvs_self_similar_event = itree_self_similar
                .find(itv.start, itv.stop)
                .filter(|itv_self_similar| {
                    let ovl_len = overlap_length(
                        itv.start,
                        itv.stop,
                        itv_self_similar.start,
                        itv_self_similar.stop,
                    ) as f64;
                    (ovl_len / itv_len) > 0.5 && *itv_self_similar.val.dv < 0.05
                })
                .flat_map(|itv_self_similar| {
                    itree_inv_events.find(itv_self_similar.start, itv_self_similar.stop)
                })
                .count();

            // One for itself
            if n_itvs_self_similar_event > 1 {
                // Store events to write to bam
                final_read_inv_events
                    .entry(itv.val.rname.clone())
                    .and_modify(|events: &mut Vec<MismatchSignal>| events.push(itv.val.clone()))
                    .or_insert_with(|| vec![itv.val.clone()]);
                writeln!(&mut outfile_inv_fh, "{}", itv.val.as_bed())?;
            }
        }
    }

    let outfile_inv_junc = output_dir.join("calls_inv_junc.bed");
    let mut outfile_inv_junc_fh = BufWriter::new(File::create(outfile_inv_junc)?);
    writeln!(&mut outfile_inv_junc_fh, "{}", SupplJuncSignal::header())?;
    for itree_inv_junc_events in itrees_chrom_inv_junc_events.values() {
        for itv in itree_inv_junc_events.iter() {
            writeln!(&mut outfile_inv_junc_fh, "{}", itv.val.as_bed())?;
        }
    }

    // Then check overlaps, to be confident require that other end also produces suppl on same side so at least 2 ovl
    let outfile_del = output_dir.join("calls_del.bed");
    let mut outfile_del_fh = BufWriter::new(File::create(outfile_del)?);
    writeln!(&mut outfile_del_fh, "{}", SupplSignal::header())?;

    for itree_del_events in itrees_chrom_del_events.values() {
        for itv in itree_del_events.iter() {
            let ovl_cnt = itree_del_events.count(itv.start, itv.stop);
            if ovl_cnt >= args.del_min_ovl_cnt {
                // Store events to write to bam
                final_read_del_events
                    .entry(itv.val.rname.clone())
                    .and_modify(|events: &mut Vec<SupplSignal>| events.push(itv.val.clone()))
                    .or_insert_with(|| vec![itv.val.clone()]);
                writeln!(&mut outfile_del_fh, "{}", itv.val.as_bed())?;
            }
        }
    }

    if let Some(out_bam) = args.output_bam {
        eprintln!("Generating tagged BAM.",);
        tag_bam(
            &args.bam,
            &out_bam,
            &final_read_inv_events,
            &final_read_del_events,
        )?;
        // bam::fs::index(out_bam)?;
    }

    eprintln!("Done!");

    Ok(())
}
