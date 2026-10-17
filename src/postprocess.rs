use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use rust_lapper::{Interval, Lapper};

use crate::{
    events::{Event, MismatchSignal, SV, SVType, Signal},
    utils::merge_intervals,
};

/// Postprocess candidate signals/events.
/// 1. Store mismatch signals by read name for secondary alignments
/// 2. Merge complete events split by type found in a single window. (suppl inversions/deletions)
/// 3. For each msignal (inversion):
///     1. Find overlaps between other mismatch signals 
///     2. For the same read msignal, get secondary msignals and find overlaps between other mismatch signals 
///     3. Check that both overlap sets (primary/secondary) share at least two reads.
///     4. Check that at least two unique unbalanced directions in both sets. Likely complentary
/// 
/// # Arguments
/// * read_events: map of chromosome and events
/// 
/// # Returns
/// * map of chromsome and SVs.
pub fn postprocess_events(read_events: HashMap<String, Vec<Event>>) -> HashMap<String, Vec<SV>> {
    let mut chrom_svs = HashMap::new();

    for (chrom, events) in read_events {
        let mut svs = vec![];
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
                            val: HashSet::from_iter(sv.reads()),
                        });
                    }
                    SVType::Inversion => {
                        inversion_events.push(Interval {
                            start: sv.itv.start,
                            stop: sv.itv.stop,
                            val: HashSet::from_iter(sv.reads()),
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
        // Expect only two.
        let merged_inversion_events = merge_intervals(
            inversion_events.into_iter(),
            1,
            |_a, _b| true,
            |a, b| a.val.extend(b.val.iter()),
        );
        let merged_deletion_events = merge_intervals(
            deletion_events.into_iter(),
            1,
            |_a, _b| true,
            |a, b| a.val.extend(b.val.iter()),
        );
        svs.extend(merged_inversion_events.into_iter().map(|itv| {
            SV {
                itv: Interval {
                    start: itv.start,
                    stop: itv.stop,
                    val: chrom.clone(),
                },
                reads: itv
                    .val
                    .into_iter()
                    .map(|read| Some(read.to_owned()))
                    .collect_array()
                    .expect("More than 2 reads in merged inversion event"),
                typ: SVType::Inversion,
            }
        }));
        svs.extend(merged_deletion_events.into_iter().map(|itv| {
            SV {
                itv: Interval {
                    start: itv.start,
                    stop: itv.stop,
                    val: chrom.clone(),
                },
                reads: itv
                    .val
                    .into_iter()
                    .map(|read| Some(read.to_owned()))
                    .collect_array()
                    .expect("More than 2 reads in merged deletion event"),
                typ: SVType::Deletion,
            }
        }));

        // TODO
        log::debug!("Supplementary signal length: {}",suppl_signals.len());

        let itree_mismatch_events = Lapper::new(mismatch_signals);
        // Ignore read if used already
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

            // Needs to share primary/secondary reads
            let shared_reads: HashSet<_> = ovl_reads.intersection(&ovl_sec_reads).collect();
            if shared_reads.len() < 2 {
                continue;
            }

            ovl.retain(|itv| shared_reads.contains(&&itv.val.read[..]));
            ovl_sec.retain(|itv| shared_reads.contains(&&itv.val.read[..]));

            // Of the shared reads, are they complementary in their balance direction?
            // Is there a read that is unbalanced on the left and one on the right?
            // TODO: Should probably be sorted and check that makes sense
            let n_dir_uniq = ovl
                .iter()
                .map(|itv| itv.val.unbalanced_dir)
                .unique()
                .count();
            let n_dir_uniq_sec = ovl_sec
                .iter()
                .map(|itv| itv.val.unbalanced_dir)
                .unique()
                .count();

            if n_dir_uniq != 2 || n_dir_uniq_sec != 2 {
                continue;
            }

            // Is there at least one primary alignment in this pairing
            let has_prim = ovl.iter().any(|itv| !itv.val.is_secondary);
            if !has_prim {
                continue;
            }

            let min_stop = std::cmp::min(ovl[0].stop, ovl_sec[0].stop);
            let max_start = std::cmp::max(ovl[0].start, ovl_sec[0].start);
            let Some(reads) = shared_reads
                .iter()
                .sorted()
                .map(|read| Some(read.to_string()))
                .collect_array()
            else {
                log::debug!("More than 2 reads: {shared_reads:?} for {itv_msignal:?}");
                continue;
            };
            let inversion = SV {
                itv: Interval {
                    start: min_stop,
                    stop: max_start,
                    val: chrom.clone(),
                },
                reads,
                typ: SVType::Inversion,
            };
            svs.push(inversion);

            eval_reads.extend(shared_reads);
        }

        chrom_svs.insert(chrom, svs);
    }
    chrom_svs
}
