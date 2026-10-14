use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

use itertools::Itertools;
use noodles::{
    bam::io::indexed_reader,
    bgzf::{self},
};
use rust_lapper::Interval;

use crate::baseline::{IndelSummaryStats, ReadIndelSummaryStats};

pub fn read_bed(bed: &Path) -> Option<HashMap<String, Vec<Interval<usize, String>>>> {
    let mut intervals: HashMap<String, Vec<Interval<usize, String>>> = HashMap::new();
    let bed_fh = File::open(bed).expect("Cannot open bedfile");
    let bed_reader = BufReader::new(bed_fh);

    for line in bed_reader.lines() {
        let Ok(line) = line else {
            log::error!("Invalid line: '{line:?}'");
            continue;
        };
        let (name, start, stop, _other_cols) =
            if let Some((name, start, stop, other_cols)) = line.splitn(4, '\t').collect_tuple() {
                (name, start, stop, other_cols)
            } else if let Some((name, start, stop)) = line.splitn(3, '\t').collect_tuple() {
                (name, start, stop, "")
            } else {
                log::error!("Invalid line: '{line}'");
                continue;
            };
        let (Ok(start), Ok(stop)) = (start.parse::<usize>(), stop.parse::<usize>()) else {
            log::error!("Cannot parse {start} or {stop} in line: '{line}'");
            continue;
        };
        let itv = Interval {
            start,
            stop,
            val: name.to_owned(),
        };
        if let Some(itvs) = intervals.get_mut(name) {
            itvs.push(itv);
        } else {
            intervals.insert(name.to_owned(), vec![itv]);
        }
    }
    Some(intervals)
}

pub fn read_indel_read_stats(tsv: &Path) -> eyre::Result<HashMap<String, ReadIndelSummaryStats>> {
    let fh_indel_read_stats = BufReader::new(File::open(tsv)?);
    let mut indel_read_stats: HashMap<String, ReadIndelSummaryStats> = HashMap::new();

    for line in fh_indel_read_stats.lines() {
        let line = line?;
        let Some(
            [
                chrom,
                mean,
                var,
                stdev,
                n,
                sec_mean,
                sec_var,
                sec_stdev,
                sec_n,
            ],
        ) = line.split('\t').collect_array()
        else {
            continue;
        };
        let mean = mean.parse::<f64>()?;
        let var = var.parse::<f64>()?;
        let stdev = stdev.parse::<f64>()?;
        let n = n.parse::<usize>()?;
        let sec_mean = sec_mean.parse::<f64>()?;
        let sec_var = sec_var.parse::<f64>()?;
        let sec_stdev = sec_stdev.parse::<f64>()?;
        let sec_n = sec_n.parse::<usize>()?;
        indel_read_stats.insert(
            chrom.to_owned(),
            ReadIndelSummaryStats {
                primary: IndelSummaryStats {
                    mean,
                    var,
                    stdev,
                    n,
                },
                secondary: IndelSummaryStats {
                    mean: sec_mean,
                    var: sec_var,
                    stdev: sec_stdev,
                    n: sec_n,
                },
            },
        );
    }

    Ok(indel_read_stats)
}

pub fn write_indel_read_stats(
    stats: &HashMap<String, ReadIndelSummaryStats>,
    outfile: &Path,
) -> eyre::Result<()> {
    let mut writer = BufWriter::new(File::create(outfile)?);
    writeln!(
        &mut writer,
        "#chrom\tmean\tvar\tstdev\tn\tsec_mean\tsec_var\tsec_stdev\tsec_n\t"
    )?;
    for (chrom, stats) in stats {
        writeln!(
            &mut writer,
            "{chrom}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            stats.primary.mean,
            stats.primary.var,
            stats.primary.stdev,
            stats.primary.n,
            stats.secondary.mean,
            stats.secondary.var,
            stats.secondary.stdev,
            stats.secondary.n,
        )?;
    }
    Ok(())
}

pub fn aligned_intervals_windows(
    fh: &mut indexed_reader::IndexedReader<bgzf::io::Reader<File>>,
    window: usize,
) -> eyre::Result<HashMap<String, Vec<Interval<usize, String>>>> {
    let header = fh.read_header()?;
    Ok(header
        .reference_sequences()
        .into_iter()
        .flat_map(move |(ctg, ref_seq)| {
            let length: usize = ref_seq.length().get();
            let (num, rem) = (length / window, length % window);
            let final_start = num * window;
            let final_itv = Interval {
                start: final_start,
                stop: final_start + rem,
                val: ctg.to_string(),
            };
            (1..num + 1)
                .map(move |i| {
                    // One-based half closed, half closed intervals
                    let start = ((i - 1) * window).clamp(1, usize::MAX);
                    let stop = (i * window).clamp(1, usize::MAX);
                    Interval {
                        start,
                        stop,
                        val: ctg.to_string(),
                    }
                })
                .chain([final_itv])
        })
        .fold(HashMap::new(), |mut acc, itv| {
            if let Some(itvs) = acc.get_mut(&itv.val) {
                itvs.push(itv);
            } else {
                acc.insert(itv.val.to_owned(), vec![itv]);
            }
            acc
        }))
}
