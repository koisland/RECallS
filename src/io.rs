use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use itertools::Itertools;
use noodles::{
    bam::io::indexed_reader,
    bgzf,
    core::{Position, Region},
};

pub fn read_bed(bed: &Path) -> Option<Vec<Region>> {
    let mut intervals = Vec::new();
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

        intervals.push(Region::new(
            name,
            Position::new(start)?..=Position::new(stop)?,
        ))
    }
    Some(intervals)
}

pub fn aligned_intervals_windows(
    fh: &mut indexed_reader::IndexedReader<bgzf::io::Reader<File>>,
    window: usize,
) -> eyre::Result<Vec<Region>> {
    let header = fh.read_header()?;
    Ok(header
        .reference_sequences()
        .into_iter()
        .flat_map(move |(ctg, ref_seq)| {
            let length: usize = ref_seq.length().get();
            let (num, rem) = (length / window, length % window);
            let final_start = num * window;
            let final_itv = Region::new(
                ctg.to_owned(),
                Position::new(final_start).unwrap()..=Position::new(final_start + rem).unwrap(),
            );
            (1..num + 1)
                .map(move |i| {
                    // One-based half closed, half closed intervals
                    let start = ((i - 1) * window).clamp(1, usize::MAX);
                    let stop = (i * window).clamp(1, usize::MAX);
                    Region::new(
                        ctg.to_owned(),
                        Position::new(start).unwrap()..=Position::new(stop).unwrap(),
                    )
                })
                .chain([final_itv])
        })
        .collect())
}
