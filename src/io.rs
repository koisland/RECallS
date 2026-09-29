use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
    str::FromStr,
};

use eyre::Context;
use itertools::Itertools;
use noodles::{
    bam::io::indexed_reader,
    bgzf::{self, io::IndexedReader},
    core::{Position, Region},
    fasta,
};
use ordered_float::OrderedFloat;
use rust_lapper::{Interval, Lapper};

use crate::{
    baseline::{IndelSummaryStats, ReadIndelSummaryStats},
    self_align::{Paf, Strand},
};

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

pub fn read_paf(paf: &Path, max_dv: f32) -> eyre::Result<HashMap<String, Lapper<usize, Paf>>> {
    let fh_paf_self_align = BufReader::new(File::open(paf)?);
    let mut pafs: HashMap<String, Vec<Interval<usize, Paf>>> = HashMap::new();
    for line in fh_paf_self_align.lines() {
        let line = line?;
        if let Some(
            [
                qchrom,
                qlen,
                qst,
                qend,
                strand,
                tchrom,
                tlen,
                tst,
                tend,
                matches,
                aln_len,
                dv,
            ],
        ) = line.split('\t').collect_array()
        {
            let qlen: usize = qlen.parse()?;
            let qst: usize = qst.parse()?;
            let qend: usize = qend.parse()?;
            let tlen: usize = tlen.parse()?;
            let tst: usize = tst.parse()?;
            let tend: usize = tend.parse()?;
            let matches: usize = matches.parse()?;
            let aln_len: usize = aln_len.parse()?;
            let dv = OrderedFloat(dv.parse()?);
            // Ignore highly divergent alignments
            if *dv > max_dv {
                continue;
            }
            let itv = Interval {
                start: qst,
                stop: qend,
                val: Paf {
                    qchrom: qchrom.to_owned(),
                    qlen,
                    qst,
                    qend,
                    strand: Strand::from_str(strand)?,
                    tchrom: tchrom.to_owned(),
                    tlen,
                    tst,
                    tend,
                    matches,
                    aln_len,
                    dv,
                },
            };
            if let Some(itvs) = pafs.get_mut(qchrom) {
                itvs.push(itv);
            } else {
                pafs.insert(qchrom.to_owned(), vec![itv]);
            }
        } else {
            eprintln!("Invalid PAF row: {line}");
            continue;
        }
    }
    Ok(pafs
        .into_iter()
        .map(|(chrom, itvs)| (chrom, Lapper::new(itvs)))
        .collect())
}

pub fn write_itvs_self_similar_paf(
    itvs: &HashMap<String, Lapper<usize, Paf>>,
    outfile: &Path,
) -> eyre::Result<()> {
    let mut writer = BufWriter::new(File::create(outfile)?);
    for itvs in itvs.values() {
        for itv in itvs.iter() {
            writeln!(
                &mut writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                itv.val.qchrom,
                itv.val.qlen,
                itv.val.qst,
                itv.val.qend,
                char::from(itv.val.strand),
                itv.val.tchrom,
                itv.val.tlen,
                itv.val.tst,
                itv.val.tend,
                itv.val.matches,
                itv.val.aln_len,
                itv.val.dv,
            )?;
        }
    }
    Ok(())
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

pub enum FastaReader {
    Bgzip(fasta::io::Reader<IndexedReader<File>>),
    Standard(fasta::io::Reader<BufReader<File>>),
}

pub struct FastaHandle {
    pub reader: FastaReader,
    pub fai: fasta::fai::Index,
}

impl FastaHandle {
    /// Create new handle.
    pub fn new(infile: impl AsRef<Path>) -> eyre::Result<Self> {
        let (fai, gzi) = Self::get_faidx(&infile)?;
        let fh = Self::read_fa(&infile, gzi.as_ref())?;
        Ok(Self { reader: fh, fai })
    }

    fn get_faidx(
        fa: &impl AsRef<Path>,
    ) -> eyre::Result<(fasta::fai::Index, Option<bgzf::gzi::Index>)> {
        // https://www.ginkgobioworks.com/2023/03/17/even-more-rapid-retrieval-from-very-large-files-with-rust/
        let fa_path = fa.as_ref().canonicalize()?;
        let is_bgzipped = fa_path.extension().and_then(|e| e.to_str()) == Some("gz");
        let mut fai_fname = fa_path.clone();
        fai_fname.as_mut_os_string().push(".fai");

        let fai = fasta::fai::fs::read(fai_fname);
        if is_bgzipped {
            let index_reader = bgzf::io::indexed_reader::Builder::default()
                .build_from_path(fa)
                .with_context(|| format!("Failed to read gzi for {fa_path:?}"))?;
            let gzi = index_reader.index().clone();

            if let Ok(fai) = fai {
                return Ok((fai, Some(gzi)));
            }
            log::debug!("No existing faidx for {fa_path:?}. Generating...");
            let mut records = Vec::new();
            let mut indexer = fasta::io::Indexer::new(index_reader);
            while let Some(record) = indexer.index_record()? {
                records.push(record);
            }

            Ok((fasta::fai::Index::from(records), Some(gzi)))
        } else {
            if let Ok(fai) = fai {
                return Ok((fai, None));
            }
            log::debug!("No existing faidx for {fa_path:?}. Generating...");
            Ok((fasta::fs::index(fa)?, None))
        }
    }

    /// Fetch coordinates. noodles use 1-based coordinates.
    pub fn fetch(
        &mut self,
        ctg_name: &str,
        start: usize,
        stop: usize,
    ) -> eyre::Result<fasta::Record> {
        let start_pos = Position::new(start.clamp(1, usize::MAX)).unwrap();
        let stop_pos = Position::new(stop.clamp(1, usize::MAX)).unwrap();
        let region = Region::new(ctg_name, start_pos..=stop_pos);
        match &mut self.reader {
            FastaReader::Bgzip(reader) => Ok(reader.query(&self.fai, &region)?),
            FastaReader::Standard(reader) => Ok(reader.query(&self.fai, &region)?),
        }
    }

    fn read_fa(
        fa: &impl AsRef<Path>,
        fa_gzi: Option<&bgzf::gzi::Index>,
    ) -> eyre::Result<FastaReader> {
        let fa_file = std::fs::File::open(fa);
        if let Some(fa_gzi) = fa_gzi {
            Ok(FastaReader::Bgzip(
                fa_file
                    .map(|file| bgzf::io::IndexedReader::new(file, fa_gzi.clone()))
                    .map(fasta::io::Reader::new)?,
            ))
        } else {
            Ok(FastaReader::Standard(
                fa_file
                    .map(std::io::BufReader::new)
                    .map(fasta::io::Reader::new)?,
            ))
        }
    }
}
