use std::{collections::HashMap, path::Path, process::Command, str::FromStr};

use eyre::bail;
use itertools::Itertools;
use noodles::fasta::{self, Record, record::Definition};
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use rust_lapper::{Interval, Lapper};

use crate::io::FastaHandle;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum Strand {
    Forward,
    Reverse,
}

impl FromStr for Strand {
    type Err = eyre::Report;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Forward),
            "-" => Ok(Strand::Reverse),
            _ => bail!("Invalid strand."),
        }
    }
}

impl From<Strand> for char {
    fn from(value: Strand) -> Self {
        match value {
            Strand::Forward => '+',
            Strand::Reverse => '-',
        }
    }
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Paf {
    pub qchrom: String,
    pub qlen: usize,
    pub qst: usize,
    pub qend: usize,
    pub strand: Strand,
    pub tchrom: String,
    pub tlen: usize,
    pub tst: usize,
    pub tend: usize,
    pub matches: usize,
    pub aln_len: usize,
    pub dv: OrderedFloat<f32>,
}

pub fn generate_contig_self_alignment(
    fa: &Path,
    seq_lens: &HashMap<String, usize>,
    max_dv: f32,
) -> eyre::Result<HashMap<String, Lapper<usize, Paf>>> {
    Ok(seq_lens
        .par_iter()
        .flat_map(|(name, ctg_len)| {
            let mut fh = FastaHandle::new(fa).expect("Cannot open fasta file");
            let rec = fh
                .fetch(name, 0, *ctg_len)
                .expect("Failed to query sequence for minimap2");

            // Remove coords in chrom name.
            let name =
                String::from_utf8(rec.name().to_vec()).expect("Invalid utf8 for sequence name.");
            let Some((_, name)) = name.rsplitn(2, ':').collect_tuple() else {
                panic!("Invalid chrom name: {name}")
            };
            let definition = Definition::new(name, None);
            let rec = Record::new(definition, rec.sequence().clone());

            // Create named tempfile and write single sequence
            let tempfile = tempfile::NamedTempFile::new()
                .expect("Unable to make tempfile for chrom fasta and minimap2");
            let tempfile_path = tempfile
                .as_ref()
                .to_str()
                .map(|s| s.to_owned())
                .expect("Unable to get tempfile name for chrom fasta and minimap2");

            let mut writer = fasta::io::Writer::new(tempfile);
            writer
                .write_record(&rec)
                .expect("Failed to write fasta record");

            let out_mm2 = Command::new("minimap2")
                .args([
                    "-PD",
                    "-k19",
                    "-w19",
                    "-m200",
                    &tempfile_path,
                    &tempfile_path,
                ])
                .output()
                .expect("Failed to spawn minimap2");

            if out_mm2.status.success() {
                let paf =
                    str::from_utf8(&out_mm2.stdout).expect("Invalid utf-8 in minimap2 output");
                let mut paf_itvs = vec![];
                for row in paf.split('\n').filter(|r| !r.is_empty()) {
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
                            _mapq,
                            _tp,
                            _cm,
                            _s1,
                            dv,
                            _rl,
                        ],
                    ) = row.split('\t').collect_array()
                    {
                        let qlen: usize = qlen.parse().unwrap();
                        let qst: usize = qst.parse().unwrap();
                        let qend: usize = qend.parse().unwrap();
                        let tlen: usize = tlen.parse().unwrap();
                        let tst: usize = tst.parse().unwrap();
                        let tend: usize = tend.parse().unwrap();
                        let matches: usize = matches.parse().unwrap();
                        let aln_len: usize = aln_len.parse().unwrap();
                        let Some((_, _, dv)) = dv.splitn(3, ':').collect_tuple() else {
                            eprintln!("Invalid dv {dv} tag for {row}");
                            continue;
                        };
                        let dv = OrderedFloat(dv.parse::<f32>().expect("Invalid dv float"));
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
                                strand: Strand::from_str(strand).unwrap(),
                                tchrom: tchrom.to_owned(),
                                tlen,
                                tst,
                                tend,
                                matches,
                                aln_len,
                                dv,
                            },
                        };
                        paf_itvs.push(itv);
                    } else {
                        eprintln!("Invalid minimap2 paf row: {row}")
                    }
                }
                Some((name.to_owned(), Lapper::new(paf_itvs)))
            } else {
                let err = str::from_utf8(&out_mm2.stderr)
                    .expect("Invalid utf-8 in minimap2 stderr output");
                eprintln!("Minimap2 failed: {err}");
                None
            }
        })
        .collect())
}
