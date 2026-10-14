use itertools::Itertools;
use rust_lapper::Interval;

use crate::{call::SupplIntervalSignals, utils::Direction};

#[derive(Debug)]
pub enum Signal {
    SupplSignal(SupplIntervalSignals),
    MismatchSignal(MismatchSignal),
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct MismatchSignal {
    pub itv: Interval<usize, String>,
    pub read: String,
    pub n_indels: usize,
    pub aln_len: usize,
    pub is_secondary: bool,
    pub unbalanced_dir: Direction,
    pub unbalanced_span: Interval<usize, ()>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SVType {
    Deletion,
    Inversion,
}

impl SVType {
    fn color(&self) -> &'static str {
        match self {
            SVType::Deletion => "255,0,0",
            SVType::Inversion => "0,0,255",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SV {
    pub itv: Interval<usize, String>,
    pub reads: [Option<String>; 2],
    pub typ: SVType,
}

impl SV {
    pub fn reads(&self) -> impl Iterator<Item = &str> {
        self.reads.iter().flatten().map(|read| &read[..])
    }

    pub fn as_bed(&self) -> String {
        let reads = if let Some((r1, r2)) = self.reads().collect_tuple() {
            (r1, r2)
        } else {
            ("", "")
        };
        format!(
            "{}\t{}\t{}\t{}~{}\t0\t.\t{}\t{}\t{}",
            self.itv.val,
            self.itv.start,
            self.itv.stop,
            reads.0,
            reads.1,
            self.itv.start,
            self.itv.stop,
            self.typ.color()
        )
    }
}

#[derive(Debug)]
pub enum Event {
    SV(SV),
    Signal(Signal),
}
