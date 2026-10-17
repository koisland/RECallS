use rust_lapper::Interval;

use crate::{call::SupplIntervalSignals, utils::Direction};

/// Intermediate signal
#[derive(Debug)]
pub enum Signal {
    /// Supplementary alignment signal
    SupplSignal(SupplIntervalSignals),
    /// Mismatch alignment signal
    MismatchSignal(MismatchSignal),
}

/// Mismatch signal
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct MismatchSignal {
    /// Interval of region
    pub itv: Interval<usize, String>,
    /// Read name
    pub read: String,
    /// Length of indels
    pub n_indels: usize,
    /// Total read alignment length
    pub aln_len: usize,
    /// If is secondary alignment
    pub is_secondary: bool,
    /// Direction of unbalanced read
    pub unbalanced_dir: Direction,
    /// Span of unbalanced read
    pub unbalanced_span: Interval<usize, ()>,
}

/// SV type detected
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
    /// Interval span of SV
    pub itv: Interval<usize, String>,
    /// Reads associated with SV
    pub reads: Vec<String>,
    /// What type of SV?
    pub typ: SVType,
}

impl SV {
    pub fn as_bed(&self) -> String {
        let reads = self.reads.join("~");
        format!(
            "{}\t{}\t{}\t{}\t0\t.\t{}\t{}\t{}",
            self.itv.val,
            self.itv.start,
            self.itv.stop,
            reads,
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
