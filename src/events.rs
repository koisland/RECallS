use rust_lapper::Interval;

use crate::call::SupplIntervalSignals;

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
    pub is_unbalanced: bool,
}
impl MismatchSignal {
    pub fn header() -> &'static str {
        "#chrom\tstart\tstop\tread_name\tn_indels\taln_len\tis_secondary\tis_unbalanced"
    }

    pub fn as_bed(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.itv.val,
            self.itv.start,
            self.itv.stop,
            self.read,
            self.n_indels,
            self.aln_len,
            self.is_secondary,
            self.is_unbalanced
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SVType {
    Deletion,
    Inversion,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SV {
    pub itv: Interval<usize, String>,
    pub read: String,
    pub typ: SVType,
}

#[derive(Debug)]
pub enum Event {
    SV(SV),
    Signal(Signal),
}
