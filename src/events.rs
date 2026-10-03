#[derive(Debug, PartialEq, Eq, Clone)]
pub struct MismatchSignal {
    pub chrom: String,
    pub start: usize,
    pub stop: usize,
    pub rname: String,
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
            self.chrom,
            self.start,
            self.stop,
            self.rname,
            self.n_indels,
            self.aln_len,
            self.is_secondary,
            self.is_unbalanced
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SupplJuncSignal {
    pub chrom: String,
    pub start_1: usize,
    pub stop_1: usize,
    pub rname_1: String,
    pub start_2: usize,
    pub stop_2: usize,
    pub rname_2: String,
}
impl SupplJuncSignal {
    pub fn header() -> &'static str {
        "#chrom\tstart\tstop\tread_name\tstart_2\tstop_2\tread_name_2"
    }

    pub fn as_bed(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.chrom,
            self.start_1,
            self.stop_1,
            self.rname_1,
            self.start_2,
            self.stop_2,
            self.rname_2
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SupplSignal {
    pub chrom: String,
    pub start: usize,
    pub stop: usize,
    pub suppl_start: usize,
    pub suppl_stop: usize,
    pub rname: String,
}
impl SupplSignal {
    pub fn header() -> &'static str {
        "#chrom\tstart\tstop\tread_name\tsuppl_start\tsuppl_stop"
    }

    pub fn as_bed(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.chrom, self.start, self.stop, self.rname, self.suppl_start, self.suppl_stop
        )
    }
}

pub enum Event {
    Deletion(SupplSignal),
    Inversion(SupplJuncSignal),
    InversionInferred(MismatchSignal),
}
