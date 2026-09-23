use std::path::PathBuf;

use clap::Parser;

/// Detect putative intrachromosomal recombination events
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    /// DSA read alignment
    #[arg(
        short = 'i',
        long,
        default_value = "test/CT22_ENA_CBCUDK010000011_CBCUDK010000011.1_6335921-6341074.bam"
    )]
    pub bam: PathBuf,

    /// Genome assembly
    #[arg(
        short,
        long,
        default_value = "test/CT22_ENA_CBCUDK010000011_CBCUDK010000011.1.fa.gz"
    )]
    pub fa: PathBuf,

    /// BED file to restrict search
    #[arg(
        short,
        long,
        default_value = "test/single_read/CT22_ENA_CBCUDK010000011_CBCUDK010000011.1.bed"
    )]
    pub bed: Option<PathBuf>,

    /// Whole genome window size if no bed file provided.
    #[arg(short, long, default_value_t = 5_000_000)]
    pub wg_window: usize,
}
