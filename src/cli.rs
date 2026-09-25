use std::path::PathBuf;

use clap::Parser;

/// Detect putative intrachromosomal recombination events from long-read sequencing of sperm data
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    /// DSA read alignment
    #[arg(short = 'i', long, default_value = "test/amplicon/chm13_chr7_sim.bam")]
    pub bam: PathBuf,

    /// Genome assembly
    #[arg(short, long, default_value = "test/amplicon/chm13_chr7.fa.gz")]
    pub fa: PathBuf,

    /// BED file to restrict search
    #[arg(short, long)]
    pub bed: Option<PathBuf>,

    /// BED file to ignore
    #[arg[short, long]]
    pub ignore_bed: Option<PathBuf>,

    /// Whole genome window size if no bed file provided.
    #[arg(short, long, default_value_t = 5_000_000)]
    pub wg_window: usize,
}
