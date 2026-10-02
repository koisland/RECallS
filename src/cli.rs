use std::path::PathBuf;

use clap::Parser;

/// Detect putative intrachromosomal recombination events from long-read sequencing of sperm data
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    /// DSA read alignment
    #[arg(short = 'i', long)]
    pub bam: PathBuf,

    /// Genome assembly
    #[arg(short, long)]
    pub fa: PathBuf,

    /// BED file to restrict search
    #[arg(short, long)]
    pub bed: Option<PathBuf>,

    #[arg(short, long, default_value = "./recalls")]
    pub output_dir: PathBuf,

    /// BED file to ignore
    #[arg[short = 'n', long]]
    pub ignore_bed: Option<PathBuf>,

    /// Number of threads
    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,

    /// Whole genome window size if no bed file provided.
    #[arg(short, long, default_value_t = 5_000_000)]
    pub wg_window: usize,

    /// Inversion indel z-score
    #[arg(long, default_value_t = 3.4)]
    pub inv_indel_zscore: f32,

    /// Inversion indel minimum aligned read length
    #[arg(long, default_value_t = 10_000)]
    pub inv_min_aln_len: usize,

    /// Inversion indel threshold unbalanced
    #[arg(long, default_value_t = 0.5)]
    pub inv_thr_unbalanced: f32,

    /// Inversion minimum number of SNVs
    #[arg(long, default_value_t = 10)]
    pub inv_min_num_snvs: usize,

    /// Deletion maximum region divergence between supplementary aligned regions.
    #[arg(long, default_value_t = 0.05)]
    pub del_max_rgn_dv: f32,

    /// Deletion minimum MAPQ of supplementary alignment.
    #[arg(long, default_value_t = 60)]
    pub del_min_mapq: u8,

    /// Deletion minimum number of overlaps with other supplementary alignment events required to call deletion event.
    #[arg(long, default_value_t = 2)]
    pub del_min_ovl_cnt: usize,
}
