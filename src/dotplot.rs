use eyre::bail;
use rammap::Strand;
use rammap::align::map::AlignFlags;
use rammap::api::{Aligner, Preset, apply_preset_str};

fn run_whole_contig_dotplot(path: &str, chrom: &str) -> eyre::Result<()> {
    let mut aligner = Aligner::from_fasta(path, Preset::MapOnt)?;
    // -PD -k19 -w19 -m200
    // https://lh3.github.io/minimap2/minimap2.html#10
    // -P - Retain all chains and don’t attempt to set primary chains.
    // -D - If query sequence name/length are identical to the target name/length, ignore diagonal anchors.
    //      This option also reduces DP-based extension along the diagonal.
    // -k - Minimizer k-mer length [15]
    // -w - Minimizer window size [10].
    //      A minimizer is the smallest k-mer in a window of w consecutive k-mers.
    // -m - Discard chains with chaining score <INT [40].
    //      Chaining score equals the approximate number of matching bases minus a concave gap penalty.
    let opt = aligner.options_mut();
    apply_preset_str(opt, &mut 19usize, &mut 19usize, &mut true, "map-ont")
        .map_err(eyre::Report::msg)?;
    // min_chain_score
    opt.chaining.min_chain_score = 200;
    // all_chains
    opt.flags.insert(AlignFlags::ALL_CHAINS);
    // no_diag
    opt.flags.insert(AlignFlags::NO_DIAG);

    let results = aligner.map_seq("read1", b"ACGTACGTACGT...");
    for m in &results.mappings {
        println!(
            "{}\t{}\t{}\t{}\tMAPQ={}",
            m.target_name,
            m.target_start,
            m.target_end,
            if m.strand == Strand::Forward {
                "+"
            } else {
                "-"
            },
            m.mapq
        );
    }

    Ok(())
}
