use std::{collections::HashMap, ops::Bound};

use eyre::bail;
use itertools::Itertools;
use noodles::{
    bam::{self, Record},
    core::{Position, Region},
    sam::alignment::record::{Flags, cigar::op::Kind},
};
use petgraph::{
    Graph,
    dot::{Config, Dot},
    prelude::NodeIndex,
};
use rust_lapper::{Interval, Lapper};

type ReadBoundItvs = (Interval<usize, NodeIndex>, Interval<usize, NodeIndex>);

/// Convert cigar string to operations.
/// * Adapted from <https://github.com/pysam-developers/pysam/blob/3e3c8b0b5ac066d692e5c720a85d293efc825200/pysam/libcalignedsegment.pyx#L2009>
pub(crate) fn get_aligned_pairs(
    cg: impl Iterator<Item = (Kind, usize)>,
    pos: usize,
) -> eyre::Result<Vec<(usize, usize, Kind)>> {
    let mut pos: usize = pos;
    let mut qpos: usize = 0;
    let mut pairs = vec![];
    // Matches only
    for (op, l) in cg {
        match op {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in pos..(pos + l) {
                    pairs.push((qpos, i, op));
                    qpos += 1
                }
                pos += l
            }
            // Track indels and softclips.
            Kind::Pad | Kind::Insertion | Kind::SoftClip => {
                qpos += l;
                continue;
            }
            Kind::Deletion => {
                for i in pos..(pos + l) {
                    pairs.push((qpos, i, op));
                }
                pos += l
            }
            Kind::HardClip => {
                continue;
            }
            Kind::Skip => pos += l,
        }
    }
    Ok(pairs)
}

#[derive(Debug, Default, PartialEq, Eq, Hash, Clone)]
struct PosNt(usize, char);

fn get_or_add_node(
    wt: PosNt,
    graph: &mut Graph<PosNt, usize>,
    node_wt_map: &mut HashMap<PosNt, (NodeIndex, bool)>,
    is_mismatch: bool,
) -> NodeIndex {
    if let Some((node_idx, _)) = node_wt_map.get(&wt) {
        *node_idx
    } else {
        let node_idx: petgraph::prelude::NodeIndex = graph.add_node(wt.clone());
        node_wt_map.insert(wt, (node_idx, is_mismatch));
        node_idx
    }
}

fn generate_mismatch_dag(aln: &str, region: Region) -> eyre::Result<()> {
    let mut indexed_reader = bam::io::indexed_reader::Builder::default().build_from_path(&aln)?;
    let header = indexed_reader.read_header()?;
    let (Bound::Included(st), Bound::Included(end)) = (
        region.start().map(|b| b.get()),
        region.end().map(|b| b.get()),
    ) else {
        bail!("Invalid st or end")
    };
    let length = end - st;
    let query = indexed_reader.query(&header, &region)?;

    // Mismatch DAG
    let mut graph: Graph<PosNt, usize> = Graph::new();
    let mut node_wt_map: HashMap<PosNt, (NodeIndex, bool)> = HashMap::new();
    let mut read_recs: HashMap<usize, Record> = HashMap::new();
    // Intervals of read edges in refpos
    let mut read_edge_bound_itvs: HashMap<usize, ReadBoundItvs> = HashMap::new();

    // Global metrics
    let mut mapq: Vec<usize> = vec![0; length + 1];
    let mut cov: Vec<usize> = vec![0; length + 1];

    for (i, rec) in query
        .records()
        .flatten()
        .filter(|aln| !aln.flags().contains(Flags::SECONDARY))
        .enumerate()
    {
        let cg: bam::record::Cigar<'_> = rec.cigar();
        let aln_pairs = get_aligned_pairs(
            cg.iter().flatten().map(|op| (op.kind(), op.len())),
            rec.alignment_start().unwrap()?.get(),
        )?;
        let seq = rec.sequence();
        let rmapq = rec.mapping_quality().unwrap_or_default().get() as usize;

        // Store nodes and add edges later
        let mut nodes = vec![];
        // Keep track of first and last ref position
        let (Some((_, min_ref_pos, _)), Some((_, max_ref_pos, _))) =
            (aln_pairs.first().cloned(), aln_pairs.last().cloned())
        else {
            continue;
        };
        for (qpos, refpos, kind) in aln_pairs
            .into_iter()
            .filter(|(_, refpos, _)| *refpos >= st && *refpos <= end)
        {
            let ipos = refpos - st;
            let nt = seq.get(qpos).map(char::from).unwrap();
            let cnt = match kind {
                Kind::Insertion | Kind::SoftClip => 0,
                Kind::SequenceMatch | Kind::Deletion => 1,
                Kind::SequenceMismatch => {
                    // Add node for refpos and the nt.
                    let wt = PosNt(refpos, nt);
                    let node_idx = get_or_add_node(wt, &mut graph, &mut node_wt_map, true);
                    nodes.push((refpos, node_idx));
                    1
                }
                _ => 1,
            };
            // Store coverage and MAPQ (total)
            cov[ipos] += cnt;
            mapq[ipos] += rmapq;
        }

        // Get the bounding intervals between the first and last mismatch.
        // Need to pick up CO events
        // | *   * |
        //  ^     ^
        let (
            Some((first_mismatch_refpos, first_mismatch_node_idx)),
            Some((last_mismatch_refpos, last_mismatch_node_idx)),
        ) = (nodes.first(), nodes.last())
        else {
            continue;
        };
        read_edge_bound_itvs.insert(
            i,
            (
                Interval {
                    start: min_ref_pos,
                    stop: *first_mismatch_refpos + 1,
                    val: *first_mismatch_node_idx,
                },
                Interval {
                    start: *last_mismatch_refpos,
                    stop: max_ref_pos + 1,
                    val: *last_mismatch_node_idx,
                },
            ),
        );

        // Add nodes and edges between.
        nodes.sort_by(|a, b| a.0.cmp(&b.0));
        graph.extend_with_edges(nodes.windows(2).map(|w| {
            let [a, b] = w else { panic!() };
            (a.1, b.1, i)
        }));

        read_recs.insert(i, rec);
    }

    // Construct intervaltree of mismatched regions
    let itree_mism = Lapper::new(
        node_wt_map
            .iter()
            .map(|(pos_nt, node_idx)| Interval {
                start: pos_nt.0,
                stop: pos_nt.0 + 1,
                val: *node_idx,
            })
            .collect(),
    );

    // TODO: Node 6339682 is incorrect.
    // Then fill in the gaps in reads at the end
    for (id, rec) in read_recs.iter() {
        let (itv_st, itv_end) = &read_edge_bound_itvs[id];
        let seq = rec.sequence();
        // For each overlap, build edge between existing mismatch nodes (^ = itv)
        // | |  | |
        //  ^
        for (mism_1, mism_2) in itree_mism
            .find(itv_st.start, itv_st.stop)
            .sorted()
            .tuple_windows()
        {
            let pos_1 = mism_1.start - st;
            let pos_2 = mism_2.start - st;
            let nt_1 = seq.get(pos_1).map(char::from).unwrap();
            let nt_2 = seq.get(pos_2).map(char::from).unwrap();
            let node_idx_mism_1 = get_or_add_node(
                PosNt(mism_1.start, nt_1),
                &mut graph,
                &mut node_wt_map,
                false,
            );
            let node_idx_mism_2 = get_or_add_node(
                PosNt(mism_2.start, nt_2),
                &mut graph,
                &mut node_wt_map,
                false,
            );

            graph.add_edge(node_idx_mism_1, node_idx_mism_2, *id);
        }
        // | |  | |
        //       ^
        for (mism_1, mism_2) in itree_mism
            .find(itv_end.start, itv_end.stop)
            .sorted()
            .tuple_windows()
        {
            let pos_1 = mism_1.start - st;
            let pos_2 = mism_2.start - st;
            let nt_1 = seq.get(pos_1).map(char::from).unwrap();
            let nt_2 = seq.get(pos_2).map(char::from).unwrap();
            let node_idx_mism_1 = get_or_add_node(
                PosNt(mism_1.start, nt_1),
                &mut graph,
                &mut node_wt_map,
                false,
            );
            let node_idx_mism_2 = get_or_add_node(
                PosNt(mism_2.start, nt_2),
                &mut graph,
                &mut node_wt_map,
                false,
            );
            graph.add_edge(node_idx_mism_1, node_idx_mism_2, *id);
        }
    }

    // Calculate MAPQ based on coverage
    for (ipos, mq) in mapq.iter_mut().enumerate() {
        let ncov = cov[ipos];
        *mq /= ncov
    }

    // TODO: Find paths

    // https://wintertee.github.io/Graphviz-Visualizer/
    let f_get_node_attr = |_g, (_node_idx, node_val)| {
        let (_, is_mismatch) = &node_wt_map[node_val];
        let color = if *is_mismatch { "red" } else { "white" };
        format!(
            r#"label = "{} ({})", style="filled", fillcolor="{}""#,
            node_val.1, node_val.0, color
        )
    };
    let dot_graph = Dot::with_attr_getters(
        &graph,
        // Omit builtin edge and node labels
        &[Config::EdgeNoLabel, Config::NodeNoLabel],
        // TODO: read as label
        // TODO: color best paths
        &|_g, edge_ref| format!(r#"label = "{}""#, edge_ref.weight()),
        // Color by mismatch
        &f_get_node_attr,
    );
    println!("{dot_graph:?}");

    Ok(())
}

fn main() -> eyre::Result<()> {
    let bam = "test/CT22_ENA_CBCUDK010000011_CBCUDK010000011.1_6335921-6341074.bam";
    // let fa = "test/CT22_ENA_CBCUDK010000011_CBCUDK010000011.1.fa.gz";

    // ENA_CBCUDK010000011_CBCUDK010000011.1:6334613-6342169
    let region = Region::new(
        "ENA_CBCUDK010000011_CBCUDK010000011.1",
        Position::new(6334613).unwrap()..=Position::new(6342169).unwrap(),
    );
    generate_mismatch_dag(bam, region)?;
    Ok(())
}
