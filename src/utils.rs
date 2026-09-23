use std::ops::Bound;

use eyre::bail;

use noodles::{core::Region, sam::alignment::record::cigar::op::Kind};

/// Convert cigar string to operations.
/// * Adapted from <https://github.com/pysam-developers/pysam/blob/3e3c8b0b5ac066d692e5c720a85d293efc825200/pysam/libcalignedsegment.pyx#L2009>
pub fn get_aligned_pairs(
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

pub fn get_coords_from_region(region: &Region) -> eyre::Result<(usize, usize)> {
    let (Bound::Included(st), Bound::Included(end)) = (
        region.start().map(|b| b.get()),
        region.end().map(|b| b.get()),
    ) else {
        bail!("Invalid st or end")
    };
    Ok((st, end))
}
