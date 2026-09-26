use integrate::prelude::*;
use kernel_density_estimation::prelude::*;

#[derive(Debug)]
pub struct UnbalancedSummary {
    pub is_unbalanced: bool,
    pub integral_left: f64,
    pub integral_right: f64,
}

/// Checks if the alignment is "unbalanced" where mismatches are localized to one end of the read.
/// * Generates a 1D KDE of mismatched positions
/// * Finds area under PDF on both ends
/// * Subtract each side to get absolute difference
/// * If greater than 50%, is considered unbalanced
///
/// # Arguments
/// * mismatch_pos: all mismatch positions on read
/// * read_len: read length
/// * min_n_events: minimum number of events to consider alignment
///
/// # Returns
/// * If is an unbalanced alignment
pub fn is_unbalanced_alignment(
    marker_qpos: &[f64],
    read_len: f64,
    min_n_events: usize,
) -> eyre::Result<Option<UnbalancedSummary>> {
    if marker_qpos.is_empty() || marker_qpos.len() < min_n_events {
        return Ok(None);
    }
    // Get midpoint of read
    let midpt = read_len / 2.0;
    let (left_lower, left_upper) = (0f64, midpt);
    let (right_lower, right_upper) = (midpt, read_len);

    if !marker_qpos.is_empty() {
        // Calculate PDF of read mismatch positions
        // https://aakinshin.net/posts/kde-bw/
        let kde = KernelDensityEstimator::new(
            marker_qpos,
            |data: &[f64]| Silverman.bandwidth(data),
            Normal,
        );
        // Find area under curve of left and right side
        let (integral_left, _) =
            gauss_kronrod_rule(|x| kde.pdf(&[x])[0], left_lower, left_upper, 7)
                .map_err(|err| eyre::Report::msg(err))?;
        let (integral_right, _) =
            gauss_kronrod_rule(|x| kde.pdf(&[x])[0], right_lower, right_upper, 7)
                .map_err(|err| eyre::Report::msg(err))?;

        // if integral_left.is_nan() || integral_right.is_nan() {
        //     eprintln!("{marker_qpos:?}")
        // }

        let abs_diff_area = (integral_left - integral_right).abs();
        Ok(Some(UnbalancedSummary {
            is_unbalanced: abs_diff_area > 0.4,
            integral_left,
            integral_right,
        }))
    } else {
        return Ok(None);
    }
}

#[cfg(test)]
mod test {
    use std::{
        collections::{HashMap, HashSet},
        fs::File,
    };

    use itertools::Itertools;
    use noodles::{
        bam::{self, io::IndexedReader},
        bgzf,
        core::{Position, Region},
        sam::alignment::{
            Record,
            record::{Flags, cigar::op::Kind},
        },
    };

    use crate::{
        unbalanced_aln::is_unbalanced_alignment,
        utils::{get_aligned_pairs, get_coords_from_region},
    };

    fn get_mismatch_pos(
        fh: &mut IndexedReader<bgzf::io::Reader<File>>,
        region: &Region,
        reads: HashSet<&str>,
    ) -> eyre::Result<HashMap<String, (Vec<f64>, f64)>> {
        let header = fh.read_header()?;
        let (st, end) = get_coords_from_region(&region)?;
        let query = fh.query(&header, &region)?;
        let mut mismatches_per_read = HashMap::new();

        let mut mismatches = vec![];

        for rec in query
            .records()
            .flatten()
            .filter(|aln| !aln.flags().contains(Flags::SECONDARY))
        {
            let name = str::from_utf8(rec.name().as_ref().unwrap())?;
            if !reads.contains(name) {
                continue;
            }
            let cg = rec.cigar();
            let aln_pairs = get_aligned_pairs(
                cg.iter().flatten().map(|op| (op.kind(), op.len())),
                rec.alignment_start().unwrap()?.get(),
            )?;
            let qscores = rec.quality_scores().as_bytes();
            let aln_len = rec.alignment_span().unwrap()? as f64;
            for (qpos, _refpos, kind) in aln_pairs
                .into_iter()
                .filter(|(_, refpos, _)| *refpos >= st && *refpos <= end)
            {
                if let Kind::SequenceMismatch = kind {
                    let qscore = qscores[qpos];
                    if qscore > 30 {
                        mismatches.push(qpos as f64);
                    }
                }
            }
            mismatches_per_read.insert(name.to_owned(), (mismatches.clone(), aln_len));
            mismatches.clear();
        }
        Ok(mismatches_per_read)
    }

    #[test]
    fn test_check_unbalanced_read() -> eyre::Result<()> {
        let reads: HashSet<&str> = HashSet::from_iter([
            "m84108_241123_231350_s2/62194242",
            "m84108_241123_231350_s2/93456693",
            "m84108_241123_231350_s2/130221242",
            "m84108_241123_231350_s2/126487957",
        ]);
        let bam = "test/single_read/CT22_ENA_CBCUDK010000011_CBCUDK010000011.1_6335921-6341074.bam";
        let mut fh = bam::io::indexed_reader::Builder::default().build_from_path(bam)?;
        let region = Region::new(
            "ENA_CBCUDK010000011_CBCUDK010000011.1",
            Position::new(6317808).unwrap()..=Position::new(6359969).unwrap(),
        );
        let res: Vec<(String, bool)> = get_mismatch_pos(&mut fh, &region, reads)?
            .into_iter()
            .map(|(rname, (mismatch_pos, read_len))| {
                (
                    rname,
                    is_unbalanced_alignment(&mismatch_pos, read_len, 5)
                        .unwrap()
                        .unwrap()
                        .is_unbalanced,
                )
            })
            .sorted_by(|a, b| a.0.cmp(&b.0))
            .collect();

        assert_eq!(
            res,
            vec![
                ("m84108_241123_231350_s2/126487957".to_owned(), false),
                ("m84108_241123_231350_s2/130221242".to_owned(), true),
                ("m84108_241123_231350_s2/62194242".to_owned(), false),
                ("m84108_241123_231350_s2/93456693".to_owned(), false)
            ]
        );
        Ok(())
    }
}
