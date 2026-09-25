use scirs2_integrate::quad::{QuadOptions, quad};
use scirs2_stats::kde::KernelDensityEstimate;

pub fn is_unbalanced_alignment(
    mismatch_pos: &[f64],
    read_len: f64,
    min_n_events: usize,
    verbose: bool,
) -> eyre::Result<bool> {
    let quad_opts = QuadOptions {
        max_evals: 10000,
        ..Default::default()
    };
    if mismatch_pos.is_empty() || mismatch_pos.len() < min_n_events {
        return Ok(false);
    }
    // Get midpoint of read
    let midpt = read_len / 2.0;
    let (left_lower, left_upper) = (0f64, midpt);
    let (right_lower, right_upper) = (midpt, read_len);

    // TODO: Also check that mismatches covers equal amount.
    if !mismatch_pos.is_empty() {
        // Calculate PDF of read mismatch positions
        let kde = KernelDensityEstimate::new(mismatch_pos, scirs2_stats::Kernel::Gaussian);
        // Find area under curve of left and right side
        let res_left = quad(
            |x| kde.evaluate(x),
            left_lower,
            left_upper,
            Some(quad_opts.clone()),
        )?;
        let res_right = quad(
            |x| kde.evaluate(x),
            right_lower,
            right_upper,
            Some(quad_opts),
        )?;
        let integral_left = res_left.value;
        let integral_right = res_right.value;
        let abs_diff_area = (integral_left - integral_right).abs();
        if verbose {
            eprintln!("{integral_left},{integral_right},{mismatch_pos:?}")
        }
        Ok(abs_diff_area > 0.5)
    } else {
        return Ok(false);
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
    fn test_check_read() -> eyre::Result<()> {
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
                    is_unbalanced_alignment(&mismatch_pos, read_len, 5, false).unwrap(),
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
