import pysam
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, integrate

# https://github.com/vollgerlab/NucFreq/blob/88ac15a819b2fdbd81806d1a5123e1d8ba03f674/NucPlot.py#L70-L83
CIGAR_OPS = {
    0: "M",
    1: "I",
    2: "D",
    3: "N",
    4: "S",
    5: "H",
    6: "P",
    7: "E",
    8: "X",
    9: "B",
    10: "NM",
}
# M    BAM_CMATCH      0
# I    BAM_CINS        1
# D    BAM_CDEL        2
# N    BAM_CREF_SKIP   3
# S    BAM_CSOFT_CLIP  4
# H    BAM_CHARD_CLIP  5
# P    BAM_CPAD        6
# =    BAM_CEQUAL      7
# X    BAM_CDIFF       8
# B    BAM_CBACK       9
# NM   NM tag          10

def main():
    pass
    # https://stackoverflow.com/questions/49293019/calculating-probability-distribution-from-time-series-data-in-python
    fh = pysam.AlignmentFile("/home/koisland/Projects/RECallS/test/single_read/CT22_ENA_CBCUDK010000011_CBCUDK010000011.1_6335921-6341074.bam")

    events = {
        "no_event": "m84108_241123_231350_s2/62194242",
        "no_event_noisy": "m84108_241123_231350_s2/93456693",
        "co_event": "m84108_241123_231350_s2/130221242",
        "other_event": "m84108_241123_231350_s2/126487957",
    }
    qry = fh.fetch(
        "ENA_CBCUDK010000011_CBCUDK010000011.1",
        6317808,
        6359969,
    )
    window = 1000
    for read in qry:
        read_name = read.query_name
        if not read_name:
            continue

        if read_name not in events.values():
            continue
        # Take only aligned length
        read_aln_len = read.query_alignment_length
        mismatch_pos = []
        qscores = read.query_alignment_qualities
        if not qscores:
            continue
        for qpos, _, op in read.get_aligned_pairs(with_cigar=True):
            op = CIGAR_OPS[op]
            if op == "X":
                qscore = qscores[qpos] if qpos else 0
                if qscore > 30:
                    mismatch_pos.append(qpos)
        
        x = np.linspace(0, read_aln_len, window)
        midpt = read_aln_len / 2
        x_left = (0, midpt)
        x_right = (midpt, read_aln_len)
        x_integral_left = np.linspace(*x_left, window)
        x_integral_right = np.linspace(*x_right, window)

        if mismatch_pos:
            kde = stats.gaussian_kde(mismatch_pos)
            integral_left, err_left = integrate.quad_vec(kde, *x_left)
            integral_left = integral_left[0]
            integral_right, err_right = integrate.quad_vec(kde, *x_right)
            integral_right = integral_right[0]
        else:
            kde = lambda x: [0] * len(x)
            integral_left = 0.0
            integral_right = 0.0

        y = kde(x)
        plt.plot(x, y, label="KDE")
        abs_diff_area = abs(integral_left - integral_right)
        if abs_diff_area > 0.5:
            status = "Split"
        else:
            status = None

        plt.fill_between(
            x_integral_left,
            0,
            kde(x_integral_left),
            alpha=0.3,
            color='b',
            label="Area: {:.3f}".format(integral_left)
        )
        plt.fill_between(
            x_integral_right,
            0,
            kde(x_integral_right),
            alpha=0.3,
            color='r',
            label="Area: {:.3f}".format(integral_right)
        )
        plt.legend()
        plt.savefig(f"{read_name.replace("/", "_")}_{status}.png")
        plt.clf()


    fh.close()

if __name__ == "__main__":
    raise SystemExit(main())
