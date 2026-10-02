import sys
from matplotlib.axes import Axes
import polars as pl
import matplotlib.pyplot as plt

def main():
    df = pl.read_csv(
        sys.argv[1],
        separator="\t",
        new_columns=["chrom", "typ", "perc"]
    )
    ax: Axes = plt.subplot()

    for color, (grp, df_grp) in zip(["red", "blue"], df.group_by(["typ"])):
        ax.hist(
            df_grp["perc"],
            color=color,
            label=grp[0],
            alpha=0.3,
            bins=100
        )
        mean = df_grp["perc"].mean()
        stdev = df_grp["perc"].std()
        mean_sd = mean + (stdev * 3.4)
        print(mean, stdev, mean_sd)
        ax.axvline(mean, color=color, linestyle="dotted")
        ax.axvline(mean_sd, color=color, linestyle="dotted")
        

    ax.set_yscale("log")
    ax.legend()
    plt.savefig("out.png")

if __name__ == "__main__":
    raise SystemExit(main())
