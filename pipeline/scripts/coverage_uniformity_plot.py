#!/usr/bin/env python3
"""
coverage_uniformity_plot.py
Bar chart comparing mean sequencing depth across the four EGFR amplicons.
Highlights the depth imbalance between amplicons and annotates each bar
with mean depth and coefficient of variation (CV).

Usage: coverage_uniformity_plot.py <annotated_depth.txt> <out.pdf>
"""
import sys
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COLORS = {
    "Amplicon_1": "#2166ac",
    "Amplicon_2": "#4dac26",
    "Amplicon_3": "#d01c8b",
    "Amplicon_4": "#f1a340",
}
ORDER = ["Amplicon_2", "Amplicon_1", "Amplicon_4", "Amplicon_3"]
EXON  = {"Amplicon_2": 18, "Amplicon_1": 19, "Amplicon_4": 20, "Amplicon_3": 21}


def parse_annotated_depth(path):
    """Parse an annotated depth file and group depth values by amplicon.

    Args:
        path: Path to the annotated depth file produced by annotate_amplicon_depth.py
              (columns: chrom, pos, depth, amplicon_name).

    Returns:
        Dict mapping amplicon_name (str) to a list of depth (int) values.
    """
    data = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            data[parts[3]].append(int(parts[2]))
    return data


def main(depth_path, out_path):
    """Generate a bar chart of mean sequencing depth per EGFR amplicon.

    Args:
        depth_path: Path to the annotated depth file (chrom, pos, depth, amplicon).
        out_path:   Destination path for the output PDF figure.
    """
    data = parse_annotated_depth(depth_path)

    means = [np.mean(data[a]) for a in ORDER]
    stds  = [np.std(data[a])  for a in ORDER]
    cvs   = [100 * s / m for s, m in zip(stds, means)]

    labels = [f"{a.replace('_', ' ')}\n(Ex{EXON[a]})" for a in ORDER]
    colors = [COLORS[a] for a in ORDER]
    x      = np.arange(len(ORDER))

    fig, ax = plt.subplots(figsize=(8, 5))

    bars = ax.bar(x, means, color=colors, alpha=0.75, edgecolor=[c for c in colors],
                  linewidth=1.2, width=0.55)

    # Error bars (±1 sd)
    ax.errorbar(x, means, yerr=stds, fmt="none", color="#444444",
                capsize=5, linewidth=1.2, capthick=1.2)

    # Annotate each bar: mean depth and CV
    for xi, (mean, cv) in enumerate(zip(means, cvs)):
        ax.text(xi, mean + max(means) * 0.03,
                f"{mean/1e3:.1f}k×\nCV={cv:.1f}%",
                ha="center", va="bottom", fontsize=8.5, color="#222222")

    # Overall CV across amplicons
    overall_cv = 100 * np.std(means) / np.mean(means)
    ax.text(0.98, 0.97,
            f"Inter-amplicon CV = {overall_cv:.1f}\\%",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=9, color="#333333",
            bbox=dict(facecolor="white", edgecolor="#cccccc", pad=4))

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Mean depth (×)", fontsize=10)
    ax.set_title("Coverage uniformity across EGFR amplicons", fontsize=11)
    ax.set_ylim(0, max(means) * 1.22)
    ax.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(
            lambda y, _: f"{y/1e3:.0f}k" if y >= 1000 else str(int(y))
        )
    )
    ax.spines[["top", "right"]].set_visible(False)

    # Print summary to stderr
    print("Coverage uniformity summary:", file=sys.stderr)
    print(f"{'Amplicon':<12} {'Exon':>5} {'Mean':>10} {'SD':>10} {'CV%':>7}",
          file=sys.stderr)
    for amp, mean, std, cv in zip(ORDER, means, stds, cvs):
        print(f"{amp:<12} {EXON[amp]:>5} {mean:>10.0f} {std:>10.0f} {cv:>7.1f}",
              file=sys.stderr)
    print(f"\nInter-amplicon CV: {overall_cv:.1f}%", file=sys.stderr)

    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight")
    print(f"Written: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <annotated_depth.txt> <out.pdf>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
