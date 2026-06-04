#!/usr/bin/env python3
"""
amplicon_depth_plot.py
Per-base depth line plot for the four BLAST-defined EGFR amplicons.

Creates a 2x2 figure, one panel per amplicon (ordered by exon number),
showing depth along the genomic position with exon boundaries marked.
Summary statistics (mean, min, max, %>=10x, %>=100x) are printed to
stderr for inclusion in the figure caption.

Usage: amplicon_depth_plot.py <annotated_depth.txt> <amplicons.bed>
                               <refFlat.txt> <out.pdf>
"""
import sys
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

CANONICAL = "NM_005228"
TARGET_EXONS = [18, 19, 20, 21]

# Amplicon colours consistent with egfr_exon_figure.py
COLORS = {
    "Amplicon_1": "#2166ac",   # Ex19 — blue
    "Amplicon_2": "#4dac26",   # Ex18 — green
    "Amplicon_3": "#d01c8b",   # Ex21 — magenta
    "Amplicon_4": "#f1a340",   # Ex20 — orange
}

# Display order: by exon number (18→19→20→21)
PANEL_ORDER = ["Amplicon_2", "Amplicon_1", "Amplicon_4", "Amplicon_3"]
EXON_LABEL  = {"Amplicon_2": 18, "Amplicon_1": 19,
                "Amplicon_4": 20, "Amplicon_3": 21}


# ── Parsers ──────────────────────────────────────────────────────────────────

def parse_annotated_depth(path):
    """Parse an annotated depth file into per-amplicon position/depth lists.

    Args:
        path: Path to the annotated depth file produced by annotate_amplicon_depth.py
              (columns: chrom, pos, depth, amplicon_name).

    Returns:
        Dict mapping amplicon_name (str) to a list of (pos, depth) int tuples.
    """
    data = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            pos, depth, amp = int(parts[1]), int(parts[2]), parts[3]
            data[amp].append((pos, depth))
    return data


def parse_amplicons_bed(path):
    """Parse a BED file into a dict of amplicon coordinates.

    Args:
        path: Path to a BED file (at least 4 columns: chrom, start, end, name).

    Returns:
        Dict mapping amplicon name (str) to (chrom, start, end) tuple.
    """
    result = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.strip().split("\t")
            result[p[3]] = (p[0], int(p[1]), int(p[2]))
    return result


def parse_refflat_exons(path, transcript):
    """Return exon intervals for a transcript from a UCSC refFlat file.

    Args:
        path:       Path to the UCSC hg19 refFlat flat-file table.
        transcript: RefSeq accession to extract (e.g. 'NM_005228').

    Returns:
        List of (exon_start, exon_end) tuples in ascending coordinate order.

    Raises:
        ValueError: If the transcript is not found in the file.
    """
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 11 or parts[1] != transcript:
                continue
            starts = [int(x) for x in parts[9].rstrip(",").split(",")]
            ends   = [int(x) for x in parts[10].rstrip(",").split(",")]
            return list(zip(starts, ends))
    raise ValueError(f"{transcript} not found")


# ── Summary statistics ────────────────────────────────────────────────────────

def summarise(depths):
    """Compute summary statistics for a list of per-base depth values.

    Args:
        depths: List of integer depth values.

    Returns:
        Dict with keys 'mean', 'min', 'max', 'pct10' (% bases >= 10×),
        and 'pct100' (% bases >= 100×).
    """
    n     = len(depths)
    mean_ = sum(depths) / n
    return {
        "mean":  mean_,
        "min":   min(depths),
        "max":   max(depths),
        "pct10": 100 * sum(d >= 10  for d in depths) / n,
        "pct100":100 * sum(d >= 100 for d in depths) / n,
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def main(depth_path, bed_path, refflat_path, out_path):
    """Generate a 2×2 per-base depth plot for the four EGFR amplicons.

    Args:
        depth_path:   Path to the annotated depth file (chrom, pos, depth, amplicon).
        bed_path:     Path to amplicons.bed (for amplicon coordinates).
        refflat_path: Path to the UCSC hg19 refFlat table (for exon boundaries).
        out_path:     Destination path for the output PDF figure.
    """
    depth_data = parse_annotated_depth(depth_path)
    amplicons  = parse_amplicons_bed(bed_path)
    exons      = parse_refflat_exons(refflat_path, CANONICAL)

    fig, axes = plt.subplots(2, 2, figsize=(13, 7))
    axes_flat = axes.flatten()

    print("Summary statistics per amplicon:", file=sys.stderr)
    print(f"{'Amplicon':<12} {'Exon':>5} {'Mean':>8} {'Min':>8} "
          f"{'Max':>8} {'>=10x':>7} {'>=100x':>8}", file=sys.stderr)

    for idx, amp_name in enumerate(PANEL_ORDER):
        ax    = axes_flat[idx]
        color = COLORS[amp_name]
        ex_num = EXON_LABEL[amp_name]

        # Depth data
        pts    = sorted(depth_data[amp_name], key=lambda x: x[0])
        pos    = [p for p, _ in pts]
        depths = [d for _, d in pts]

        # Summary stats
        stats = summarise(depths)
        print(
            f"{amp_name:<12} {ex_num:>5} {stats['mean']:>8.0f} "
            f"{stats['min']:>8} {stats['max']:>8} "
            f"{stats['pct10']:>6.1f}% {stats['pct100']:>7.1f}%",
            file=sys.stderr
        )

        # ── Fill + line ──────────────────────────────────────────────────────
        ax.fill_between(pos, depths, alpha=0.25, color=color)
        ax.plot(pos, depths, color=color, linewidth=1.0)

        # ── Exon boundary lines ──────────────────────────────────────────────
        es, ee = exons[ex_num - 1]
        for x, label in ((es, "exon start"), (ee, "exon end")):
            ax.axvline(x, color="#333333", linewidth=1.0,
                       linestyle="--", alpha=0.7)
        # Shade the exon region
        ax.axvspan(es, ee, color="#cccccc", alpha=0.25, zorder=0)

        # ── Threshold lines (10x and 100x) ───────────────────────────────────
        ymax = max(depths)
        for thresh, label in ((10, "10×"), (100, "100×")):
            ax.axhline(thresh, color="#999999", linewidth=0.7,
                       linestyle=":", alpha=0.8)
            # Only label if threshold is visible (>1% of y range)
            if thresh > ymax * 0.01:
                ax.text(pos[-1], thresh, f" {label}",
                        va="center", fontsize=6.5, color="#888888")

        # ── Axes formatting ──────────────────────────────────────────────────
        ax.set_title(f"{amp_name.replace('_', ' ')}  —  Exon {ex_num}",
                     fontsize=10, fontweight="bold", color=color, pad=5)
        ax.set_xlabel("Genomic position (hg19, chr7)", fontsize=8)
        ax.set_ylabel("Depth", fontsize=8)
        ax.tick_params(labelsize=7)

        # Format x-axis as Mb with enough precision for ~200 bp amplicons
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f"{x/1e6:.5f} Mb")
        )
        ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=4))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=15, ha="right")

        # Format y-axis with k suffix
        ax.yaxis.set_major_formatter(
            mticker.FuncFormatter(lambda y, _: f"{y/1e3:.0f}k" if y >= 1000 else str(int(y)))
        )

        # Annotate mean depth inside panel
        ax.text(0.98, 0.96, f"mean = {stats['mean']/1e3:.1f}k×",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=8, color=color,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=2))

        ax.set_xlim(pos[0], pos[-1])
        ax.set_ylim(bottom=0)
        ax.spines[["top", "right"]].set_visible(False)

    plt.suptitle("Per-base sequencing depth within EGFR amplicons (hg19)",
                 fontsize=11, y=1.01)
    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight")
    print(f"Written: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} <annotated_depth.txt> <amplicons.bed> "
              f"<refFlat.txt> <out.pdf>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
