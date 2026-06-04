#!/usr/bin/env python3
"""
insert_size_plot.py
Per-amplicon insert size distribution plot for paired-end amplicon sequencing.

For each amplicon defined in amplicons.bed, reads that overlap the amplicon
region are extracted from the BAM via `samtools view`, and the TLEN field
(template/insert length) of properly paired, mapped reads is collected.
A 2x2 multi-panel histogram is produced, one panel per amplicon.

The read length (150 bp by default) is marked as a vertical dashed line to
show which amplicons are shorter than the read length and therefore produce
read-through into the adaptor sequence.

Usage:
    insert_size_plot.py <aln.bam> <amplicons.bed> <out.pdf> [--read-length N]
"""
import subprocess
import sys
import collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


READ_LENGTH = 150  # default Illumina MiSeq read length


def parse_amplicons(bed_path):
    """Return amplicon records from a BED file.

    Args:
        bed_path: Path to a BED file with at least 4 columns
                  (chrom, start, end, name).

    Returns:
        List of (chrom, start, end, name) tuples with 0-based coordinates.
    """
    amplicons = []
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            amplicons.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))
    return amplicons


def collect_tlens(bam, chrom, start, end):
    """
    Run samtools view over the region and return a list of positive TLEN values
    for properly paired, mapped reads (flags: -f 2 = proper pair, -F 4 = mapped,
    -F 256 = primary alignment only).
    Only positive TLEN is kept to count each fragment once.
    The 'chr' prefix is stripped from the chromosome name to match the BAM's
    NCBI-style naming convention (e.g. '7' instead of 'chr7').
    """
    chrom_nochr = chrom.lstrip("chr")
    region = f"{chrom_nochr}:{start+1}-{end}"  # samtools uses 1-based regions
    cmd = ["samtools", "view", "-f", "2", "-F", "2308", bam, region]
    result = subprocess.run(cmd, capture_output=True, text=True)
    tlens = []
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 9:
            continue
        try:
            tlen = int(fields[8])
        except ValueError:
            continue
        if tlen > 0:  # keep one per pair (positive strand mate)
            tlens.append(tlen)
    return tlens


def main(bam, bed, out_pdf, read_length=READ_LENGTH):
    """Generate a per-amplicon insert size distribution figure.

    Args:
        bam:         Path to the aligned BAM file.
        bed:         Path to amplicons.bed defining amplicon regions.
        out_pdf:     Destination path for the output PDF figure.
        read_length: Nominal read length in bp (default 150). Marked as a
                     vertical reference line in each panel.
    """
    amplicons = parse_amplicons(bed)
    # Sort by amplicon name for consistent panel order
    amplicons.sort(key=lambda x: x[3])

    n = len(amplicons)
    ncols = 2
    nrows = (n + 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(10, 4 * nrows))
    axes = axes.flatten()

    for i, (chrom, start, end, name) in enumerate(amplicons):
        ax = axes[i]
        tlens = collect_tlens(bam, chrom, start, end)

        if tlens:
            # Bin range: clip at [50, 500] to keep the plot readable
            clipped = [t for t in tlens if 50 <= t <= 500]
            ax.hist(clipped, bins=60, color="steelblue", edgecolor="white",
                    linewidth=0.3, alpha=0.85)
            median = sorted(tlens)[len(tlens) // 2]
            mean = sum(tlens) / len(tlens)
            ax.axvline(median, color="darkorange", linewidth=1.5,
                       linestyle="--", label=f"Median {median} bp")
            ax.axvline(read_length, color="crimson", linewidth=1.2,
                       linestyle=":", label=f"Read length {read_length} bp")
            ax.set_yscale("log")
            ax.legend(fontsize=7)
            print(f"{name}: n={len(tlens):,}  mean={mean:.0f}  median={median}",
                  file=sys.stderr)
        else:
            ax.text(0.5, 0.5, "No paired reads", transform=ax.transAxes,
                    ha="center", va="center", color="gray")
            print(f"{name}: no paired reads found", file=sys.stderr)

        expected_len = end - start
        ax.set_title(f"{name}  (expected ≈{expected_len} bp)", fontsize=10)
        ax.set_xlabel("Insert size (bp)", fontsize=9)
        ax.set_ylabel("Read pairs", fontsize=9)
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(
            lambda x, _: f"{int(x):,}"))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(
            lambda x, _: f"{int(x):,}"))

    # Hide unused panels
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("Per-amplicon insert size distribution\n"
                 "(red dotted line = read length; orange dashed = median)",
                 fontsize=11, y=1.01)
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"Saved {out_pdf}", file=sys.stderr)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("bam")
    parser.add_argument("bed")
    parser.add_argument("out_pdf")
    parser.add_argument("--read-length", type=int, default=READ_LENGTH)
    args = parser.parse_args()
    main(args.bam, args.bed, args.out_pdf, args.read_length)
