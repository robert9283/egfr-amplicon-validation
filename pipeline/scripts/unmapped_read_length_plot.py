#!/usr/bin/env python3
"""
unmapped_read_length_plot.py
Read-length histogram of unmapped reads from a FASTQ file.

Plots the distribution of read lengths to determine whether unmapped reads
are short (e.g. adapter dimers) or full-length. A vertical line marks the
nominal read length (default 150 bp) for reference.

Usage:
    unmapped_read_length_plot.py <unmapped.fastq.gz> <out.pdf> [--read-length N]
"""
import sys
import gzip
import collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


READ_LENGTH = 150


def read_lengths(fastq_gz):
    """Yield the length of each read sequence in a gzipped FASTQ file.

    Args:
        fastq_gz: Path to a gzipped FASTQ file (.fastq.gz).

    Yields:
        Integer length of each read sequence.
    """
    with gzip.open(fastq_gz, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip("\n")
            fh.readline()  # '+'
            fh.readline()  # quality
            yield len(seq)


def main(fastq_gz, out_pdf, read_length=READ_LENGTH):
    """Generate a read-length histogram of unmapped reads.

    Args:
        fastq_gz:    Path to the unmapped reads gzipped FASTQ file.
        out_pdf:     Destination path for the output PDF figure.
        read_length: Nominal read length in bp (default 150). Marked as a
                     vertical reference line in the histogram.
    """
    lengths = list(read_lengths(fastq_gz))
    total = len(lengths)
    counter = collections.Counter(lengths)

    print(f"Total unmapped reads: {total:,}", file=sys.stderr)
    for length, count in sorted(counter.items()):
        print(f"  {length} bp: {count:,} ({100*count/total:.1f}%)", file=sys.stderr)

    fig, ax = plt.subplots(figsize=(7, 4))

    bins = range(min(counter) - 1, max(counter) + 3)
    ax.hist(lengths, bins=list(bins), color="steelblue",
            edgecolor="white", linewidth=0.4, alpha=0.85)
    ax.axvline(read_length, color="crimson", linewidth=1.5,
               linestyle=":", label=f"Nominal read length ({read_length} bp)")

    ax.set_xlabel("Read length (bp)", fontsize=10)
    ax.set_ylabel("Number of reads", fontsize=10)
    ax.set_title(f"Read-length distribution of unmapped reads (n = {total:,})",
                 fontsize=11)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(
        lambda x, _: f"{int(x):,}"))
    ax.legend(fontsize=9)

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"Saved {out_pdf}", file=sys.stderr)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq_gz")
    parser.add_argument("out_pdf")
    parser.add_argument("--read-length", type=int, default=READ_LENGTH)
    args = parser.parse_args()
    main(args.fastq_gz, args.out_pdf, args.read_length)
