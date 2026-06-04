#!/usr/bin/env python3
"""
annotate_amplicon_depth.py
Annotate a samtools depth file with the amplicon name for each position.

Reads amplicons.bed to get amplicon intervals, then for each row in the
depth file finds which amplicon the position falls in and appends the
amplicon name as a fourth column.

Chromosome names are normalised before matching (chr prefix stripped from
the BED file to match the NCBI-style names used in the BAM / depth file).

Usage: annotate_amplicon_depth.py <amplicons.bed> <amplicon_depth.txt> <out.txt>

Output columns: chrom, pos (1-based), depth, amplicon_name
"""
import sys


def parse_bed(path):
    """Return amplicon intervals from a BED file with chr prefix stripped.

    Args:
        path: Path to the BED file (at least 4 columns: chrom, start, end, name).

    Returns:
        List of (chrom_nochr, start, end, name) tuples where chrom_nochr has
        any 'chr' prefix removed to match NCBI-style BAM chromosome names.
    """
    intervals = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            chrom = parts[0].lstrip("chr")   # normalise to no-chr
            intervals.append((chrom, int(parts[1]), int(parts[2]), parts[3]))
    return intervals


def find_amplicon(chrom, pos, intervals):
    """Return amplicon name if pos falls within any interval, else None.
    BED intervals are 0-based half-open; samtools depth pos is 1-based.
    """
    for ivl_chrom, start, end, name in intervals:
        if chrom == ivl_chrom and start < pos <= end:
            return name
    return None


def main(bed_path, depth_path, out_path):
    """Annotate a samtools depth file with amplicon names and write output.

    Args:
        bed_path:   Path to amplicons.bed defining amplicon intervals.
        depth_path: Path to the samtools depth output file (chrom, pos, depth).
        out_path:   Destination path for the annotated depth file
                    (chrom, pos, depth, amplicon_name).
    """
    intervals = parse_bed(bed_path)

    unmatched = 0
    with open(depth_path) as depth_fh, open(out_path, "w") as out_fh:
        for line in depth_fh:
            if not line.strip():
                continue
            chrom, pos, depth = line.strip().split("\t")
            pos = int(pos)
            amp = find_amplicon(chrom, pos, intervals)
            if amp is None:
                unmatched += 1
                continue
            out_fh.write(f"{chrom}\t{pos}\t{depth}\t{amp}\n")

    if unmatched:
        print(f"Warning: {unmatched} positions not matched to any amplicon",
              file=sys.stderr)
    print(f"Written: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <amplicons.bed> <amplicon_depth.txt> <out.txt>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
