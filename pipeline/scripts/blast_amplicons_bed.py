#!/usr/bin/env python3
"""
blast_amplicons_bed.py
Generate a BED file of amplicon coordinates from BLAST primer results.

Uses the strict search as the primary source; falls back to the permissive
search for any primer missing from the strict results (e.g. amplicon3_rev).

Output BED6: chrom, start, end, name, score, strand
Coordinates are 0-based half-open (BED convention).

Usage: blast_amplicons_bed.py <strict.txt> <permissive.txt> <out.bed>
"""
import sys
from collections import defaultdict

COLS = ("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sstrand")


def parse_best_hits(filepath):
    """Parse a BLAST tabular output file and return the best hit per primer.

    Args:
        filepath: Path to a blastn -outfmt 6 output file with the columns
                  defined in COLS (including sstrand as column 13).

    Returns:
        Dict mapping primer name (qseqid) to its best-hit row dict.
        Only the first hit per primer is retained (blastn orders by score).
    """
    hits = {}
    with open(filepath) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            d = dict(zip(COLS, parts))
            if d["qseqid"] not in hits:
                hits[d["qseqid"]] = d
    return hits


def main(strict_path, permissive_path, out_path):
    """Generate a BED file of amplicon coordinates from BLAST primer results.

    Combines strict and permissive BLAST hits: the strict hit is preferred for
    each primer; the permissive hit is used as fallback when a primer has no
    strict hit (e.g. amplicon3_rev whose E-value just exceeded the strict threshold).

    Args:
        strict_path:     Path to the strict BLAST results (high-stringency hits).
        permissive_path: Path to the permissive BLAST results (fallback hits).
        out_path:        Destination path for the output BED6 file.
    """
    strict     = parse_best_hits(strict_path)
    permissive = parse_best_hits(permissive_path)

    # Determine amplicon IDs in order
    all_primers = list(dict.fromkeys(
        list(strict.keys()) + list(permissive.keys())
    ))
    amp_ids = sorted({p.rsplit("_", 1)[0] for p in all_primers})

    with open(out_path, "w") as fh:
        for amp_id in amp_ids:
            fwd = strict.get(f"{amp_id}_fwd") or permissive.get(f"{amp_id}_fwd")
            rev = strict.get(f"{amp_id}_rev") or permissive.get(f"{amp_id}_rev")

            if not fwd or not rev:
                continue  # skip if either primer has no hit at all

            chrom = fwd["sseqid"]
            # Amplicon span: min of all primer coordinates to max
            coords = [int(fwd["sstart"]), int(fwd["send"]),
                      int(rev["sstart"]), int(rev["send"])]
            start = min(coords) - 1   # convert to 0-based
            end   = max(coords)       # half-open

            amp_label = amp_id.replace("amplicon", "Amplicon_")
            fh.write(f"{chrom}\t{start}\t{end}\t{amp_label}\t.\t+\n")

    print(f"Written: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <strict.txt> <permissive.txt> <out.bed>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
