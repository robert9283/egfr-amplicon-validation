#!/usr/bin/env python3
"""
Shared utilities for parsing blastn tabular output (-outfmt 6).
Imported by blast_primer_table.py and blast_primer_table_relaxed.py.
"""


def parse_blast_hits(path):
    """Parse a blastn tabular output and return the best hit per query name.

    Handles both the standard 12-column format and the extended 13-column
    format that includes an sstrand column appended at the end.

    Args:
        path: Path to a blastn -outfmt 6 output file with columns:
              qseqid, sseqid, pident, length, mismatch, gapopen,
              qstart, qend, sstart, send, evalue, bitscore[, sstrand].

    Returns:
        Dict mapping query name (str) to a hit dict with keys:
        'chrom' (str), 'left' (int), 'right' (int), 'pident' (float),
        'evalue' (float), 'strand' (str or None).
        Only the first (best-scoring) hit per query is kept.
    """
    hits = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            qseqid = parts[0]
            if qseqid in hits:
                continue  # keep only best (first) hit per query
            sstart, send = int(parts[8]), int(parts[9])
            hits[qseqid] = {
                "chrom":  parts[1],
                "left":   min(sstart, send),
                "right":  max(sstart, send),
                "pident": float(parts[2]),
                "evalue": float(parts[10]),
                "strand": parts[12] if len(parts) > 12 else None,
            }
    return hits
