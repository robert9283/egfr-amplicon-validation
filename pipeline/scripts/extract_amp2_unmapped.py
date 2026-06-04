"""Extract Amplicon 2 unmapped read pairs and subsample N pairs to FASTA.

Usage:
    python3 extract_amp2_unmapped.py R1.fastq.gz R2.fastq.gz out.fa [--n 200] [--seed 42]

Logic:
  - Reads the unmapped R1/R2 FASTQ files (produced by unmapped_split_reads).
  - Assigns each read pair to an amplicon by primer matching of the first 20 bp
    of R1 against forward primers (≤2 mismatches), with fallback to R2 vs
    reverse primers.
  - Retains only pairs assigned to Amplicon 2.
  - Randomly subsamples --n pairs (default 200) and writes them as FASTA with
    headers  >pair_NNNN_R1  and  >pair_NNNN_R2.
"""

import gzip
import random
import sys
import argparse
from pathlib import Path

# Primers (from config.yaml, hard-coded here to avoid YAML dependency in Docker)
FORWARD_PRIMERS = [
    "TGCCAGCACATGAGCTACAC",  # Amplicon 1
    "CTGCATGATGAGCTGCACGA",  # Amplicon 2
    "GCAAGCTCCAAGGACATCGA",  # Amplicon 3
    "GGAGCCCAGCACTTTGATCT",  # Amplicon 4
]
REVERSE_PRIMERS = [
    "TGGCAGAGGTGGAAATCAGG",  # Amplicon 1
    "AGTGCCTCAGGAGCCTGTAG",  # Amplicon 2
    "GCACGGTGGAGGTGAGGTTG",  # Amplicon 3
    "ACAGCAGAGCCCAGCAAGTT",  # Amplicon 4
]


def hamming(a: str, b: str) -> int:
    """Count mismatches between two strings up to the length of the shorter one.

    Args:
        a: First sequence string.
        b: Second sequence string.

    Returns:
        Number of positions where the characters differ.
    """
    return sum(x != y for x, y in zip(a, b))


def assign_amplicon(r1_seq: str, r2_seq: str, max_mm: int = 2, trim: int = 20) -> int:
    """Return the 1-based amplicon index for a read pair, or 0 if unassigned.

    Matches the first `trim` bases of R1 against forward primers, then falls
    back to matching R2 against reverse primers.

    Args:
        r1_seq: R1 read sequence string.
        r2_seq: R2 read sequence string.
        max_mm: Maximum allowed mismatches for a primer match (default 2).
        trim:   Number of bases from the read start to compare (default 20).

    Returns:
        1-based amplicon index (1–4), or 0 if no primer matched.
    """
    r1_prefix = r1_seq[:trim].upper()
    for i, fwd in enumerate(FORWARD_PRIMERS):
        if hamming(r1_prefix, fwd[:trim]) <= max_mm:
            return i + 1
    r2_prefix = r2_seq[:trim].upper()
    for i, rev in enumerate(REVERSE_PRIMERS):
        if hamming(r2_prefix, rev[:trim]) <= max_mm:
            return i + 1
    return 0


def read_fastq_pairs(r1_path: str, r2_path: str):
    """Yield paired read records from two FASTQ (or gzipped FASTQ) files.

    Args:
        r1_path: Path to the R1 FASTQ or .fastq.gz file.
        r2_path: Path to the R2 FASTQ or .fastq.gz file.

    Yields:
        Tuples of (header1, seq1, header2, seq2) strings for each read pair.
    """
    open_fn = gzip.open if r1_path.endswith(".gz") else open
    with open_fn(r1_path, "rt") as f1, \
         (gzip.open(r2_path, "rt") if r2_path.endswith(".gz") else open(r2_path)) as f2:
        while True:
            h1 = f1.readline().rstrip("\n")
            if not h1:
                break
            s1 = f1.readline().rstrip("\n")
            f1.readline()  # +
            f1.readline()  # qual
            h2 = f2.readline().rstrip("\n")
            s2 = f2.readline().rstrip("\n")
            f2.readline()
            f2.readline()
            yield h1, s1, h2, s2


def main():
    """Extract Amplicon 2 unmapped read pairs and write a subsampled FASTA.

    Parses command-line arguments for R1/R2 FASTQ paths, output FASTA path,
    subsample size (--n, default 200), and random seed (--seed, default 42).
    Only pairs assigned to Amplicon 2 by primer matching are retained.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("r1")
    parser.add_argument("r2")
    parser.add_argument("out_fasta")
    parser.add_argument("--n", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    amp2_pairs = []
    for h1, s1, h2, s2 in read_fastq_pairs(args.r1, args.r2):
        if assign_amplicon(s1, s2) == 2:
            amp2_pairs.append((s1, s2))

    print(f"Amplicon 2 unmapped pairs found: {len(amp2_pairs)}", file=sys.stderr)

    random.seed(args.seed)
    sample = random.sample(amp2_pairs, min(args.n, len(amp2_pairs)))

    with open(args.out_fasta, "w") as fh:
        for i, (s1, s2) in enumerate(sample):
            fh.write(f">pair_{i:04d}_R1\n{s1}\n")
            fh.write(f">pair_{i:04d}_R2\n{s2}\n")

    print(f"Wrote {len(sample)} pairs ({2*len(sample)} records) to {args.out_fasta}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
