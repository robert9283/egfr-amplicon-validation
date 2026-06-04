#!/usr/bin/env python3
"""
check_adaptors.py — Verify which adaptor sequences are present in R1 and R2 reads.

Background
----------
In paired-end Illumina sequencing, adaptors appear at the 3' end of reads only
when a read is long enough to sequence through the entire insert (read-through).
The question is: which adaptor appears at the 3' end of R1, and which at R2?

Two possible scenarios:
  Scenario A — Same adaptor on both ends:
    5'-[Adp1]-[Fwd Primer]-[Insert]-[Rev Primer RC]-[Adp1 RC]-3'
    → R1 reads through → hits Adp1 RC → trim Adp1 from R1 (what we did)

  Scenario B — Different adaptors on each end (standard paired-end):
    5'-[Adp1]-[Fwd Primer]-[Insert]-[Rev Primer RC]-[Adp2 RC]-3'
    → R1 reads through → hits Adp2 RC → should trim Adp2 RC from R1

This script counts occurrences of both adaptors and their reverse complements
in R1 and R2 reads to determine which scenario applies.

Usage
-----
python3 check_adaptors.py <R1.fastq.gz> <R2.fastq.gz>
"""
import gzip
import sys


ADAPTOR1 = "AAGACTCGGCAGCATCTCCA"
ADAPTOR2 = "GCGATCGTCACTGTTCTCCA"


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence.

    Args:
        seq: DNA sequence string (A/C/G/T, case-insensitive).

    Returns:
        Reverse complement string in the same case convention.
    """
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def count_occurrences(fastq_gz, sequences):
    """Count reads containing each of the given sequences in a gzipped FASTQ file.

    Args:
        fastq_gz:  Path to a gzipped FASTQ file.
        sequences: List of DNA sequence strings to search for (exact substring match).

    Returns:
        Tuple of (total_reads, counts) where total_reads is the number of reads
        scanned and counts is a dict mapping each sequence to its hit count.
    """
    counts = {seq: 0 for seq in sequences}
    total = 0
    with gzip.open(fastq_gz, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:          # sequence line only
                total += 1
                line = line.strip()
                for seq in sequences:
                    if seq in line:
                        counts[seq] += 1
    return total, counts


def pct(n, total):
    """Format a count as a percentage string to three decimal places.

    Args:
        n:     Numerator count.
        total: Denominator count.

    Returns:
        Percentage string (e.g. '2.345%'), or 'N/A' if total is zero.
    """
    return f"{100 * n / total:.3f}%" if total else "N/A"


def main(r1_path, r2_path):
    """Check which adaptor sequences are present in R1 and R2 reads.

    Scans both FASTQ files for all four adaptor orientations (Adp1, Adp1 RC,
    Adp2, Adp2 RC) and prints a summary table to stdout, followed by an
    interpretation of Scenario A vs Scenario B.

    Args:
        r1_path: Path to R1 gzipped FASTQ file.
        r2_path: Path to R2 gzipped FASTQ file.
    """
    adp1    = ADAPTOR1
    adp1_rc = reverse_complement(ADAPTOR1)
    adp2    = ADAPTOR2
    adp2_rc = reverse_complement(ADAPTOR2)

    sequences = {
        "Adp1    (fwd)": adp1,
        "Adp1 RC (rev)": adp1_rc,
        "Adp2    (fwd)": adp2,
        "Adp2 RC (rev)": adp2_rc,
    }

    print("Sequences searched:")
    for label, seq in sequences.items():
        print(f"  {label}: {seq}")
    print()

    print(f"Scanning R1: {r1_path} ...")
    total_r1, counts_r1 = count_occurrences(r1_path, list(sequences.values()))

    print(f"Scanning R2: {r2_path} ...")
    total_r2, counts_r2 = count_occurrences(r2_path, list(sequences.values()))
    print()

    col_w = 16
    print(f"{'Sequence':<20} {'R1 count':>10} {'R1 %':>8}   {'R2 count':>10} {'R2 %':>8}")
    print("-" * 62)
    for label, seq in sequences.items():
        r1_n = counts_r1[seq]
        r2_n = counts_r2[seq]
        print(f"{label:<20} {r1_n:>10,} {pct(r1_n, total_r1):>8}   {r2_n:>10,} {pct(r2_n, total_r2):>8}")

    print("-" * 62)
    print(f"{'Total reads':<20} {total_r1:>10,} {'':>8}   {total_r2:>10,}")

    print()
    print("Interpretation:")
    print("  Scenario A (same adaptor both ends):")
    print("    → Adp1 RC dominant in R1, Adp2 RC dominant in R2")
    print("    → Cutadapt command was correct as-is")
    print()
    print("  Scenario B (different adaptors each end, standard paired-end):")
    print("    → Adp2 RC dominant in R1, Adp1 RC dominant in R2")
    print("    → Cutadapt should use -a RC(Adp2) for R1 and -A RC(Adp1) for R2")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <R1.fastq.gz> <R2.fastq.gz>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
