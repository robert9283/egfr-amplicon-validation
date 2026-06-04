#!/usr/bin/env python3
"""
Compare primer sequences to the reference genome at BLAST-determined positions.

For each primer in primers.fa, the reference sequence is extracted at the
BLAST hit coordinates using samtools faidx. Forward primers are compared
directly; reverse primers are compared to the reverse complement of the
reference forward strand. Mismatches are reported per position.
A LaTeX table summarising the comparison is written to the output path.

Usage: primer_vs_ref.py <primers.fa> <blast_hits.txt> <ref.fa> <output.tex>
"""

import re
import subprocess
import sys


def read_primers(fa_path):
    """Read a FASTA file of primer sequences into a dictionary.

    Args:
        fa_path: Path to the FASTA file containing primer sequences.

    Returns:
        Dict mapping primer name (str) to sequence (str, uppercase).
    """
    primers = {}
    name = None
    with open(fa_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
            elif name is not None:
                primers[name] = line.upper()
    return primers


def read_best_blast_hits(blast_path):
    """Read a BLAST tabular output and return the best (first) hit per query.

    Expects columns: qseqid sseqid pident length mismatch gapopen
    qstart qend sstart send evalue bitscore [sstrand].

    Args:
        blast_path: Path to the BLAST tabular output file (12 or 13 columns).

    Returns:
        Dict mapping query name (str) to a dict with keys:
            chrom (str), sstart (int), send (int), sstrand (str).
    """
    hits = {}
    with open(blast_path) as fh:
        for line in fh:
            fields = line.strip().split("\t")
            qseqid = fields[0]
            if qseqid in hits:
                continue
            sstart = int(fields[8])
            send = int(fields[9])
            if len(fields) >= 13:
                sstrand = fields[12]
            else:
                sstrand = "plus" if sstart <= send else "minus"
            hits[qseqid] = {
                "chrom": fields[1],
                "sstart": sstart,
                "send": send,
                "sstrand": sstrand,
            }
    return hits


def revcomp(seq):
    """Return the reverse complement of a DNA sequence.

    Args:
        seq: DNA sequence string (A/C/G/T, any case).

    Returns:
        Reverse complement string in uppercase.
    """
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1]


def fetch_ref(ref_fa, chrom, start, end):
    """Extract a sequence from a FASTA file using samtools faidx.

    Coordinates are 1-based inclusive (samtools convention). The chr prefix
    is stripped from chrom to match NCBI-style FASTA headers.

    Args:
        ref_fa:  Path to the indexed reference FASTA file.
        chrom:   Chromosome name (with or without 'chr' prefix).
        start:   1-based start coordinate (int).
        end:     1-based end coordinate (int).

    Returns:
        Reference sequence string (uppercase), empty string on failure.
    """
    chrom_name = re.sub(r"^chr", "", chrom)
    coord = f"{chrom_name}:{start}-{end}"
    result = subprocess.run(
        ["samtools", "faidx", ref_fa, coord],
        capture_output=True, text=True,
    )
    seq = "".join(
        line.strip()
        for line in result.stdout.splitlines()
        if not line.startswith(">")
    )
    return seq.upper()


def compare_sequences(primer, ref_fwd, sstrand):
    """Compare a primer sequence to the reference at its binding site.

    For minus-strand primers the reference forward sequence is reverse-
    complemented before comparison so both sequences read 5'→3'.

    Args:
        primer:   Primer sequence string (5'→3', uppercase).
        ref_fwd:  Reference sequence on the forward strand (uppercase).
        sstrand:  BLAST strand of the hit ('plus' or 'minus').

    Returns:
        Tuple (ref_oriented, mismatches) where ref_oriented (str) is the
        reference sequence oriented 5'→3' relative to the primer, and
        mismatches is a list of (position, primer_base, ref_base) tuples
        (1-based positions).
    """
    ref_oriented = ref_fwd if sstrand == "plus" else revcomp(ref_fwd)
    mismatches = [
        (i + 1, p, r)
        for i, (p, r) in enumerate(zip(primer, ref_oriented))
        if p != r
    ]
    return ref_oriented, mismatches


def write_table(comparisons, out_path):
    """Write a LaTeX table comparing primer sequences to the reference genome.

    Args:
        comparisons: List of dicts, each with keys: name (str), primer (str),
                     ref (str), sstrand (str), mismatches (list of tuples).
        out_path:    Destination path for the generated .tex table file.

    Returns:
        None.
    """
    # Output only the tabular environment (no float wrapper) so this file
    # can be \input-ed directly inside a tcolorbox or figure environment.
    lines = [
        "% Auto-generated by pipeline/scripts/primer_vs_ref.py -- do not edit manually.",
        r"\centering\small",
        r"\begin{tabular}{llccp{4.5cm}}",
        r"\toprule",
        r"\textbf{Primer} & \textbf{Strand} & \textbf{Length} "
        r"& \textbf{Mismatches} & \textbf{Detail} \\",
        r"\midrule",
    ]

    for c in comparisons:
        name_tex = c["name"].replace("_", r"\_")
        strand = "+" if c["sstrand"] == "plus" else "$-$"
        length = len(c["primer"])
        n_mm = len(c["mismatches"])
        if n_mm == 0:
            detail = r"\textit{exact match}"
        else:
            parts = [
                f"pos~{pos}: primer~{p}$\\to$ref~{r}"
                for pos, p, r in c["mismatches"]
            ]
            detail = "; ".join(parts)
        mm_cell = str(n_mm) if n_mm > 0 else "0"
        lines.append(
            f"\\texttt{{{name_tex}}} & {strand} & {length} & {mm_cell} & {detail} \\\\"
        )

    lines += [
        r"\bottomrule",
        r"\end{tabular}",
    ]

    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main(primers_fa, blast_path, ref_fa, out_path):
    """Run the primer-vs-reference comparison and write a LaTeX table.

    Args:
        primers_fa: Path to the primers FASTA file.
        blast_path: Path to the BLAST tabular output (best hits).
        ref_fa:     Path to the indexed reference FASTA.
        out_path:   Destination path for the LaTeX table.

    Returns:
        None.
    """
    primers = read_primers(primers_fa)
    hits = read_best_blast_hits(blast_path)

    comparisons = []
    for name, primer_seq in primers.items():
        if name not in hits:
            print(f"WARNING: no BLAST hit for {name}, skipping", file=sys.stderr)
            continue
        hit = hits[name]
        start = min(hit["sstart"], hit["send"])
        end = max(hit["sstart"], hit["send"])
        ref_fwd = fetch_ref(ref_fa, hit["chrom"], start, end)
        if not ref_fwd:
            print(f"WARNING: could not fetch reference for {name}", file=sys.stderr)
            continue
        ref_oriented, mismatches = compare_sequences(
            primer_seq, ref_fwd, hit["sstrand"]
        )
        comparisons.append({
            "name": name,
            "primer": primer_seq,
            "ref": ref_oriented,
            "sstrand": hit["sstrand"],
            "mismatches": mismatches,
        })
        status = "MATCH" if not mismatches else f"{len(mismatches)} MISMATCH(ES)"
        print(f"{name}: {status}", file=sys.stderr)

    write_table(comparisons, out_path)
    print(f"Written: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} <primers.fa> <blast_hits.txt> <ref.fa> <output.tex>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
