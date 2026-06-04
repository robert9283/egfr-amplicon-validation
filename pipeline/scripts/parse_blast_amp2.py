"""Parse the local phiX174 BLAST TSV and render a LaTeX summary table.

Input TSV columns (blastn outfmt 6 custom):
  qseqid  stitle  pident  length  evalue  bitscore

Queries are labelled pair_NNNN_R1 (R1 reads only — they carry the Amp2
forward primer and are therefore the informative end).  All hits are
against the phiX174 database, so every hit counts as phiX174.  Reads
absent from the TSV returned no hit.
"""

import sys


def main(blast_tsv: str, out_tex: str, n_reads: int = 200):
    """Parse a phiX174 BLAST TSV and write a LaTeX summary table.

    Counts R1 reads with a phiX174 hit vs. no hit and renders a two-row
    LaTeX table with counts and percentages.

    Args:
        blast_tsv: Path to the blastn output TSV (outfmt 6 custom:
                   qseqid stitle pident length evalue bitscore).
        out_tex:   Destination path for the generated .tex table file.
        n_reads:   Total number of reads in the subsample (default 200).
    """
    hit_ids = set()

    with open(blast_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 1:
                continue
            qid = parts[0]
            if qid.endswith("_R1"):
                hit_ids.add(qid)

    n_phix = len(hit_ids)
    n_nohit = n_reads - n_phix

    rows = (
        f"  \\textbf{{$\\phi$X174 (Sinsheimervirus phiX174)}} & "
        f"\\textbf{{{n_phix}}} & \\textbf{{{100 * n_phix // n_reads}\\%}} \\\\\n"
        f"  no hit & {n_nohit} & {100 * n_nohit // n_reads}\\% \\\\"
    )

    tex = r"""\begin{table}[H]
  \centering
  \caption{BLAST results for 200 subsampled Amplicon~2 unmapped R1 reads
    queried against the local $\phi$X174 genome (accession LC786485,
    \textit{Sinsheimervirus phiX174}; \texttt{blastn}, \(E\)-value
    \(\leq 10^{-5}\)).  R1 reads carry the Amp2 forward primer and are
    the informative end.  Every read with a hit matches $\phi$X174,
    confirming off-target priming on the Illumina spike-in control.}
  \label{tab:blast-amp2}
  \begin{tabular}{lrr}
    \toprule
    \textbf{Category} & \textbf{R1 reads} & \textbf{\%} \\
    \midrule
""" + rows + r"""
    \bottomrule
  \end{tabular}
\end{table}
"""

    with open(out_tex, "w") as fh:
        fh.write(tex)

    print(f"Written {out_tex}")
    print(f"phiX174 hits: {n_phix}/{n_reads}, no hit: {n_nohit}/{n_reads}")


if __name__ == "__main__":
    blast_tsv = sys.argv[1]
    out_tex = sys.argv[2]
    n_reads = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    main(blast_tsv, out_tex, n_reads)
