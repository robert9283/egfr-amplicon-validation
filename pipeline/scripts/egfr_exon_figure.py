#!/usr/bin/env python3
"""
egfr_exon_figure.py
Plot EGFR exon structure (hg19) with the four amplicons overlaid.

Reads the UCSC hg19 refFlat table to extract exon coordinates for the
canonical EGFR transcript (NM_005228), then draws a schematic showing
exon boxes with exon numbers, intron lines with direction chevrons, and
the four amplicon spans below the gene track.

Usage: egfr_exon_figure.py <refFlat.txt> <amplicons.bed> <out.pdf>
"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

CANONICAL = "NM_005228"          # canonical EGFR RefSeq transcript (hg19)
AMPLICON_COLORS = ["#2166ac", "#4dac26", "#d01c8b", "#f1a340"]
PADDING = 1500                   # bp to show either side of the amplicon span


# ── Parsers ──────────────────────────────────────────────────────────────────

def parse_refflat(path, transcript):
    """Return strand and exon intervals for a transcript from a UCSC refFlat file.

    Args:
        path:       Path to the UCSC hg19 refFlat flat-file table.
        transcript: RefSeq accession to extract (e.g. 'NM_005228').

    Returns:
        Tuple of (strand, exon_list) where strand is '+' or '-' and
        exon_list is a list of (exon_start, exon_end) int tuples in
        ascending coordinate order (0-based half-open, BED convention).

    Raises:
        ValueError: If the transcript is not found in the file.
    """
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 11 or parts[1] != transcript:
                continue
            strand      = parts[3]
            exon_starts = [int(x) for x in parts[9].rstrip(",").split(",")]
            exon_ends   = [int(x) for x in parts[10].rstrip(",").split(",")]
            return strand, list(zip(exon_starts, exon_ends))
    raise ValueError(f"Transcript {transcript} not found in {path}")


def parse_bed(path):
    """Return amplicon records from a BED file, sorted by start position.

    Args:
        path: Path to a BED file (at least 4 columns: chrom, start, end, name).

    Returns:
        List of (chrom, start, end, name) tuples sorted by start coordinate.
    """
    records = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            records.append((parts[0], int(parts[1]), int(parts[2]), parts[3]))
    return sorted(records, key=lambda r: r[1])   # sort by start


def exon_number(idx, total, strand):
    """Return the 1-based exon number in transcript order.

    Exon coordinates are stored in ascending genomic order. For minus-strand
    genes the highest-coordinate exon is exon 1.

    Args:
        idx:    0-based index into the ascending-coordinate exon list.
        total:  Total number of exons in the transcript.
        strand: Strand character ('+' or '-').

    Returns:
        1-based exon number in transcript order.
    """
    if strand == "-":
        # minus strand: highest coord exon is exon 1
        return total - idx
    return idx + 1


def assign_label_rows(amplicons, char_width_bp):
    """
    Assign each amplicon to a label row (0 = bottom, 1 = lower) so that
    labels do not overlap.  Uses greedy left-to-right placement.
    Returns list of row indices parallel to amplicons.
    """
    # Sort by centre, keeping original index
    order = sorted(range(len(amplicons)), key=lambda i: (amplicons[i][1] + amplicons[i][2]) / 2)
    rows  = [0] * len(amplicons)
    # row_end tracks the rightmost x used by each row
    row_end = {}
    for i in order:
        _, s, e, name = amplicons[i]
        label_len = len(name) * char_width_bp
        cx        = (s + e) / 2
        lx        = cx - label_len / 2
        rx        = cx + label_len / 2
        placed    = False
        for row in range(4):           # try up to 4 rows
            if row_end.get(row, -1e12) < lx - 50:
                rows[i]    = row
                row_end[row] = rx
                placed      = True
                break
        if not placed:
            rows[i]    = 3
            row_end[3] = max(row_end.get(3, -1e12), rx)
    return rows


# ── Main ──────────────────────────────────────────────────────────────────────

def main(refflat_path, bed_path, out_path):
    """Generate a schematic EGFR exon structure figure with amplicons overlaid.

    Args:
        refflat_path: Path to the UCSC hg19 refFlat table.
        bed_path:     Path to amplicons.bed (BLAST-derived amplicon coordinates).
        out_path:     Destination path for the output PDF figure.
    """
    strand, exons = parse_refflat(refflat_path, CANONICAL)
    amplicons     = parse_bed(bed_path)

    amp_min    = min(s for _, s, _, _ in amplicons)
    amp_max    = max(e for _, _, e, _ in amplicons)
    view_start = amp_min - PADDING
    view_end   = amp_max + PADDING
    view_span  = view_end - view_start

    total_exons = len(exons)
    visible     = [(s, e, i) for i, (s, e) in enumerate(exons)
                   if s < view_end and e > view_start]

    # Approximate character width in bp for label overlap detection
    char_width_bp = view_span * 0.008

    label_rows    = assign_label_rows(amplicons, char_width_bp)
    row_y_offset  = 0.32   # vertical step per label row

    # ── Canvas ──────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(13, 4.2))
    ax.set_xlim(view_start, view_end)
    ax.set_ylim(-2.2, 2.4)
    ax.axis("off")

    # ── Gene backbone ────────────────────────────────────────────────────────
    ax.hlines(1.0, view_start, view_end, color="#555555", linewidth=1.4, zorder=1)

    # ── Direction chevrons on introns ────────────────────────────────────────
    chevron_step = 1500
    chevron_dir  = 1 if strand == "+" else -1
    exon_ranges  = [(s, e) for s, e, _ in visible]

    for cx in range(view_start + chevron_step, view_end, chevron_step):
        if any(s - 80 <= cx <= e + 80 for s, e in exon_ranges):
            continue
        ax.annotate("",
                    xy=(cx + chevron_dir * 220, 1.0),
                    xytext=(cx, 1.0),
                    arrowprops=dict(arrowstyle="-|>", color="#aaaaaa",
                                    lw=0.9, mutation_scale=9))

    # ── Exon boxes ───────────────────────────────────────────────────────────
    EH = 0.50
    for s, e, idx in visible:
        xs = max(s, view_start)
        xe = min(e, view_end)
        ax.add_patch(mpatches.FancyBboxPatch(
            (xs, 1.0 - EH / 2), xe - xs, EH,
            boxstyle="square,pad=0", linewidth=0.9,
            edgecolor="#1a6496", facecolor="#4393c3", zorder=3
        ))
        ex_num = exon_number(idx, total_exons, strand)
        ax.text((xs + xe) / 2, 1.0 + EH / 2 + 0.13, f"Ex{ex_num}",
                ha="center", va="bottom", fontsize=7.5, color="#1a3a5c",
                fontweight="bold")

    # ── Amplicon bars ────────────────────────────────────────────────────────
    AY = 0.10
    AH = 0.32
    for (_, s, e, name), color, row in zip(amplicons, AMPLICON_COLORS, label_rows):
        ax.add_patch(mpatches.FancyBboxPatch(
            (s, AY - AH / 2), e - s, AH,
            boxstyle="square,pad=0", linewidth=1.0,
            edgecolor=color, facecolor=color, alpha=0.35, zorder=3
        ))
        ax.add_patch(mpatches.FancyBboxPatch(
            (s, AY - AH / 2), e - s, AH,
            boxstyle="square,pad=0", linewidth=1.0,
            edgecolor=color, facecolor="none", zorder=4
        ))
        label  = name.replace("_", "\u2009")
        label_y = AY - AH / 2 - 0.18 - row * row_y_offset
        ax.text((s + e) / 2, label_y, label,
                ha="center", va="top", fontsize=8, color=color,
                fontweight="bold")

    # ── Dashed guides: amplicon edges → gene track ───────────────────────────
    for _, s, e, _ in amplicons:
        for x in (s, e):
            ax.vlines(x, AY + AH / 2, 1.0 - EH / 2,
                      colors="#cccccc", linewidth=0.65,
                      linestyles="dashed", zorder=2)

    # ── Gene label ───────────────────────────────────────────────────────────
    strand_label = "plus" if strand == "+" else "minus"
    ax.text(view_start + 150, 1.68,
            f"EGFR  (NM_005228, chr7, {strand_label} strand, hg19)",
            ha="left", va="bottom", fontsize=9.5,
            fontstyle="italic", color="#222222")

    # ── Genomic coordinate axis ──────────────────────────────────────────────
    axis_y = -1.35
    ax.hlines(axis_y, view_start, view_end, color="#888888", linewidth=0.7)

    step   = 5000
    start_ = (view_start // step + 1) * step
    for tp in range(start_, view_end, step):
        ax.vlines(tp, axis_y - 0.07, axis_y + 0.07,
                  colors="#888888", linewidth=0.7)
        ax.text(tp, axis_y - 0.18, f"{tp / 1e6:.3f} Mb",
                ha="center", va="top", fontsize=6.5, color="#666666")

    plt.tight_layout(pad=0.4)
    plt.savefig(out_path, bbox_inches="tight")
    print(f"Written: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <refFlat.txt> <amplicons.bed> <out.pdf>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
