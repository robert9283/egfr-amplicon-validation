#!/usr/bin/env python3
"""
variant_table.py
Parse PASS variants from Mutect2 (VCF.gz) and VarDict (VCF), annotate each
with the amplicon and exon it falls in (from amplicons.bed + exon mapping),
and render a LaTeX table.

Usage:
    variant_table.py <mutect2_filtered.vcf.gz> <vardict.vcf> \
                     <amplicons.bed> <amplicon_exon_overlap.tex> <out.tex>

amplicon_exon_overlap.tex is parsed to extract the exon number for each
amplicon name (e.g. Amplicon_1 -> exon 19).
"""
import gzip
import re
import sys


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def open_vcf(path):
    """Open a VCF or gzipped VCF file for reading.

    Args:
        path: Path to a .vcf or .vcf.gz file.

    Returns:
        File-like object opened in text mode.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_amplicons(bed_path):
    """Parse amplicon intervals from a BED file with chr prefix stripped.

    Args:
        bed_path: Path to a BED file (at least 4 columns: chrom, start, end, name).

    Returns:
        List of (chrom_nochr, start, end, name) tuples where chrom_nochr has
        any 'chr' prefix removed to match NCBI-style VCF chromosome names.
    """
    regions = []
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            chrom = parts[0].lstrip("chr")
            regions.append((chrom, int(parts[1]), int(parts[2]), parts[3]))
    return regions


def assign_amplicon(chrom, pos, regions):
    """Return the amplicon name containing a 1-based genomic position, or None.

    Args:
        chrom:   Chromosome string (with or without 'chr' prefix).
        pos:     1-based genomic position (int).
        regions: List of (chrom_nochr, start, end, name) tuples from parse_amplicons().

    Returns:
        Amplicon name string if pos falls within any region, or None.
    """
    chrom = str(chrom).lstrip("chr")
    for (c, start, end, name) in regions:
        if c == chrom and start <= pos <= end:
            return name
    return None


def parse_exon_map(tex_path):
    """
    Extract amplicon->exon mapping from amplicon_exon_overlap_table.tex.
    Looks for rows like: Amplicon\_1 & 19 & ...
    Falls back to a hardcoded map if parsing fails.
    """
    exon_map = {}
    try:
        with open(tex_path) as fh:
            for line in fh:
                # Match lines like: Amplicon\_1 & 19 &
                m = re.search(r"Amplicon\\_(\d+)\s*&\s*(\d+)", line)
                if m:
                    amp_name = f"Amplicon_{m.group(1)}"
                    exon_map[amp_name] = int(m.group(2))
    except FileNotFoundError:
        pass

    if not exon_map:
        # Fallback from our blast/refFlat analysis
        exon_map = {
            "Amplicon_1": 19,
            "Amplicon_2": 18,
            "Amplicon_3": 21,
            "Amplicon_4": 20,
        }
    return exon_map


def parse_ann(info):
    """Extract the first HGVSp protein annotation from a SnpEff ANN INFO field.

    Args:
        info: VCF INFO field string (semicolon-separated key=value pairs).

    Returns:
        HGVSp annotation string from the first ANN entry (e.g. 'p.Glu746_Ala750del'),
        or None if the ANN field is absent or the annotation is empty.
    """
    for field in info.split(";"):
        if field.startswith("ANN="):
            first_ann = field[4:].split(",")[0]
            parts = first_ann.split("|")
            if len(parts) > 10 and parts[10]:
                return parts[10]  # e.g. p.Thr790Met
    return None


def af_from_info(info, fmt, sample):
    """Extract allele frequency (AF) from FORMAT/SAMPLE or INFO fields.

    Tries the FORMAT AF field first (as used by Mutect2), then falls back
    to the INFO AF= field (as used by VarDict).

    Args:
        info:   VCF INFO field string.
        fmt:    VCF FORMAT field string (colon-separated keys).
        sample: VCF SAMPLE field string (colon-separated values).

    Returns:
        AF as a float (0.0–1.0), or None if not found or not parseable.
    """
    # Try FORMAT AF field first (Mutect2)
    keys = fmt.split(":")
    vals = sample.split(":")
    fmt_dict = dict(zip(keys, vals))
    if "AF" in fmt_dict:
        try:
            return float(fmt_dict["AF"].split(",")[0])
        except ValueError:
            pass
    # Try INFO AF field (VarDict)
    for field in info.split(";"):
        if field.startswith("AF="):
            try:
                return float(field[3:])
            except ValueError:
                pass
    return None


def dp_from_format(fmt, sample):
    """Extract read depth (DP) from FORMAT/SAMPLE fields.

    Args:
        fmt:    VCF FORMAT field string (colon-separated keys).
        sample: VCF SAMPLE field string (colon-separated values).

    Returns:
        DP as an integer, or None if not found or not parseable.
    """
    keys = fmt.split(":")
    vals = sample.split(":")
    fmt_dict = dict(zip(keys, vals))
    if "DP" in fmt_dict:
        try:
            return int(fmt_dict["DP"])
        except ValueError:
            pass
    return None


def parse_vcf(path, caller, regions, exon_map):
    """Parse PASS-filtered variants from a VCF or gzipped VCF file.

    Args:
        path:     Path to a .vcf or .vcf.gz file (SnpEff-annotated).
        caller:   Caller label string (e.g. 'Mutect2' or 'VarDict').
        regions:  List of amplicon intervals from parse_amplicons().
        exon_map: Dict mapping amplicon name to exon number from parse_exon_map().

    Returns:
        List of variant dicts with keys: chrom, pos, ref, alt, af, dp,
        type, amp, exon, caller, mutation.
    """
    variants = []
    with open_vcf(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            chrom, pos, _, ref, alt, _, filt = parts[:7]
            if filt != "PASS":
                continue
            pos = int(pos)
            info = parts[7] if len(parts) > 7 else "."
            fmt  = parts[8] if len(parts) > 8 else "."
            samp = parts[9] if len(parts) > 9 else "."

            af = af_from_info(info, fmt, samp)
            dp = dp_from_format(fmt, samp)
            mutation = parse_ann(info)

            amp = assign_amplicon(chrom, pos, regions)
            exon = exon_map.get(amp) if amp else None
            amp_display = amp.replace("_", "~") if amp else "---"

            # Determine variant type
            if len(ref) == len(alt) == 1:
                vtype = "SNV"
            elif len(ref) > len(alt):
                vtype = "Deletion"
            else:
                vtype = "Insertion"

            variants.append({
                "chrom":    chrom.lstrip("chr"),
                "pos":      pos,
                "ref":      ref,
                "alt":      alt,
                "af":       af,
                "dp":       dp,
                "type":     vtype,
                "amp":      amp_display,
                "exon":     exon,
                "caller":   caller,
                "mutation": mutation,
            })
    return variants


def fmt_ref_alt(seq, maxlen=16):
    """Truncate a long REF or ALT sequence for display in a table.

    Args:
        seq:    DNA sequence string.
        maxlen: Maximum display length before truncation (default 16).

    Returns:
        Original string if len <= maxlen, else first 13 characters + '...'.
    """
    if len(seq) > maxlen:
        return seq[:13] + "..."
    return seq


def fmt_num(n):
    """Format an integer with thousands separators, or '---' for None.

    Args:
        n: Integer to format, or None.

    Returns:
        String with comma-separated thousands groups, or '---'.
    """
    if n is None:
        return "---"
    return f"{n:,}"


def fmt_af(af):
    """Format an allele frequency as a LaTeX percentage string, or '---' for None.

    Args:
        af: Allele frequency as a float (0.0–1.0), or None.

    Returns:
        Percentage string (e.g. '48.0\\%'), or '---'.
    """
    if af is None:
        return "---"
    return f"{100*af:.1f}\\%"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(mutect2_vcf, vardict_vcf, bed_path, exon_tex, out_path):
    """Parse PASS variants from both callers and write a LaTeX variant table.

    Args:
        mutect2_vcf: Path to the SnpEff-annotated Mutect2 VCF (.vcf.gz).
        vardict_vcf: Path to the SnpEff-annotated VarDict VCF (.vcf).
        bed_path:    Path to amplicons.bed (for amplicon assignment).
        exon_tex:    Path to amplicon_exon_overlap_table.tex (for exon mapping).
        out_path:    Destination path for the generated .tex table file.
    """
    regions  = parse_amplicons(bed_path)
    exon_map = parse_exon_map(exon_tex)

    variants  = parse_vcf(mutect2_vcf, "Mutect2",   regions, exon_map)
    variants += parse_vcf(vardict_vcf,  "VarDict",  regions, exon_map)

    # Sort by position, then caller
    variants.sort(key=lambda v: (v["pos"], v["caller"]))

    lines = [
        r"% Auto-generated by pipeline/scripts/variant_table.py",
        r"% via rule variant_table.",
        r"\begin{table}[H]",
        r"  \centering",
        r"  \footnotesize",
        r"  \resizebox{\textwidth}{!}{%",
        r"  \begin{tabular}{lllllrrlll}",
        r"    \toprule",
        r"    \textbf{Chr} & \textbf{Position} & \textbf{Ref} & \textbf{Alt}"
        r" & \textbf{Mutation} & \textbf{VAF} & \textbf{Depth}"
        r" & \textbf{Amplicon} & \textbf{Exon} & \textbf{Caller} \\",
        r"    \midrule",
    ]

    for v in variants:
        ref = fmt_ref_alt(v["ref"])
        alt = fmt_ref_alt(v["alt"])
        exon_str = str(v["exon"]) if v["exon"] is not None else "---"
        mut_str = (v["mutation"].replace("_", "\\_") if v["mutation"] else "---")
        lines.append(
            f"    {v['chrom']} & {fmt_num(v['pos'])} & "
            f"\\texttt{{{ref}}} & \\texttt{{{alt}}} & "
            f"\\texttt{{{mut_str}}} & "
            f"{fmt_af(v['af'])} & {fmt_num(v['dp'])} & "
            f"{v['amp']} & {exon_str} & {v['caller']} \\\\"
        )

    lines += [
        r"    \bottomrule",
        r"  \end{tabular}%",
        r"  }% end resizebox",
        r"  \caption{%",
        r"    PASS-filtered variants identified by Mutect2 and VarDictJava,",
        r"    annotated with SnpEff (hg19 database).",
        r"    Mutation: protein-level consequence in HGVS notation (e.g.\ "
        r"\texttt{p.Thr790Met} = T790M); \texttt{---} if no coding effect.",
        r"    Amplicon and exon assignments are based on the BLAST-derived",
        r"    amplicon coordinates (\texttt{amplicons.bed}) and the RefFlat",
        r"    exon mapping (Table~\ref{tab:amplicon-exon-overlap}).",
        r"    VAF: variant allele frequency (fraction of reads carrying the",
        r"    alternative allele). Depth reflects the effective pileup depth",
        r"    used by each caller.",
        r"  }",
        r"  \label{tab:variants}",
        r"\end{table}",
    ]

    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"Written {len(variants)} PASS variants to {out_path}",
          file=sys.stderr)
    for v in variants:
        print(f"  {v['chrom']}:{v['pos']} {v['ref']}->{v['alt']} "
              f"VAF={fmt_af(v['af'])} {v['amp']} ex{v['exon']} [{v['caller']}]",
              file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 6:
        print(f"Usage: {sys.argv[0]} <mutect2.vcf.gz> <vardict.vcf> "
              f"<amplicons.bed> <amplicon_exon_overlap.tex> <out.tex>")
        sys.exit(1)
    main(*sys.argv[1:])
