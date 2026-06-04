# -------------------------------------------------------
# BLAST rules: primer localisation against hg19
# Goal: map each primer pair to exact hg19 coordinates,
#       linking amplicon IDs (Table 5) to genomic regions
#       (Table 9).
# -------------------------------------------------------

# (DATA, RESULTS, RESOURCES defined in Snakefile)

HG19_FA      = f"{RESOURCES}/hg19.fa.gz"
BLAST_DB     = f"{RESOURCES}/blast_db/hg19"
HG19_REFFLAT = f"{RESOURCES}/hg19_refFlat.txt"


rule download_hg19_fasta:
    """
    Download the hg19 reference genome FASTA from UCSC (~1 GB compressed).
    Run once; Snakemake skips this rule if the output already exists.
    """
    output:
        HG19_FA,
    log:
        f"{RESULTS}/logs/download_hg19_fasta.log",
    shell:
        """
        wget -qO {output} \
            http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz \
            2> {log}
        """


rule build_blast_db:
    """
    Build a BLAST nucleotide database from the hg19 FASTA.
    makeblastdb requires uncompressed input, so we decompress on the fly.
    """
    input:
        HG19_FA,
    output:
        multiext(BLAST_DB, ".nhr", ".nin", ".nsq"),
    log:
        f"{RESULTS}/logs/build_blast_db.log",
    shell:
        """
        mkdir -p $(dirname {BLAST_DB})
        gunzip -c {input} | makeblastdb \
            -in - \
            -dbtype nucl \
            -title hg19 \
            -out {BLAST_DB} \
            2> {log}
        """


rule blast_primers:
    """
    Align all forward and reverse primer sequences against the hg19 BLAST
    database to determine the exact chromosomal coordinates of each amplicon.

    The 8 primers (4 forward + 4 reverse) are written to a FASTA file and
    queried with blastn using parameters tuned for short sequences:
      -task blastn-short : preset optimised for sequences < 50 bp
      -word_size 7       : smaller seed for short queries
      -evalue 1e-3       : permissive E-value (primers are ~20-30 bp)
      -perc_identity 90  : require near-perfect match
      -max_target_seqs 1 : report only the best hit per primer

    Output is tabular (-outfmt 6): qseqid sseqid pident length mismatch
    gapopen qstart qend sstart send evalue bitscore.
    """
    input:
        db = multiext(BLAST_DB, ".nhr", ".nin", ".nsq"),
    output:
        primers_fa = f"{RESULTS}/blast/primers.fa",
        blast_out  = f"{RESULTS}/blast/primers_vs_hg19.txt",
    log:
        f"{RESULTS}/logs/blast_primers.log",
    params:
        db        = BLAST_DB,
        primers   = config["primers"],
    run:
        # Write primer FASTA
        fwd = params.primers["forward"]
        rev = params.primers["reverse"]
        with open(output.primers_fa, "w") as fh:
            for i, seq in enumerate(fwd, 1):
                fh.write(f">amplicon{i}_fwd\n{seq}\n")
            for i, seq in enumerate(rev, 1):
                fh.write(f">amplicon{i}_rev\n{seq}\n")

        # Run BLAST
        shell(
            "blastn"
            " -query {output.primers_fa}"
            " -db {params.db}"
            " -task blastn-short"
            " -word_size 7"
            " -evalue 1e-3"
            " -perc_identity 90"
            " -max_target_seqs 1"
            " -outfmt '6 qseqid sseqid pident length mismatch gapopen"
            " qstart qend sstart send evalue bitscore sstrand'"
            " -out {output.blast_out}"
            " 2> {log}"
        )


rule blast_primers_offtarget:
    """
    Search for secondary (off-target) primer binding sites genome-wide.

    Uses permissive parameters compared to the primary search:
      -max_target_seqs 10  : report up to 10 hits per primer (not just the best)
      -perc_identity 70    : allow up to 30% mismatches (partial complementarity
                             sufficient for non-specific PCR priming)
      -evalue 1            : very permissive; short primers with partial matches
                             will have high E-values
    Hits on the primary chr7 target are expected; hits elsewhere (e.g. chr4)
    explain off-target reads observed in the BAM coverage analysis.
    """
    input:
        primers_fa = f"{RESULTS}/blast/primers.fa",
        db         = multiext(BLAST_DB, ".nhr", ".nin", ".nsq"),
    output:
        f"{RESULTS}/blast/primers_vs_hg19_offtarget.txt",
    params:
        db = BLAST_DB,
    log:
        f"{RESULTS}/logs/blast_primers_offtarget.log",
    shell:
        """
        blastn \
            -query {input.primers_fa} \
            -db {params.db} \
            -task blastn-short \
            -word_size 7 \
            -evalue 1 \
            -perc_identity 70 \
            -max_target_seqs 10 \
            -outfmt '6 qseqid sseqid pident length mismatch gapopen \
                     qstart qend sstart send evalue bitscore sstrand' \
            -out {output} \
            2> {log}
        """


rule blast_amplicons_bed:
    """
    Generate a BED file of the four amplicon spans from the BLAST primer
    results (strict search, falling back to permissive for missing primers).
    """
    input:
        strict     = f"{RESULTS}/blast/primers_vs_hg19.txt",
        permissive = f"{RESULTS}/blast/primers_vs_hg19_offtarget.txt",
    output:
        bed = f"{RESULTS}/blast/amplicons.bed",
    log:
        f"{RESULTS}/logs/blast_amplicons_bed.log",
    shell:
        """
        python3 /pipeline/scripts/blast_amplicons_bed.py \
            {input.strict} {input.permissive} {output.bed} 2> {log}
        """


rule blast_gene_intersect:
    """
    Intersect BLAST-defined amplicon coordinates with hg19 RefFlat gene
    annotation and render a LaTeX table of gene overlaps (bp and %).
    """
    input:
        bed      = f"{RESULTS}/blast/amplicons.bed",
        genes    = f"{RESOURCES}/hg19_genes.bed",
    output:
        intersect = f"{RESULTS}/blast/amplicons_vs_genes.txt",
        tex       = f"{RESULTS}/tables/blast_gene_intersect_table.tex",
    log:
        f"{RESULTS}/logs/blast_gene_intersect.log",
    shell:
        """
        mkdir -p $(dirname {output.tex})
        bedtools intersect -a {input.bed} -b {input.genes} -wo \
            | sort -k1,1 -k2,2n \
            > {output.intersect} 2> {log}
        python3 /pipeline/scripts/blast_gene_intersect_table.py \
            {output.intersect} {output.tex} 2>> {log}
        """


rule amplicon_length_table:
    """
    Derive amplicon, primer, insert, and R1/R2 overlap lengths from the
    relaxed primer coordinate table and the primer sequences in config.yaml.
    config.yaml is a static mount at /pipeline/config.yaml; it is passed as
    a param rather than an input to avoid Snakemake DAG resolution issues
    with absolute paths outside the results directory.
    """
    input:
        relaxed_tex = f"{RESULTS}/tables/blast_primer_table_relaxed.tex",
    output:
        tex = f"{RESULTS}/tables/amplicon_length_table.tex",
    params:
        config_file = "/pipeline/config.yaml",
    log:
        f"{RESULTS}/logs/amplicon_length_table.log",
    shell:
        """
        mkdir -p $(dirname {output.tex})
        python3 /pipeline/scripts/amplicon_length_table.py \
            {params.config_file} {input.relaxed_tex} {output.tex} 2> {log}
        """


rule blast_offtarget_table:
    """
    Parse the permissive off-target BLAST output and render a LaTeX table
    listing all secondary primer hits (primary EGFR hit excluded per primer).
    """
    input:
        permissive = f"{RESULTS}/blast/primers_vs_hg19_offtarget.txt",
        strict     = f"{RESULTS}/blast/primers_vs_hg19.txt",
    output:
        tex = f"{RESULTS}/tables/blast_offtarget_table.tex",
    log:
        f"{RESULTS}/logs/blast_offtarget_table.log",
    shell:
        """
        mkdir -p $(dirname {output.tex})
        python3 /pipeline/scripts/blast_offtarget_table.py \
            {input.permissive} {input.strict} {output.tex} 2> {log}
        """


rule blast_primer_table_relaxed:
    """
    Build a complete primer coordinate table by combining strict and permissive
    BLAST results. Primers missing from the strict search (e.g. amplicon3_rev)
    are filled from the permissive search and marked with a dagger footnote.
    """
    input:
        strict     = f"{RESULTS}/blast/primers_vs_hg19.txt",
        permissive = f"{RESULTS}/blast/primers_vs_hg19_offtarget.txt",
    output:
        tex = f"{RESULTS}/tables/blast_primer_table_relaxed.tex",
    log:
        f"{RESULTS}/logs/blast_primer_table_relaxed.log",
    shell:
        """
        mkdir -p $(dirname {output.tex})
        python3 /pipeline/scripts/blast_primer_table_relaxed.py \
            {input.strict} {input.permissive} {output.tex} 2> {log}
        """


rule blast_primer_table:
    """
    Parse blastn output and render a LaTeX table of primer genomic coordinates.
    Links each amplicon ID to its hg19 span (fwd primer start → rev primer end).
    """
    input:
        blast_out = f"{RESULTS}/blast/primers_vs_hg19.txt",
    output:
        tex = f"{RESULTS}/tables/blast_primer_table.tex",
    log:
        f"{RESULTS}/logs/blast_primer_table.log",
    shell:
        """
        mkdir -p $(dirname {output.tex})
        python3 /pipeline/scripts/blast_primer_table.py \
            {input.blast_out} {output.tex} 2> {log}
        """


rule download_hg19_refflat:
    """
    Download the UCSC hg19 refFlat gene annotation table.
    Unlike hg19_genes.bed (gene-level intervals), refFlat contains per-transcript
    exon coordinates required for the EGFR exon structure figure.
    Run once; Snakemake skips if the output already exists.
    """
    output:
        HG19_REFFLAT,
    log:
        f"{RESULTS}/logs/download_hg19_refflat.log",
    shell:
        """
        wget -qO - \
            http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz \
            2> {log} \
            | gunzip -c > {output}
        """


rule primer_properties_table:
    """
    Compute GC content and Tm (nearest-neighbour, SantaLucia 1998) for all
    8 primers and render a LaTeX table. Used to assess whether Amplicon 2's
    primers are outliers relative to the other pairs (primer efficiency
    hypothesis for the observed depth imbalance).
    """
    output:
        tex = f"{RESULTS}/tables/primer_properties_table.tex",
    params:
        config_file = "/pipeline/config.yaml",
    log:
        f"{RESULTS}/logs/primer_properties_table.log",
    shell:
        """
        mkdir -p $(dirname {output.tex})
        python3 /pipeline/scripts/primer_properties_table.py \
            {params.config_file} {output.tex} 2> {log}
        """


rule amplicon_exon_overlap_table:
    """
    Generate a LaTeX table comparing each amplicon span with its target EGFR
    exon, showing amplicon coordinates, exon coordinates, and the flanking
    intronic regions on each side.
    """
    input:
        refflat = HG19_REFFLAT,
        bed     = f"{RESULTS}/blast/amplicons.bed",
    output:
        tex = f"{RESULTS}/tables/amplicon_exon_overlap_table.tex",
    log:
        f"{RESULTS}/logs/amplicon_exon_overlap_table.log",
    shell:
        """
        mkdir -p $(dirname {output.tex})
        python3 /pipeline/scripts/amplicon_exon_overlap_table.py \
            {input.refflat} {input.bed} {output.tex} 2> {log}
        """


rule egfr_exon_table:
    """
    Generate a LaTeX table of EGFR kinase-domain exons (18-21) with hg19
    coordinates from the UCSC refFlat table and known somatic mutation
    annotations for each exon.
    """
    input:
        refflat = HG19_REFFLAT,
    output:
        tex = f"{RESULTS}/tables/egfr_exon_table.tex",
    log:
        f"{RESULTS}/logs/egfr_exon_table.log",
    shell:
        """
        mkdir -p $(dirname {output.tex})
        python3 /pipeline/scripts/egfr_exon_table.py \
            {input.refflat} {output.tex} 2> {log}
        """


rule egfr_exon_figure:
    """
    Generate a PDF figure showing the EGFR exon structure (hg19, NM_005228)
    with the four amplicon spans overlaid.  Uses the UCSC refFlat table for
    exon coordinates and the BLAST-derived amplicons.bed for amplicon positions.
    """
    input:
        refflat  = HG19_REFFLAT,
        bed      = f"{RESULTS}/blast/amplicons.bed",
    output:
        pdf = f"{RESULTS}/figures/egfr_exon_figure.pdf",
    log:
        f"{RESULTS}/logs/egfr_exon_figure.log",
    shell:
        """
        mkdir -p $(dirname {output.pdf})
        python3 /pipeline/scripts/egfr_exon_figure.py \
            {input.refflat} {input.bed} {output.pdf} 2> {log}
        """


rule extract_amp2_unmapped:
    """
    Extract Amplicon 2 unmapped read pairs and subsample 200 pairs to FASTA.

    Reads the split unmapped R1/R2 FASTQs (produced by unmapped_split_reads),
    assigns each pair to an amplicon by primer matching (same logic as
    primer_assignment_table.py: first 20 bp of R1 vs forward primers, ≤2 mm),
    retains only Amplicon 2 pairs, randomly subsamples 200, and writes FASTA
    with headers >pair_NNNN_R1 / >pair_NNNN_R2 for BLAST submission.
    """
    input:
        r1 = f"{RESULTS}/mapping/unmapped_R1.fastq.gz",
        r2 = f"{RESULTS}/mapping/unmapped_R2.fastq.gz",
    output:
        fa = f"{RESULTS}/blast/amp2_unmapped_sample.fa",
    log:
        f"{RESULTS}/logs/extract_amp2_unmapped.log",
    shell:
        """
        mkdir -p $(dirname {output.fa})
        python3 /pipeline/scripts/extract_amp2_unmapped.py \
            {input.r1} {input.r2} {output.fa} --n 200 --seed 42 2> {log}
        """


rule blast_amp2_vs_phix174:
    """
    BLAST 200 subsampled Amplicon 2 unmapped R1 reads against the local
    phiX174 database (accession LC786485, built by build_phix174_blast_db).

    R1 reads carry the Amp2 forward primer at their 5' end and are the
    informative end; R2 reads are excluded.  The local database keeps the
    pipeline self-contained (no internet required) and runs inside Docker.

    E-value 1e-5 catches full-length phiX174 matches; reads that return no
    hit at this threshold are of unknown origin.
    """
    input:
        fa = f"{RESULTS}/blast/amp2_unmapped_sample.fa",
        db = multiext(f"{RESOURCES}/phix174_db/phix174", ".nhr", ".nin", ".nsq"),
    output:
        tsv = f"{RESULTS}/blast/amp2_vs_phix174.tsv",
    params:
        db = f"{RESOURCES}/phix174_db/phix174",
    log:
        f"{RESULTS}/logs/blast_amp2_vs_phix174.log",
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        # Extract R1 reads only into a temp file
        python3 -c "
import sys
write = False
for line in open('{input.fa}'):
    if line.startswith('>'):
        write = '_R1' in line
    if write:
        sys.stdout.write(line)
" > /tmp/amp2_r1.fa

        blastn \
            -query /tmp/amp2_r1.fa \
            -db {params.db} \
            -task blastn \
            -evalue 1e-5 \
            -max_target_seqs 1 \
            -outfmt "6 qseqid stitle pident length evalue bitscore" \
            -out {output.tsv} \
            2> {log}
        echo "Hits: $(wc -l < {output.tsv}) / 200 R1 reads" >> {log}
        """


rule parse_blast_amp2:
    """
    Parse the local phiX174 BLAST TSV and render a LaTeX table showing how
    many of the 200 subsampled Amplicon 2 R1 reads match phiX174 vs have
    no hit.  R1 reads are used because they carry the Amp2 forward primer.
    """
    input:
        tsv = f"{RESULTS}/blast/amp2_vs_phix174.tsv",
    output:
        tex = f"{RESULTS}/tables/blast_amp2_summary.tex",
    log:
        f"{RESULTS}/logs/parse_blast_amp2.log",
    shell:
        """
        mkdir -p $(dirname {output.tex})
        python3 /pipeline/scripts/parse_blast_amp2.py \
            {input.tsv} {output.tex} 200 2> {log}
        """


rule download_phix174:
    """
    Download the bacteriophage phiX174 genome (accession LC786485) from NCBI
    in FASTA format.  This is the same entry that matched the Amplicon 2
    unmapped reads in the remote BLAST search.  The genome is only 5,386 bp
    so the download is fast.  Run once; Snakemake skips if output exists.
    """
    output:
        fa = f"{RESOURCES}/phix174_LC786485.fa",
    log:
        f"{RESULTS}/logs/download_phix174.log",
    shell:
        """
        mkdir -p $(dirname {output.fa})
        curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi\
?db=nucleotide&id=LC786485&rettype=fasta&retmode=text" \
            > {output.fa} 2> {log}
        """


rule build_phix174_blast_db:
    """
    Build a local BLAST nucleotide database from the phiX174 genome.
    Running BLAST locally against this tiny 5 kb database is instant and
    does not require internet access, allowing a permissive e-value without
    flooding the output with spurious hits from unrelated organisms.
    """
    input:
        fa = f"{RESOURCES}/phix174_LC786485.fa",
    output:
        multiext(f"{RESOURCES}/phix174_db/phix174", ".nhr", ".nin", ".nsq"),
    log:
        f"{RESULTS}/logs/build_phix174_blast_db.log",
    shell:
        """
        mkdir -p {RESOURCES}/phix174_db
        makeblastdb \
            -in {input.fa} \
            -dbtype nucl \
            -title phiX174 \
            -out {RESOURCES}/phix174_db/phix174 \
            2> {log}
        """


