# -------------------------------------------------------
# Mapping rules: flagstat, idxstats, per-base depth
# Uses the pre-aligned BAM file (hg19) from the task.
# -------------------------------------------------------

# (DATA, RESULTS, SAMPLE, BAM, R1, R2 defined in Snakefile)


rule bam_header:
    """
    Extract the full BAM header and write it to a plain-text file.
    The header encodes key metadata about the dataset:
      @HD  - file format version and sort order
      @SQ  - reference sequence dictionary (chromosome names and lengths)
      @RG  - read group (sample name, sequencing platform)
      @PG  - processing history (which programs were run and in which order)
    Inspecting the @PG lines reveals the upstream alignment pipeline and
    whether any preprocessing steps (e.g. adaptor trimming) were applied
    before we received the data.
    """
    input:
        bam = f"{DATA}/{BAM}",
    output:
        f"{RESULTS}/qc/bam_header.txt",
    log:
        f"{RESULTS}/logs/bam_header.log",
    shell:
        """
        mkdir -p $(dirname {output})
        samtools view -H {input.bam} > {output} 2> {log}
        """


rule task3_flagstat:
    """
    Overall mapping statistics: total reads, mapped, paired, duplicates, etc.
    Answers Task 3: percentage of reads mapped to the human genome.
    """
    input:
        bam = f"{DATA}/{BAM}",
    output:
        f"{RESULTS}/mapping/flagstat.txt",
    log:
        f"{RESULTS}/logs/flagstat.log",
    shell:
        """
        mkdir -p $(dirname {output})
        samtools flagstat {input.bam} > {output} 2> {log}
        """


rule download_hg19_genes:
    """
    Download and convert the UCSC hg19 refFlat table to BED6 format.
    Run once; Snakemake will skip this rule if the output already exists.
    """
    output:
        "/pipeline/resources/hg19_genes.bed",
    shell:
        """
        wget -qO- http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz \
            | gunzip \
            | awk 'BEGIN{{OFS="\\t"}} {{print $3, $5, $6, $1, ".", $4}}' \
            | sort -k1,1 -k2,2n \
            > {output}
        """


rule task2_gene_region:
    """
    Identify which gene regions the reads map to (Task 2).

    Pipeline:
      1. bedtools genomecov -bg  : per-base coverage BED (only covered bases)
      2. bedtools merge          : collapse adjacent covered bases into amplicon intervals
      3. bedtools intersect      : annotate each interval with overlapping hg19 RefFlat genes
      4. gene_region_table_bedtools.py : render LaTeX table
    """
    input:
        bam        = f"{DATA}/{BAM}",
        annotation = "/pipeline/resources/hg19_genes.bed",
    output:
        covered    = f"{RESULTS}/gene_region/covered.bed",
        merged     = f"{RESULTS}/gene_region/amplicons.bed",
        annotated  = f"{RESULTS}/gene_region/annotated.bed",
        tex        = f"{RESULTS}/tables/gene_region_table.tex",
    params:
        min_depth  = config["bedtools"]["min_amplicon_depth"],
    log:
        f"{RESULTS}/logs/gene_region_table.log",
    shell:
        """
        mkdir -p $(dirname {output.tex}) $(dirname {output.covered})

        # 1. Per-base coverage in BED graph format (-bg: only output covered bases)
        #    Filter to positions with depth >= min_amplicon_depth to exclude scattered noise reads.
        bedtools genomecov -ibam {input.bam} -bg \
            | awk '$4 >= {params.min_depth}' \
            > {output.covered}

        # 2. Merge contiguous high-coverage intervals into amplicon-level regions.
        #    Add chr prefix to match UCSC-style names in hg19_genes.bed.
        bedtools merge -i {output.covered} \
            | awk 'BEGIN{{OFS="\\t"}} !/^chr/{{$1="chr"$1}} {{print}}' \
            > {output.merged}

        # 3. Annotate each amplicon interval with overlapping gene names
        #    (hg19_genes.bed uses UCSC chr-prefixed names; the BAM uses NCBI names
        #     without the prefix — normalised in step 2 above)
        #    -wa: write the original amplicon interval
        #    -wb: append the matching gene annotation columns
        bedtools intersect \
            -a {output.merged} \
            -b {input.annotation} \
            -wa -wb \
            | sort -k1,1 -k2,2n \
            > {output.annotated}

        # 4. Render LaTeX table
        python3 /pipeline/scripts/gene_region_table_bedtools.py \
            {output.annotated} {output.tex} 2>> {log}
        """


rule task2_idxstats:
    """
    Per-chromosome read counts.
    Helps identify which chromosomes/regions the reads map to (Task 2).
    """
    input:
        bam = f"{DATA}/{BAM}",
        bai = f"{DATA}/{BAM}.bai",
    output:
        f"{RESULTS}/mapping/idxstats.txt",
    log:
        f"{RESULTS}/logs/idxstats.log",
    shell:
        """
        samtools idxstats {input.bam} > {output} 2> {log}
        """


rule task4_depth:
    """
    Per-base sequencing depth across all mapped positions.
    Answers Task 4: coverage of each amplicon.
    -a: output all positions including zero-coverage bases.
    -Q: minimum base quality threshold.
    """
    input:
        bam = f"{DATA}/{BAM}",
    output:
        f"{RESULTS}/mapping/depth.txt",
    params:
        min_base_quality  = config["samtools"]["min_base_quality"],
    log:
        f"{RESULTS}/logs/depth.log",
    shell:
        """
        samtools depth \
            -a \
            -Q {params.min_base_quality} \
            {input.bam} \
            > {output} \
            2> {log}
        """


rule amplicon_depth:
    """
    Per-base sequencing depth restricted to the four BLAST-defined amplicon
    spans (amplicons.bed).  Used for panel validation: checks that all four
    amplicons were sequenced at sufficient depth.

    Flags:
      -a  : report every position, including zero-depth bases (needed for
            accurate min-depth and breadth-of-coverage calculations)
      -b  : restrict output to positions within the given BED file
      -Q  : minimum base quality (same threshold as task4_depth)

    Output columns: chrom, pos (1-based), depth
    """
    input:
        bam = f"{DATA}/{BAM}",
        bed = f"{RESULTS}/blast/amplicons.bed",
    output:
        f"{RESULTS}/mapping/amplicon_depth.txt",
    params:
        min_base_quality = config["samtools"]["min_base_quality"],
    log:
        f"{RESULTS}/logs/amplicon_depth.log",
    shell:
        """
        mkdir -p $(dirname {output})
        samtools depth \
            -a \
            -b <(sed 's/^chr//' {input.bed}) \
            -Q {params.min_base_quality} \
            {input.bam} \
            > {output} \
            2> {log}
        """


rule annotate_amplicon_depth:
    """
    Annotate samtools depth output with the amplicon name for each position.
    Adds a fourth column (amplicon name) by intersecting positions with
    amplicons.bed.  Required as input for the per-amplicon depth summary.
    """
    input:
        bed   = f"{RESULTS}/blast/amplicons.bed",
        depth = f"{RESULTS}/mapping/amplicon_depth.txt",
    output:
        f"{RESULTS}/mapping/amplicon_depth_annotated.txt",
    log:
        f"{RESULTS}/logs/annotate_amplicon_depth.log",
    shell:
        """
        python3 /pipeline/scripts/annotate_amplicon_depth.py \
            {input.bed} {input.depth} {output} 2> {log}
        """


rule amplicon_depth_plot:
    """
    Per-base depth line plot for the four EGFR amplicons (2x2 panels,
    ordered by exon number).  Exon boundaries are marked with dashed
    vertical lines; the exon body is shaded.  Summary statistics are
    printed to the log for inclusion in the figure caption.
    """
    input:
        depth   = f"{RESULTS}/mapping/amplicon_depth_annotated.txt",
        bed     = f"{RESULTS}/blast/amplicons.bed",
        refflat = "/pipeline/resources/hg19_refFlat.txt",
    output:
        pdf = f"{RESULTS}/figures/amplicon_depth_plot.pdf",
    log:
        f"{RESULTS}/logs/amplicon_depth_plot.log",
    shell:
        """
        mkdir -p $(dirname {output.pdf})
        python3 /pipeline/scripts/amplicon_depth_plot.py \
            {input.depth} {input.bed} {input.refflat} {output.pdf} 2> {log}
        """


rule coverage_uniformity_plot:
    """
    Bar chart comparing mean depth across the four EGFR amplicons.
    Annotates each bar with mean depth and within-amplicon CV.
    Reports inter-amplicon CV to the log.
    """
    input:
        depth = f"{RESULTS}/mapping/amplicon_depth_annotated.txt",
    output:
        pdf = f"{RESULTS}/figures/coverage_uniformity_plot.pdf",
    log:
        f"{RESULTS}/logs/coverage_uniformity_plot.log",
    shell:
        """
        mkdir -p $(dirname {output.pdf})
        python3 /pipeline/scripts/coverage_uniformity_plot.py \
            {input.depth} {output.pdf} 2> {log}
        """


rule on_target_rate:
    """
    Compute the fraction of mapped reads that fall within the four
    BLAST-defined amplicon spans.

    Steps:
      1. Count total mapped reads from flagstat output.
      2. Use bedtools intersect -u to extract reads overlapping any
         amplicon, then count with samtools view -c.
         (-u: output each read at most once, even if it overlaps multiple
          intervals; strip chr prefix from BED to match BAM naming)
      3. Write a small LaTeX table with total, on-target, and off-target
         counts and percentages.
    """
    input:
        bam      = f"{DATA}/{BAM}",
        bed      = f"{RESULTS}/blast/amplicons.bed",
        flagstat = f"{RESULTS}/mapping/flagstat.txt",
    output:
        tex = f"{RESULTS}/tables/on_target_rate.tex",
    log:
        f"{RESULTS}/logs/on_target_rate.log",
    shell:
        """
        mkdir -p $(dirname {output.tex})
        python3 /pipeline/scripts/on_target_rate.py \
            {input.bam} {input.bed} {input.flagstat} {output.tex} 2> {log}
        """


rule task3_flagstat_table:
    """
    Parse samtools flagstat output and render a LaTeX table of key mapping
    metrics (total, mapped, unmapped, properly paired, singletons, duplicates).
    Answers Task 3: percentage of reads mapped to the human genome.
    """
    input:
        f"{RESULTS}/mapping/flagstat.txt",
    output:
        f"{RESULTS}/tables/flagstat_table.tex",
    log:
        f"{RESULTS}/logs/flagstat_table.log",
    shell:
        """
        mkdir -p $(dirname {output})
        python3 /pipeline/scripts/flagstat_table.py {input} {output} 2> {log}
        """


rule task3_unmapped_reads:
    """
    Extract unmapped reads into a separate FASTQ for downstream analysis
    (e.g. BLAST to identify contamination or adapter dimers). Task 3.
    """
    input:
        bam = f"{DATA}/{BAM}",
    output:
        f"{RESULTS}/mapping/unmapped.fastq.gz",
    log:
        f"{RESULTS}/logs/unmapped.log",
    shell:
        """
        samtools view -f 4 -b {input.bam} \
            | samtools fastq - \
            | gzip > {output} \
            2> {log}
        """


rule unmapped_split_reads:
    """
    Split the interleaved unmapped FASTQ (from task3_unmapped_reads) into
    separate R1 and R2 files.  samtools fastq appends /1 or /2 to each
    read name to indicate read-end; this rule splits on that suffix so
    downstream rules can address each end independently.
    """
    input:
        fastq = f"{RESULTS}/mapping/unmapped.fastq.gz",
    output:
        r1 = f"{RESULTS}/mapping/unmapped_R1.fastq.gz",
        r2 = f"{RESULTS}/mapping/unmapped_R2.fastq.gz",
    log:
        f"{RESULTS}/logs/unmapped_split_reads.log",
    shell:
        """
        mkdir -p $(dirname {output.r1})
        python3 /pipeline/scripts/split_unmapped_fastq.py \
            {input.fastq} {output.r1} {output.r2} 2> {log}
        """


rule unmapped_primer_assignment:
    """
    Assign each unmapped read pair to an amplicon by primer matching,
    reusing primer_assignment_table.py. Shows which amplicons contribute
    disproportionately to the unmapped fraction.
    """
    input:
        r1 = f"{RESULTS}/mapping/unmapped_R1.fastq.gz",
        r2 = f"{RESULTS}/mapping/unmapped_R2.fastq.gz",
    output:
        f"{RESULTS}/tables/unmapped_primer_assignment.tex",
    log:
        f"{RESULTS}/logs/unmapped_primer_assignment.log",
    shell:
        """
        mkdir -p $(dirname {output})
        python3 /pipeline/scripts/primer_assignment_table.py \
            {input.r1} {input.r2} {output} \
            --caption "Primer assignment of the 47{{,}}663 unmapped read pairs
(from 97{{,}}021 individual unmapped reads; two reads per pair).
The same primer-matching approach as Table~\\ref{{tab:primer-assignment}}:
first 20\,bp of R1 matched against forward primers (up to 2 mismatches);
R2 against reverse primers if R1 unmatched.
The \\% column shows each amplicon's share of all 47{{,}}663 unmapped read pairs." \
            --label "tab:unmapped-primer-assignment" \
            2> {log}
        """


rule unmapped_read_length_plot:
    """
    Read-length histogram of unmapped reads.
    Shows whether unmapped reads are short (adapter dimers) or full-length,
    helping to identify the dominant cause of mapping failure.
    """
    input:
        f"{RESULTS}/mapping/unmapped.fastq.gz",
    output:
        f"{RESULTS}/figures/unmapped_read_length_plot.pdf",
    log:
        f"{RESULTS}/logs/unmapped_read_length_plot.log",
    shell:
        """
        mkdir -p $(dirname {output})
        python3 /pipeline/scripts/unmapped_read_length_plot.py \
            {input} {output} --read-length 150 \
            2> {log}
        """


rule amplicon_unmapped_rate:
    """
    Per-amplicon unmapped rate table.
    Cross-references primer assignment of all read pairs (to get totals per
    amplicon) with primer assignment of unmapped read pairs, and computes
    the fraction of each amplicon's reads that failed to align to hg19.
    The Unassigned row shows the same metric for reads with no recognisable
    primer, which are expected to have a higher unmapped rate.
    """
    input:
        all_r1      = f"{DATA}/{R1}",
        all_r2      = f"{DATA}/{R2}",
        unmapped_r1 = f"{RESULTS}/mapping/unmapped_R1.fastq.gz",
        unmapped_r2 = f"{RESULTS}/mapping/unmapped_R2.fastq.gz",
    output:
        f"{RESULTS}/tables/amplicon_unmapped_rate.tex",
    log:
        f"{RESULTS}/logs/amplicon_unmapped_rate.log",
    shell:
        """
        mkdir -p $(dirname {output})
        python3 /pipeline/scripts/amplicon_unmapped_rate.py \
            {input.all_r1} {input.all_r2} \
            {input.unmapped_r1} {input.unmapped_r2} \
            {output} 2> {log}
        """
