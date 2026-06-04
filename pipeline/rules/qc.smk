# -------------------------------------------------------
# QC rules: FastQC on raw reads, Cutadapt trimming,
# FastQC on trimmed reads, insert-size plot, primer vs ref
# (DATA, RESULTS, SAMPLE, R1, R2, BAM defined in Snakefile)
# -------------------------------------------------------


rule task1_fastqc_raw:
    """Run FastQC on raw paired-end FASTQ files."""
    input:
        r1 = f"{DATA}/{R1}",
        r2 = f"{DATA}/{R2}",
    output:
        html_r1 = f"{RESULTS}/qc/raw/{SAMPLE}_R1_fastqc.html",
        html_r2 = f"{RESULTS}/qc/raw/{SAMPLE}_R2_fastqc.html",
        zip_r1  = f"{RESULTS}/qc/raw/{SAMPLE}_R1_fastqc.zip",
        zip_r2  = f"{RESULTS}/qc/raw/{SAMPLE}_R2_fastqc.zip",
    params:
        outdir = f"{RESULTS}/qc/raw",
    log:
        f"{RESULTS}/logs/fastqc_raw.log",
    shell:
        """
        mkdir -p {params.outdir}
        fastqc --outdir {params.outdir} {input.r1} {input.r2} 2> {log}
        # Rename outputs to match expected names
        mv {params.outdir}/$(basename {input.r1} .fastq.gz)_fastqc.html {output.html_r1}
        mv {params.outdir}/$(basename {input.r2} .fastq.gz)_fastqc.html {output.html_r2}
        mv {params.outdir}/$(basename {input.r1} .fastq.gz)_fastqc.zip  {output.zip_r1}
        mv {params.outdir}/$(basename {input.r2} .fastq.gz)_fastqc.zip  {output.zip_r2}
        """


rule task1_cutadapt:
    """
    Trim Illumina adaptors from paired-end reads.
    -a / -A: adaptor sequences on R1 and R2 respectively.
    Also trims low-quality 3' bases and discards very short reads.
    """
    input:
        r1 = f"{DATA}/{R1}",
        r2 = f"{DATA}/{R2}",
    output:
        r1 = f"{RESULTS}/trimmed/{SAMPLE}_R1_trimmed.fastq.gz",
        r2 = f"{RESULTS}/trimmed/{SAMPLE}_R2_trimmed.fastq.gz",
        report = f"{RESULTS}/trimmed/{SAMPLE}_cutadapt_report.txt",
    params:
        trim_r1       = config["trim_r1"],
        trim_r2       = config["trim_r2"],
        min_length    = config["cutadapt"]["min_length"],
        quality_cutoff = config["cutadapt"]["quality_cutoff"],
    log:
        f"{RESULTS}/logs/cutadapt.log",
    shell:
        """
        mkdir -p $(dirname {output.r1})
        cutadapt \
            -a {params.trim_r1} \
            -A {params.trim_r2} \
            --minimum-length {params.min_length} \
            --quality-cutoff {params.quality_cutoff} \
            --json {output.report} \
            -o {output.r1} \
            -p {output.r2} \
            {input.r1} {input.r2} \
            2> {log}
        """


rule task1_fastqc_trimmed:
    """Run FastQC on adapter-trimmed reads to confirm clean output."""
    input:
        r1 = rules.task1_cutadapt.output.r1,
        r2 = rules.task1_cutadapt.output.r2,
    output:
        html_r1 = f"{RESULTS}/qc/trimmed/{SAMPLE}_R1_trimmed_fastqc.html",
        html_r2 = f"{RESULTS}/qc/trimmed/{SAMPLE}_R2_trimmed_fastqc.html",
        zip_r1  = f"{RESULTS}/qc/trimmed/{SAMPLE}_R1_trimmed_fastqc.zip",
        zip_r2  = f"{RESULTS}/qc/trimmed/{SAMPLE}_R2_trimmed_fastqc.zip",
    params:
        outdir = f"{RESULTS}/qc/trimmed",
    log:
        f"{RESULTS}/logs/fastqc_trimmed.log",
    shell:
        """
        mkdir -p {params.outdir}
        fastqc --outdir {params.outdir} {input.r1} {input.r2} 2> {log}
        """




rule task1_check_adaptors:
    """
    Verify which adaptor sequences are present in R1 and R2 reads by searching
    for both adaptors and their reverse complements. Determines whether
    Scenario A (same adaptor both ends) or Scenario B (standard paired-end)
    applies to this library.
    """
    input:
        r1 = f"{DATA}/{R1}",
        r2 = f"{DATA}/{R2}",
    output:
        f"{RESULTS}/qc/adaptor_check.txt",
    log:
        f"{RESULTS}/logs/adaptor_check.log",
    shell:
        """
        python3 /pipeline/scripts/check_adaptors.py {input.r1} {input.r2} \
            > {output} 2> {log}
        """


rule task1_adaptor_check_table:
    """Generate a LaTeX table from the adaptor check results."""
    input:
        rules.task1_check_adaptors.output,
    output:
        f"{RESULTS}/tables/adaptor_check_table.tex",
    log:
        f"{RESULTS}/logs/adaptor_check_table.log",
    shell:
        """
        mkdir -p $(dirname {output})
        python3 /pipeline/scripts/adaptor_check_table.py {input} {output} 2> {log}
        """


rule task1_per_amplicon_readthrough:
    """Per-amplicon read-through statistics: assigns reads to amplicons by primer
    matching and checks for far-adapter RC (read-through signal)."""
    input:
        r1 = f"{DATA}/{R1}",
        r2 = f"{DATA}/{R2}",
    output:
        f"{RESULTS}/tables/per_amplicon_readthrough.tex",
    log:
        f"{RESULTS}/logs/per_amplicon_readthrough.log",
    shell:
        """
        mkdir -p $(dirname {output})
        python3 /pipeline/scripts/per_amplicon_readthrough.py \
            {input.r1} {input.r2} {output} 2> {log}
        """


rule task1_primer_assignment_table:
    """
    Primer assignment table: assigns each read pair to one of the four amplicons
    by primer matching and reports the assignment rate per amplicon.
    Confirms that primer sequences are present at the expected 5' position of reads.
    """
    input:
        r1 = f"{DATA}/{R1}",
        r2 = f"{DATA}/{R2}",
    output:
        f"{RESULTS}/tables/primer_assignment_table.tex",
    log:
        f"{RESULTS}/logs/primer_assignment_table.log",
    shell:
        """
        mkdir -p $(dirname {output})
        python3 /pipeline/scripts/primer_assignment_table.py \
            {input.r1} {input.r2} {output} 2> {log}
        """


rule task1_fastqc_table:
    """Generate a LaTeX table from FastQC results for raw and trimmed reads."""
    input:
        raw_r1  = rules.task1_fastqc_raw.output.zip_r1,
        raw_r2  = rules.task1_fastqc_raw.output.zip_r2,
        trim_r1 = rules.task1_fastqc_trimmed.output.zip_r1,
        trim_r2 = rules.task1_fastqc_trimmed.output.zip_r2,
    output:
        f"{RESULTS}/tables/fastqc_table.tex",
    log:
        f"{RESULTS}/logs/fastqc_table.log",
    shell:
        """
        mkdir -p $(dirname {output})
        python3 /pipeline/scripts/fastqc_table.py \
            {input.raw_r1} {input.raw_r2} {input.trim_r1} {input.trim_r2} \
            {output} 2> {log}
        """


rule task1_cutadapt_table:
    """Generate a LaTeX table from the Cutadapt JSON report."""
    input:
        rules.task1_cutadapt.output.report,
    output:
        f"{RESULTS}/tables/cutadapt_table.tex",
    log:
        f"{RESULTS}/logs/cutadapt_table.log",
    shell:
        """
        mkdir -p $(dirname {output})
        python3 /pipeline/scripts/cutadapt_table.py {input} {output} 2> {log}
        """


rule insert_size_plot:
    """
    Per-amplicon insert size distribution.

    For each EGFR amplicon, extracts the TLEN (template length) field from
    properly paired mapped reads via samtools view and plots a histogram.
    A vertical line marks the read length (150 bp) to indicate which amplicons
    are shorter than the read and therefore generate read-through artefacts.
    """
    input:
        bam = f"{DATA}/{BAM}",
        bai = f"{DATA}/{BAM}.bai",
        bed = f"{RESULTS}/blast/amplicons.bed",
    output:
        f"{RESULTS}/figures/insert_size_plot.pdf",
    log:
        f"{RESULTS}/logs/insert_size_plot.log",
    shell:
        """
        mkdir -p $(dirname {output})
        python3 /pipeline/scripts/insert_size_plot.py \
            {input.bam} {input.bed} {output} \
            --read-length 150 \
            2> {log}
        """


rule primer_vs_ref:
    """Compare each primer sequence to the hg19 reference at its BLAST binding site.

    Extracts the reference sequence at each primer's BLAST-determined coordinates
    and compares it base-by-base to the primer sequence (reverse-complementing
    minus-strand hits before comparison). Outputs a LaTeX table of matches and
    any mismatches, which is used to confirm whether primer-end artefacts are
    caused by primer/reference divergence or by other mechanisms (e.g. polymerase
    slippage at the primer-insert junction).
    """
    input:
        primers  = f"{RESULTS}/blast/primers.fa",
        blast    = f"{RESULTS}/blast/primers_vs_hg19.txt",
        ref      = f"{RESOURCES}/hg19_chr7.fa",
    output:
        f"{RESULTS}/tables/primer_vs_ref.tex",
    log:
        f"{RESULTS}/logs/primer_vs_ref.log",
    shell:
        """
        mkdir -p $(dirname {output})
        python3 /pipeline/scripts/primer_vs_ref.py \
            {input.primers} {input.blast} {input.ref} {output} 2> {log}
        """
