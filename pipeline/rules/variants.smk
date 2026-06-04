# -------------------------------------------------------
# Variant calling rules: Mutect2 (GATK4) + VarDictJava
# Answers Task 5: identify variants and evaluate quality.
#
# Reference note: the BAM uses NCBI chromosome names ("7")
# while the downloaded hg19.fa.gz uses UCSC names ("chr7").
# prepare_reference extracts chr7 and strips the prefix so
# the reference matches the BAM header exactly.
# -------------------------------------------------------

# (DATA, RESULTS, RESOURCES, SAMPLE, BAM defined in Snakefile)

HG19_FA_GZ   = f"{RESOURCES}/hg19.fa.gz"
HG19_CHR7_FA = f"{RESOURCES}/hg19_chr7.fa"


rule prepare_reference:
    """
    Extract chr7 from the downloaded hg19 FASTA, strip the 'chr' prefix
    to match the BAM's NCBI-style chromosome names, then index with
    samtools faidx and create a GATK sequence dictionary.

    Only chr7 is needed since all four EGFR amplicons lie on chr7.
    Extracting a single chromosome keeps disk usage and runtime low
    compared to decompressing and renaming the full ~3 GB reference.
    """
    input:
        HG19_FA_GZ,
    output:
        fa   = HG19_CHR7_FA,
        fai  = HG19_CHR7_FA + ".fai",
        dic  = HG19_CHR7_FA.replace(".fa", ".dict"),
    log:
        f"{RESULTS}/logs/prepare_reference.log",
    shell:
        """
        # hg19.fa.gz is standard gzip (not bgzip), so samtools faidx
        # cannot do random access on it.  Decompress the full file and
        # extract chr7 with awk, then strip the chr prefix.
        gunzip -c {input} \
            | awk 'BEGIN{{f=0}} /^>/{{f=($0==">chr7")}} f' \
            | sed 's/^>chr7/>7/' \
            > {output.fa} 2> {log}
        samtools faidx {output.fa}           2>> {log}
        gatk CreateSequenceDictionary \
            -R {output.fa}                   2>> {log}
        """


rule mutect2_call:
    """
    Call somatic SNVs and indels with GATK Mutect2 in tumour-only mode.

    Flags:
      -L amplicons.bed  : restrict to the four EGFR amplicon spans
      --tumor-sample    : sample name read from the BAM read-group (SM tag)
      --dont-use-soft-clipped-bases : ignore soft-clipped bases at primer
                          boundaries to reduce false positives at amplicon ends

    The chr prefix is stripped from amplicons.bed to match the BAM and
    reference naming convention (process substitution via bash -c).
    """
    input:
        bam  = f"{DATA}/{BAM}",
        bai  = f"{DATA}/{BAM}.bai",
        fa   = HG19_CHR7_FA,
        fai  = HG19_CHR7_FA + ".fai",
        dic  = HG19_CHR7_FA.replace(".fa", ".dict"),
        bed  = f"{RESULTS}/blast/amplicons.bed",
    output:
        vcf  = f"{RESULTS}/variants/mutect2_raw.vcf.gz",
        tbi  = f"{RESULTS}/variants/mutect2_raw.vcf.gz.tbi",
        stat = f"{RESULTS}/variants/mutect2_raw.vcf.gz.stats",
    log:
        f"{RESULTS}/logs/mutect2_call.log",
    shell:
        """
        mkdir -p $(dirname {output.vcf})

        # Get sample name from BAM read group
        sample=$(samtools view -H {input.bam} \
                    | grep '^@RG' \
                    | grep -oP 'SM:\\K[^\\t]+' \
                    | head -1)
        [ -z "$sample" ] && sample="{SAMPLE}"

        # GATK -L does not support process substitution; write to a temp file
        bed_nochr=$(mktemp /tmp/amplicons_nochr_XXXXXX.bed)
        sed 's/^chr//' {input.bed} > "$bed_nochr"

        gatk Mutect2 \
            -R {input.fa} \
            -I {input.bam} \
            -O {output.vcf} \
            -L "$bed_nochr" \
            --tumor-sample "$sample" \
            --dont-use-soft-clipped-bases \
            2>> {log}

        rm -f "$bed_nochr"
        """


rule mutect2_filter:
    """
    Apply GATK's learned somatic filtering model to the raw Mutect2 calls.
    FilterMutectCalls uses orientation bias, strand artefact, and contamination
    models to assign PASS/FAIL to each variant.  PASS variants are retained
    for downstream analysis.
    """
    input:
        vcf  = f"{RESULTS}/variants/mutect2_raw.vcf.gz",
        stat = f"{RESULTS}/variants/mutect2_raw.vcf.gz.stats",
        fa   = HG19_CHR7_FA,
    output:
        vcf  = f"{RESULTS}/variants/mutect2_filtered.vcf.gz",
        tbi  = f"{RESULTS}/variants/mutect2_filtered.vcf.gz.tbi",
    log:
        f"{RESULTS}/logs/mutect2_filter.log",
    shell:
        """
        gatk FilterMutectCalls \
            -R {input.fa} \
            -V {input.vcf} \
            -O {output.vcf} \
            2> {log}
        """


rule vardict_call:
    """
    Call somatic variants with VarDictJava, optimised for amplicon sequencing.

    Pipeline:
      1. vardict-java  : call variants in each amplicon interval
      2. teststrandbias.R : apply strand-bias Fisher test
      3. var2vcf_valid.pl : convert to VCF, keeping only PASS variants

    Key parameters:
      -f 0.01  : minimum allele frequency 1%
      -c 1 -S 2 -E 3 -g 4 : BED column indices (chr, start, end, gene)
      -q 20    : minimum base quality
      -Q 30    : minimum mapping quality
      -N       : sample name

    The chr prefix is stripped from amplicons.bed (process substitution)
    to match the BAM's NCBI chromosome naming.
    """
    input:
        bam = f"{DATA}/{BAM}",
        bai = f"{DATA}/{BAM}.bai",
        fa  = HG19_CHR7_FA,
        fai = HG19_CHR7_FA + ".fai",
        bed = f"{RESULTS}/blast/amplicons.bed",
    output:
        vcf = f"{RESULTS}/variants/vardict.vcf",
    log:
        f"{RESULTS}/logs/vardict_call.log",
    shell:
        """
        mkdir -p $(dirname {output.vcf})

        vardict-java \
            -G {input.fa} \
            -f 0.01 \
            -N {SAMPLE} \
            -b {input.bam} \
            -q 20 \
            -Q 30 \
            -c 1 -S 2 -E 3 -g 4 \
            <(sed 's/^chr//' {input.bed}) \
            2>> {log} \
        | teststrandbias.R 2>> {log} \
        | var2vcf_valid.pl \
            -N {SAMPLE} \
            -E \
            -f 0.01 \
            > {output.vcf} 2>> {log}
        """


rule download_snpeff_db:
    """
    Download the SnpEff hg19 database (UCSC chromosome names, chr-prefixed)
    into the resources directory.  Runs once; Snakemake skips if the sentinel
    file already exists.
    """
    output:
        sentinel = f"{RESOURCES}/snpeff_data/hg19/snpEffectPredictor.bin",
    log:
        f"{RESULTS}/logs/download_snpeff_db.log",
    shell:
        """
        mkdir -p {RESOURCES}/snpeff_data
        snpEff download \
            -dataDir {RESOURCES}/snpeff_data \
            hg19 \
            2> {log}
        """


rule snpeff_annotate_mutect2:
    """
    Annotate Mutect2 PASS variants with SnpEff (hg19 database).
    Adds ANN INFO field with gene, effect, and HGVSp (protein change).
    Mutect2 VCF already uses chr-prefixed chromosome names.
    """
    input:
        vcf      = f"{RESULTS}/variants/mutect2_filtered.vcf.gz",
        sentinel = f"{RESOURCES}/snpeff_data/hg19/snpEffectPredictor.bin",
    output:
        vcf = f"{RESULTS}/variants/mutect2_annotated.vcf.gz",
        tbi = f"{RESULTS}/variants/mutect2_annotated.vcf.gz.tbi",
    log:
        f"{RESULTS}/logs/snpeff_annotate_mutect2.log",
    shell:
        """
        snpEff \
            -dataDir {RESOURCES}/snpeff_data \
            -noStats \
            hg19 \
            {input.vcf} \
            2> {log} \
        | bgzip > {output.vcf}
        tabix {output.vcf}
        """


rule snpeff_annotate_vardict:
    """
    Annotate VarDictJava PASS variants with SnpEff (hg19 database).
    VarDictJava strips the chr prefix; add it back before annotation
    so chromosome names match the hg19 SnpEff database.
    """
    input:
        vcf      = f"{RESULTS}/variants/vardict.vcf",
        sentinel = f"{RESOURCES}/snpeff_data/hg19/snpEffectPredictor.bin",
    output:
        vcf = f"{RESULTS}/variants/vardict_annotated.vcf",
    log:
        f"{RESULTS}/logs/snpeff_annotate_vardict.log",
    shell:
        """
        sed '/^#/!s/^/chr/' {input.vcf} \
        | snpEff \
            -dataDir {RESOURCES}/snpeff_data \
            -noStats \
            hg19 \
            - \
            2> {log} \
        > {output.vcf}
        """


rule variant_table:
    """
    Parse PASS variants from SnpEff-annotated Mutect2 and VarDictJava VCFs,
    annotate each with amplicon, exon, and protein change (from SnpEff ANN
    field), and render a LaTeX summary table.
    """
    input:
        mutect2  = f"{RESULTS}/variants/mutect2_annotated.vcf.gz",
        vardict  = f"{RESULTS}/variants/vardict_annotated.vcf",
        bed      = f"{RESULTS}/blast/amplicons.bed",
        exon_tex = f"{RESULTS}/tables/amplicon_exon_overlap_table.tex",
    output:
        f"{RESULTS}/tables/variant_table.tex",
    log:
        f"{RESULTS}/logs/variant_table.log",
    shell:
        """
        mkdir -p $(dirname {output})
        python3 /pipeline/scripts/variant_table.py \
            {input.mutect2} {input.vardict} \
            {input.bed} {input.exon_tex} \
            {output} 2> {log}
        """


# -------------------------------------------------------
# Legacy VarScan2 rules (kept for reference, not used in
# the primary variant calling workflow)
# -------------------------------------------------------

rule task5_mpileup:
    """
    Generate a pileup of base calls at each position.
    Used as input to VarScan2 for variant calling.
    """
    input:
        bam = f"{DATA}/{BAM}",
    output:
        f"{RESULTS}/variants/mpileup.txt",
    log:
        f"{RESULTS}/logs/mpileup.log",
    shell:
        """
        mkdir -p $(dirname {output})
        samtools mpileup \
            -B \
            -q 20 \
            -Q 20 \
            {input.bam} \
            > {output} \
            2> {log}
        """


rule task5_varscan_snp:
    """Call SNPs with VarScan2 (legacy; use mutect2_filter instead)."""
    input:
        pileup = rules.task5_mpileup.output,
    output:
        f"{RESULTS}/variants/{SAMPLE}_snps.vcf",
    log:
        f"{RESULTS}/logs/varscan_snp.log",
    shell:
        """
        varscan mpileup2snp {input.pileup} \
            --min-coverage 10 \
            --min-var-freq 0.01 \
            --min-avg-qual 20 \
            --p-value 0.05 \
            --output-vcf 1 \
            > {output} \
            2> {log}
        """


rule task5_varscan_indel:
    """Call small indels with VarScan2 (legacy; use mutect2_filter instead)."""
    input:
        pileup = rules.task5_mpileup.output,
    output:
        f"{RESULTS}/variants/{SAMPLE}_indels.vcf",
    log:
        f"{RESULTS}/logs/varscan_indel.log",
    shell:
        """
        varscan mpileup2indel {input.pileup} \
            --min-coverage 10 \
            --min-var-freq 0.01 \
            --min-avg-qual 20 \
            --p-value 0.05 \
            --output-vcf 1 \
            > {output} \
            2> {log}
        """
