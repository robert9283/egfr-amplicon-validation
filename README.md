# EGFR Amplicon Sequencing Validation

A reproducible bioinformatics pipeline for validating a targeted amplicon sequencing panel
designed to detect clinically actionable mutations in the *EGFR* gene relevant to
non-small-cell lung cancer (NSCLC).

## Overview

The panel covers four amplicons spanning exons 18–21, including the exon 19 in-frame
deletions (e.g. E746\_A750del) and the exon 21 L858R substitution that confer sensitivity
to tyrosine kinase inhibitors (TKIs) such as gefitinib, erlotinib, and osimertinib.

The analysis addresses five validation questions:

1. Read constitution and adaptor content
2. Target gene identification via BLAST-based primer alignment
3. Overall mapping rate and read quality
4. Per-amplicon coverage depth and uniformity
5. Somatic variant calling with orthogonal callers (Mutect2, VarDictJava)

**Key findings:**
- All four amplicons map on-target to EGFR exon sequences on chromosome 7
- Mapping rates exceed 99% overall
- Three amplicons show high, uniform coverage (>1,000× per base); Amplicon 2 shows substantially lower coverage (~9,000 vs ~200,000 read pairs), with 89% of its read pairs unmapped — BLAST identifies 87% of these as phiX174 bacteriophage (Illumina spike-in), indicating the Amplicon 2 forward primer cross-reacts with the phage genome during PCR
- One confirmed somatic EGFR driver mutation: exon 19 in-frame deletion E746\_A750del (VAF 48%, COSM6223)
- T790M and L858R are absent, confirmed by >200,000× coverage at both hotspot positions
- Two Mutect2 PASS calls at amplicon boundaries are primer-end artefacts caused by missing primer trimming (`ivar trim`)

The full analysis report is available at [`report/report.pdf`](report/report.pdf).

## Requirements

- [Docker](https://docs.docker.com/get-docker/)
- Input: paired-end FASTQ files and a pre-aligned BAM file

## Reproducing the analysis

**1. Clone the repository**

```bash
git clone git@github.com:robert9283/egfr-amplicon-validation.git
cd egfr-amplicon-validation
```

**2. Place your input data**

Put your FASTQ and BAM files in `sequenceFile/` and update `pipeline/config.yaml` if needed:

```bash
mkdir sequenceFile
cp /path/to/aln1.fastq.gz /path/to/aln2.fastq.gz /path/to/aln.bam sequenceFile/
```

**3. Build the Docker image**

```bash
docker build -t amplicon-pipeline pipeline/
```

This installs all tools via conda inside the container (~5–10 min on first build).

**4. Run the pipeline**

```bash
bash pipeline/run.sh --cores 4
```

This mounts the input data, results, pipeline code, and resources into the container
and runs the full Snakemake workflow. All outputs are written to `pipeline/results/`.

Alternatively, run Docker directly:

```bash
docker run --rm \
    -v "$(pwd)/sequenceFile:/data:ro" \
    -v "$(pwd)/pipeline/results:/results" \
    -v "$(pwd)/pipeline:/pipeline:ro" \
    -v "$(pwd)/pipeline/resources:/pipeline/resources" \
    amplicon-pipeline \
    --snakefile /pipeline/Snakefile \
    --configfile /pipeline/config.yaml \
    --directory /results \
    --cores 4
```

**5. Compile the report**

```bash
bash report/compile.sh
```

Copies pipeline outputs into `report/tables/` and `report/figures/`, then compiles
`report/report.tex` to `report/report.pdf` using pdflatex + bibtex.

**Useful Snakemake options**

```bash
# Preview which rules will run without executing them
bash pipeline/run.sh --cores 4 --dry-run

# Run up to a specific rule
bash pipeline/run.sh --cores 4 --until task3_flagstat

# Force re-run of all rules
bash pipeline/run.sh --cores 4 --forceall
```

## Repository structure

```
├── pipeline/
│   ├── Dockerfile          # single image with all tools
│   ├── environment.yaml    # pinned conda tool versions
│   ├── config.yaml         # paths, sequences, parameters
│   ├── Snakefile           # top-level workflow
│   ├── run.sh              # Docker build + run wrapper
│   ├── rules/
│   │   ├── qc.smk          # FastQC, Cutadapt, adaptor checking, insert size
│   │   ├── mapping.smk     # samtools flagstat/idxstats/depth, coverage plots
│   │   ├── variants.smk    # Mutect2, VarDictJava, VarScan2, SnpEff annotation
│   │   └── blast.smk       # BLAST primer alignment, gene/exon intersection
│   └── scripts/            # Python scripts for table and figure generation
└── report/
    ├── report.tex           # main LaTeX report
    ├── compile.sh           # report compilation script
    ├── figures/             # TikZ figures
    └── tables/              # auto-generated LaTeX tables
```

## Tools

| Tool | Version | Purpose |
|---|---|---|
| FastQC | 0.12.1 | Raw and trimmed read QC |
| Cutadapt | 4.6 | Adaptor trimming |
| samtools | 1.19.2 | BAM processing, mapping stats, depth |
| bowtie2 | 2.5.3 | Read alignment |
| BLAST+ | 2.14.1 | Primer-to-genome alignment |
| bedtools | 2.31.1 | Genomic interval operations |
| Mutect2 | 4.4.0.0 | Somatic variant calling |
| VarDictJava | 1.8.3 | Somatic variant calling |
| VarScan2 | 2.4.6 | Somatic variant calling |
| SnpEff | 5.1 | Variant annotation |
| Snakemake | 7.32.4 | Workflow management |
