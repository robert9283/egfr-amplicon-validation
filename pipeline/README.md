# Amplicon Sequencing Analysis Pipeline

Snakemake pipeline for analysing NGS amplicon sequencing data (Illumina MiSeq, paired-end).
Runs inside a single Docker container — no local tool installation required beyond Docker.

## Requirements

- [Docker](https://docs.docker.com/get-docker/)

## Quick start

```bash
bash run.sh --cores 4
```

This will:
1. Build the Docker image (first run only, ~5–10 min)
2. Mount the sequence data and results directory
3. Run the full Snakemake pipeline

## Directory structure

```
pipeline/
├── Dockerfile          # single image with all tools
├── environment.yaml    # pinned conda tool versions
├── config.yaml         # paths, sequences, parameters
├── Snakefile           # top-level workflow
├── run.sh              # build + run wrapper
├── rules/
│   ├── qc.smk          # FastQC, Cutadapt, MultiQC
│   ├── mapping.smk     # samtools flagstat, idxstats, depth
│   └── variants.smk    # mpileup + VarScan2
└── results/            # output directory (created on first run)
```

## Pipeline steps

| Rule | Tool | Output | Task |
|---|---|---|---|
| `task1_fastqc_raw` | FastQC | `qc/raw/` | 1 |
| `task1_cutadapt` | Cutadapt | `trimmed/` | 1 |
| `task1_fastqc_trimmed` | FastQC | `qc/trimmed/` | 1 |
| `task1_multiqc` | MultiQC | `qc/multiqc_report.html` | 1 |
| `task2_idxstats` | samtools | `mapping/idxstats.txt` | 2 |
| `task3_flagstat` | samtools | `mapping/flagstat.txt` | 3 |
| `task3_unmapped_reads` | samtools | `mapping/unmapped.fastq.gz` | 3 |
| `task4_depth` | samtools | `mapping/depth.txt` | 4 |
| `task5_mpileup` | samtools | `variants/mpileup.txt` | 5 |
| `task5_varscan_snp` | VarScan2 | `variants/aln_snps.vcf` | 5 |
| `task5_varscan_indel` | VarScan2 | `variants/aln_indels.vcf` | 5 |

## Useful options

```bash
# Dry run — preview which rules will execute without running them
bash run.sh --cores 4 --dry-run

# Run only a specific rule
bash run.sh --cores 4 --until task3_flagstat

# Force re-run of all rules
bash run.sh --cores 4 --forceall

# Pass any other Snakemake flag the same way
bash run.sh --cores 4 --verbose
```

## Configuration

Edit `config.yaml` to change:
- Input file names
- Adaptor / primer sequences
- Cutadapt quality and length thresholds
- samtools base quality threshold

## Tools and versions

| Tool | Version | Purpose |
|---|---|---|
| FastQC | 0.12.1 | Raw and trimmed read QC |
| Cutadapt | 4.6 | Adaptor trimming |
| MultiQC | 1.19 | QC report aggregation |
| samtools | 1.19.2 | BAM processing, coverage, mapping stats |
| VarScan2 | 2.4.6 | Variant calling |
| Snakemake | 7.32.4 | Workflow management |
