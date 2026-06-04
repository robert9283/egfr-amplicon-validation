#!/usr/bin/env bash
# -------------------------------------------------------
# Copy auto-generated files from the pipeline results
# into the report folder before compiling.
# Add new entries here as more tables/figures are generated
# by the pipeline.
# -------------------------------------------------------
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_RESULTS="${SCRIPT_DIR}/../pipeline/results"
REPORT_TABLES="${SCRIPT_DIR}/tables"

mkdir -p "${REPORT_TABLES}"

echo "Copying generated files into report/..."

# --- Tables ---
cp "${PIPELINE_RESULTS}/tables/flagstat_table.tex"       "${REPORT_TABLES}/flagstat_table.tex"
echo "  [ok] tables/flagstat_table.tex"

cp "${PIPELINE_RESULTS}/tables/cutadapt_table.tex"       "${REPORT_TABLES}/cutadapt_table.tex"
echo "  [ok] tables/cutadapt_table.tex"

cp "${PIPELINE_RESULTS}/tables/adaptor_check_table.tex"  "${REPORT_TABLES}/adaptor_check_table.tex"
echo "  [ok] tables/adaptor_check_table.tex"

cp "${PIPELINE_RESULTS}/tables/fastqc_table.tex"                 "${REPORT_TABLES}/fastqc_table.tex"
echo "  [ok] tables/fastqc_table.tex"

cp "${PIPELINE_RESULTS}/tables/per_amplicon_readthrough.tex"     "${REPORT_TABLES}/per_amplicon_readthrough.tex"
echo "  [ok] tables/per_amplicon_readthrough.tex"

cp "${PIPELINE_RESULTS}/tables/primer_assignment_table.tex"      "${REPORT_TABLES}/primer_assignment_table.tex"
echo "  [ok] tables/primer_assignment_table.tex"

cp "${PIPELINE_RESULTS}/tables/primer_vs_ref.tex"               "${REPORT_TABLES}/primer_vs_ref.tex"
echo "  [ok] tables/primer_vs_ref.tex"

cp "${PIPELINE_RESULTS}/tables/gene_region_table.tex"            "${REPORT_TABLES}/gene_region_table.tex"
echo "  [ok] tables/gene_region_table.tex"

cp "${PIPELINE_RESULTS}/tables/blast_primer_table.tex"           "${REPORT_TABLES}/blast_primer_table.tex"
echo "  [ok] tables/blast_primer_table.tex"

cp "${PIPELINE_RESULTS}/tables/blast_primer_table_relaxed.tex"   "${REPORT_TABLES}/blast_primer_table_relaxed.tex"
echo "  [ok] tables/blast_primer_table_relaxed.tex"

cp "${PIPELINE_RESULTS}/tables/blast_offtarget_table.tex"        "${REPORT_TABLES}/blast_offtarget_table.tex"
echo "  [ok] tables/blast_offtarget_table.tex"

cp "${PIPELINE_RESULTS}/tables/amplicon_length_table.tex"        "${REPORT_TABLES}/amplicon_length_table.tex"
echo "  [ok] tables/amplicon_length_table.tex"

cp "${PIPELINE_RESULTS}/tables/blast_gene_intersect_table.tex"   "${REPORT_TABLES}/blast_gene_intersect_table.tex"
echo "  [ok] tables/blast_gene_intersect_table.tex"

cp "${PIPELINE_RESULTS}/tables/egfr_exon_table.tex"              "${REPORT_TABLES}/egfr_exon_table.tex"
echo "  [ok] tables/egfr_exon_table.tex"

cp "${PIPELINE_RESULTS}/tables/amplicon_exon_overlap_table.tex"  "${REPORT_TABLES}/amplicon_exon_overlap_table.tex"
echo "  [ok] tables/amplicon_exon_overlap_table.tex"

cp "${PIPELINE_RESULTS}/tables/primer_properties_table.tex"      "${REPORT_TABLES}/primer_properties_table.tex"
echo "  [ok] tables/primer_properties_table.tex"

cp "${PIPELINE_RESULTS}/tables/unmapped_primer_assignment.tex"   "${REPORT_TABLES}/unmapped_primer_assignment.tex"
echo "  [ok] tables/unmapped_primer_assignment.tex"

cp "${PIPELINE_RESULTS}/tables/on_target_rate.tex"               "${REPORT_TABLES}/on_target_rate.tex"
echo "  [ok] tables/on_target_rate.tex"

cp "${PIPELINE_RESULTS}/tables/amplicon_unmapped_rate.tex"       "${REPORT_TABLES}/amplicon_unmapped_rate.tex"
echo "  [ok] tables/amplicon_unmapped_rate.tex"

cp "${PIPELINE_RESULTS}/tables/blast_amp2_summary.tex"           "${REPORT_TABLES}/blast_amp2_summary.tex"
echo "  [ok] tables/blast_amp2_summary.tex"

cp "${PIPELINE_RESULTS}/tables/variant_table.tex"                "${REPORT_TABLES}/variant_table.tex"
echo "  [ok] tables/variant_table.tex"

# --- Figures ---
REPORT_FIGURES="${SCRIPT_DIR}/figures"
mkdir -p "${REPORT_FIGURES}"

cp "${PIPELINE_RESULTS}/figures/egfr_exon_figure.pdf"            "${REPORT_FIGURES}/egfr_exon_figure.pdf"
echo "  [ok] figures/egfr_exon_figure.pdf"

cp "${PIPELINE_RESULTS}/figures/amplicon_depth_plot.pdf"         "${REPORT_FIGURES}/amplicon_depth_plot.pdf"
echo "  [ok] figures/amplicon_depth_plot.pdf"

cp "${PIPELINE_RESULTS}/figures/coverage_uniformity_plot.pdf"    "${REPORT_FIGURES}/coverage_uniformity_plot.pdf"
echo "  [ok] figures/coverage_uniformity_plot.pdf"

cp "${PIPELINE_RESULTS}/figures/insert_size_plot.pdf"           "${REPORT_FIGURES}/insert_size_plot.pdf"
echo "  [ok] figures/insert_size_plot.pdf"

cp "${PIPELINE_RESULTS}/figures/unmapped_read_length_plot.pdf"  "${REPORT_FIGURES}/unmapped_read_length_plot.pdf"
echo "  [ok] figures/unmapped_read_length_plot.pdf"

echo "Done."
