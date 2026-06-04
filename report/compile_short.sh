#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "Copying generated files..."
bash "${SCRIPT_DIR}/copy_generated_files.sh"

echo "Cleaning auxiliary files..."
rm -f report_long.aux report_long.toc report_long.lof report_long.lot report_long.out report_long.bbl report_long.blg

echo "Pass 1: pdflatex..."
pdflatex -interaction=nonstopmode report_long.tex || true

echo "Running bibtex..."
bibtex report_long || true

echo "Pass 2: pdflatex..."
pdflatex -interaction=nonstopmode report_long.tex || true

echo "Pass 3: pdflatex (resolve references)..."
pdflatex -interaction=nonstopmode report_long.tex || true

# Fail loudly if no PDF was produced
test -f report_long.pdf || { echo "ERROR: report_long.pdf was not produced"; exit 1; }

echo "Done. Output: $SCRIPT_DIR/report_long.pdf"
