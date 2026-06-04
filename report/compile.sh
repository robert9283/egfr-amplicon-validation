#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "Copying generated files..."
bash "${SCRIPT_DIR}/copy_generated_files.sh"

echo "Cleaning auxiliary files..."
rm -f report.aux report.toc report.lof report.lot report.out report.bbl report.blg

echo "Pass 1: pdflatex..."
pdflatex -interaction=nonstopmode report.tex || true

echo "Running bibtex..."
bibtex report || true

echo "Pass 2: pdflatex..."
pdflatex -interaction=nonstopmode report.tex || true

echo "Pass 3: pdflatex (resolve references)..."
pdflatex -interaction=nonstopmode report.tex || true

# Fail loudly if no PDF was produced
test -f report.pdf || { echo "ERROR: report.pdf was not produced"; exit 1; }

echo "Done. Output: $SCRIPT_DIR/report.pdf"
