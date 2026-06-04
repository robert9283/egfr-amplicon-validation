#!/usr/bin/env bash
# -------------------------------------------------------
# Build the Docker image and run the Snakemake pipeline.
# Usage: bash run.sh [snakemake options]
# Example: bash run.sh --cores 4
#          bash run.sh --cores 4 --dry-run
# -------------------------------------------------------
set -e

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$(cd "${PIPELINE_DIR}/../sequenceFile" && pwd)"
RESULTS_DIR="${PIPELINE_DIR}/results"
RESOURCES_DIR="${PIPELINE_DIR}/resources"
IMAGE_NAME="amplicon-pipeline"

mkdir -p "${RESULTS_DIR}" "${RESOURCES_DIR}"

echo "Building Docker image: ${IMAGE_NAME}..."
docker build -t "${IMAGE_NAME}" "${PIPELINE_DIR}"

echo "Running pipeline..."
echo "  Data:    ${DATA_DIR}"
echo "  Results: ${RESULTS_DIR}"
echo ""

docker run --rm \
    -v "${DATA_DIR}:/data:ro" \
    -v "${RESULTS_DIR}:/results" \
    -v "${PIPELINE_DIR}:/pipeline:ro" \
    -v "${RESOURCES_DIR}:/pipeline/resources" \
    "${IMAGE_NAME}" \
    --snakefile /pipeline/Snakefile \
    --configfile /pipeline/config.yaml \
    --directory /results \
    --cores "${CORES:-4}" \
    "$@"
