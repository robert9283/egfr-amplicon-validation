FROM continuumio/miniconda3:24.1.2-0

SHELL ["/bin/bash", "-c"]

# Install mamba for faster dependency resolution
RUN conda install -n base -c conda-forge mamba -y && \
    conda clean -afy

# Install all bioinformatics tools via conda environment file
COPY environment.yaml /tmp/environment.yaml
RUN mamba env update -n base -f /tmp/environment.yaml && \
    conda clean -afy

# Install build tools needed to compile cutadapt on linux-aarch64, plus curl and BLAST
# (bioconda blast package has no linux-aarch64 build; ncbi-blast+ from Debian apt works on ARM)
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential curl ncbi-blast+ && \
    rm -rf /var/lib/apt/lists/*

# Install cutadapt via pip (bioconda conda package unavailable for linux-aarch64)
# Install Python visualisation packages separately so this layer can be updated
# without invalidating the heavy conda layer above.
RUN pip install cutadapt==4.6 matplotlib==3.9 numpy

WORKDIR /pipeline

# Default: run snakemake with the pipeline config
ENTRYPOINT ["snakemake"]
CMD ["--help"]
