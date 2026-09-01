# syntax=docker/dockerfile:1.7

# Dockerfile for running the PBI Snakemake pipeline
FROM condaforge/mambaforge:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    curl \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy environment file and create conda environment
COPY workflow/envs/base_environment.yaml /tmp/base_environment.yaml

# Install snakemake and core dependencies using conda env
# NOTE: SSL verification is temporarily disabled due to self-signed certificates in
# the GitHub Actions CI environment. In production, ensure proper CA certificates are
# configured or use a trusted mirror. Re-enable with: conda config --set ssl_verify true
# Update the base environment instead of creating a new one
RUN --mount=type=cache,target=/opt/conda/pkgs \
    conda config --set ssl_verify false && \
    conda env update -n base -f /tmp/base_environment.yaml && \
    conda config --set ssl_verify true

# Copy workflow files after dependency installation so source changes do not
# invalidate the cached conda environment layer.
COPY workflow /app/workflow
COPY src /app/src
COPY setup.py /app/

# Normalize source permissions in image to avoid host-umask-dependent failures
# when running container as a non-root UID/GID from docker-compose.
RUN chmod -R a+rX /app/workflow /app/src && chmod a+r /app/setup.py

# Create data directory structure.
# chmod 777 ensures that when Docker seeds the named volumes from these image
# directories on first run, the non-root host UID (set via `user:` in
# docker-compose.yml) can write inside them without permission errors.
RUN mkdir -p /data/processed/databases /data/processed/sequences /data/processed/gff3 /data/processed/blast_db /data/intermediate /data/raw /cache \
    && mkdir -p /pipeline-logs/reports /pipeline-logs/logs \
    && mkdir -p /app/workflow/.snakemake \
    && chmod -R 777 /data /cache /pipeline-logs /app/workflow/.snakemake

# Set the working directory to workflow for Snakemake
WORKDIR /app/workflow

# Create a wrapper entrypoint that patches /etc/passwd for arbitrary UID/GID.
# When docker-compose runs as host UID that has no entry in /etc/passwd,
# getpass.getuser() (called by Snakemake) raises OSError.
# This wrapper adds a temporary passwd entry before exec-ing Snakemake.
RUN printf '#!/bin/sh\n\
    CURRENT_UID=$(id -u)\n\
    CURRENT_GID=$(id -g)\n\
    export USER="${USER:-pipeline}"\n\
    export LOGNAME="${LOGNAME:-${USER}}"\n\
    if ! getent passwd "${CURRENT_UID}" > /dev/null 2>&1 && [ -w /etc/passwd ]; then\n\
    echo "pipeline:x:${CURRENT_UID}:${CURRENT_GID}:Pipeline user:${HOME:-/app}:/bin/sh" >> /etc/passwd 2>/dev/null || true\n\
    fi\n\
    exec "$@"\n' > /entrypoint.sh && chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]

# Default command to run the pipeline 
CMD ["snakemake", "--cores", "all", "--use-conda", "--printshellcmds", "--cache", "--conda-prefix", "/cache/conda", "--directory", "/app/workflow", "--snakefile", "/app/workflow/Snakefile"]
