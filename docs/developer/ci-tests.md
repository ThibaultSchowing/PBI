# CI Pipeline Tests

PBI-Scope uses GitHub Actions to automatically build the full Snakemake pipeline and validate the `pbi` Python package against real data on every pull request and push to `main`.

---

## Overview

| Property | Value |
|----------|-------|
| Runner | `ubuntu-latest` |
| Timeout | 120 minutes |
| Triggers | Pull requests to `main`, pushes to `main`, manual dispatch |
| Expected runtime | 20-40 minutes (with Buildx cache) |

The CI pipeline performs two jobs in sequence:

1. **Run the Snakemake pipeline** against a reduced subset of PhageScope data (2 of 26 sources).
2. **Run smoke tests** that exercise the `pbi` Python package against the database and FASTA files produced by the pipeline.

If the pipeline step fails, the smoke tests still run (they use `if: always()`), so you can see whether the package works even when upstream data sources are temporarily unavailable.

---

## Workflow Steps

### 1. Checkout and Buildx Setup

```yaml
- uses: actions/checkout@v4
- uses: docker/setup-buildx-action@v3
```

The repository is checked out and Docker Buildx is initialized. Buildx enables advanced build features, most importantly **layer caching** (see [Docker Buildx Caching](#docker-buildx-caching) below).

### 2. Build Docker Image

```yaml
- uses: docker/build-push-action@v6
  with:
    context: .
    file: workflow/Dockerfile
    tags: pbiscope-ci:latest
    load: true
    cache-from: type=gha,scope=pbiscope-pipeline
    cache-to: type=gha,mode=max,scope=pbiscope-pipeline
```

The Docker image is built from `workflow/Dockerfile`. It installs Snakemake, DuckDB, pyfaidx, and all other dependencies via conda. The image is tagged `pbiscope-ci:latest` and loaded into the local Docker daemon so subsequent steps can use it.

### 3. Run Pipeline (CI Subset)

```bash
docker run --rm \
  -v ${{ github.workspace }}/data:/data \
  snakemake --cores 2 --use-conda \
    --snakefile /app/workflow/Snakefile.ci
```

The Snakemake pipeline runs inside the built container using `Snakefile.ci`, which loads the full default config and then overrides it with `config.ci.yaml`. This produces a DuckDB database, indexed FASTA files, and merged CSV outputs under `/data/processed/`.

### 4. Run Smoke Tests

```bash
docker run --rm \
  -v ${{ github.workspace }}/tests:/app/tests \
  -e DATA_PATH="/data/processed" \
  bash -c "pip install -e /app/. pytest && \
           python -m pytest /app/tests/test_pbi_ci_smoke.py -v"
```

A second container mounts the `tests/` directory and the pipeline output, installs the `pbi` package in editable mode, then runs the smoke test suite with pytest.

---

## CI Subset Configuration

The CI pipeline uses a reduced configuration to keep runtime under 40 minutes. The full production pipeline downloads from 26 sources; CI uses only 2.

| Setting | Production | CI |
|---------|-----------|-----|
| Phage metadata sources | 26 | 2 (RefSeq, PhagesDB) |
| Metadata features | 9 | All 9 included |
| FASTA sources | 26 | 2 (RefSeq, PhagesDB) |
| GFF3 sources | 26 | 2 (RefSeq, PhagesDB) |
| Host genome downloads | Configurable | Enabled (`metadata_only_mode: false`) |

The CI config lives at `workflow/config/config.ci.yaml`. The CI Snakefile (`workflow/Snakefile.ci`) loads the full default config first, then replaces top-level keys with CI overrides.

---

## Smoke Tests

The smoke tests live in `tests/test_pbi_ci_smoke.py` and run inside the Docker container against the real pipeline output. They verify that the `pbi` package works end-to-end without asserting exact data values.

| TestClass | Tests | What It Validates |
|------------|-------|-------------------|
| TestConnection | 2 | Database file exists, `quick_connect()` returns a `SequenceRetriever` |
| TestMetadataQueries | 4 | Phage metadata, protein metadata, phage-host pairs, structured filters all return results |
| TestFastaAccess | 3 | FASTA files exist, `get_phage_sequence()` returns a valid string |
| TestPhageGenomeStreaming | 6 | Genome retrieval in all modes (concat, first, list, dict), gap handling, missing phage error |
| TestHostGenomeStreaming | 8 | Host genome retrieval in all modes, genome stats, gap handling, graceful skip when host data unavailable |
| TestPhageHostPairStreaming | 2 | Phage-host pairs with sequences, concat mode |
| TestErrorHandling | 1 | Nonexistent phage ID returns `None` |

**Total: 26 smoke tests**

The CI also runs `test_private_data_ingestion.py` (17 tests) which validates the private data ingestion pipeline independently of the Snakemake output. These tests use temporary fixtures and do not require the pipeline to have completed.

### Running Tests Locally

After a local pipeline run, you can run both test suites inside the Docker container:

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/tests:/app/tests \
  -e DATA_PATH="/data/processed" \
  pbiscope-ci:latest \
  bash -c "pip install -e /app/. pytest && \
           python -m pytest /app/tests/test_pbi_ci_smoke.py /app/tests/test_private_data_ingestion.py -v"
```

The private data ingestion tests (`test_private_data_ingestion.py`) can also run without the pipeline, since they use temporary fixtures. Run them directly on your host:

```bash
python -m pytest tests/test_private_data_ingestion.py -v
```

---

## Docker Buildx Caching

The CI workflow uses Docker Buildx with GitHub Actions cache to speed up image builds.

### What is Buildx?

Buildx is an extended build backend for Docker that supports advanced features like multi-platform builds, build secrets, and **build cache backends**. The `docker/setup-buildx-action` GitHub Action installs and configures it automatically.

### How the Cache Works

```yaml
cache-from: type=gha,scope=pbiscope-pipeline
cache-to: type=gha,mode=max,scope=pbiscope-pipeline
```

- **`type=gha`**: Uses the GitHub Actions cache backend. Cache layers are stored as GitHub Actions cache entries, scoped to the repository.
- **`scope=pbiscope-pipeline`**: Isolates this workflow's cache from other workflows in the same repository.
- **`mode=max`**: Exports all layers (not just the final image layers) to the cache. This maximizes cache hits for intermediate layers.

### Why It Matters

The Docker image installs conda environments, which can take 5-10 minutes to resolve and download. With Buildx caching:

- **First run**: Full build, layers are pushed to cache (~10-15 min).
- **Subsequent runs** (when only code changes): Docker reuses cached layers for unchanged steps (conda install, system packages). Only the `COPY` layers for modified files are rebuilt (~1-3 min).

The cache is automatically evicted by GitHub after 7 days of inactivity.

---

## Artifacts

After each CI run, two artifacts are uploaded:

| Artifact | Contents | Retention |
|----------|----------|-----------|
| `pipeline-logs` | Snakemake logs, reports, and run metadata from `/pipeline-logs/` | 7 days |
| `merged-data` | Merged CSV outputs from `data/intermediate/csv/merged/` | 7 days |

These are available for download from the GitHub Actions run summary page.

---

## See Also

- [Code Structure](code-structure.md) — project layout and test overview
- `workflow/config/config.ci.yaml` — CI subset configuration
- `tests/test_pbi_ci_smoke.py` — smoke test source code
- `workflow/Snakefile.ci` — CI Snakefile with config override logic
