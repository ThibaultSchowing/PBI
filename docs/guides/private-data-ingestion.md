# Private Data Ingestion

> **⚠️ Important:** After adding, removing, or modifying private data sources, you **must** delete the manifest file to force the pipeline to re-scan your data:
>
> ```bash
> docker exec <container_name> rm /private-data/private_manifest.json
> ```
>
> Without this step, the pipeline will reuse the cached manifest and your changes will be ignored.

PBI-Scope can ingest private sources from `private_data/` in addition to public PhageScope data.

## Directory structure

Each private source is a subdirectory under the private data root. There are two modes:

### Phage + Host mode (full)

```text
private_data/
  <Source_DB>/
    metadata.csv
    phage.fasta
    hosts/
      <Host_ID>.fna
```

### Phage-only mode

```text
private_data/
  <Source_DB>/
    metadata.csv
    phage.fasta
```

When you do not have host genomes, use **phage-only mode**. Set `Host_ID` and `Host_name` to `unknown` in `metadata.csv` and omit the `hosts/` directory entirely.

## metadata.csv format

### Required columns

| Column | Description |
|---|---|
| `Phage_ID` | Unique phage identifier (must exist in `phage.fasta`) |
| `Host_ID` | Host identifier (`unknown` when no host is available) |
| `Host_name` | Host species name (`unknown` when no host is available) |
| `Source_DB` | Must match the source directory name exactly |
| `interaction` | Either `temperate` or `virulent` |

Any additional columns are stored as entity attributes.

### Example: Phage + Host

```csv
Phage_ID,Host_ID,Host_name,Source_DB,interaction
MyPhage_1,Escherichia_coli_K12,Escherichia coli,MyPrivateDB,virulent
MyPhage_2,Escherichia_coli_K12,Escherichia coli,MyPrivateDB,temperate
```

### Example: Phage-only

```csv
Phage_ID,Host_ID,Host_name,Source_DB,interaction
MyPhage_1,unknown,unknown,MyPrivateDB,virulent
MyPhage_2,unknown,unknown,MyPrivateDB,temperate
```

## Mandatory rules

- `metadata.csv` is required
- `phage.fasta` is required
- `hosts/` directory is optional — required only when `Host_ID` values are real identifiers
- When `Host_ID` and `Host_name` are both `unknown`, no `hosts/` directory is needed
- Every `Phage_ID` in metadata must exist in `phage.fasta`

## Validate before pipeline

```bash
pbi validate-private
```

## Runtime behavior

- Valid private sources are ingested and linked with `source_type=private`
- Invalid sources are skipped (public pipeline still completes)
- Re-running pipeline synchronizes removals/additions
- `Source_DB` in `metadata.csv` must match the source folder name exactly
- Phage-only sources (no host genomes) are ingested without host mapping entries in `dim_hosts`

## Validate what was ingested

Use DuckDB (or `SequenceRetriever`) to inspect available source labels:

```sql
SELECT Source_DB, source_type, COUNT(*) AS phage_count
FROM fact_phages
GROUP BY Source_DB, source_type
ORDER BY source_type, Source_DB;
```

If you filter `Source_DB = 'test_private'` and get `0`, first check this query to confirm the exact source name currently present (for example `test_private_2`).

## Output mappings

- `private_phage_mapping.json` routes private phage retrieval
- `host_fasta_mapping.json` includes host paths (public + private)

## Logs

In Docker runs, logs/reports are available in `./pipeline_logs/`.

Private-source validation details are written to:

- `private_data/private_manifest.json` (host path)
- `/private-data/private_manifest.json` (inside container)

This manifest explicitly lists:

- `is_valid` per source
- validation `errors`
- skipped/ingested source counts

For provenance/version-pinning details and public-source diagnostics, see:

- [How It Works](how-it-works.md#data-provenance)
