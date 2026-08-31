from __future__ import annotations

import hashlib
import json
import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Iterable, List, Optional, Set, Tuple

import pandas as pd

if TYPE_CHECKING:
    import duckdb


MANDATORY_COLUMNS = ["Phage_ID", "Host_ID", "Host_name", "Source_DB", "interaction"]
MAX_ERROR_EXAMPLES = 20
# Extension precedence for per-host files. ``.fna`` comes first because it is
# the canonical extension emitted by the host-download workflow; when multiple
# files exist for the same Host_ID, the first extension in this tuple wins.
HOST_FASTA_EXTENSIONS = (".fna", ".fasta", ".fa")

logger = logging.getLogger(__name__)


@dataclass
class PrivateSourceValidation:
    source_dir: Path
    source_db: str
    errors: List[str]
    warnings: List[str]
    stats: Dict[str, int]
    metadata_df: Optional[pd.DataFrame] = None
    extra_columns: Optional[List[str]] = None

    @property
    def is_valid(self) -> bool:
        return len(self.errors) == 0


def _extract_fasta_identifier(header_line: str) -> str:
    if not header_line.startswith(">"):
        raise ValueError("FASTA header must start with '>'")
    header = header_line[1:].strip()
    # Empty identifiers are ignored by callers to keep error messages actionable elsewhere.
    return header.split()[0] if header else ""


def parse_fasta_ids(fasta_path: Path) -> tuple[Set[str], Set[str]]:
    """
    Parse FASTA identifiers using the first token of each header line.

    Args:
        fasta_path: Path to FASTA file.

    Returns:
        Tuple (all_ids, duplicate_ids), where both items are sets of IDs.
    """
    ids: Set[str] = set()
    duplicates: Set[str] = set()

    with fasta_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            seq_id = _extract_fasta_identifier(line)
            if not seq_id:
                continue
            if seq_id in ids:
                duplicates.add(seq_id)
            ids.add(seq_id)

    return ids, duplicates


def _required_files(source_dir: Path) -> Dict[str, Path]:
    return {
        "metadata.csv": source_dir / "metadata.csv",
        "phage.fasta": source_dir / "phage.fasta",
    }


def _host_candidates(host_dir: Path, host_id: str) -> List[Path]:
    """Return all candidate per-host FASTA paths for supported extensions."""
    return [host_dir / f"{host_id}{ext}" for ext in HOST_FASTA_EXTENSIONS]


def _resolve_host_file(host_dir: Path, host_id: str) -> Optional[Path]:
    """Return the first existing per-host FASTA path, or None when absent."""
    for candidate in _host_candidates(host_dir, host_id):
        if candidate.exists():
            return candidate
    return None


def validate_private_source(source_dir: Path, include_dataframe: bool = False) -> PrivateSourceValidation:
    source_dir = Path(source_dir)
    source_db = source_dir.name
    errors: List[str] = []
    warnings: List[str] = []
    stats: Dict[str, int] = {"rows": 0, "unique_phages": 0, "unique_hosts": 0, "extra_attributes": 0}

    required = _required_files(source_dir)
    missing = [name for name, path in required.items() if not path.exists()]
    if missing:
        errors.append(f"Missing required files: {', '.join(missing)}")
        return PrivateSourceValidation(source_dir, source_db, errors, warnings, stats)

    try:
        df = pd.read_csv(required["metadata.csv"], dtype=str, keep_default_na=False)
    except Exception as exc:
        errors.append(f"metadata.csv cannot be parsed: {exc}")
        return PrivateSourceValidation(source_dir, source_db, errors, warnings, stats)

    missing_columns = [col for col in MANDATORY_COLUMNS if col not in df.columns]
    if missing_columns:
        errors.append(f"metadata.csv missing required columns: {missing_columns}")
        return PrivateSourceValidation(source_dir, source_db, errors, warnings, stats)

    df = df.copy()
    for col in df.select_dtypes(include=["object", "string"]).columns:
        df[col] = df[col].str.strip()

    stats["rows"] = len(df)
    if len(df) == 0:
        errors.append("metadata.csv is empty")

    for col in MANDATORY_COLUMNS:
        invalid_rows = df.index[df[col].astype(str).str.len() == 0].tolist()
        if invalid_rows:
            errors.append(
                f"metadata.csv column '{col}' has empty values at rows "
                f"{[row + 2 for row in invalid_rows[:MAX_ERROR_EXAMPLES]]}"
            )

    unique_sources = sorted(set(v for v in df["Source_DB"].astype(str).tolist() if v))
    if unique_sources != [source_db]:
        errors.append(
            f"Source_DB mismatch: expected only '{source_db}', found {unique_sources or ['<empty>']}"
        )

    normalized_interactions = df["interaction"].astype(str).str.lower()
    df["interaction"] = normalized_interactions

    phage_ids, phage_duplicates = parse_fasta_ids(required["phage.fasta"])
    if phage_duplicates:
        errors.append(
            f"Duplicate FASTA identifiers in phage.fasta: "
            f"{sorted(phage_duplicates)[:MAX_ERROR_EXAMPLES]}"
        )

    csv_phage_ids = set(df["Phage_ID"].tolist())
    csv_host_ids = set(df["Host_ID"].tolist())

    missing_phages = sorted(csv_phage_ids - phage_ids)
    if missing_phages:
        errors.append(f"Phage_ID not found in phage.fasta: {missing_phages[:MAX_ERROR_EXAMPLES]}")

    host_dir = source_dir / "hosts"
    if host_dir.exists() and host_dir.is_dir():
        missing_hosts = []
        invalid_host_files = []
        for host_id in sorted(csv_host_ids):
            host_file = _resolve_host_file(host_dir, host_id)
            if host_file is None:
                missing_hosts.append(host_id)
                continue
            file_host_ids, host_duplicates = parse_fasta_ids(host_file)
            if host_duplicates:
                errors.append(
                    f"Duplicate FASTA identifiers in hosts/{host_file.name}: "
                    f"{sorted(host_duplicates)[:MAX_ERROR_EXAMPLES]}"
                )
            if host_id not in file_host_ids:
                invalid_host_files.append(f"{host_id} -> {host_file.name}")
        if missing_hosts:
            errors.append(
                "Host_ID missing corresponding FASTA file in hosts/: "
                f"{missing_hosts[:MAX_ERROR_EXAMPLES]}"
            )
        if invalid_host_files:
            errors.append(
                "Host FASTA files in hosts/ do not contain expected Host_ID headers "
                "(expected Host_ID -> file): "
                f"{invalid_host_files[:MAX_ERROR_EXAMPLES]}"
            )
    else:
        # hosts/ directory absent — only valid if all Host_ID values are "unknown"
        unknown_only = all(hid == "unknown" for hid in csv_host_ids)
        if not unknown_only:
            errors.append(
                "Missing required host sequences: provide a hosts/<Host_ID>.fna file for each Host_ID"
            )
        else:
            warnings.append(
                "No hosts/ directory found; all Host_ID values are 'unknown'. "
                "Host sequences will not be available for these phages."
            )

    duplicated_rows = int(df.duplicated(subset=["Phage_ID", "Host_ID", "Source_DB"]).sum())
    if duplicated_rows:
        warnings.append(
            f"Found {duplicated_rows} duplicated (Phage_ID, Host_ID, Source_DB) rows; they will be deduplicated"
        )

    stats["unique_phages"] = int(df["Phage_ID"].nunique())
    stats["unique_hosts"] = int(df["Host_ID"].nunique())
    extra_columns = [col for col in df.columns if col not in MANDATORY_COLUMNS]
    stats["extra_attributes"] = int((df[extra_columns] != "").sum().sum()) if extra_columns else 0

    return PrivateSourceValidation(
        source_dir=source_dir,
        source_db=source_db,
        errors=errors,
        warnings=warnings,
        stats=stats,
        metadata_df=df if include_dataframe else None,
        extra_columns=extra_columns,
    )


def _normalize_excluded_source_dirs(excluded_source_dirs: Optional[Iterable[str | Path]]) -> Set[Path]:
    if not excluded_source_dirs:
        return set()
    normalized: Set[Path] = set()
    for excluded in excluded_source_dirs:
        try:
            normalized.add(Path(excluded).expanduser().resolve())
        except OSError:
            continue
    return normalized


def discover_private_sources(
    roots: Iterable[str],
    excluded_source_dirs: Optional[Iterable[str | Path]] = None,
) -> List[Path]:
    discovered: List[Path] = []
    seen: Set[Path] = set()
    excluded_dirs = _normalize_excluded_source_dirs(excluded_source_dirs)
    for root in roots:
        root_path = Path(root).expanduser()
        if not root_path.exists() or not root_path.is_dir():
            continue
        for child in sorted(root_path.iterdir()):
            if not child.is_dir():
                continue
            resolved_child = child.resolve()
            if resolved_child in excluded_dirs:
                continue
            if resolved_child not in seen:
                discovered.append(child)
                seen.add(resolved_child)
    return discovered


def validate_private_roots(
    roots: Iterable[str],
    include_dataframe: bool = False,
    excluded_source_dirs: Optional[Iterable[str | Path]] = None,
) -> Dict:
    sources = discover_private_sources(roots, excluded_source_dirs=excluded_source_dirs)
    duplicate_names = {}
    for src in sources:
        duplicate_names[src.name] = duplicate_names.get(src.name, 0) + 1

    validations: List[PrivateSourceValidation] = []
    for src in sources:
        result = validate_private_source(src, include_dataframe=include_dataframe)
        if duplicate_names[src.name] > 1:
            result.errors.append(
                f"Duplicate source directory name '{src.name}' across roots is not allowed"
            )
        validations.append(result)

    return {
        "roots": [str(Path(r).expanduser()) for r in roots],
        "sources_found": len(sources),
        "sources_valid": sum(1 for v in validations if v.is_valid),
        "sources_invalid": sum(1 for v in validations if not v.is_valid),
        "sources": validations,
    }


def _source_fingerprint(source_dir: Path) -> str:
    hasher = hashlib.sha256()
    for file_name in ("metadata.csv", "phage.fasta"):
        path = source_dir / file_name
        if not path.exists():
            continue
        hasher.update(file_name.encode("utf-8"))
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                hasher.update(chunk)
    hosts_dir = source_dir / "hosts"
    if hosts_dir.is_dir():
        for host_file in sorted(hosts_dir.iterdir()):
            if host_file.is_file():
                hasher.update(host_file.name.encode("utf-8"))
                with host_file.open("rb") as handle:
                    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                        hasher.update(chunk)
    return hasher.hexdigest()


def build_private_manifest(
    roots: Iterable[str],
    excluded_source_dirs: Optional[Iterable[str | Path]] = None,
) -> Dict:
    summary = validate_private_roots(
        roots,
        include_dataframe=False,
        excluded_source_dirs=excluded_source_dirs,
    )
    manifest_sources = []
    for src in summary["sources"]:
        manifest_sources.append(
            {
                "source_db": src.source_db,
                "source_dir": str(src.source_dir),
                "is_valid": src.is_valid,
                "errors": src.errors,
                "warnings": src.warnings,
                "stats": src.stats,
                "fingerprint": _source_fingerprint(src.source_dir),
            }
        )

    return {
        "roots": summary["roots"],
        "sources_found": summary["sources_found"],
        "sources_valid": summary["sources_valid"],
        "sources_invalid": summary["sources_invalid"],
        "sources": manifest_sources,
    }


def write_private_manifest(manifest: Dict, output_path: Path) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")


def _ensure_duplicate_columns(conn: "duckdb.DuckDBPyConnection") -> None:
    """Add duplicate detection columns to fact_phages if they don't exist."""
    try:
        columns = {
            row[1] for row in conn.execute("PRAGMA table_info(fact_phages)").fetchall()
        }
        cols_to_add = []
        if "is_duplicate_of_public" not in columns:
            cols_to_add.append(("is_duplicate_of_public", "BOOLEAN DEFAULT FALSE"))
        if "duplicate_public_id" not in columns:
            cols_to_add.append(("duplicate_public_id", "VARCHAR"))
        if "duplicate_pident" not in columns:
            cols_to_add.append(("duplicate_pident", "DOUBLE"))
        if "duplicate_qcovs" not in columns:
            cols_to_add.append(("duplicate_qcovs", "DOUBLE"))

        for col_name, col_type in cols_to_add:
            conn.execute(f"ALTER TABLE fact_phages ADD COLUMN {col_name} {col_type}")
    except Exception:
        # Columns may already exist or table not yet created - non-critical
        pass


def ingest_private_sources_into_db(
    conn: duckdb.DuckDBPyConnection,
    source_dirs: Iterable[str],
    blast_db_dir: Path | str | None = None,
    private_phage_mapping_path: Path | str | None = None,
    pident_threshold: float = 99.0,
    qcovs_threshold: float = 90.0,
) -> Dict:
    ingested = []
    skipped = []
    allowed_source_tables = {"fact_phages", "private_interactions", "private_entity_attributes"}

    def _delete_private_rows_for_sources(table_name: str, source_dbs: List[str]) -> None:
        if table_name not in allowed_source_tables:
            raise ValueError(f"Unsupported table for private-source cleanup: {table_name}")
        if not source_dbs:
            return
        placeholders = ", ".join(["?"] * len(source_dbs))
        conn.execute(
            f"""
            DELETE FROM {table_name}
            WHERE source_type = 'private'
              AND Source_DB IN ({placeholders})
            """,
            source_dbs,
        )

    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS private_entity_attributes (
            Source_DB VARCHAR,
            source_type VARCHAR,
            entity_type VARCHAR,
            entity_id VARCHAR,
            attribute_key VARCHAR,
            attribute_value VARCHAR,
            value_type VARCHAR
        )
        """
    )
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS private_interactions (
            Phage_ID VARCHAR,
            Host_ID VARCHAR,
            Host_name VARCHAR,
            interaction VARCHAR,
            Source_DB VARCHAR,
            source_type VARCHAR
        )
        """
    )
    conn.execute("DROP TABLE IF EXISTS private_phage_host_associations")
    conn.execute(
        """
        -- Materialized helper table used by create_duckdb views and downstream
        -- joins that need source-aware private phage-host links.
        CREATE TABLE private_phage_host_associations (
            Phage_ID VARCHAR,
            Host_ID VARCHAR,
            Source_DB VARCHAR,
            source_type VARCHAR
        )
        """
    )

    validations: List[PrivateSourceValidation] = []
    for source_dir_str in source_dirs:
        validation = validate_private_source(Path(source_dir_str), include_dataframe=True)
        if not validation.is_valid or validation.metadata_df is None:
            skipped.append(
                {
                    "source_db": validation.source_db,
                    "source_dir": str(validation.source_dir),
                    "errors": validation.errors,
                }
            )
            continue
        validations.append(validation)

    current_source_dbs = sorted({validation.source_db for validation in validations})
    _delete_private_rows_for_sources("fact_phages", current_source_dbs)
    _delete_private_rows_for_sources("private_interactions", current_source_dbs)
    _delete_private_rows_for_sources("private_entity_attributes", current_source_dbs)

    # Ensure duplicate detection columns exist on fact_phages
    _ensure_duplicate_columns(conn)

    existing_private_sources = {
        row[0]
        for row in conn.execute(
            "SELECT DISTINCT Source_DB FROM fact_phages WHERE source_type = 'private'"
        ).fetchall()
        if row[0]
    }
    stale_source_dbs = sorted(existing_private_sources - set(current_source_dbs))
    _delete_private_rows_for_sources("fact_phages", stale_source_dbs)
    _delete_private_rows_for_sources("private_interactions", stale_source_dbs)
    _delete_private_rows_for_sources("private_entity_attributes", stale_source_dbs)

    for validation in validations:
        df = validation.metadata_df.drop_duplicates(subset=["Phage_ID", "Host_ID", "Source_DB"]).copy()
        df["source_type"] = "private"

        private_phages = (
            df[["Phage_ID", "Source_DB", "Host_name", "interaction", "source_type"]]
            .drop_duplicates(subset=["Phage_ID", "Source_DB"])
            .rename(columns={"Host_name": "Host", "interaction": "Lifestyle"})
        )
        private_phages["Length"] = pd.NA
        private_phages["GC_content"] = pd.NA
        private_phages["Taxonomy"] = pd.NA
        private_phages["Completeness"] = pd.NA
        private_phages["Cluster"] = pd.NA
        private_phages["Subcluster"] = pd.NA
        private_phages["is_duplicate_of_public"] = False
        private_phages["duplicate_public_id"] = pd.NA
        private_phages["duplicate_pident"] = pd.NA
        private_phages["duplicate_qcovs"] = pd.NA
        private_phages = private_phages[
            [
                "Phage_ID",
                "Source_DB",
                "Length",
                "GC_content",
                "Taxonomy",
                "Completeness",
                "Host",
                "Lifestyle",
                "Cluster",
                "Subcluster",
                "source_type",
                "is_duplicate_of_public",
                "duplicate_public_id",
                "duplicate_pident",
                "duplicate_qcovs",
            ]
        ]
        conn.register("private_phages_df", private_phages)
        conn.execute(
            """
            INSERT INTO fact_phages BY NAME
            SELECT * FROM private_phages_df
            """
        )

        private_interactions = df[
            ["Phage_ID", "Host_ID", "Host_name", "interaction", "Source_DB", "source_type"]
        ].copy()
        conn.register("private_interactions_df", private_interactions)
        conn.execute(
            """
            INSERT INTO private_interactions
            SELECT * FROM private_interactions_df
            """
        )

        extra_columns = validation.extra_columns or []
        if extra_columns:
            attrs = df[["Phage_ID", "Host_ID", "Source_DB"] + extra_columns].copy()
            attrs = attrs.melt(
                id_vars=["Phage_ID", "Host_ID", "Source_DB"],
                value_vars=extra_columns,
                var_name="attribute_key",
                value_name="attribute_value",
            )
            attrs = attrs[attrs["attribute_value"].astype(str).str.len() > 0].copy()
            attrs["source_type"] = "private"
            attrs["entity_type"] = "interaction"
            attrs["entity_id"] = attrs["Phage_ID"] + "|" + attrs["Host_ID"]
            attrs["value_type"] = "string"
            attrs = attrs[
                [
                    "Source_DB",
                    "source_type",
                    "entity_type",
                    "entity_id",
                    "attribute_key",
                    "attribute_value",
                    "value_type",
                ]
            ]
            if not attrs.empty:
                conn.register("private_attrs_df", attrs)
                conn.execute(
                    """
                    INSERT INTO private_entity_attributes
                    SELECT * FROM private_attrs_df
                    """
                )

        ingested.append(
            {
                "source_db": validation.source_db,
                "source_dir": str(validation.source_dir),
                "rows": validation.stats["rows"],
                "unique_phages": validation.stats["unique_phages"],
                "unique_hosts": validation.stats["unique_hosts"],
            }
        )

    # Recompute the association table from canonical private_interactions to keep
    # link rows exactly in sync after source additions/removals/updates.
    conn.execute(
        """
        INSERT INTO private_phage_host_associations
        SELECT DISTINCT Phage_ID, Host_ID, Source_DB, source_type
        FROM private_interactions
        WHERE source_type = 'private'
        """
    )

    if "dim_hosts" not in set(conn.execute("SHOW TABLES").fetchnumpy()["name"]):
        conn.execute(
            """
            CREATE TABLE dim_hosts (
                Host_ID VARCHAR,
                Species_Name VARCHAR,
                Strain_Name VARCHAR,
                Assembly_Accession VARCHAR,
                Assembly_Name VARCHAR,
                Assembly_Level VARCHAR,
                Genome_Length BIGINT,
                GC_Content DOUBLE,
                RefSeq_Category VARCHAR,
                Download_Date VARCHAR,
                Source VARCHAR,
                source_type VARCHAR
            )
            """
        )

    # Private host rows are regenerated from current private_interactions so stale
    # hosts disappear automatically when sources are removed.
    conn.execute("DELETE FROM dim_hosts WHERE source_type = 'private'")
    conn.execute(
        """
        INSERT INTO dim_hosts
        SELECT
            Host_ID,
            Host_name AS Species_Name,
            NULL AS Strain_Name,
            NULL AS Assembly_Accession,
            NULL AS Assembly_Name,
            NULL AS Assembly_Level,
            NULL AS Genome_Length,
            NULL AS GC_Content,
            NULL AS RefSeq_Category,
            NULL AS Download_Date,
            'private' AS Source,
            'private' AS source_type
        FROM (
            SELECT DISTINCT Host_ID, Host_name
            FROM private_interactions
            WHERE source_type = 'private'
              AND Host_ID != 'unknown'
        )
        """
    )

    # Run duplicate detection against public BLAST DB if available
    duplicate_stats: Dict = {"checked": 0, "duplicates_found": 0}
    if blast_db_dir and private_phage_mapping_path:
        try:
            duplicate_stats = check_private_duplicates_against_public(
                conn=conn,
                blast_db_dir=blast_db_dir,
                private_phage_mapping_path=private_phage_mapping_path,
                pident_threshold=pident_threshold,
                qcovs_threshold=qcovs_threshold,
            )
        except Exception as exc:
            logger.warning("Duplicate detection failed (non-blocking): %s", exc)

    return {"ingested": ingested, "skipped": skipped, "duplicate_check": duplicate_stats}


def _iter_fasta_records(fasta_path: Path) -> Iterable[Tuple[str, str, str]]:
    """Yield FASTA records as (id_token, full_header_without_>, sequence)."""
    header: Optional[str] = None
    seq_parts: List[str] = []

    with fasta_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_parts)
                    record_id = header.split()[0] if header else ""
                    if record_id:
                        yield record_id, header, seq
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)

    if header is not None:
        seq = "".join(seq_parts)
        record_id = header.split()[0] if header else ""
        if record_id:
            yield record_id, header, seq


def _write_fasta_record(handle, header: str, sequence: str, line_width: int = 80) -> None:
    handle.write(f">{header}\n")
    for i in range(0, len(sequence), line_width):
        handle.write(sequence[i:(i + line_width)] + "\n")


def _hash_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    """Compute a SHA-256 digest for a file using chunked reads."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def prepare_private_sequence_artifacts(
    manifest: Dict,
    private_phage_dir: Path,
    private_phage_mapping_path: Path,
    private_host_mapping_path: Path,
) -> Dict:
    """
    Build sequence artifacts used by SequenceRetriever for private datasets.

    For each valid private source this function:
      1. Copies the source ``phage.fasta`` (filtering to wanted Phage_IDs) into a
         per-source writable directory (``private_phage_dir/<source_db>/phage.fasta``)
         and indexes it with pyfaidx so that SequenceRetriever can open it on-demand.
      2. Registers host sequences from the source ``hosts/`` directory directly,
         mapping each ``Host_ID`` to ``<source_dir>/hosts/<Host_ID>.fna`` (no copy).

    Outputs written to disk:
      - ``private_phage_mapping_path``  (JSON) – maps ``source_db`` to the path of
        the copied phage FASTA.  SequenceRetriever uses this to route private-phage
        sequence lookups without mixing private data into ``all_phages.fasta``.
      - ``private_host_mapping_path``  (JSON) – maps ``Host_ID`` to the path of the
        per-host ``.fna`` file inside the source ``hosts/`` directory.  Merged into
        the global ``host_fasta_mapping.json`` by the ``create_host_mapping`` rule.

    Args:
        manifest: Parsed private manifest JSON (``{"sources": [...]}``)
        private_phage_dir: Writable base directory for per-source phage FASTAs.
            Each source gets its own sub-directory:
            ``private_phage_dir/<source_db>/phage.fasta``.
        private_phage_mapping_path: Output path for the phage mapping JSON.
        private_host_mapping_path: Output path for the host mapping JSON.

    Returns:
        Dict with counts of processed sources and written/skipped sequences.
    """
    import pyfaidx  # imported here to keep it out of the top-level import chain for tests

    private_phage_dir = Path(private_phage_dir)
    private_phage_mapping_path = Path(private_phage_mapping_path)
    private_host_mapping_path = Path(private_host_mapping_path)

    private_phage_dir.mkdir(parents=True, exist_ok=True)
    private_phage_mapping_path.parent.mkdir(parents=True, exist_ok=True)
    private_host_mapping_path.parent.mkdir(parents=True, exist_ok=True)

    valid_sources = [
        src for src in manifest.get("sources", [])
        if src.get("is_valid", False) and src.get("source_dir")
    ]

    # source_db → path to per-source copied+indexed phage.fasta
    phage_mapping: Dict[str, str] = {}
    host_mapping: Dict[str, str] = {}
    host_hashes: Dict[str, str] = {}
    stats: Dict = {
        "sources_processed": 0,
        "phages_written": 0,
        "hosts_written": 0,
        "host_duplicates_identical": 0,
        "host_duplicates_conflicting": 0,
        "missing_phage_ids": 0,
        "missing_host_ids": 0,
    }

    for src in valid_sources:
        source_dir = Path(src["source_dir"])
        source_db = src.get("source_db") or source_dir.name
        metadata_path = source_dir / "metadata.csv"
        phage_fasta_path = source_dir / "phage.fasta"
        host_dir_path = source_dir / "hosts"

        if not (metadata_path.exists() and phage_fasta_path.exists()):
            continue

        stats["sources_processed"] += 1
        df = pd.read_csv(metadata_path, dtype=str, keep_default_na=False)
        for col in df.columns:
            if pd.api.types.is_string_dtype(df[col]):
                df[col] = df[col].str.strip()

        wanted_phages = set(df["Phage_ID"].tolist())
        wanted_hosts = set(df["Host_ID"].tolist())

        # ── Phage sequences ────────────────────────────────────────────────────
        # Copy the filtered phage sequences to a writable per-source directory so
        # that pyfaidx can create the .fai index file next to the FASTA.
        source_output_dir = private_phage_dir / source_db
        source_output_dir.mkdir(parents=True, exist_ok=True)
        dest_phage_fasta = source_output_dir / "phage.fasta"

        found_phages: Set[str] = set()
        with dest_phage_fasta.open("w", encoding="utf-8") as phage_out:
            for rec_id, header, sequence in _iter_fasta_records(phage_fasta_path):
                if rec_id not in wanted_phages:
                    continue
                found_phages.add(rec_id)
                _write_fasta_record(phage_out, header, sequence)
                stats["phages_written"] += 1

        # Index the copied FASTA with pyfaidx using the same split_char as the
        # main phage index so that SequenceRetriever's key_function works identically.
        try:
            pyfaidx.Fasta(
                str(dest_phage_fasta),
                split_char="\x00",
                rebuild=True,
                read_long_names=True,
            )
            logger.info("Indexed private phage FASTA for source '%s': %s", source_db, dest_phage_fasta)
        except Exception as exc:
            logger.warning("Failed to index phage FASTA for source '%s': %s", source_db, exc)

        phage_mapping[source_db] = str(dest_phage_fasta)
        stats["missing_phage_ids"] += len(wanted_phages - found_phages)

        # ── Host sequences ─────────────────────────────────────────────────────
        # Host FASTA files live in the source's hosts/ directory.
        # Map each Host_ID directly to its file — no copying needed.
        if host_dir_path.exists() and host_dir_path.is_dir():
            logger.info("Using hosts/ directory for private source '%s'", source_db)
            found_hosts: Set[str] = set()
            for host_id in sorted(wanted_hosts):
                source_host_file = _resolve_host_file(host_dir_path, host_id)
                if source_host_file is None:
                    continue

                found_hosts.add(host_id)
                seq_hash = _hash_file(source_host_file)

                if host_id in host_mapping:
                    existing_seq_hash = host_hashes.get(host_id)
                    if existing_seq_hash == seq_hash:
                        stats["host_duplicates_identical"] += 1
                    else:
                        stats["host_duplicates_conflicting"] += 1
                    continue

                # Point directly at the source file — no copy under any intermediate dir.
                host_mapping[host_id] = str(source_host_file)
                host_hashes[host_id] = seq_hash
                stats["hosts_written"] += 1

            stats["missing_host_ids"] += len(wanted_hosts - found_hosts)
        else:
            # No hosts/ directory — valid for phage-only sources where all
            # Host_ID values are "unknown".  Skip host mapping entirely.
            logger.info(
                "No hosts/ directory for private source '%s'; skipping host sequence mapping",
                source_db,
            )

    private_phage_mapping_path.write_text(
        json.dumps(dict(sorted(phage_mapping.items())), indent=2),
        encoding="utf-8",
    )

    private_host_mapping_path.write_text(
        json.dumps(dict(sorted(host_mapping.items())), indent=2),
        encoding="utf-8",
    )

    return stats


def check_private_duplicates_against_public(
    conn: "duckdb.DuckDBPyConnection",
    blast_db_dir: Path | str,
    private_phage_mapping_path: Path | str,
    pident_threshold: float = 99.0,
    qcovs_threshold: float = 90.0,
) -> Dict:
    """Search private phage sequences against the public BLAST DB to detect duplicates.

    Uses the private phage mapping (source_db -> phage.fasta path) to locate
    private sequences, then searches each against the public ``phages`` BLAST
    database. Matches with identity above *pident_threshold* and query coverage
    above *qcovs_threshold* are flagged as duplicates in ``fact_phages``.

    Args:
        conn: Active DuckDB connection (must have ``fact_phages`` table).
        blast_db_dir: Path to the BLAST database directory containing the
            ``phages/all_phages`` public database.
        private_phage_mapping_path: Path to the ``private_phage_mapping.json``
            produced by ``prepare_private_sequences``.
        pident_threshold: Minimum percent identity to flag as duplicate (default 99.0).
        qcovs_threshold: Minimum query coverage % to flag as duplicate (default 90.0).

    Returns:
        Dict with keys ``checked``, ``duplicates_found``, ``sources``, and
        ``duplicates`` (mapping Phage_ID -> match info).
    """
    blast_db_dir = Path(blast_db_dir)
    private_phage_mapping_path = Path(private_phage_mapping_path)

    stats: Dict = {
        "checked": 0,
        "duplicates_found": 0,
        "sources": [],
        "duplicates": {},
    }

    # If no private phage mapping or BLAST DB, skip gracefully
    if not private_phage_mapping_path.exists():
        logger.warning(
            "Private phage mapping not found at %s; skipping duplicate check",
            private_phage_mapping_path,
        )
        return stats

    if not (blast_db_dir / "phages" / "all_phages").exists():
        logger.warning(
            "Public BLAST DB not found at %s/phages/all_phages; skipping duplicate check",
            blast_db_dir,
        )
        return stats

    with open(private_phage_mapping_path, encoding="utf-8") as f:
        phage_mapping: Dict[str, str] = json.load(f)

    if not phage_mapping:
        logger.info("No private phage sources to check for duplicates")
        return stats

    # Import BlastSearcher here to avoid circular imports at module level
    from pbi.blast_search import BlastSearcher

    searcher = BlastSearcher(blast_db_dir)

    # Concatenate all private phage FASTA files into a temporary file
    import tempfile

    tmp_fasta = tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False, encoding="utf-8"
    )
    seq_source: Dict[str, str] = {}  # seq_id -> source_db

    try:
        for source_db, fasta_path in sorted(phage_mapping.items()):
            fasta_path = Path(fasta_path)
            if not fasta_path.exists():
                logger.warning(
                    "Private FASTA not found for source '%s': %s",
                    source_db, fasta_path,
                )
                continue
            stats["sources"].append(source_db)
            with open(fasta_path, encoding="utf-8") as in_f:
                for line in in_f:
                    tmp_fasta.write(line)
                    if line.startswith(">"):
                        seq_id = line[1:].strip().split()[0]
                        seq_source[seq_id] = source_db
        tmp_fasta.close()

        if not seq_source:
            logger.info("No private sequences found to check for duplicates")
            return stats

        stats["checked"] = len(seq_source)

        # Search private sequences against the public phages DB
        results = searcher.search_fasta(
            Path(tmp_fasta.name),
            program="blastn",
            db="phages",
            max_hits=1,
            evalue=1e-5,
        )

        if results.empty:
            logger.info("No BLAST hits found against public data")
            return stats

        # Filter by identity threshold
        matches = results[results["pident"] >= pident_threshold].copy()

        if matches.empty:
            logger.info(
                "No private phages with >= %.1f%% identity to public data",
                pident_threshold,
            )
            return stats

        # Build duplicate mapping: private Phage_ID -> best public match
        duplicates: Dict[str, Dict] = {}
        for _, row in matches.iterrows():
            query_id = row["qseqid"]
            subject_id = row["sseqid"]
            pident = float(row["pident"])

            # Keep only the best hit per query
            if query_id not in duplicates or pident > duplicates[query_id]["pident"]:
                duplicates[query_id] = {
                    "public_phage_id": subject_id,
                    "pident": pident,
                    "evalue": float(row.get("evalue", 0)),
                    "source_db": seq_source.get(query_id, "unknown"),
                }

        stats["duplicates_found"] = len(duplicates)
        stats["duplicates"] = duplicates

        # Update fact_phages with duplicate annotations
        for phage_id, dup_info in duplicates.items():
            source_db = dup_info["source_db"]
            conn.execute(
                """
                UPDATE fact_phages
                SET is_duplicate_of_public = TRUE,
                    duplicate_public_id = ?,
                    duplicate_pident = ?,
                    duplicate_qcovs = NULL
                WHERE Phage_ID = ?
                  AND Source_DB = ?
                  AND source_type = 'private'
                """,
                [dup_info["public_phage_id"], dup_info["pident"], phage_id, source_db],
            )

        logger.info(
            "Duplicate check complete: %d/%d private phages match public data "
            "(pident>=%.1f%%)",
            stats["duplicates_found"], stats["checked"], pident_threshold,
        )

    except FileNotFoundError as exc:
        logger.warning("BLAST not available or DB missing; skipping duplicate check: %s", exc)
    except Exception as exc:
        logger.warning("Duplicate check failed (non-blocking): %s", exc)
    finally:
        try:
            os.unlink(tmp_fasta.name)
        except OSError:
            pass

    return stats
