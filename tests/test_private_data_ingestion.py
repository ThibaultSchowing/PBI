#!/usr/bin/env python3

from pathlib import Path

import duckdb
import json
import pandas as pd

from pbi.private_data import (
    build_private_manifest,
    ingest_private_sources_into_db,
    prepare_private_sequence_artifacts,
    validate_private_source,
)


def _write_text(path: Path, content: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _create_private_source(
    root: Path,
    name: str,
    metadata_rows: list[dict],
    phage_ids: list[str],
    host_ids: list[str],
):
    source_dir = root / name
    source_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(metadata_rows).to_csv(source_dir / "metadata.csv", index=False)
    _write_text(source_dir / "phage.fasta", "".join([f">{pid} desc\nATGC\n" for pid in phage_ids]))
    hosts_dir = source_dir / "hosts"
    hosts_dir.mkdir(parents=True, exist_ok=True)
    for hid in host_ids:
        _write_text(hosts_dir / f"{hid}.fna", f">{hid} desc\nATGC\n")
    return source_dir


def _prepare_minimal_db(conn: duckdb.DuckDBPyConnection):
    conn.execute(
        """
        CREATE TABLE fact_phages (
            Phage_ID VARCHAR,
            Source_DB VARCHAR,
            Length INTEGER,
            GC_content DOUBLE,
            Taxonomy VARCHAR,
            Completeness VARCHAR,
            Host VARCHAR,
            Lifestyle VARCHAR,
            Cluster VARCHAR,
            Subcluster VARCHAR,
            Provider_Name VARCHAR,
            Provider_Release VARCHAR,
            Provider_Snapshot_Date VARCHAR,
            Provider_Schema_Profile VARCHAR,
            Input_Source_Key VARCHAR,
            Input_File VARCHAR,
            Input_Retrieved_At VARCHAR,
            source_type VARCHAR,
            is_duplicate_of_public BOOLEAN DEFAULT FALSE,
            duplicate_public_id VARCHAR,
            duplicate_pident DOUBLE,
            duplicate_qcovs DOUBLE
        )
        """
    )


def test_validate_private_source_success(tmp_path):
    source = _create_private_source(
        tmp_path,
        "Project_A",
        [
            {
                "Phage_ID": "P1",
                "Host_ID": "H1",
                "Host_name": "Host One",
                "Source_DB": "Project_A",
                "interaction": "Virulent",
                "Lab_Condition": "anaerobic",
            }
        ],
        ["P1"],
        ["H1"],
    )
    result = validate_private_source(source)
    assert result.is_valid
    assert result.stats["rows"] == 1


def test_validate_private_source_missing_required_file(tmp_path):
    source = tmp_path / "Project_A"
    source.mkdir()
    # No phage.fasta — this is still a required file
    _write_text(source / "metadata.csv", "Phage_ID,Host_ID,Host_name,Source_DB,interaction\n")
    result = validate_private_source(source)
    assert not result.is_valid
    assert any("Missing required files" in e for e in result.errors)


def test_validate_private_source_missing_required_column(tmp_path):
    source = _create_private_source(
        tmp_path,
        "Project_A",
        [{"Phage_ID": "P1", "Host_ID": "H1", "Source_DB": "Project_A", "interaction": "virulent"}],
        ["P1"],
        ["H1"],
    )
    result = validate_private_source(source)
    assert not result.is_valid
    assert any("missing required columns" in e for e in result.errors)


def test_validate_private_source_empty_interaction(tmp_path):
    source = _create_private_source(
        tmp_path,
        "Project_A",
        [
            {
                "Phage_ID": "P1",
                "Host_ID": "H1",
                "Host_name": "Host One",
                "Source_DB": "Project_A",
                "interaction": "",
            }
        ],
        ["P1"],
        ["H1"],
    )
    result = validate_private_source(source)
    assert not result.is_valid
    assert any("empty values" in e for e in result.errors)


def test_validate_private_source_custom_interaction_types(tmp_path):
    """Custom interaction types (e.g. 'interaction', 'no interaction') should be accepted."""
    source = _create_private_source(
        tmp_path,
        "Project_Custom",
        [
            {
                "Phage_ID": "P1",
                "Host_ID": "H1",
                "Host_name": "Host One",
                "Source_DB": "Project_Custom",
                "interaction": "interaction",
            },
            {
                "Phage_ID": "P2",
                "Host_ID": "H2",
                "Host_name": "Host Two",
                "Source_DB": "Project_Custom",
                "interaction": "no interaction",
            },
        ],
        ["P1", "P2"],
        ["H1", "H2"],
    )
    result = validate_private_source(source, include_dataframe=True)
    assert result.is_valid
    assert result.stats["rows"] == 2
    # Values should be normalized to lowercase
    assert list(result.metadata_df["interaction"]) == ["interaction", "no interaction"]


def test_validate_private_source_any_interaction_value(tmp_path):
    """Any non-empty interaction string should be accepted."""
    source = _create_private_source(
        tmp_path,
        "Project_Any",
        [
            {
                "Phage_ID": "P1",
                "Host_ID": "H1",
                "Host_name": "Host One",
                "Source_DB": "Project_Any",
                "interaction": "Latent",
            }
        ],
        ["P1"],
        ["H1"],
    )
    result = validate_private_source(source, include_dataframe=True)
    assert result.is_valid
    # Should be normalized to lowercase
    assert result.metadata_df.iloc[0]["interaction"] == "latent"


def test_validate_private_source_id_mismatch(tmp_path):
    source = _create_private_source(
        tmp_path,
        "Project_A",
        [
            {
                "Phage_ID": "P1",
                "Host_ID": "H1",
                "Host_name": "Host One",
                "Source_DB": "Project_A",
                "interaction": "virulent",
            }
        ],
        ["P2"],
        ["H2"],
    )
    result = validate_private_source(source)
    assert not result.is_valid
    assert any("Phage_ID not found" in e for e in result.errors)
    assert any("Host_ID missing corresponding FASTA file in hosts/" in e for e in result.errors)


def test_validate_private_source_missing_host_sequences_is_invalid(tmp_path):
    """Private sources must provide a hosts/<Host_ID>.fna file for each Host_ID."""
    source = tmp_path / "Project_B"
    source.mkdir()
    pd.DataFrame(
        [
            {
                "Phage_ID": "P1",
                "Host_ID": "H1",
                "Host_name": "Bacteroides dorei",
                "Source_DB": "Project_B",
                "interaction": "virulent",
            }
        ]
    ).to_csv(source / "metadata.csv", index=False)
    _write_text(source / "phage.fasta", ">P1 desc\nATGC\n")
    # Intentionally no hosts/ directory

    result = validate_private_source(source)
    assert not result.is_valid
    assert any("Missing required host sequences" in e for e in result.errors)


def test_prepare_private_sequence_artifacts_requires_host_directory(tmp_path):
    """Sources without a hosts/ directory should be skipped via manifest validity."""
    source = tmp_path / "Project_NoHost"
    source.mkdir()
    pd.DataFrame(
        [
            {
                "Phage_ID": "P1",
                "Host_ID": "H1",
                "Host_name": "Bacteroides dorei",
                "Source_DB": "Project_NoHost",
                "interaction": "virulent",
            }
        ]
    ).to_csv(source / "metadata.csv", index=False)
    _write_text(source / "phage.fasta", ">P1 desc\nATGC\n")
    # No hosts/ directory

    manifest = {"sources": [{"source_db": "Project_NoHost", "source_dir": str(source), "is_valid": False}]}

    private_phage_dir = tmp_path / "private" / "phages"
    private_phage_mapping = tmp_path / "private" / "private_phage_mapping.json"
    private_host_mapping = tmp_path / "private" / "private_host_mapping.json"

    stats = prepare_private_sequence_artifacts(
        manifest=manifest,
        private_phage_dir=private_phage_dir,
        private_phage_mapping_path=private_phage_mapping,
        private_host_mapping_path=private_host_mapping,
    )

    assert stats["sources_processed"] == 0
    assert stats["phages_written"] == 0
    assert stats["hosts_written"] == 0

    # Phage mapping and host mapping should both be empty JSON objects
    assert json.loads(private_phage_mapping.read_text(encoding="utf-8")) == {}
    assert json.loads(private_host_mapping.read_text(encoding="utf-8")) == {}
    valid_source = _create_private_source(
        tmp_path,
        "Valid_Source",
        [
            {
                "Phage_ID": "P1",
                "Host_ID": "H1",
                "Host_name": "Host One",
                "Source_DB": "Valid_Source",
                "interaction": "temperate",
                "Batch": "B1",
            }
        ],
        ["P1"],
        ["H1"],
    )
    invalid_source = _create_private_source(
        tmp_path,
        "Invalid_Source",
        [
            {
                "Phage_ID": "PX",
                "Host_ID": "HX",
                "Host_name": "Host X",
                "Source_DB": "Wrong_Source",
                "interaction": "virulent",
            }
        ],
        ["PX"],
        ["HX"],
    )

    conn = duckdb.connect(":memory:")
    _prepare_minimal_db(conn)

    summary = ingest_private_sources_into_db(conn, [str(valid_source), str(invalid_source)])
    assert len(summary["ingested"]) == 1
    assert len(summary["skipped"]) == 1

    phage_count = conn.execute("SELECT COUNT(*) FROM fact_phages WHERE source_type = 'private'").fetchone()[0]
    assert phage_count == 1

    attrs_count = conn.execute("SELECT COUNT(*) FROM private_entity_attributes").fetchone()[0]
    assert attrs_count == 1


def test_prepare_private_sequence_artifacts_generates_private_fasta_and_mapping(tmp_path):
    source = _create_private_source(
        tmp_path,
        "Project_A",
        [
            {
                "Phage_ID": "P1",
                "Host_ID": "H1",
                "Host_name": "Host One",
                "Source_DB": "Project_A",
                "interaction": "virulent",
            }
        ],
        ["P1"],
        ["H1"],
    )
    manifest = {
        "sources": [
            {
                "source_db": "Project_A",
                "source_dir": str(source),
                "is_valid": True,
            }
        ]
    }

    private_phage_dir = tmp_path / "private" / "phages"
    private_phage_mapping = tmp_path / "private" / "private_phage_mapping.json"
    private_host_mapping = tmp_path / "private" / "private_host_mapping.json"

    stats = prepare_private_sequence_artifacts(
        manifest=manifest,
        private_phage_dir=private_phage_dir,
        private_phage_mapping_path=private_phage_mapping,
        private_host_mapping_path=private_host_mapping,
    )

    assert stats["sources_processed"] == 1
    assert stats["phages_written"] == 1
    assert stats["hosts_written"] == 1

    # The phage mapping JSON should point to the per-source phage FASTA copy
    with private_phage_mapping.open("r", encoding="utf-8") as handle:
        pmapping = json.load(handle)
    assert "Project_A" in pmapping
    phage_fasta_path = Path(pmapping["Project_A"])
    assert phage_fasta_path.exists()
    assert ">P1 desc" in phage_fasta_path.read_text(encoding="utf-8")
    # Index file must exist next to the phage FASTA
    assert Path(str(phage_fasta_path) + ".fai").exists()

    with private_host_mapping.open("r", encoding="utf-8") as handle:
        hmapping = json.load(handle)
    assert "H1" in hmapping
    mapped_path = Path(hmapping["H1"])
    # Mapping points directly to the source's hosts/ file — no intermediate copy.
    assert mapped_path == source / "hosts" / "H1.fna"
    assert mapped_path.exists()
    assert ">H1 desc" in mapped_path.read_text(encoding="utf-8")


def test_prepare_private_sequence_artifacts_skips_invalid_sources(tmp_path):
    valid_source = _create_private_source(
        tmp_path,
        "Valid_Source",
        [
            {
                "Phage_ID": "P1",
                "Host_ID": "H1",
                "Host_name": "Host One",
                "Source_DB": "Valid_Source",
                "interaction": "temperate",
            }
        ],
        ["P1"],
        ["H1"],
    )
    invalid_source = _create_private_source(
        tmp_path,
        "Invalid_Source",
        [
            {
                "Phage_ID": "PX",
                "Host_ID": "HX",
                "Host_name": "Host X",
                "Source_DB": "Invalid_Source",
                "interaction": "virulent",
            }
        ],
        ["PX"],
        ["HX"],
    )

    manifest = {
        "sources": [
            {"source_db": "Valid_Source", "source_dir": str(valid_source), "is_valid": True},
            {"source_db": "Invalid_Source", "source_dir": str(invalid_source), "is_valid": False},
        ]
    }

    private_phage_dir = tmp_path / "private" / "phages"
    private_phage_mapping = tmp_path / "private" / "private_phage_mapping.json"
    private_host_mapping = tmp_path / "private" / "private_host_mapping.json"

    stats = prepare_private_sequence_artifacts(
        manifest=manifest,
        private_phage_dir=private_phage_dir,
        private_phage_mapping_path=private_phage_mapping,
        private_host_mapping_path=private_host_mapping,
    )

    assert stats["sources_processed"] == 1
    with private_phage_mapping.open("r", encoding="utf-8") as handle:
        pmapping = json.load(handle)
    assert "Valid_Source" in pmapping
    phage_fasta_path = Path(pmapping["Valid_Source"])
    assert ">P1 desc" in phage_fasta_path.read_text(encoding="utf-8")
    assert "Invalid_Source" not in pmapping


def test_build_private_manifest_excludes_generated_phage_directory(tmp_path):
    valid_source = _create_private_source(
        tmp_path,
        "Project_A",
        [
            {
                "Phage_ID": "P1",
                "Host_ID": "H1",
                "Host_name": "Host One",
                "Source_DB": "Project_A",
                "interaction": "virulent",
            }
        ],
        ["P1"],
        ["H1"],
    )

    generated_phage_dir = tmp_path / "phages"
    generated_phage_dir.mkdir()

    manifest = build_private_manifest(
        [str(tmp_path)],
        excluded_source_dirs=[str(generated_phage_dir)],
    )

    assert manifest["sources_found"] == 1
    assert manifest["sources_valid"] == 1
    assert manifest["sources_invalid"] == 0
    assert len(manifest["sources"]) == 1
    assert manifest["sources"][0]["source_db"] == valid_source.name


def test_ingest_private_sources_resyncs_deleted_sources_and_prevents_duplicates(tmp_path):
    source_a = _create_private_source(
        tmp_path,
        "Source_A",
        [
            {
                "Phage_ID": "PA",
                "Host_ID": "HA",
                "Host_name": "Host A",
                "Source_DB": "Source_A",
                "interaction": "virulent",
            }
        ],
        ["PA"],
        ["HA"],
    )
    source_b = _create_private_source(
        tmp_path,
        "Source_B",
        [
            {
                "Phage_ID": "PB",
                "Host_ID": "HB",
                "Host_name": "Host B",
                "Source_DB": "Source_B",
                "interaction": "temperate",
            }
        ],
        ["PB"],
        ["HB"],
    )

    conn = duckdb.connect(":memory:")
    _prepare_minimal_db(conn)

    first_run = ingest_private_sources_into_db(conn, [str(source_a), str(source_b)])
    assert len(first_run["ingested"]) == 2
    assert conn.execute("SELECT COUNT(*) FROM fact_phages WHERE source_type = 'private'").fetchone()[0] == 2

    # Rerun with only source_a (simulates source_b deleted from private_data folder)
    second_run = ingest_private_sources_into_db(conn, [str(source_a)])
    assert len(second_run["ingested"]) == 1

    private_rows = conn.execute(
        "SELECT Phage_ID, Source_DB FROM fact_phages WHERE source_type = 'private' ORDER BY Phage_ID"
    ).fetchall()
    assert private_rows == [("PA", "Source_A")]

    interaction_rows = conn.execute(
        "SELECT Phage_ID, Source_DB FROM private_interactions ORDER BY Phage_ID"
    ).fetchall()
    assert interaction_rows == [("PA", "Source_A")]

    host_rows = conn.execute(
        "SELECT Host_ID, source_type FROM dim_hosts WHERE source_type = 'private' ORDER BY Host_ID"
    ).fetchall()
    assert host_rows == [("HA", "private")]


# ── Phage-only source tests ─────────────────────────────────────────────────


def test_validate_private_source_phage_only_without_hosts_dir_is_valid(tmp_path):
    """Phage-only sources with Host_ID='unknown' and no hosts/ should be valid."""
    source = tmp_path / "PhageOnly_Project"
    source.mkdir()
    pd.DataFrame(
        [
            {
                "Phage_ID": "P1",
                "Host_ID": "unknown",
                "Host_name": "unknown",
                "Source_DB": "PhageOnly_Project",
                "interaction": "virulent",
            },
            {
                "Phage_ID": "P2",
                "Host_ID": "unknown",
                "Host_name": "unknown",
                "Source_DB": "PhageOnly_Project",
                "interaction": "temperate",
            },
        ]
    ).to_csv(source / "metadata.csv", index=False)
    _write_text(source / "phage.fasta", ">P1 desc\nATGC\n>P2 desc\nGCTA\n")
    # No hosts/ directory — this should be fine for phage-only sources

    result = validate_private_source(source)
    assert result.is_valid
    assert result.stats["unique_phages"] == 2
    assert result.stats["unique_hosts"] == 1  # single "unknown" host


def test_validate_private_source_mixed_host_ids_without_hosts_dir_is_invalid(tmp_path):
    """Mixing real Host_IDs with 'unknown' without hosts/ should be invalid."""
    source = tmp_path / "Mixed_Project"
    source.mkdir()
    pd.DataFrame(
        [
            {
                "Phage_ID": "P1",
                "Host_ID": "H1",
                "Host_name": "Real Host",
                "Source_DB": "Mixed_Project",
                "interaction": "virulent",
            },
            {
                "Phage_ID": "P2",
                "Host_ID": "unknown",
                "Host_name": "unknown",
                "Source_DB": "Mixed_Project",
                "interaction": "temperate",
            },
        ]
    ).to_csv(source / "metadata.csv", index=False)
    _write_text(source / "phage.fasta", ">P1 desc\nATGC\n>P2 desc\nGCTA\n")
    # No hosts/ directory — invalid because H1 is a real host ID

    result = validate_private_source(source)
    assert not result.is_valid
    assert any("Missing required host sequences" in e for e in result.errors)


def test_prepare_private_sequence_artifacts_phage_only(tmp_path):
    """Phage-only sources should produce phage mapping but no host mapping."""
    source = tmp_path / "PhageOnly_Source"
    source.mkdir()
    pd.DataFrame(
        [
            {
                "Phage_ID": "P1",
                "Host_ID": "unknown",
                "Host_name": "unknown",
                "Source_DB": "PhageOnly_Source",
                "interaction": "virulent",
            }
        ]
    ).to_csv(source / "metadata.csv", index=False)
    _write_text(source / "phage.fasta", ">P1 desc\nATGC\n")

    manifest = {
        "sources": [
            {
                "source_db": "PhageOnly_Source",
                "source_dir": str(source),
                "is_valid": True,
            }
        ]
    }

    private_phage_dir = tmp_path / "private" / "phages"
    private_phage_mapping = tmp_path / "private" / "private_phage_mapping.json"
    private_host_mapping = tmp_path / "private" / "private_host_mapping.json"

    stats = prepare_private_sequence_artifacts(
        manifest=manifest,
        private_phage_dir=private_phage_dir,
        private_phage_mapping_path=private_phage_mapping,
        private_host_mapping_path=private_host_mapping,
    )

    assert stats["sources_processed"] == 1
    assert stats["phages_written"] == 1
    assert stats["hosts_written"] == 0

    # Phage mapping should exist and point to the copied FASTA
    with private_phage_mapping.open("r", encoding="utf-8") as handle:
        pmapping = json.load(handle)
    assert "PhageOnly_Source" in pmapping
    phage_fasta_path = Path(pmapping["PhageOnly_Source"])
    assert phage_fasta_path.exists()
    assert ">P1 desc" in phage_fasta_path.read_text(encoding="utf-8")

    # Host mapping should be empty
    assert json.loads(private_host_mapping.read_text(encoding="utf-8")) == {}


def test_ingest_private_sources_phage_only(tmp_path):
    """Phage-only sources should ingest phages but not create dim_hosts entries for 'unknown'."""
    source = tmp_path / "PhageOnly_Ingest"
    source.mkdir()
    pd.DataFrame(
        [
            {
                "Phage_ID": "P1",
                "Host_ID": "unknown",
                "Host_name": "unknown",
                "Source_DB": "PhageOnly_Ingest",
                "interaction": "virulent",
            }
        ]
    ).to_csv(source / "metadata.csv", index=False)
    _write_text(source / "phage.fasta", ">P1 desc\nATGC\n")

    conn = duckdb.connect(":memory:")
    _prepare_minimal_db(conn)

    summary = ingest_private_sources_into_db(conn, [str(source)])
    assert len(summary["ingested"]) == 1
    assert summary["ingested"][0]["unique_phages"] == 1

    # Phage should be in fact_phages
    phage_count = conn.execute("SELECT COUNT(*) FROM fact_phages WHERE source_type = 'private'").fetchone()[0]
    assert phage_count == 1

    # Interaction should exist in private_interactions
    interaction_rows = conn.execute(
        "SELECT Phage_ID, Host_ID FROM private_interactions ORDER BY Phage_ID"
    ).fetchall()
    assert interaction_rows == [("P1", "unknown")]

    # dim_hosts should NOT contain an "unknown" entry
    host_rows = conn.execute(
        "SELECT Host_ID FROM dim_hosts WHERE source_type = 'private'"
    ).fetchall()
    assert host_rows == []


# ── Duplicate detection tests ─────────────────────────────────────────────────


def test_ensure_duplicate_columns_added():
    """Ingestion should add duplicate detection columns to fact_phages."""
    source = _create_private_source(
        Path("/tmp/_test_dup_cols"),
        "DupCols",
        [{"Phage_ID": "P1", "Host_ID": "H1", "Host_name": "Host",
          "Source_DB": "DupCols", "interaction": "virulent"}],
        ["P1"],
        ["H1"],
    )
    conn = duckdb.connect(":memory:")
    _prepare_minimal_db(conn)

    ingest_private_sources_into_db(conn, [str(source)])

    columns = {row[1] for row in conn.execute("PRAGMA table_info(fact_phages)").fetchall()}
    assert "is_duplicate_of_public" in columns
    assert "duplicate_public_id" in columns
    assert "duplicate_pident" in columns
    assert "duplicate_qcovs" in columns


def test_private_phages_default_not_duplicate():
    """Private phages should default to is_duplicate_of_public = FALSE."""
    source = _create_private_source(
        Path("/tmp/_test_dup_default"),
        "DupDefault",
        [{"Phage_ID": "P1", "Host_ID": "H1", "Host_name": "Host",
          "Source_DB": "DupDefault", "interaction": "virulent"}],
        ["P1"],
        ["H1"],
    )
    conn = duckdb.connect(":memory:")
    _prepare_minimal_db(conn)

    ingest_private_sources_into_db(conn, [str(source)])

    row = conn.execute(
        "SELECT is_duplicate_of_public, duplicate_public_id FROM fact_phages "
        "WHERE Phage_ID = 'P1'"
    ).fetchone()
    assert row[0] is False
    assert row[1] is None


def test_duplicate_check_with_no_blast_db(tmp_path):
    """When blast_db_dir is None, duplicate check is skipped gracefully."""
    source = _create_private_source(
        tmp_path / "src",
        "NoBlast",
        [{"Phage_ID": "P1", "Host_ID": "H1", "Host_name": "Host",
          "Source_DB": "NoBlast", "interaction": "virulent"}],
        ["P1"],
        ["H1"],
    )
    conn = duckdb.connect(":memory:")
    _prepare_minimal_db(conn)

    summary = ingest_private_sources_into_db(
        conn, [str(source)], blast_db_dir=None, private_phage_mapping_path=None
    )
    assert len(summary["ingested"]) == 1
    dup_stats = summary.get("duplicate_check", {})
    assert dup_stats.get("checked", 0) == 0


def test_duplicate_check_with_missing_blast_db(tmp_path):
    """When blast_db_dir points to non-existent path, duplicate check is skipped."""
    source = _create_private_source(
        tmp_path / "src",
        "MissingBlast",
        [{"Phage_ID": "P1", "Host_ID": "H1", "Host_name": "Host",
          "Source_DB": "MissingBlast", "interaction": "virulent"}],
        ["P1"],
        ["H1"],
    )
    conn = duckdb.connect(":memory:")
    _prepare_minimal_db(conn)

    summary = ingest_private_sources_into_db(
        conn,
        [str(source)],
        blast_db_dir=str(tmp_path / "nonexistent_blast_db"),
        private_phage_mapping_path=str(tmp_path / "nonexistent_mapping.json"),
    )
    assert len(summary["ingested"]) == 1
    dup_stats = summary.get("duplicate_check", {})
    assert dup_stats.get("checked", 0) == 0


def test_duplicate_check_empty_private_data(tmp_path):
    """When no private phage mapping exists, duplicate check returns empty stats."""
    conn = duckdb.connect(":memory:")
    _prepare_minimal_db(conn)

    from pbi.private_data import check_private_duplicates_against_public

    stats = check_private_duplicates_against_public(
        conn=conn,
        blast_db_dir=tmp_path / "blast_db",
        private_phage_mapping_path=tmp_path / "nonexistent_mapping.json",
    )
    assert stats["checked"] == 0
    assert stats["duplicates_found"] == 0


def test_blast_searcher_lists_private_and_combined():
    """BlastSearcher should list private and combined as valid databases."""
    from pbi.blast_search import BlastSearcher

    # Create a minimal BlastSearcher with a temp dir
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        searcher = BlastSearcher(tmpdir)
        dbs = searcher.list_databases()
        assert "private" in dbs
        assert "combined" in dbs
        assert dbs["private"]["type"] == "nucl"
        assert dbs["combined"]["type"] == "nucl"


def test_blast_searcher_get_db_prefix_private_combined():
    """BlastSearcher.get_db_prefix should work for private and combined."""
    from pbi.blast_search import BlastSearcher

    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        searcher = BlastSearcher(tmpdir)
        # Create done markers so get_db_prefix finds them
        for db_name in ("private", "combined"):
            db_dir = Path(tmpdir) / db_name
            db_dir.mkdir(parents=True, exist_ok=True)
            (db_dir / f"makeblastdb_{db_name}.done").touch()
        prefix = searcher.get_db_prefix("private")
        assert "all_private" in prefix
        prefix = searcher.get_db_prefix("combined")
        assert "all_combined" in prefix
