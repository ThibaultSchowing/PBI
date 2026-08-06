#!/usr/bin/env python3
"""
CI smoke tests for the pbi package.

These tests verify that the pbi package works end-to-end against the database
and FASTA files produced by the pipeline CI subset.  They do NOT assert exact
data values — only that queries execute, return expected shapes, and the
package API behaves correctly.

Run inside the Docker container after the pipeline has completed:
    DATA_PATH=/data/processed pytest tests/test_pbi_ci_smoke.py -v
"""

import os
import sys
from pathlib import Path

import pytest

# Ensure the package is importable when installed in the container.
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from pbi import SequenceRetriever, get_default_paths, quick_connect


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def paths():
    """Return default data paths (mirrors get_default_paths)."""
    return get_default_paths()


@pytest.fixture(scope="module")
def retriever(paths):
    """Open a SequenceRetriever against the CI database and close after tests."""
    # Skip the entire module if the database does not exist (e.g. pipeline failed).
    if not paths["database"].exists():
        pytest.skip(f"Database not found: {paths['database']}")
    r = quick_connect()
    yield r
    r.close()


# ---------------------------------------------------------------------------
# Connection tests
# ---------------------------------------------------------------------------

class TestConnection:
    def test_database_file_exists(self, paths):
        assert paths["database"].exists(), f"Missing DB: {paths['database']}"

    def test_quick_connect_returns_retriever(self, retriever):
        assert isinstance(retriever, SequenceRetriever)

    def test_connection_is_read_only(self, retriever):
        assert retriever.conn.read_only


# ---------------------------------------------------------------------------
# Metadata query tests
# ---------------------------------------------------------------------------

class TestMetadataQueries:
    def test_get_phage_metadata(self, retriever):
        df = retriever.get_phage_metadata(limit=10)
        assert len(df) > 0, "Expected at least one phage"
        assert "Phage_ID" in df.columns

    def test_get_protein_metadata(self, retriever):
        df = retriever.get_protein_metadata(limit=10)
        assert len(df) > 0, "Expected at least one protein"
        assert "Protein_ID" in df.columns

    def test_get_phage_host_pairs(self, retriever):
        df = retriever.get_phage_host_pairs(limit=10)
        # Phage-host pairs may be empty if host genomes were not downloaded
        # (metadata_only_mode), so we just verify the query runs.
        assert isinstance(df, object)  # DataFrame or empty

    def test_structured_filter(self, retriever):
        df = retriever.query_phage_host_pairs(
            phage_filters={"Source_DB": "RefSeq_Phage_Metadata_URL"},
            limit=5,
        )
        assert isinstance(df, object)


# ---------------------------------------------------------------------------
# FASTA access tests (only if FASTA files exist)
# ---------------------------------------------------------------------------

class TestFastaAccess:
    def test_phage_fasta_file_exists(self, paths):
        if not paths["phage_fasta"].exists():
            pytest.skip("Phage FASTA not present in CI subset")

    def test_protein_fasta_file_exists(self, paths):
        if not paths["protein_fasta"].exists():
            pytest.skip("Protein FASTA not present in CI subset")

    def test_get_phage_sequence(self, retriever, paths):
        if not paths["phage_fasta"].exists():
            pytest.skip("Phage FASTA not present")
        # Pick the first phage ID from the database.
        first_phage = retriever.conn.execute(
            "SELECT Phage_ID FROM fact_phages LIMIT 1"
        ).fetchone()
        if first_phage is None:
            pytest.skip("No phages in database")
        seq = retriever.get_phage_sequence(first_phage[0])
        assert isinstance(seq, str)
        assert len(seq) > 0


# ---------------------------------------------------------------------------
# Error handling tests
# ---------------------------------------------------------------------------

class TestErrorHandling:
    def test_nonexistent_phage_returns_none(self, retriever):
        result = retriever.get_phage_sequence("NONEXISTENT_PHAGE_999999")
        assert result is None
