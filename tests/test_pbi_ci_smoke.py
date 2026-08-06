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
        assert isinstance(df, object)

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
        first_phage = retriever.conn.execute(
            "SELECT Phage_ID FROM fact_phages LIMIT 1"
        ).fetchone()
        if first_phage is None:
            pytest.skip("No phages in database")
        seq = retriever.get_phage_sequence(first_phage[0])
        assert isinstance(seq, str)
        assert len(seq) > 0


# ---------------------------------------------------------------------------
# Phage genome streaming tests
# ---------------------------------------------------------------------------

class TestPhageGenomeStreaming:
    def _get_any_phage_id(self, retriever):
        row = retriever.conn.execute(
            "SELECT Phage_ID FROM fact_phages LIMIT 1"
        ).fetchone()
        return row[0] if row else None

    def test_get_phage_genome_concat(self, retriever, paths):
        if not paths["phage_fasta"].exists():
            pytest.skip("Phage FASTA not present")
        phage_id = self._get_any_phage_id(retriever)
        if phage_id is None:
            pytest.skip("No phages in database")
        seq = retriever.get_phage_genome(phage_id, mode="concat")
        assert isinstance(seq, str)
        assert len(seq) > 0

    def test_get_phage_genome_first(self, retriever, paths):
        if not paths["phage_fasta"].exists():
            pytest.skip("Phage FASTA not present")
        phage_id = self._get_any_phage_id(retriever)
        if phage_id is None:
            pytest.skip("No phages in database")
        seq = retriever.get_phage_genome(phage_id, mode="first")
        assert isinstance(seq, str)
        assert len(seq) > 0

    def test_get_phage_genome_list(self, retriever, paths):
        if not paths["phage_fasta"].exists():
            pytest.skip("Phage FASTA not present")
        phage_id = self._get_any_phage_id(retriever)
        if phage_id is None:
            pytest.skip("No phages in database")
        result = retriever.get_phage_genome(phage_id, mode="list")
        assert isinstance(result, list)
        assert len(result) > 0
        assert all(isinstance(s, str) for s in result)

    def test_get_phage_genome_dict(self, retriever, paths):
        if not paths["phage_fasta"].exists():
            pytest.skip("Phage FASTA not present")
        phage_id = self._get_any_phage_id(retriever)
        if phage_id is None:
            pytest.skip("No phages in database")
        result = retriever.get_phage_genome(phage_id, mode="dict")
        assert isinstance(result, dict)
        assert len(result) > 0
        assert all(isinstance(k, str) and isinstance(v, str) for k, v in result.items())

    def test_get_phage_genome_concat_with_gap(self, retriever, paths):
        if not paths["phage_fasta"].exists():
            pytest.skip("Phage FASTA not present")
        phage_id = self._get_any_phage_id(retriever)
        if phage_id is None:
            pytest.skip("No phages in database")
        seq_no_gap = retriever.get_phage_genome(phage_id, mode="concat", gap=0)
        seq_with_gap = retriever.get_phage_genome(phage_id, mode="concat", gap=100)
        # For single-contig phages, gap has no effect; lengths should be equal.
        # For multi-contig, gap version should be longer.
        assert len(seq_with_gap) >= len(seq_no_gap)

    def test_nonexistent_phage_raises_key_error(self, retriever, paths):
        if not paths["phage_fasta"].exists():
            pytest.skip("Phage FASTA not present")
        with pytest.raises(KeyError):
            retriever.get_phage_genome("NONEXISTENT_PHAGE_999999")


# ---------------------------------------------------------------------------
# Host genome streaming tests (skipped if host data not available)
# ---------------------------------------------------------------------------

class TestHostGenomeStreaming:
    def _get_any_host_id(self, retriever):
        if not retriever._has_host_data:
            return None
        row = retriever.conn.execute(
            "SELECT Host_ID FROM dim_hosts LIMIT 1"
        ).fetchone()
        return row[0] if row else None

    def test_has_host_data(self, retriever):
        # Just verify the flag is accessible; skip if no host data.
        if not retriever._has_host_data:
            pytest.skip("Host FASTA not configured (metadata_only_mode or no host genomes)")

    def test_get_host_genome_concat(self, retriever):
        host_id = self._get_any_host_id(retriever)
        if host_id is None:
            pytest.skip("No host genomes available")
        seq = retriever.get_host_genome(host_id, mode="concat")
        assert isinstance(seq, str)
        assert len(seq) > 0

    def test_get_host_genome_first(self, retriever):
        host_id = self._get_any_host_id(retriever)
        if host_id is None:
            pytest.skip("No host genomes available")
        seq = retriever.get_host_genome(host_id, mode="first")
        assert isinstance(seq, str)
        assert len(seq) > 0

    def test_get_host_genome_list(self, retriever):
        host_id = self._get_any_host_id(retriever)
        if host_id is None:
            pytest.skip("No host genomes available")
        result = retriever.get_host_genome(host_id, mode="list")
        assert isinstance(result, list)
        assert len(result) > 0
        assert all(isinstance(s, str) for s in result)

    def test_get_host_genome_dict(self, retriever):
        host_id = self._get_any_host_id(retriever)
        if host_id is None:
            pytest.skip("No host genomes available")
        result = retriever.get_host_genome(host_id, mode="dict")
        assert isinstance(result, dict)
        assert len(result) > 0

    def test_get_host_genome_stats(self, retriever):
        host_id = self._get_any_host_id(retriever)
        if host_id is None:
            pytest.skip("No host genomes available")
        stats = retriever.get_host_genome_stats(host_id)
        assert isinstance(stats, dict)
        assert "contig_count" in stats
        assert "total_length" in stats
        assert "lengths" in stats
        assert stats["contig_count"] > 0
        assert stats["total_length"] > 0
        assert len(stats["lengths"]) == stats["contig_count"]

    def test_get_host_genome_concat_with_gap(self, retriever):
        host_id = self._get_any_host_id(retriever)
        if host_id is None:
            pytest.skip("No host genomes available")
        seq_no_gap = retriever.get_host_genome(host_id, mode="concat", gap=0)
        seq_with_gap = retriever.get_host_genome(host_id, mode="concat", gap=100)
        assert len(seq_with_gap) >= len(seq_no_gap)

    def test_nonexistent_host_raises_key_error(self, retriever):
        if not retriever._has_host_data:
            pytest.skip("Host FASTA not configured")
        with pytest.raises(KeyError):
            retriever.get_host_genome("NONEXISTENT_HOST_999999")


# ---------------------------------------------------------------------------
# Phage-host pair streaming tests (with sequences)
# ---------------------------------------------------------------------------

class TestPhageHostPairStreaming:
    def test_get_phage_host_pairs_with_sequences(self, retriever):
        df = retriever.get_phage_host_pairs(limit=5)
        if len(df) == 0:
            pytest.skip("No phage-host pairs (host genomes may not be downloaded)")
        assert "Phage_Sequence" in df.columns or "Host_Sequence" in df.columns

    def test_get_phage_host_pairs_concat_mode(self, retriever):
        df = retriever.get_phage_host_pairs(limit=5, host_contig_mode="concat")
        if len(df) == 0:
            pytest.skip("No phage-host pairs")


# ---------------------------------------------------------------------------
# Error handling tests
# ---------------------------------------------------------------------------

class TestErrorHandling:
    def test_nonexistent_phage_returns_none(self, retriever):
        result = retriever.get_phage_sequence("NONEXISTENT_PHAGE_999999")
        assert result is None
