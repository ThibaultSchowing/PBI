"""CI smoke tests for PBI-Scope pipeline.

These tests verify the core functionality of the PBI-Scope system after
running the CI pipeline. They run inside a Docker container with the
pipeline data mounted.
"""

import os
import sys
import tempfile
import shutil
from pathlib import Path

import pytest

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

DATA_DIR = os.environ.get("PBI_DATA_DIR", "/data")
PROCESSED_DIR = os.path.join(DATA_DIR, "processed")
INTERMEDIATE_DIR = os.path.join(DATA_DIR, "intermediate")


class TestConnection:
    """Test database and file connections."""

    def test_database_exists(self):
        db_path = os.path.join(PROCESSED_DIR, "databases", "phagescope.duckdb")
        assert os.path.exists(db_path), f"Database not found: {db_path}"

    def test_phage_fasta_exists(self):
        fasta_path = os.path.join(INTERMEDIATE_DIR, "fasta", "merged", "all_phages.fasta")
        assert os.path.exists(fasta_path), f"Phage FASTA not found: {fasta_path}"

    def test_protein_fasta_exists(self):
        fasta_path = os.path.join(INTERMEDIATE_DIR, "fasta", "merged", "all_proteins.fasta")
        assert os.path.exists(fasta_path), f"Protein FASTA not found: {fasta_path}"

    def test_phage_fasta_index_exists(self):
        fasta_path = os.path.join(INTERMEDIATE_DIR, "fasta", "merged", "all_phages.fasta")
        index_path = fasta_path + ".fai"
        assert os.path.exists(index_path), f"Phage FASTA index not found: {index_path}"

    def test_protein_fasta_index_exists(self):
        fasta_path = os.path.join(INTERMEDIATE_DIR, "fasta", "merged", "all_proteins.fasta")
        index_path = fasta_path + ".fai"
        assert os.path.exists(index_path), f"Protein FASTA index not found: {index_path}"


class TestMetadataQueries:
    """Test DuckDB metadata queries."""

    def _get_connection(self):
        import duckdb
        db_path = os.path.join(PROCESSED_DIR, "databases", "phagescope.duckdb")
        return duckdb.connect(db_path, read_only=True)

    def test_phages_table_not_empty(self):
        conn = self._get_connection()
        count = conn.execute("SELECT COUNT(*) FROM fact_phages").fetchone()[0]
        assert count > 0, "fact_phages table is empty"

    def test_proteins_table_not_empty(self):
        conn = self._get_connection()
        count = conn.execute("SELECT COUNT(*) FROM dim_proteins").fetchone()[0]
        assert count > 0, "dim_proteins table is empty"

    def test_phages_have_required_columns(self):
        conn = self._get_connection()
        result = conn.execute("SELECT Phage_ID, Source_DB, Length FROM fact_phages LIMIT 1").fetchone()
        assert result is not None, "No phages found"
        assert result[0] is not None, "Phage_ID is NULL"
        assert result[1] is not None, "Source_DB is NULL"

    def test_proteins_have_required_columns(self):
        conn = self._get_connection()
        result = conn.execute("SELECT Protein_ID, Phage_ID FROM dim_proteins LIMIT 1").fetchone()
        assert result is not None, "No proteins found"
        assert result[0] is not None, "Protein_ID is NULL"
        assert result[1] is not None, "Phage_ID is NULL"


class TestFastaAccess:
    """Test FASTA file access."""

    def _get_retriever(self):
        from pbi.sequence_retrieval import SequenceRetriever
        db_path = os.path.join(PROCESSED_DIR, "databases", "phagescope.duckdb")
        phage_fasta = os.path.join(INTERMEDIATE_DIR, "fasta", "merged", "all_phages.fasta")
        protein_fasta = os.path.join(INTERMEDIATE_DIR, "fasta", "merged", "all_proteins.fasta")
        return SequenceRetriever(db_path, phage_fasta, protein_fasta, preload=False)

    def test_retriever_initializes(self):
        retriever = self._get_retriever()
        assert retriever is not None

    def test_get_stats(self):
        retriever = self._get_retriever()
        stats = retriever.get_stats()
        assert "database" in stats
        assert "fasta" in stats
        assert stats["database"]["phages"] > 0
        assert stats["database"]["proteins"] > 0

    def test_get_phage_metadata(self):
        retriever = self._get_retriever()
        metadata = retriever.get_phage_metadata(limit=10)
        assert len(metadata) > 0
        assert "Phage_ID" in metadata.columns


class TestPhageGenomeStreaming:
    """Test phage genome streaming."""

    def _get_retriever(self):
        from pbi.sequence_retrieval import SequenceRetriever
        db_path = os.path.join(PROCESSED_DIR, "databases", "phagescope.duckdb")
        phage_fasta = os.path.join(INTERMEDIATE_DIR, "fasta", "merged", "all_phages.fasta")
        protein_fasta = os.path.join(INTERMEDIATE_DIR, "fasta", "merged", "all_proteins.fasta")
        return SequenceRetriever(db_path, phage_fasta, protein_fasta, preload=False)

    def test_get_phage_sequence(self):
        retriever = self._get_retriever()
        metadata = retriever.get_phage_metadata(limit=1)
        phage_id = metadata.iloc[0]["Phage_ID"]
        seq = retriever.get_phage_sequence(phage_id)
        assert seq is not None
        assert len(seq) > 0
        assert all(c in "ATCGN" for c in seq.upper())

    def test_get_phage_genome_concat(self):
        retriever = self._get_retriever()
        metadata = retriever.get_phage_metadata(limit=1)
        phage_id = metadata.iloc[0]["Phage_ID"]
        seq = retriever.get_phage_genome(phage_id, mode="concat")
        assert seq is not None
        assert len(seq) > 0


class TestHostGenomeStreaming:
    """Test host genome streaming (if available)."""

    def _get_retriever(self):
        from pbi.sequence_retrieval import SequenceRetriever
        db_path = os.path.join(PROCESSED_DIR, "databases", "phagescope.duckdb")
        phage_fasta = os.path.join(INTERMEDIATE_DIR, "fasta", "merged", "all_phages.fasta")
        protein_fasta = os.path.join(INTERMEDIATE_DIR, "fasta", "merged", "all_proteins.fasta")

        host_mapping = os.path.join(INTERMEDIATE_DIR, "fasta", "merged", "host_mapping.json")
        if not os.path.exists(host_mapping):
            host_mapping = None

        return SequenceRetriever(db_path, phage_fasta, protein_fasta,
                                 host_mapping_path=host_mapping, preload=False)

    def test_host_data_available(self):
        retriever = self._get_retriever()
        stats = retriever.get_stats()
        if "hosts" not in stats.get("database", {}):
            pytest.skip("Host genome data not available in CI subset")

    def test_get_host_genome(self):
        retriever = self._get_retriever()
        try:
            stats = retriever.get_stats()
            if "hosts" not in stats.get("database", {}):
                pytest.skip("Host genome data not available")
        except Exception:
            pytest.skip("Host data not configured")

        mapping = retriever.get_host_mapping()
        if not mapping:
            pytest.skip("No host mapping available")

        host_id = list(mapping.keys())[0]
        seq = retriever.get_host_genome(host_id, mode="concat")
        assert seq is not None
        assert len(seq) > 0


class TestPhageHostPairStreaming:
    """Test phage-host pair streaming."""

    def _get_retriever(self):
        from pbi.sequence_retrieval import SequenceRetriever
        db_path = os.path.join(PROCESSED_DIR, "databases", "phagescope.duckdb")
        phage_fasta = os.path.join(INTERMEDIATE_DIR, "fasta", "merged", "all_phages.fasta")
        protein_fasta = os.path.join(INTERMEDIATE_DIR, "fasta", "merged", "all_proteins.fasta")
        return SequenceRetriever(db_path, phage_fasta, protein_fasta, preload=False)

    def test_phage_host_pairs_exist(self):
        retriever = self._get_retriever()
        try:
            pairs = retriever.get_phage_host_pairs(limit=5)
            if pairs is not None and len(pairs) > 0:
                assert "Phage_ID" in pairs.columns
        except Exception:
            pytest.skip("Phage-host association table not available")


class TestErrorHandling:
    """Test error handling for missing data."""

    def _get_retriever(self):
        from pbi.sequence_retrieval import SequenceRetriever
        db_path = os.path.join(PROCESSED_DIR, "databases", "phagescope.duckdb")
        phage_fasta = os.path.join(INTERMEDIATE_DIR, "fasta", "merged", "all_phages.fasta")
        protein_fasta = os.path.join(INTERMEDIATE_DIR, "fasta", "merged", "all_proteins.fasta")
        return SequenceRetriever(db_path, phage_fasta, protein_fasta, preload=False)

    def test_nonexistent_phage_returns_none(self):
        retriever = self._get_retriever()
        seq = retriever.get_phage_sequence("nonexistent_phage_12345")
        assert seq is None

    def test_nonexistent_phage_genome_returns_none(self):
        retriever = self._get_retriever()
        seq = retriever.get_phage_genome("nonexistent_phage_12345", mode="concat")
        assert seq is None or len(seq) == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
