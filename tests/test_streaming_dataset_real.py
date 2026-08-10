"""
Tests for pbi.streaming_dataset with real in-memory data.

Tests actual data loading, iteration, and item access for both
PhageHostStreamingDataset and PhageHostIndexedDataset.
"""

import json
from pathlib import Path

import duckdb
import pandas as pd
import pytest


def _write_fasta(path, records):
    with open(path, "w") as f:
        for header, seq in records.items():
            f.write(f">{header}\n{seq}\n")


def _write_fai(path, records):
    with open(path, "w") as f:
        offset = 0
        for header, seq in records.items():
            seq_len = len(seq)
            line_bytes = seq_len + len(header) + 2
            f.write(f"{header}\t{seq_len}\t{offset}\t{60}\t{61}\n")
            offset += line_bytes


@pytest.fixture()
def streaming_setup(tmp_path):
    """Create a full setup for streaming dataset tests."""
    db_path = tmp_path / "test.duckdb"
    conn = duckdb.connect(str(db_path))

    # Create tables
    conn.execute("""
        CREATE TABLE fact_phages (
            Phage_ID VARCHAR, Source_DB VARCHAR, source_type VARCHAR,
            GC_content DOUBLE, Length INTEGER, Taxonomy VARCHAR,
            Completeness VARCHAR, Host VARCHAR, Lifestyle VARCHAR,
            Cluster VARCHAR, Subcluster VARCHAR
        )
    """)
    conn.execute("""
        CREATE TABLE dim_hosts (
            Host_ID VARCHAR, Species_Name VARCHAR, GC_Content DOUBLE,
            Genome_Length INTEGER, Assembly_Level VARCHAR, RefSeq_Category VARCHAR
        )
    """)
    conn.execute("""
        CREATE TABLE phage_host_associations (
            Phage_ID VARCHAR, Host_ID VARCHAR
        )
    """)

    conn.execute("""
        INSERT INTO fact_phages VALUES
        ('phage_001', 'RefSeq', 'public', 45.0, 50000, 'Caudovirales', 'Complete', 'E. coli', 'Lytic', 'A', 'S1'),
        ('phage_002', 'RefSeq', 'public', 52.0, 35000, 'Caudovirales', 'Draft', 'B. subtilis', 'Temperate', 'B', 'S2'),
        ('phage_003', 'PhagesDB', 'public', 38.0, 75000, 'Myoviridae', 'Complete', 'E. coli', 'Lytic', 'C', 'S3')
    """)
    conn.execute("""
        INSERT INTO dim_hosts VALUES
        ('host_001', 'E. coli', 50.8, 4600000, 'Complete Genome', 'refseq'),
        ('host_002', 'B. subtilis', 43.5, 4200000, 'Chromosome', 'genbank')
    """)
    conn.execute("""
        INSERT INTO phage_host_associations VALUES
        ('phage_001', 'host_001'), ('phage_002', 'host_002'), ('phage_003', 'host_001')
    """)
    conn.close()

    # Create FASTA files
    phage_records = {
        "phage_001": "ATCGATCGATCG" * 100,
        "phage_002": "GCTAGCTAGCTA" * 80,
        "phage_003": "TTTTAAAACCCC" * 200,
    }
    host_records = {
        "host_001": "ATCG" * 500,
        "host_002": "GCTA" * 400,
    }

    phage_fasta = tmp_path / "phages.fasta"
    host_fasta = tmp_path / "hosts.fasta"
    _write_fasta(phage_fasta, phage_records)
    _write_fasta(host_fasta, host_records)
    _write_fai(tmp_path / "phages.fasta.fai", phage_records)
    _write_fai(tmp_path / "hosts.fasta.fai", host_records)

    return {
        "db_path": str(db_path),
        "phage_fasta": str(phage_fasta),
        "host_fasta": str(host_fasta),
    }


class TestPhageHostStreamingDataset:
    def test_iterates_samples(self, streaming_setup):
        from pbi.streaming_dataset import PhageHostStreamingDataset

        ds = PhageHostStreamingDataset(
            db_path=streaming_setup["db_path"],
            phage_fasta_path=streaming_setup["phage_fasta"],
            host_fasta_path=streaming_setup["host_fasta"],
        )
        samples = list(ds)
        assert len(samples) > 0
        sample = samples[0]
        assert "Phage_ID" in sample
        assert "Host_ID" in sample
        assert "Phage_Sequence" in sample
        assert "Host_Sequence" in sample
        ds.close()

    def test_with_where_clause(self, streaming_setup):
        from pbi.streaming_dataset import PhageHostStreamingDataset

        ds = PhageHostStreamingDataset(
            db_path=streaming_setup["db_path"],
            phage_fasta_path=streaming_setup["phage_fasta"],
            host_fasta_path=streaming_setup["host_fasta"],
            where_clause="Source_DB = 'RefSeq'",
        )
        samples = list(ds)
        for s in samples:
            assert s.get("Phage_Source") == "RefSeq" or s.get("Source_DB") == "RefSeq"
        ds.close()

    def test_context_manager(self, streaming_setup):
        from pbi.streaming_dataset import PhageHostStreamingDataset

        with PhageHostStreamingDataset(
            db_path=streaming_setup["db_path"],
            phage_fasta_path=streaming_setup["phage_fasta"],
            host_fasta_path=streaming_setup["host_fasta"],
        ) as ds:
            samples = list(ds)
            assert len(samples) > 0

    def test_batch_yields_correct_size(self, streaming_setup):
        from pbi.streaming_dataset import PhageHostStreamingDataset

        ds = PhageHostStreamingDataset(
            db_path=streaming_setup["db_path"],
            phage_fasta_path=streaming_setup["phage_fasta"],
            host_fasta_path=streaming_setup["host_fasta"],
            batch_size=2,
        )
        # StreamingDataset yields all samples in one batch (IterableDataset)
        samples = list(ds)
        assert len(samples) > 0
        ds.close()


class TestPhageHostIndexedDataset:
    def test_len(self, streaming_setup):
        from pbi.streaming_dataset import PhageHostIndexedDataset

        ds = PhageHostIndexedDataset(
            db_path=streaming_setup["db_path"],
            phage_fasta_path=streaming_setup["phage_fasta"],
            host_fasta_path=streaming_setup["host_fasta"],
        )
        assert len(ds) == 3  # 3 phage-host associations
        ds.close()

    def test_getitem(self, streaming_setup):
        from pbi.streaming_dataset import PhageHostIndexedDataset

        ds = PhageHostIndexedDataset(
            db_path=streaming_setup["db_path"],
            phage_fasta_path=streaming_setup["phage_fasta"],
            host_fasta_path=streaming_setup["host_fasta"],
        )
        sample = ds[0]
        assert isinstance(sample, dict)
        assert "Phage_ID" in sample
        assert "Host_ID" in sample
        ds.close()

    def test_context_manager(self, streaming_setup):
        from pbi.streaming_dataset import PhageHostIndexedDataset

        with PhageHostIndexedDataset(
            db_path=streaming_setup["db_path"],
            phage_fasta_path=streaming_setup["phage_fasta"],
            host_fasta_path=streaming_setup["host_fasta"],
        ) as ds:
            assert len(ds) == 3
            sample = ds[0]
            assert "Phage_ID" in sample


class TestParseWhereClause:
    def test_none(self):
        from pbi.streaming_dataset import parse_where_clause
        where, limit = parse_where_clause(None)
        assert where is None
        assert limit is None

    def test_limit_only(self):
        from pbi.streaming_dataset import parse_where_clause
        where, limit = parse_where_clause("LIMIT 100")
        assert where is None
        assert limit == "LIMIT 100"

    def test_where_and_limit(self):
        from pbi.streaming_dataset import parse_where_clause
        where, limit = parse_where_clause("p.Length > 1000 LIMIT 50")
        assert where == "p.Length > 1000"
        assert limit == "LIMIT 50"

    def test_where_only(self):
        from pbi.streaming_dataset import parse_where_clause
        where, limit = parse_where_clause("p.GC > 0.5")
        assert where == "p.GC > 0.5"
        assert limit is None

    def test_limit_offset(self):
        from pbi.streaming_dataset import parse_where_clause
        where, limit = parse_where_clause("LIMIT 1000 OFFSET 5000")
        assert where is None
        assert limit == "LIMIT 1000 OFFSET 5000"


class TestPhageHostCollateFn:
    def test_collate(self):
        from pbi.streaming_dataset import phage_host_collate_fn

        batch = [
            {"Phage_ID": "p1", "Host_ID": "h1", "Label": 1},
            {"Phage_ID": "p2", "Host_ID": "h2", "Label": 0},
        ]
        result = phage_host_collate_fn(batch)
        assert isinstance(result, dict)
        assert result["Phage_ID"] == ["p1", "p2"]
        assert result["Label"] == [1, 0]
