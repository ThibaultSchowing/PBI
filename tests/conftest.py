"""
Shared pytest fixtures for the PBI test suite.

Provides in-memory DuckDB databases, temporary FASTA files, and
pre-configured SequenceRetriever instances for unit and integration tests.
"""

import json
import textwrap
from pathlib import Path

import duckdb
import pytest


# ---------------------------------------------------------------------------
# DuckDB helpers
# ---------------------------------------------------------------------------

def _create_test_db(conn: duckdb.DuckDBPyConnection) -> None:
    """Populate an in-memory DuckDB with minimal test data."""
    conn.execute("""
        CREATE TABLE fact_phages (
            Phage_ID VARCHAR,
            Source_DB VARCHAR,
            source_type VARCHAR,
            GC_content DOUBLE,
            Length INTEGER,
            Taxonomy VARCHAR,
            Completeness VARCHAR,
            Host VARCHAR,
            Lifestyle VARCHAR,
            Cluster VARCHAR,
            Subcluster VARCHAR
        )
    """)
    conn.execute("""
        CREATE TABLE dim_hosts (
            Host_ID VARCHAR,
            Species_Name VARCHAR,
            GC_Content DOUBLE,
            Genome_Length INTEGER,
            Assembly_Level VARCHAR,
            RefSeq_Category VARCHAR
        )
    """)
    conn.execute("""
        CREATE TABLE phage_host_associations (
            Phage_ID VARCHAR,
            Host_ID VARCHAR
        )
    """)
    conn.execute("""
        CREATE TABLE dim_proteins (
            Protein_ID VARCHAR,
            Phage_ID VARCHAR,
            Source_DB VARCHAR,
            Protein_Length INTEGER,
            Annotation VARCHAR
        )
    """)

    # Insert test phages
    conn.execute("""
        INSERT INTO fact_phages VALUES
        ('phage_001', 'RefSeq', 'public', 45.0, 50000, 'Caudovirales', 'Complete', 'Escherichia coli', 'Lytic', 'ClusterA', 'Sub1'),
        ('phage_002', 'RefSeq', 'public', 52.0, 35000, 'Caudovirales', 'Draft', 'Bacillus subtilis', 'Temperate', 'ClusterB', 'Sub2'),
        ('phage_003', 'PhagesDB', 'public', 38.0, 75000, 'Myoviridae', 'Complete', 'Escherichia coli', 'Lytic', 'ClusterC', 'Sub3')
    """)

    # Insert test hosts
    conn.execute("""
        INSERT INTO dim_hosts VALUES
        ('host_001', 'Escherichia coli', 50.8, 4600000, 'Complete Genome', 'refseq'),
        ('host_002', 'Bacillus subtilis', 43.5, 4200000, 'Chromosome', 'genbank'),
        ('host_003', 'Pseudomonas aeruginosa', 66.6, 6300000, 'Complete Genome', 'refseq')
    """)

    # Insert phage-host associations
    conn.execute("""
        INSERT INTO phage_host_associations VALUES
        ('phage_001', 'host_001'),
        ('phage_001', 'host_002'),
        ('phage_002', 'host_003'),
        ('phage_003', 'host_001')
    """)

    # Insert test proteins
    conn.execute("""
        INSERT INTO dim_proteins VALUES
        ('prot_001', 'phage_001', 'RefSeq', 300, 'hypothetical protein'),
        ('prot_002', 'phage_001', 'RefSeq', 450, 'tail fiber protein'),
        ('prot_003', 'phage_002', 'RefSeq', 200, 'capsid protein')
    """)


# ---------------------------------------------------------------------------
# FASTA file helpers
# ---------------------------------------------------------------------------

def _write_fasta(path: Path, records: dict) -> None:
    """Write a simple FASTA file.  records = {header: sequence}."""
    with open(path, "w") as f:
        for header, seq in records.items():
            f.write(f">{header}\n")
            # Wrap at 80 chars
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")


def _write_fai(path: Path, records: dict) -> None:
    """Write a minimal FASTA index (.fai) file matching *records*."""
    with open(path, "w") as f:
        offset = 0
        for header, seq in records.items():
            seq_len = len(seq)
            line_bytes = seq_len + len(header) + 2  # >header\n + seq\n
            f.write(f"{header}\t{seq_len}\t{offset}\t{60}\t{61}\n")
            offset += line_bytes


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def tmp_db(tmp_path):
    """Create a temporary DuckDB database file with test data."""
    db_path = tmp_path / "test.duckdb"
    conn = duckdb.connect(str(db_path))
    _create_test_db(conn)
    conn.close()
    return db_path


@pytest.fixture()
def in_memory_db():
    """Return an in-memory DuckDB connection populated with test data."""
    conn = duckdb.connect(":memory:")
    _create_test_db(conn)
    yield conn
    conn.close()


@pytest.fixture()
def tmp_fasta_files(tmp_path):
    """Create temporary phage, protein, and host FASTA files with indices.

    Returns a dict of Path objects.
    """
    phage_records = {
        "phage_001": "ATCGATCGATCG" * 100,   # 1300 bp
        "phage_002": "GCTAGCTAGCTA" * 80,     # 960 bp
        "phage_003": "TTTTAAAACCCC" * 200,    # 2400 bp
    }
    protein_records = {
        "prot_001": "MKTILINTLLSI" * 25,      # 300 aa
        "prot_002": "AAACCCGGGTTT" * 37 + "A", # 444+1 aa
        "prot_003": "WWKKYYEEQQR" * 18 + "W", # 198+1 aa
    }
    host_records = {
        "host_001": "ATCG" * 500,             # 2000 bp
        "host_002": "GCTA" * 400,             # 1600 bp
        "host_003": "AAAA" * 300,             # 1200 bp
    }

    phage_fasta = tmp_path / "phages.fasta"
    protein_fasta = tmp_path / "proteins.fasta"
    host_fasta = tmp_path / "hosts.fasta"

    _write_fasta(phage_fasta, phage_records)
    _write_fasta(protein_fasta, protein_records)
    _write_fasta(host_fasta, host_records)

    _write_fai(phage_fasta.with_suffix(".fasta.fai"), phage_records)
    _write_fai(protein_fasta.with_suffix(".fasta.fai"), protein_records)
    _write_fai(host_fasta.with_suffix(".fasta.fai"), host_records)

    return {
        "phage_fasta": phage_fasta,
        "protein_fasta": protein_fasta,
        "host_fasta": host_fasta,
        "phage_records": phage_records,
        "protein_records": protein_records,
        "host_records": host_records,
    }


@pytest.fixture()
def retriever(tmp_db, tmp_fasta_files):
    """Return a SequenceRetriever connected to the test DB and FASTA files."""
    from pbi.sequence_retrieval import SequenceRetriever

    r = SequenceRetriever(
        db_path=str(tmp_db),
        phage_fasta_path=str(tmp_fasta_files["phage_fasta"]),
        protein_fasta_path=str(tmp_fasta_files["protein_fasta"]),
        host_fasta_path=str(tmp_fasta_files["host_fasta"]),
        preload=False,
    )
    yield r
    r.close()


@pytest.fixture()
def retriever_with_mapping(tmp_db, tmp_path, tmp_fasta_files):
    """Return a SequenceRetriever using host mapping mode."""
    from pbi.sequence_retrieval import SequenceRetriever

    # Create individual host FASTA files
    host_dir = tmp_path / "host_fastas"
    host_dir.mkdir()

    mapping = {}
    for host_id, seq in tmp_fasta_files["host_records"].items():
        host_file = host_dir / f"{host_id}.fasta"
        _write_fasta(host_file, {host_id: seq})
        _write_fai(host_file.with_suffix(".fasta.fai"), {host_id: seq})
        mapping[host_id] = str(host_file)

    mapping_path = tmp_path / "host_mapping.json"
    mapping_path.write_text(json.dumps(mapping))

    r = SequenceRetriever(
        db_path=str(tmp_db),
        phage_fasta_path=str(tmp_fasta_files["phage_fasta"]),
        protein_fasta_path=str(tmp_fasta_files["protein_fasta"]),
        host_mapping_path=str(mapping_path),
        preload=False,
    )
    yield r
    r.close()


@pytest.fixture()
def gff3_index(tmp_path):
    """Create a temporary GFF3 index and file for testing GFF3Retriever."""
    gff3_dir = tmp_path / "gff3"
    gff3_dir.mkdir()

    # Write a GFF3 file
    gff3_content = textwrap.dedent("""\
        ##gff-version 3
        phage_001\tRefSeq\tgene\t1\t300\t.\t+\t0\tID=gene001;Name=orf001
        phage_001\tRefSeq\tCDS\t1\t300\t.\t+\t0\tID=cds001;Parent=gene001
        phage_002\tRefSeq\tgene\t1\t200\t.\t+\t0\tID=gene002;Name=orf002
        phage_002\tRefSeq\tCDS\t1\t200\t.\t+\t0\tID=cds002;Parent=gene002
    """)
    gff3_file = gff3_dir / "RefSeq.gff3"
    gff3_file.write_text(gff3_content)

    # Build index
    index = {}
    lines = gff3_content.split("\n")
    byte_offset = 0
    for line in lines:
        if line.startswith("##") or not line.strip():
            byte_offset += len(line.encode("utf-8")) + 1
            continue
        parts = line.split("\t")
        if len(parts) >= 1:
            phage_id = parts[0]
            if phage_id not in index:
                index[phage_id] = {
                    "source_db": "RefSeq",
                    "file_name": "RefSeq.gff3",
                    "byte_offset": byte_offset,
                    "byte_length": len(line.encode("utf-8")) + 1,
                }
        byte_offset += len(line.encode("utf-8")) + 1

    index_path = gff3_dir / "gff3_index.json"
    index_path.write_text(json.dumps(index))

    return {"gff3_dir": gff3_dir, "index_path": index_path, "index": index}


@pytest.fixture()
def api_client_factory(retriever):
    """Return a factory that creates a FastAPI TestClient with the test retriever."""
    from unittest.mock import patch

    def _factory():
        from fastapi.testclient import TestClient
        import api.app as app_module

        # Patch the global retriever
        with patch.object(app_module, "retriever", retriever):
            client = TestClient(app_module.app)
            yield client
            client.close()

    return _factory
