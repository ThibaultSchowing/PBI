"""
Tests for the BLAST search wrapper (pbi/blast_search.py).

Uses mocked subprocess calls to test BLAST integration without requiring
actual BLAST+ installation or database files.
"""

import os
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch, mock_open

import pandas as pd
import pytest


# ---------------------------------------------------------------------------
# BlastSearcher initialization tests
# ---------------------------------------------------------------------------

class TestBlastSearcherInit:
    def test_init_with_explicit_path(self, tmp_path):
        from pbi.blast_search import BlastSearcher
        searcher = BlastSearcher(tmp_path)
        assert searcher.blast_db_dir == tmp_path

    def test_init_with_string_path(self, tmp_path):
        from pbi.blast_search import BlastSearcher
        searcher = BlastSearcher(str(tmp_path))
        assert searcher.blast_db_dir == tmp_path

    def test_init_default_path_no_env(self, tmp_path):
        from pbi.blast_search import BlastSearcher
        with patch.dict(os.environ, {}, clear=True):
            with patch("pbi.blast_search.Path") as mock_path:
                # Mock the parent chain to get a valid path
                mock_path.return_value.parent.parent.parent = tmp_path
                searcher = BlastSearcher()
                assert searcher.blast_db_dir is not None

    def test_init_with_data_path_env(self, tmp_path):
        from pbi.blast_search import BlastSearcher
        blast_dir = tmp_path / "blast_db"
        with patch.dict(os.environ, {"DATA_PATH": str(tmp_path)}):
            searcher = BlastSearcher()
            assert searcher.blast_db_dir == tmp_path / "blast_db"


# ---------------------------------------------------------------------------
# Database listing tests
# ---------------------------------------------------------------------------

class TestListDatabases:
    def test_list_databases_empty(self, tmp_path):
        from pbi.blast_search import BlastSearcher
        searcher = BlastSearcher(tmp_path)
        dbs = searcher.list_databases()
        assert "phages" in dbs
        assert "proteins" in dbs
        assert "hosts" in dbs
        for name, info in dbs.items():
            assert info["exists"] is False
            assert info["type"] in ("nucl", "prot")

    def test_list_databases_with_done_marker(self, tmp_path):
        from pbi.blast_search import BlastSearcher
        # Create done marker
        db_dir = tmp_path / "phages"
        db_dir.mkdir()
        (db_dir / "makeblastdb_phages.done").touch()

        searcher = BlastSearcher(tmp_path)
        dbs = searcher.list_databases()
        assert dbs["phages"]["exists"] is True
        assert dbs["proteins"]["exists"] is False

    def test_list_databases_with_nsq_file(self, tmp_path):
        from pbi.blast_search import BlastSearcher
        # Create .nsq file (nucleotide database)
        db_dir = tmp_path / "phages"
        db_dir.mkdir()
        (db_dir / "all_phages.nsq").touch()

        searcher = BlastSearcher(tmp_path)
        dbs = searcher.list_databases()
        assert dbs["phages"]["exists"] is True

    def test_list_databases_with_psq_file(self, tmp_path):
        from pbi.blast_search import BlastSearcher
        # Create .psq file (protein database)
        db_dir = tmp_path / "proteins"
        db_dir.mkdir()
        (db_dir / "all_proteins.psq").touch()

        searcher = BlastSearcher(tmp_path)
        dbs = searcher.list_databases()
        assert dbs["proteins"]["exists"] is True


# ---------------------------------------------------------------------------
# Database prefix tests
# ---------------------------------------------------------------------------

class TestGetDbPrefix:
    def test_valid_db(self, tmp_path):
        from pbi.blast_search import BlastSearcher
        # Create done marker
        db_dir = tmp_path / "phages"
        db_dir.mkdir()
        (db_dir / "makeblastdb_phages.done").touch()

        searcher = BlastSearcher(tmp_path)
        prefix = searcher.get_db_prefix("phages")
        assert prefix == str(tmp_path / "phages" / "all_phages")

    def test_invalid_db_name(self, tmp_path):
        from pbi.blast_search import BlastSearcher
        searcher = BlastSearcher(tmp_path)
        with pytest.raises(ValueError, match="Invalid database"):
            searcher.get_db_prefix("invalid_db")

    def test_missing_db_raises(self, tmp_path):
        from pbi.blast_search import BlastSearcher
        searcher = BlastSearcher(tmp_path)
        with pytest.raises(FileNotFoundError, match="BLAST database"):
            searcher.get_db_prefix("phages")


# ---------------------------------------------------------------------------
# Sequence search tests (mocked subprocess)
# ---------------------------------------------------------------------------

class TestSearchSequence:
    @patch("pbi.blast_search.BlastSearcher._run_blast")
    def test_search_sequence_calls_run_blast(self, mock_run, tmp_path):
        from pbi.blast_search import BlastSearcher
        mock_run.return_value = pd.DataFrame(columns=[
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore",
            "qlen", "slen", "qcovs", "scovhsp",
        ])
        # Create done marker
        db_dir = tmp_path / "phages"
        db_dir.mkdir()
        (db_dir / "makeblastdb_phages.done").touch()

        searcher = BlastSearcher(tmp_path)
        result = searcher.search_sequence("ATGCGT", program="blastn", db="phages")
        assert isinstance(result, pd.DataFrame)
        mock_run.assert_called_once()

    def test_search_sequence_invalid_program(self, tmp_path):
        from pbi.blast_search import BlastSearcher
        searcher = BlastSearcher(tmp_path)
        with pytest.raises(ValueError, match="Invalid program"):
            searcher.search_sequence("ATGCGT", program="invalid")

    def test_search_sequence_invalid_db(self, tmp_path):
        from pbi.blast_search import BlastSearcher
        searcher = BlastSearcher(tmp_path)
        with pytest.raises(ValueError, match="Invalid database"):
            searcher.search_sequence("ATGCGT", program="blastn", db="invalid")


# ---------------------------------------------------------------------------
# FASTA search tests (mocked subprocess)
# ---------------------------------------------------------------------------

class TestSearchFasta:
    @patch("pbi.blast_search.BlastSearcher._run_blast")
    def test_search_fasta_calls_run_blast(self, mock_run, tmp_path):
        from pbi.blast_search import BlastSearcher
        mock_run.return_value = pd.DataFrame(columns=[
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore",
            "qlen", "slen", "qcovs", "scovhsp",
        ])
        # Create done marker
        db_dir = tmp_path / "phages"
        db_dir.mkdir()
        (db_dir / "makeblastdb_phages.done").touch()

        # Create query FASTA
        query_file = tmp_path / "query.fasta"
        query_file.write_text(">query\nATGCGT\n")

        searcher = BlastSearcher(tmp_path)
        result = searcher.search_fasta(query_file, program="blastn", db="phages")
        assert isinstance(result, pd.DataFrame)

    def test_search_fasta_missing_file(self, tmp_path):
        from pbi.blast_search import BlastSearcher
        searcher = BlastSearcher(tmp_path)
        with pytest.raises(FileNotFoundError, match="Query FASTA not found"):
            searcher.search_fasta(tmp_path / "nonexistent.fasta")

    def test_search_fasta_invalid_program(self, tmp_path):
        from pbi.blast_search import BlastSearcher
        searcher = BlastSearcher(tmp_path)
        query_file = tmp_path / "query.fasta"
        query_file.write_text(">query\nATGCGT\n")
        with pytest.raises(ValueError, match="Invalid program"):
            searcher.search_fasta(query_file, program="bad")


# ---------------------------------------------------------------------------
# Query hash tests
# ---------------------------------------------------------------------------

class TestQueryHash:
    def test_same_input_same_hash(self):
        from pbi.blast_search import BlastSearcher
        h1 = BlastSearcher.query_hash("ATGCGT", "blastn", "phages", 1e-5)
        h2 = BlastSearcher.query_hash("ATGCGT", "blastn", "phages", 1e-5)
        assert h1 == h2

    def test_different_input_different_hash(self):
        from pbi.blast_search import BlastSearcher
        h1 = BlastSearcher.query_hash("ATGCGT", "blastn", "phages", 1e-5)
        h2 = BlastSearcher.query_hash("TTTAAA", "blastn", "phages", 1e-5)
        assert h1 != h2

    def test_hash_is_string(self):
        from pbi.blast_search import BlastSearcher
        h = BlastSearcher.query_hash("ATGCGT", "blastn", "phages", 1e-5)
        assert isinstance(h, str)
        assert len(h) == 16


# ---------------------------------------------------------------------------
# Source DB extraction tests
# ---------------------------------------------------------------------------

class TestExtractSourceDb:
    def test_known_source(self):
        from pbi.blast_search import _extract_source_db
        assert _extract_source_db("refseq_phage_123") == "RefSeq"
        assert _extract_source_db("genbank_phage_456") == "GenBank"
        assert _extract_source_db("PhagesDB_phage_789") == "PhagesDB"

    def test_unknown_source(self):
        from pbi.blast_search import _extract_source_db
        assert _extract_source_db("XYZ_unknown_123") == "Unknown"

    def test_empty_string(self):
        from pbi.blast_search import _extract_source_db
        assert _extract_source_db("") == "Unknown"


# ---------------------------------------------------------------------------
# _run_blast tests (mocked subprocess)
# ---------------------------------------------------------------------------

class TestRunBlast:
    @patch("subprocess.run")
    def test_run_blast_success(self, mock_run, tmp_path):
        from pbi.blast_search import BlastSearcher
        # Create done marker
        db_dir = tmp_path / "phages"
        db_dir.mkdir()
        (db_dir / "makeblastdb_phages.done").touch()

        # Create query FASTA
        query = tmp_path / "query.fasta"
        query.write_text(">query\nATGCGT\n")

        # Mock subprocess result
        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = "query\trefseq_1\t95.5\t6\t0\t0\t1\t6\t100\t105\t1e-10\t120\n"
        mock_run.return_value = mock_result

        searcher = BlastSearcher(tmp_path)
        with patch.object(searcher, "_blast_bin", return_value="/usr/bin/blastn"):
            df = searcher._run_blast(
                str(query), "blastn", "phages",
                max_hits=10, evalue=1e-5, outfmt=6,
            )
        assert len(df) == 1
        assert df.iloc[0]["qseqid"] == "query"

    @patch("subprocess.run")
    def test_run_blast_empty_output(self, mock_run, tmp_path):
        from pbi.blast_search import BlastSearcher
        # Create done marker
        db_dir = tmp_path / "phages"
        db_dir.mkdir()
        (db_dir / "makeblastdb_phages.done").touch()

        query = tmp_path / "query.fasta"
        query.write_text(">query\nATGCGT\n")

        mock_result = MagicMock()
        mock_result.returncode = 0
        mock_result.stdout = ""
        mock_run.return_value = mock_result

        searcher = BlastSearcher(tmp_path)
        with patch.object(searcher, "_blast_bin", return_value="/usr/bin/blastn"):
            df = searcher._run_blast(
                str(query), "blastn", "phages",
                max_hits=10, evalue=1e-5, outfmt=6,
            )
        assert len(df) == 0

    @patch("subprocess.run")
    def test_run_blast_failure_raises(self, mock_run, tmp_path):
        from pbi.blast_search import BlastSearcher
        db_dir = tmp_path / "phages"
        db_dir.mkdir()
        (db_dir / "makeblastdb_phages.done").touch()

        query = tmp_path / "query.fasta"
        query.write_text(">query\nATGCGT\n")

        mock_result = MagicMock()
        mock_result.returncode = 1
        mock_result.stderr = "Error: database not found"
        mock_result.stdout = ""
        mock_run.return_value = mock_result

        searcher = BlastSearcher(tmp_path)
        with patch.object(searcher, "_blast_bin", return_value="/usr/bin/blastn"):
            with pytest.raises(RuntimeError, match="BLAST failed"):
                searcher._run_blast(
                    str(query), "blastn", "phages",
                    max_hits=10, evalue=1e-5, outfmt=6,
                )


# ---------------------------------------------------------------------------
# Lazy import tests
# ---------------------------------------------------------------------------

class TestLazyImport:
    def test_import_blast_searcher(self):
        from pbi import BlastSearcher
        assert BlastSearcher is not None

    def test_getattr_blast_searcher(self):
        from pbi import __getattr__
        cls = __getattr__("BlastSearcher")
        assert cls is not None
