"""
Tests for pbi.gff3_retrieval.GFF3Retriever.

Uses tmp_path to create mock GFF3 files and indices.
"""

import json
import textwrap
from pathlib import Path

import pytest


@pytest.fixture()
def gff3_setup(tmp_path):
    """Create a full GFF3Retriever test setup with two source databases."""
    gff3_dir = tmp_path / "gff3"
    gff3_dir.mkdir()

    # Write RefSeq GFF3
    refseq_content = textwrap.dedent("""\
        ##gff-version 3
        phage_A\tRefSeq\tgene\t1\t300\t.\t+\t0\tID=geneA1
        phage_A\tRefSeq\tCDS\t1\t300\t.\t+\t0\tID=cdsA1
        phage_B\tRefSeq\tgene\t1\t200\t.\t+\t0\tID=geneB1
    """)
    (gff3_dir / "RefSeq.gff3").write_text(refseq_content)

    # Write PhagesDB GFF3
    phagesdb_content = textwrap.dedent("""\
        ##gff-version 3
        phage_C\tPhagesDB\tgene\t1\t150\t.\t+\t0\tID=geneC1
    """)
    (gff3_dir / "PhagesDB.gff3").write_text(phagesdb_content)

    # Build index (simple: entire file per phage)
    index = {}
    byte_offset = 0
    for fname, source in [("RefSeq.gff3", "RefSeq"), ("PhagesDB.gff3", "PhagesDB")]:
        content = (gff3_dir / fname).read_text()
        for line in content.split("\n"):
            if line.startswith("##") or not line.strip():
                byte_offset += len(line.encode("utf-8")) + 1
                continue
            phage_id = line.split("\t")[0]
            if phage_id not in index:
                index[phage_id] = {
                    "source_db": source,
                    "file_name": fname,
                    "byte_offset": byte_offset,
                    "byte_length": len(line.encode("utf-8")) + 1,
                }
            byte_offset += len(line.encode("utf-8")) + 1

    index_path = gff3_dir / "gff3_index.json"
    index_path.write_text(json.dumps(index))

    from pbi.gff3_retrieval import GFF3Retriever
    retriever = GFF3Retriever(str(gff3_dir), str(index_path))
    return {"retriever": retriever, "gff3_dir": gff3_dir, "index": index}


class TestGFF3RetrieverInit:
    def test_loads_index_lazily(self, gff3_setup):
        r = gff3_setup["retriever"]
        assert r._index is None  # not loaded yet
        _ = r.index
        assert r._index is not None

    def test_missing_index_raises(self, tmp_path):
        from pbi.gff3_retrieval import GFF3Retriever
        r = GFF3Retriever(str(tmp_path), str(tmp_path / "nonexistent.json"))
        with pytest.raises(FileNotFoundError):
            _ = r.index


class TestGetGFF3:
    def test_returns_content(self, gff3_setup):
        r = gff3_setup["retriever"]
        content = r.get_gff3("phage_A")
        assert isinstance(content, str)
        assert len(content) > 0

    def test_missing_phage_returns_empty(self, gff3_setup):
        r = gff3_setup["retriever"]
        content = r.get_gff3("NONEXISTENT_PHAGE")
        assert content == ""


class TestGetGFF3Lines:
    def test_yields_lines(self, gff3_setup):
        r = gff3_setup["retriever"]
        lines = list(r.get_gff3_lines("phage_A"))
        assert len(lines) > 0
        assert all(isinstance(line, str) for line in lines)

    def test_missing_phage_yields_nothing(self, gff3_setup):
        r = gff3_setup["retriever"]
        lines = list(r.get_gff3_lines("NONEXISTENT_PHAGE"))
        assert lines == []


class TestListPhages:
    def test_all_phages(self, gff3_setup):
        r = gff3_setup["retriever"]
        phages = r.list_phages()
        assert "phage_A" in phages
        assert "phage_B" in phages
        assert "phage_C" in phages

    def test_filter_by_source(self, gff3_setup):
        r = gff3_setup["retriever"]
        refseq = r.list_phages(source_db="RefSeq")
        assert "phage_A" in refseq
        assert "phage_B" in refseq
        assert "phage_C" not in refseq

        phagesdb = r.list_phages(source_db="PhagesDB")
        assert "phage_C" in phagesdb
        assert "phage_A" not in phagesdb


class TestGetSourceDB:
    def test_existing_phage(self, gff3_setup):
        r = gff3_setup["retriever"]
        assert r.get_source_db("phage_A") == "RefSeq"
        assert r.get_source_db("phage_C") == "PhagesDB"

    def test_missing_phage(self, gff3_setup):
        r = gff3_setup["retriever"]
        assert r.get_source_db("NONEXISTENT") is None


class TestHasPhage:
    def test_existing(self, gff3_setup):
        r = gff3_setup["retriever"]
        assert r.has_phage("phage_A") is True

    def test_missing(self, gff3_setup):
        r = gff3_setup["retriever"]
        assert r.has_phage("NONEXISTENT") is False


class TestStats:
    def test_stats_structure(self, gff3_setup):
        r = gff3_setup["retriever"]
        stats = r.stats()
        assert "total_phages" in stats
        assert "sources" in stats
        assert "index_path" in stats
        assert "gff3_dir" in stats
        assert stats["total_phages"] == 3
        assert stats["sources"]["RefSeq"] == 2
        assert stats["sources"]["PhagesDB"] == 1
