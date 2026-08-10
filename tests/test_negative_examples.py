"""
Tests for pbi.negative_examples.NegativeExampleGenerator.

Uses in-memory DuckDB with test data and host FASTA files.
"""

import json
from pathlib import Path

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
def neg_gen(tmp_path, tmp_db):
    """Create a NegativeExampleGenerator with test data."""
    from pbi.sequence_retrieval import SequenceRetriever
    from pbi.negative_examples import NegativeExampleGenerator

    # Create phage and protein FASTA files
    phage_records = {
        "phage_001": "ATCG" * 325,    # 1300 bp
        "phage_002": "GCTA" * 240,    # 960 bp
        "phage_003": "TTTT" * 600,    # 2400 bp
    }
    protein_records = {
        "prot_001": "ACGT" * 75,      # 300 aa
        "prot_002": "TGCA" * 111,     # 444 aa
        "prot_003": "AAAA" * 50,      # 200 aa
    }
    host_records = {
        "host_001": "ATCG" * 500,     # 2000 bp
        "host_002": "GCTA" * 400,     # 1600 bp
        "host_003": "CCCC" * 300,     # 1200 bp
    }

    phage_fasta = tmp_path / "phages.fasta"
    protein_fasta = tmp_path / "proteins.fasta"
    host_fasta = tmp_path / "hosts.fasta"

    _write_fasta(phage_fasta, phage_records)
    _write_fasta(protein_fasta, protein_records)
    _write_fasta(host_fasta, host_records)
    _write_fai(phage_fasta.parent / "phages.fasta.fai", phage_records)
    _write_fai(protein_fasta.parent / "proteins.fasta.fai", protein_records)
    _write_fai(host_fasta.parent / "hosts.fasta.fai", host_records)

    retriever = SequenceRetriever(
        db_path=str(tmp_db),
        phage_fasta_path=str(phage_fasta),
        protein_fasta_path=str(protein_fasta),
        host_fasta_path=str(host_fasta),
        preload=False,
    )

    gen = NegativeExampleGenerator(retriever)
    yield gen
    retriever.close()


class TestInit:
    def test_requires_host_data(self, tmp_db, tmp_path):
        from pbi.sequence_retrieval import SequenceRetriever
        from pbi.negative_examples import NegativeExampleGenerator

        phage_fasta = tmp_path / "phages.fasta"
        protein_fasta = tmp_path / "proteins.fasta"
        _write_fasta(phage_fasta, {"p1": "ATCG"})
        _write_fasta(protein_fasta, {"prot1": "ACGT"})
        _write_fai(phage_fasta.parent / "phages.fasta.fai", {"p1": "ATCG"})
        _write_fai(protein_fasta.parent / "proteins.fasta.fai", {"prot1": "ACGT"})

        retriever = SequenceRetriever(
            db_path=str(tmp_db),
            phage_fasta_path=str(phage_fasta),
            protein_fasta_path=str(protein_fasta),
            preload=False,
        )
        with pytest.raises(ValueError, match="host data"):
            NegativeExampleGenerator(retriever)
        retriever.close()

    def test_caches_phages_and_hosts(self, neg_gen):
        assert len(neg_gen.all_phages) > 0
        assert len(neg_gen.all_hosts) > 0


class TestGenerateRandomNegatives:
    def test_correct_count(self, neg_gen):
        positives = pd.DataFrame({"Phage_ID": ["phage_001"], "Host_ID": ["host_001"]})
        negatives = neg_gen.generate_random_negatives(positives, ratio=1.0)
        assert len(negatives) > 0
        assert len(negatives) <= 1  # ratio=1.0 with 1 positive

    def test_no_overlap_with_positives(self, neg_gen):
        positives = pd.DataFrame({
            "Phage_ID": ["phage_001", "phage_002"],
            "Host_ID": ["host_001", "host_002"],
        })
        negatives = neg_gen.generate_random_negatives(positives, ratio=1.0)
        pos_set = set(zip(positives["Phage_ID"], positives["Host_ID"]))
        neg_set = set(zip(negatives["Phage_ID"], negatives["Host_ID"]))
        assert len(pos_set & neg_set) == 0

    def test_label_column(self, neg_gen):
        positives = pd.DataFrame({"Phage_ID": ["phage_001"], "Host_ID": ["host_001"]})
        negatives = neg_gen.generate_random_negatives(positives, ratio=0.5)
        if len(negatives) > 0:
            assert (negatives["Label"] == 0).all()


class TestGenerateGCBasedNegatives:
    def test_returns_dataframe(self, neg_gen):
        positives = pd.DataFrame({"Phage_ID": ["phage_001"], "Host_ID": ["host_001"]})
        negatives = neg_gen.generate_gc_based_negatives(positives, ratio=0.5, min_gc_difference=5.0)
        assert isinstance(negatives, pd.DataFrame)

    def test_gc_difference_threshold(self, neg_gen):
        positives = pd.DataFrame({"Phage_ID": ["phage_001"], "Host_ID": ["host_001"]})
        negatives = neg_gen.generate_gc_based_negatives(
            positives, ratio=0.5, min_gc_difference=50.0
        )
        if len(negatives) > 0 and "GC_Difference" in negatives.columns:
            assert (negatives["GC_Difference"] >= 50.0).all()


class TestGenerateTaxonomyBasedNegatives:
    def test_returns_dataframe(self, neg_gen):
        positives = pd.DataFrame({"Phage_ID": ["phage_001"], "Host_ID": ["host_001"]})
        negatives = neg_gen.generate_taxonomy_based_negatives(positives, ratio=0.5)
        assert isinstance(negatives, pd.DataFrame)

    def test_with_exclude_species(self, neg_gen):
        positives = pd.DataFrame({"Phage_ID": ["phage_001"], "Host_ID": ["host_001"]})
        negatives = neg_gen.generate_taxonomy_based_negatives(
            positives, ratio=0.5, exclude_species=["Escherichia coli"]
        )
        assert isinstance(negatives, pd.DataFrame)


class TestGenerateBalancedDataset:
    def test_random_strategy(self, neg_gen):
        positives = pd.DataFrame({
            "Phage_ID": ["phage_001", "phage_002"],
            "Host_ID": ["host_001", "host_002"],
        })
        dataset = neg_gen.generate_balanced_dataset(
            positive_pairs=positives, strategy="random", positive_ratio=0.5
        )
        assert "Label" in dataset.columns
        assert len(dataset) > 0

    def test_mixed_strategy(self, neg_gen):
        positives = pd.DataFrame({
            "Phage_ID": ["phage_001", "phage_002", "phage_003"],
            "Host_ID": ["host_001", "host_002", "host_003"],
        })
        dataset = neg_gen.generate_balanced_dataset(
            positive_pairs=positives, strategy="mixed", positive_ratio=0.5
        )
        assert "Label" in dataset.columns
        pos_count = (dataset["Label"] == 1).sum()
        neg_count = (dataset["Label"] == 0).sum()
        assert pos_count > 0
        assert neg_count > 0

    def test_unknown_strategy_raises(self, neg_gen):
        positives = pd.DataFrame({"Phage_ID": ["phage_001"], "Host_ID": ["host_001"]})
        with pytest.raises(ValueError, match="Unknown strategy"):
            neg_gen.generate_balanced_dataset(
                positive_pairs=positives, strategy="invalid"
            )
