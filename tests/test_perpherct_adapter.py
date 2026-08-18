"""
Tests for PERPHECT PBIAdapter module.
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "PERPHECT"))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_retriever(phage_seqs=None, host_seqs=None):
    """Create a mock SequenceRetriever with controllable sequence returns."""
    retriever = MagicMock()
    retriever.get_phage_sequence = MagicMock(
        side_effect=lambda pid: phage_seqs.get(pid) if phage_seqs else "ATCG" * 100
    )
    retriever.get_host_sequence = MagicMock(
        side_effect=lambda hid, contig_mode="concat": host_seqs.get(hid) if host_seqs else "ATCG" * 500
    )
    return retriever


def _make_positive_pairs(n=5):
    """Create a minimal positive pairs DataFrame."""
    return pd.DataFrame({
        "Phage_ID": [f"phage_{i}" for i in range(n)],
        "Host_ID": [f"host_{i}" for i in range(n)],
        "Phage_Sequence": ["ATCG" * 100] * n,
        "Host_Sequence": ["ATCG" * 500] * n,
    })


def _make_negative_pairs(n=5):
    """Create a minimal negative pairs DataFrame."""
    return pd.DataFrame({
        "Phage_ID": [f"phage_{i+n}" for i in range(n)],
        "Host_ID": [f"host_{i+n}" for i in range(n)],
        "Phage_Sequence": ["GCTA" * 100] * n,
        "Host_Sequence": ["TTTT" * 500] * n,
    })


# ---------------------------------------------------------------------------
# ID Mapping Tests
# ---------------------------------------------------------------------------

class TestPBIAdapterIDMapping:
    """Tests for string-to-integer ID mapping."""

    def test_host_id_mapping(self):
        from pbi_adapter import PBIAdapter

        retriever = _make_retriever()
        adapter = PBIAdapter(retriever)

        id1 = adapter._map_host_id("host_A")
        id2 = adapter._map_host_id("host_B")
        id3 = adapter._map_host_id("host_A")  # duplicate

        assert id1 == 0
        assert id2 == 1
        assert id3 == 0  # same as id1

    def test_phage_id_mapping(self):
        from pbi_adapter import PBIAdapter

        retriever = _make_retriever()
        adapter = PBIAdapter(retriever)

        id1 = adapter._map_phage_id("phage_X")
        id2 = adapter._map_phage_id("phage_Y")

        assert id1 == 0
        assert id2 == 1

    def test_get_id_maps(self):
        from pbi_adapter import PBIAdapter

        retriever = _make_retriever()
        adapter = PBIAdapter(retriever)

        adapter._map_host_id("host_A")
        adapter._map_phage_id("phage_X")

        host_map, phage_map = adapter.get_id_maps()
        assert host_map == {"host_A": 0}
        assert phage_map == {"phage_X": 0}


# ---------------------------------------------------------------------------
# Pad and Encode Tests
# ---------------------------------------------------------------------------

class TestPBIAdapterPadAndEncode:
    """Tests for sequence padding and one-hot encoding."""

    def test_short_sequence_padded(self):
        from pbi_adapter import PBIAdapter

        retriever = _make_retriever()
        adapter = PBIAdapter(retriever, phage_threshold=100)

        result = adapter._pad_and_encode("ATCG", 100)
        assert result.shape == (100, 4)
        # First 4 rows should be one-hot, rest zeros
        # ATCG: A=col0, T=col1, C=col3, G=col2
        assert result[0, 0] == 1  # A
        assert result[1, 1] == 1  # T
        assert result[2, 3] == 1  # C
        assert result[3, 2] == 1  # G
        assert np.sum(result[4:]) == 0  # padding

    def test_long_sequence_truncated(self):
        from pbi_adapter import PBIAdapter

        retriever = _make_retriever()
        adapter = PBIAdapter(retriever, phage_threshold=4)

        result = adapter._pad_and_encode("ATCGATCG", 4)
        assert result.shape == (4, 4)
        # Only first 4 bases kept: A, T, C, G
        np.testing.assert_array_equal(result[0], [1, 0, 0, 0])  # A
        np.testing.assert_array_equal(result[1], [0, 1, 0, 0])  # T
        np.testing.assert_array_equal(result[2], [0, 0, 0, 1])  # C
        np.testing.assert_array_equal(result[3], [0, 0, 1, 0])  # G

    def test_exact_length(self):
        from pbi_adapter import PBIAdapter

        retriever = _make_retriever()
        adapter = PBIAdapter(retriever, phage_threshold=4)

        result = adapter._pad_and_encode("ATCG", 4)
        assert result.shape == (4, 4)


# ---------------------------------------------------------------------------
# DataFrame Conversion Tests
# ---------------------------------------------------------------------------

class TestPBIAdapterDataframes:
    """Tests for to_perphect_dataframes."""

    def test_basic_conversion(self):
        from pbi_adapter import PBIAdapter

        retriever = _make_retriever()
        adapter = PBIAdapter(
            retriever,
            bacterium_threshold=100,
            phage_threshold=50,
            bacterium_min_length=10,
            phage_min_length=10,
        )

        positives = _make_positive_pairs(3)
        couples_df, bacteria_df, phages_df = adapter.to_perphect_dataframes(positives)

        # Check couples_df structure
        assert list(couples_df.columns) == ["id", "bacterium_id", "phage_id", "interaction_type"]
        assert len(couples_df) == 3
        assert all(couples_df["interaction_type"] == 1)

        # Check bacteria_df structure
        assert "bacterium_sequence" in bacteria_df.columns
        assert bacteria_df.index.name == "bacterium_id"

        # Check phages_df structure
        assert "phage_sequence" in phages_df.columns
        assert phages_df.index.name == "phage_id"

    def test_with_negatives(self):
        from pbi_adapter import PBIAdapter

        retriever = _make_retriever()
        adapter = PBIAdapter(
            retriever,
            bacterium_threshold=100,
            phage_threshold=50,
            bacterium_min_length=10,
            phage_min_length=10,
        )

        positives = _make_positive_pairs(2)
        negatives = _make_negative_pairs(2)
        couples_df, _, _ = adapter.to_perphect_dataframes(positives, negatives)

        assert len(couples_df) == 4
        assert sum(couples_df["interaction_type"] == 1) == 2
        assert sum(couples_df["interaction_type"] == 0) == 2

    def test_drops_too_short_sequences(self):
        from pbi_adapter import PBIAdapter

        retriever = _make_retriever()
        adapter = PBIAdapter(
            retriever,
            bacterium_min_length=10000,
            phage_min_length=10000,
        )

        positives = _make_positive_pairs(3)
        couples_df, _, _ = adapter.to_perphect_dataframes(positives)

        # All should be dropped (test sequences are too short)
        assert len(couples_df) == 0


# ---------------------------------------------------------------------------
# Generator Tests
# ---------------------------------------------------------------------------

class TestPBIAdapterGenerator:
    """Tests for create_tf_generator."""

    def test_generator_yields_batches(self):
        from pbi_adapter import PBIAdapter

        retriever = _make_retriever()
        adapter = PBIAdapter(
            retriever,
            bacterium_threshold=20,
            phage_threshold=10,
            bacterium_min_length=10,
            phage_min_length=10,
        )

        # Manually set up ID maps and sequences
        adapter._map_host_id("host_0")
        adapter._map_host_id("host_1")
        adapter._map_phage_id("phage_0")
        adapter._map_phage_id("phage_1")
        adapter._host_sequences = {"host_0": "ATCG" * 5, "host_1": "GCTA" * 5}
        adapter._phage_sequences = {"phage_0": "AAAA" * 5, "phage_1": "TTTT" * 5}

        couples = np.array([[0, 0], [1, 1], [0, 1], [1, 0]])
        labels = np.array([1.0, 1.0, 0.0, 0.0])

        gen = adapter.create_tf_generator(couples, labels, batch_size=2, shuffle=False)
        batch_inputs, batch_targets = next(gen)

        assert isinstance(batch_inputs, tuple)
        assert len(batch_inputs) == 2
        assert batch_inputs[0].shape == (2, 20, 4)  # bacterium
        assert batch_inputs[1].shape == (2, 10, 4)  # phage
        assert batch_targets.shape == (2,)


# ---------------------------------------------------------------------------
# Prepare Training Data Tests
# ---------------------------------------------------------------------------

class TestPrepareTrainingData:
    """Tests for prepare_training_data."""

    def test_returns_couples_and_labels(self):
        from pbi_adapter import PBIAdapter

        retriever = _make_retriever()
        adapter = PBIAdapter(
            retriever,
            bacterium_min_length=10,
            phage_min_length=10,
        )

        positives = _make_positive_pairs(3)
        couples, labels = adapter.prepare_training_data(positives)

        assert couples.shape == (3, 2)
        assert labels.shape == (3,)
        assert all(labels == 1.0)

    def test_with_mixed_pairs(self):
        from pbi_adapter import PBIAdapter

        retriever = _make_retriever()
        adapter = PBIAdapter(
            retriever,
            bacterium_min_length=10,
            phage_min_length=10,
        )

        positives = _make_positive_pairs(2)
        negatives = _make_negative_pairs(2)
        couples, labels = adapter.prepare_training_data(positives, negatives)

        assert len(couples) == 4
        assert sum(labels == 1.0) == 2
        assert sum(labels == 0.0) == 2

    def test_raises_on_empty_data(self):
        from pbi_adapter import PBIAdapter

        retriever = _make_retriever()
        adapter = PBIAdapter(
            retriever,
            bacterium_min_length=100000,
            phage_min_length=100000,
        )

        positives = _make_positive_pairs(3)
        with pytest.raises(ValueError, match="No valid pairs found"):
            adapter.prepare_training_data(positives)


# ---------------------------------------------------------------------------
# Classify Pairs by Interaction Type Tests
# ---------------------------------------------------------------------------

class TestClassifyPairsByInteraction:
    """Tests for classify_pairs_by_interaction."""

    def _make_retriever_with_interactions(self, interactions_df):
        """Create a mock retriever with a conn that returns interaction types."""
        retriever = _make_retriever()
        retriever.conn = MagicMock()
        retriever.conn.execute = MagicMock(return_value=MagicMock(fetchdf=MagicMock(return_value=interactions_df)))
        return retriever

    def test_classifies_no_interaction_as_negative(self):
        from pbi_adapter import PBIAdapter

        interactions_df = pd.DataFrame({
            "Phage_ID": ["phage_0"],
            "Host_ID": ["host_0"],
            "interaction": ["no interaction"],
        })
        retriever = self._make_retriever_with_interactions(interactions_df)
        adapter = PBIAdapter(retriever)

        all_pairs = _make_positive_pairs(3)
        pos_df, neg_df = adapter.classify_pairs_by_interaction(all_pairs)

        # phage_0/host_0 should be negative, others positive
        assert len(neg_df) == 1
        assert neg_df.iloc[0]["Phage_ID"] == "phage_0"
        assert neg_df.iloc[0]["negative_source"] == "private_data"
        assert len(pos_df) == 2

    def test_classifies_virulent_as_positive(self):
        from pbi_adapter import PBIAdapter

        interactions_df = pd.DataFrame({
            "Phage_ID": ["phage_0"],
            "Host_ID": ["host_0"],
            "interaction": ["virulent"],
        })
        retriever = self._make_retriever_with_interactions(interactions_df)
        adapter = PBIAdapter(retriever)

        all_pairs = _make_positive_pairs(3)
        pos_df, neg_df = adapter.classify_pairs_by_interaction(all_pairs)

        # All should be positive
        assert len(pos_df) == 3
        assert len(neg_df) == 0

    def test_handles_empty_interactions_table(self):
        from pbi_adapter import PBIAdapter

        interactions_df = pd.DataFrame(columns=["Phage_ID", "Host_ID", "interaction"])
        retriever = self._make_retriever_with_interactions(interactions_df)
        adapter = PBIAdapter(retriever)

        all_pairs = _make_positive_pairs(3)
        pos_df, neg_df = adapter.classify_pairs_by_interaction(all_pairs)

        # All should be positive (no interactions to classify)
        assert len(pos_df) == 3
        assert len(neg_df) == 0

    def test_handles_missing_interaction_table(self):
        from pbi_adapter import PBIAdapter

        retriever = _make_retriever()
        retriever.conn = MagicMock()
        retriever.conn.execute = MagicMock(side_effect=Exception("table not found"))
        adapter = PBIAdapter(retriever)

        all_pairs = _make_positive_pairs(3)
        pos_df, neg_df = adapter.classify_pairs_by_interaction(all_pairs)

        # Should gracefully fall back to all positive
        assert len(pos_df) == 3
        assert len(neg_df) == 0

    def test_multiple_negative_interactions(self):
        from pbi_adapter import PBIAdapter

        interactions_df = pd.DataFrame({
            "Phage_ID": ["phage_0", "phage_1", "phage_2"],
            "Host_ID": ["host_0", "host_1", "host_2"],
            "interaction": ["no interaction", "none", "virulent"],
        })
        retriever = self._make_retriever_with_interactions(interactions_df)
        adapter = PBIAdapter(retriever)

        all_pairs = _make_positive_pairs(3)
        pos_df, neg_df = adapter.classify_pairs_by_interaction(all_pairs)

        # phage_0 and phage_1 are negative, phage_2 is positive
        assert len(neg_df) == 2
        assert len(pos_df) == 1
        assert set(neg_df["Phage_ID"]) == {"phage_0", "phage_1"}
