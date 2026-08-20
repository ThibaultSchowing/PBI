#!/usr/bin/env python3
"""Tests for SequenceRetriever.get_pair_ids_only()."""

import pytest


class TestGetPairIdsOnly:
    """Unit tests against the in-memory test database."""

    def test_returns_all_pairs(self, retriever):
        df = retriever.get_pair_ids_only()
        assert set(df.columns) == {"Phage_ID", "Host_ID"}
        # conftest inserts 4 associations
        assert len(df) == 4

    def test_limit(self, retriever):
        df = retriever.get_pair_ids_only(limit=2)
        assert len(df) == 2

    def test_shuffle_returns_same_columns(self, retriever):
        df = retriever.get_pair_ids_only(shuffle=True)
        assert set(df.columns) == {"Phage_ID", "Host_ID"}
        assert len(df) == 4

    def test_shuffle_is_deterministic(self, retriever):
        df1 = retriever.get_pair_ids_only(shuffle=True)
        df2 = retriever.get_pair_ids_only(shuffle=True)
        assert df1.values.tolist() == df2.values.tolist()

    def test_shuffle_differs_from_default(self, retriever):
        default = retriever.get_pair_ids_only(shuffle=False)
        shuffled = retriever.get_pair_ids_only(shuffle=True)
        # The order should differ (MD5 vs natural).  We allow the pathological
        # case where they happen to coincide, but it is astronomically unlikely
        # with 4 distinct pairs.
        assert default.values.tolist() != shuffled.values.tolist()

    def test_seed_has_no_effect(self, retriever):
        df_a = retriever.get_pair_ids_only(shuffle=True, seed=42)
        df_b = retriever.get_pair_ids_only(shuffle=True, seed=99)
        assert df_a.values.tolist() == df_b.values.tolist()

    def test_limit_with_shuffle(self, retriever):
        df = retriever.get_pair_ids_only(shuffle=True, limit=1)
        assert len(df) == 1
        assert set(df.columns) == {"Phage_ID", "Host_ID"}


class TestGetPairIdsOnlyWithMapping:
    """Tests using the host-mapping fixture (same DB, different retriever)."""

    def test_returns_all_pairs(self, retriever_with_mapping):
        df = retriever_with_mapping.get_pair_ids_only()
        assert len(df) == 4

    def test_shuffle_deterministic(self, retriever_with_mapping):
        df1 = retriever_with_mapping.get_pair_ids_only(shuffle=True)
        df2 = retriever_with_mapping.get_pair_ids_only(shuffle=True)
        assert df1.values.tolist() == df2.values.tolist()
