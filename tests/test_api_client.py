"""
Tests for pbi.api_client.APIClient.

Uses unittest.mock to mock HTTP responses, so no running server is needed.
"""

import json
from unittest.mock import patch, MagicMock, PropertyMock

import pandas as pd
import pytest
import requests


@pytest.fixture()
def client():
    """Return an APIClient with a mocked session."""
    from pbi.api_client import APIClient
    c = APIClient(base_url="http://mock-server:8000", timeout=5)
    yield c
    c.close()


def _mock_response(json_data=None, text="", status_code=200):
    """Create a mock requests.Response."""
    resp = MagicMock(spec=requests.Response)
    resp.status_code = status_code
    resp.json.return_value = json_data or {}
    resp.text = text
    resp.raise_for_status = MagicMock()
    if status_code >= 400:
        error = requests.HTTPError(response=resp)
        resp.raise_for_status.side_effect = error
    return resp


class TestAPIClientInit:
    def test_creates_session(self, client):
        assert client._session is not None
        assert client._base_url == "http://mock-server:8000"

    def test_strips_trailing_slash(self):
        from pbi.api_client import APIClient
        c = APIClient("http://example.com/")
        assert c._base_url == "http://example.com"
        c.close()


class TestContextManager:
    def test_enter_exit(self):
        from pbi.api_client import APIClient
        with APIClient("http://example.com") as c:
            session = c._session
        # Session close should have been called (can't assert .called on real method)
        # Just verify no exception was raised


class TestHealth:
    def test_health(self, client):
        with patch.object(client, "_get", return_value={"status": "healthy"}):
            result = client.health()
            assert result["status"] == "healthy"
            client._get.assert_called_once_with("/health")


class TestGetStats:
    def test_get_stats(self, client):
        mock_data = {"database": {"phages": 100}}
        with patch.object(client, "_get", return_value=mock_data):
            result = client.get_stats()
            assert result == mock_data


class TestListTables:
    def test_list_tables(self, client):
        mock_data = {"data": [{"name": "fact_phages", "type": "table"}]}
        with patch.object(client, "_get", return_value=mock_data):
            df = client.list_tables()
            assert isinstance(df, pd.DataFrame)
            assert len(df) == 1
            assert df.iloc[0]["name"] == "fact_phages"


class TestQuery:
    def test_query(self, client):
        mock_data = {"data": [{"x": 1}, {"x": 2}]}
        with patch.object(client, "_post", return_value=mock_data):
            df = client.query("SELECT 1 AS x")
            assert isinstance(df, pd.DataFrame)
            assert len(df) == 2


class TestGetPhageMetadata:
    def test_no_filter(self, client):
        mock_data = {"data": [{"Phage_ID": "p1"}, {"Phage_ID": "p2"}]}
        with patch.object(client, "_get", return_value=mock_data):
            df = client.get_phage_metadata()
            assert len(df) == 2

    def test_with_filter(self, client):
        mock_data = {"data": [{"Phage_ID": "p1"}]}
        with patch.object(client, "_get", return_value=mock_data) as m:
            client.get_phage_metadata(where_clause="Source_DB = 'RefSeq'", limit=100)
            # _get is called with (path, params_dict)
            call_args = m.call_args
            params = call_args[0][1]  # second positional arg
            assert params["where"] == "Source_DB = 'RefSeq'"
            assert params["limit"] == 100


class TestGetHostMetadata:
    def test_returns_dataframe(self, client):
        mock_data = {"data": [{"Host_ID": "h1"}]}
        with patch.object(client, "_get", return_value=mock_data):
            df = client.get_host_metadata()
            assert isinstance(df, pd.DataFrame)
            assert len(df) == 1


class TestGetPhageHostPairs:
    def test_returns_dataframe(self, client):
        mock_data = {"data": [{"Phage_ID": "p1", "Host_ID": "h1"}]}
        with patch.object(client, "_get", return_value=mock_data):
            df = client.get_phage_host_pairs()
            assert len(df) == 1


class TestGetProteinMetadata:
    def test_returns_dataframe(self, client):
        mock_data = {"data": [{"Protein_ID": "prot1"}]}
        with patch.object(client, "_get", return_value=mock_data):
            df = client.get_protein_metadata()
            assert len(df) == 1


class TestGetPhageSequence:
    def test_existing(self, client):
        mock_data = {"sequence": "ATCGATCG"}
        with patch.object(client, "_get", return_value=mock_data):
            seq = client.get_phage_sequence("phage_001")
            assert seq == "ATCGATCG"

    def test_not_found_returns_none(self, client):
        resp = _mock_response(status_code=404)
        with patch.object(client, "_get", side_effect=requests.HTTPError(response=resp)):
            seq = client.get_phage_sequence("NONEXISTENT")
            assert seq is None


class TestGetPhageGenome:
    def test_concat_mode(self, client):
        mock_data = {"sequence": "ATCGATCG", "mode": "concat"}
        with patch.object(client, "_get", return_value=mock_data):
            result = client.get_phage_genome("phage_001", mode="concat")
            assert result == "ATCGATCG"

    def test_list_mode(self, client):
        mock_data = {"contigs": ["ATCG", "GCTA"], "mode": "list"}
        with patch.object(client, "_get", return_value=mock_data):
            result = client.get_phage_genome("phage_001", mode="list")
            assert isinstance(result, list)
            assert len(result) == 2

    def test_not_found_returns_none(self, client):
        resp = _mock_response(status_code=404)
        with patch.object(client, "_get", side_effect=requests.HTTPError(response=resp)):
            result = client.get_phage_genome("NONEXISTENT")
            assert result is None


class TestGetHostGenome:
    def test_concat_mode(self, client):
        mock_data = {"sequence": "GCTAGCTA", "mode": "concat"}
        with patch.object(client, "_get", return_value=mock_data):
            result = client.get_host_genome("host_001", mode="concat")
            assert result == "GCTAGCTA"

    def test_not_found_returns_none(self, client):
        resp = _mock_response(status_code=404)
        with patch.object(client, "_get", side_effect=requests.HTTPError(response=resp)):
            result = client.get_host_genome("NONEXISTENT")
            assert result is None


class TestGetHostGenomeStats:
    def test_existing(self, client):
        mock_data = {"stats": {"contig_count": 1, "total_length": 4600000}}
        with patch.object(client, "_get", return_value=mock_data):
            stats = client.get_host_genome_stats("host_001")
            assert stats["contig_count"] == 1

    def test_not_found_returns_none(self, client):
        resp = _mock_response(status_code=404)
        with patch.object(client, "_get", side_effect=requests.HTTPError(response=resp)):
            stats = client.get_host_genome_stats("NONEXISTENT")
            assert stats is None


class TestGetPhageHostPairsIterator:
    def test_yields_batches(self, client):
        # Mock get_phage_host_pairs to return a DataFrame
        df = pd.DataFrame({
            "Phage_ID": [f"p{i}" for i in range(25)],
            "Host_ID": [f"h{i}" for i in range(25)],
        })
        with patch.object(client, "get_phage_host_pairs", return_value=df):
            batches = list(client.get_phage_host_pairs_iterator(batch_size=10))
            assert len(batches) == 3  # 10 + 10 + 5
            assert len(batches[0]) == 10
            assert len(batches[2]) == 5


class TestClose:
    def test_close(self, client):
        # Verify close doesn't raise
        client.close()
