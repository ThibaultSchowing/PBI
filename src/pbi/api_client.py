"""
HTTP client that mirrors the SequenceRetriever interface via the PBI API.

Usage::

    from pbi.api_client import APIClient

    client = APIClient("http://pbi-api:8000")

    # Metadata queries (fast, no local FASTA needed)
    phages = client.get_phage_metadata(where="Source_DB = 'RefSeq'", limit=100)
    hosts = client.get_host_metadata(limit=50)

    # Single sequence retrieval
    seq = client.get_phage_sequence("NC_001330.1")

    # Genome retrieval
    genome = client.get_phage_genome("NC_001330.1", mode="concat")

    # Arbitrary SQL
    result = client.query("SELECT Source_DB, COUNT(*) FROM fact_phages GROUP BY Source_DB")

    client.close()
"""

from __future__ import annotations

import logging
from typing import Optional, Dict, List, Any, Union

import pandas as pd
import requests

logger = logging.getLogger(__name__)


class APIClient:
    """Client that provides a SequenceRetriever-like interface over HTTP.

    Parameters
    ----------
    base_url:
        Base URL of the PBI API (e.g. ``"http://pbi-api:8000"``).
    timeout:
        HTTP request timeout in seconds.
    """

    def __init__(self, base_url: str = "http://localhost:8000", timeout: int = 120):
        self._base_url = base_url.rstrip("/")
        self._timeout = timeout
        self._session = requests.Session()
        self._session.headers.update({"Content-Type": "application/json"})

    # ── Internal helpers ──────────────────────────────────────────────────

    def _get(self, path: str, params: Optional[Dict] = None) -> Dict:
        """Make a GET request and return JSON."""
        url = f"{self._base_url}{path}"
        resp = self._session.get(url, params=params, timeout=self._timeout)
        resp.raise_for_status()
        return resp.json()

    def _post(self, path: str, json_data: Optional[Dict] = None) -> Dict:
        """Make a POST request and return JSON."""
        url = f"{self._base_url}{path}"
        resp = self._session.post(url, json=json_data, timeout=self._timeout)
        resp.raise_for_status()
        return resp.json()

    def _get_text(self, path: str, params: Optional[Dict] = None) -> str:
        """Make a GET request and return plain text."""
        url = f"{self._base_url}{path}"
        resp = self._session.get(url, params=params, timeout=self._timeout)
        resp.raise_for_status()
        return resp.text

    def _records_to_df(self, data: List[Dict]) -> pd.DataFrame:
        """Convert list of dicts to DataFrame."""
        if not data:
            return pd.DataFrame()
        return pd.DataFrame(data)

    # ── Utility methods ───────────────────────────────────────────────────

    def health(self) -> Dict:
        """Check API health."""
        return self._get("/health")

    def get_stats(self) -> Dict:
        """Get database statistics."""
        return self._get("/stats")

    def list_tables(self) -> pd.DataFrame:
        """List all tables and views in the database."""
        result = self._get("/tables")
        return self._records_to_df(result.get("data", []))

    def query(self, sql: str, limit: Optional[int] = 1000) -> pd.DataFrame:
        """Execute a read-only SQL query and return results as a DataFrame.

        Parameters
        ----------
        sql:
            A SELECT statement to execute.
        limit:
            Maximum number of rows to return.
        """
        result = self._post("/query", {"query": sql, "limit": limit})
        return self._records_to_df(result.get("data", []))

    # ── Metadata methods ──────────────────────────────────────────────────

    def get_phage_metadata(
        self,
        where_clause: Optional[str] = None,
        limit: Optional[int] = 1000,
    ) -> pd.DataFrame:
        """Get phage metadata with optional SQL WHERE clause.

        Parameters
        ----------
        where_clause:
            SQL WHERE fragment (without the WHERE keyword).
            Example: ``"Source_DB = 'RefSeq' AND Length > 10000"``
        limit:
            Maximum rows to return.
        """
        params: Dict[str, Any] = {"limit": limit}
        if where_clause:
            params["where"] = where_clause
        result = self._get("/phage-metadata", params)
        return self._records_to_df(result.get("data", []))

    def get_host_metadata(
        self,
        where_clause: Optional[str] = None,
        limit: Optional[int] = 1000,
    ) -> pd.DataFrame:
        """Get host metadata with optional SQL WHERE clause.

        Parameters
        ----------
        where_clause:
            SQL WHERE fragment (without the WHERE keyword).
        limit:
            Maximum rows to return.
        """
        params: Dict[str, Any] = {"limit": limit}
        if where_clause:
            params["where"] = where_clause
        result = self._get("/host-metadata", params)
        return self._records_to_df(result.get("data", []))

    def get_phage_host_metadata(
        self,
        where_clause: Optional[str] = None,
        limit: Optional[int] = 1000,
    ) -> pd.DataFrame:
        """Get combined phage-host metadata with optional filtering.

        Parameters
        ----------
        where_clause:
            SQL WHERE fragment (without the WHERE keyword).
        limit:
            Maximum rows to return.
        """
        params: Dict[str, Any] = {"limit": limit}
        if where_clause:
            params["where"] = where_clause
        result = self._get("/phage-host-metadata", params)
        return self._records_to_df(result.get("data", []))

    def get_phage_host_pairs(
        self,
        where_clause: Optional[str] = None,
        limit: Optional[int] = 1000,
        host_contig_mode: str = "concat",
        phage_contig_mode: str = "first",
    ) -> pd.DataFrame:
        """Get phage-host pairs with sequences and metadata.

        Parameters
        ----------
        where_clause:
            SQL WHERE fragment.
        limit:
            Maximum rows to return.
        host_contig_mode:
            How to return host sequences: ``"first"``, ``"concat"``, ``"list"``, ``"dict"``.
        phage_contig_mode:
            How to return phage sequences: ``"first"``, ``"concat"``, ``"list"``, ``"dict"``.
        """
        params: Dict[str, Any] = {
            "limit": limit,
            "host_contig_mode": host_contig_mode,
            "phage_contig_mode": phage_contig_mode,
        }
        if where_clause:
            params["where"] = where_clause
        result = self._get("/phage-host-pairs", params)
        return self._records_to_df(result.get("data", []))

    def get_protein_metadata(
        self,
        where_clause: Optional[str] = None,
        limit: Optional[int] = 1000,
    ) -> pd.DataFrame:
        """Get protein metadata with optional filtering.

        Parameters
        ----------
        where_clause:
            SQL WHERE fragment.
        limit:
            Maximum rows to return.
        """
        params: Dict[str, Any] = {"limit": limit}
        if where_clause:
            params["where"] = where_clause
        result = self._get("/protein-metadata", params)
        return self._records_to_df(result.get("data", []))

    # ── Sequence methods ──────────────────────────────────────────────────

    def get_phage_sequence(self, phage_id: str) -> Optional[str]:
        """Get the DNA sequence for a single phage.

        Parameters
        ----------
        phage_id:
            Phage identifier (e.g. ``"NC_001330.1"``).

        Returns
        -------
        str or None
            The DNA sequence, or None if not found.
        """
        try:
            result = self._get(f"/phage/{phage_id}/sequence")
            return result.get("sequence")
        except requests.HTTPError as e:
            if e.response is not None and e.response.status_code == 404:
                return None
            raise

    def get_phage_sequences(
        self,
        query: str,
        limit: Optional[int] = None,
    ) -> pd.DataFrame:
        """Get phage sequences by SQL query.

        Parameters
        ----------
        query:
            SQL query that returns a ``Phage_ID`` column.
        limit:
            Maximum rows.
        """
        result = self._post("/phages", {"query": query, "limit": limit})
        return self._records_to_df(result.get("sequences", []))

    def get_protein_sequences(
        self,
        query: str,
        limit: Optional[int] = None,
    ) -> pd.DataFrame:
        """Get protein sequences by SQL query.

        Parameters
        ----------
        query:
            SQL query that returns a ``Protein_ID`` column.
        limit:
            Maximum rows.
        """
        result = self._post("/proteins", {"query": query, "limit": limit})
        return self._records_to_df(result.get("sequences", []))

    # ── Genome methods ────────────────────────────────────────────────────

    def get_phage_genome(
        self,
        phage_id: str,
        mode: str = "concat",
        gap: int = 0,
        order: str = "length_desc",
    ) -> Optional[Union[str, List[str], Dict[str, str]]]:
        """Get full phage genome.

        Parameters
        ----------
        phage_id:
            Phage identifier.
        mode:
            ``"concat"`` — single string with contigs joined by *gap*.
            ``"first"`` — first contig only.
            ``"list"`` — list of contig sequences.
            ``"dict"`` — dict mapping contig name to sequence.
        gap:
            Number of N characters between contigs in concat mode.
        order:
            Contig ordering: ``"length_desc"``, ``"length_asc"``, ``"natural"``.

        Returns
        -------
        str, list, dict, or None
            Genome data in the requested format, or None if not found.
        """
        params = {"mode": mode, "gap": gap, "order": order}
        try:
            result = self._get(f"/phage/{phage_id}/genome", params)
            if mode in ("concat", "first"):
                return result.get("sequence")
            return result.get("contigs")
        except requests.HTTPError as e:
            if e.response is not None and e.response.status_code == 404:
                return None
            raise

    def get_host_genome(
        self,
        host_id: str,
        mode: str = "concat",
        gap: int = 0,
        order: str = "length_desc",
    ) -> Optional[Union[str, List[str], Dict[str, str]]]:
        """Get full host genome.

        Parameters
        ----------
        host_id:
            Host identifier (e.g. ``"GCF_000005845.1"``).
        mode:
            ``"concat"`` — single string with contigs joined by *gap*.
            ``"first"`` — first contig only.
            ``"list"`` — list of contig sequences.
            ``"dict"`` — dict mapping contig name to sequence.
        gap:
            Number of N characters between contigs in concat mode.
        order:
            Contig ordering: ``"length_desc"``, ``"length_asc"``, ``"natural"``.

        Returns
        -------
        str, list, dict, or None
            Genome data in the requested format, or None if not found.
        """
        params = {"mode": mode, "gap": gap, "order": order}
        try:
            result = self._get(f"/host/{host_id}/genome", params)
            if mode in ("concat", "first"):
                return result.get("sequence")
            return result.get("contigs")
        except requests.HTTPError as e:
            if e.response is not None and e.response.status_code == 404:
                return None
            raise

    def get_host_genome_stats(self, host_id: str) -> Optional[Dict]:
        """Get contig statistics for a host genome.

        Parameters
        ----------
        host_id:
            Host identifier.

        Returns
        -------
        dict or None
            Statistics dict with contig count, lengths, total, etc.
        """
        try:
            result = self._get(f"/host/{host_id}/genome-stats")
            return result.get("stats")
        except requests.HTTPError as e:
            if e.response is not None and e.response.status_code == 404:
                return None
            raise

    # ── Compatibility aliases ─────────────────────────────────────────────

    def get_phage_host_pairs_iterator(
        self,
        where_clause: Optional[str] = None,
        batch_size: int = 1000,
        host_contig_mode: str = "first",
        phage_contig_mode: str = "first",
    ):
        """Yield batches of phage-host pairs (mimics SequenceRetriever iterator).

        Note: Unlike the direct retriever, this fetches the full result and
        yields chunks. For very large result sets, consider using
        ``get_phage_host_pairs()`` with a LIMIT/OFFSET where clause.
        """
        df = self.get_phage_host_pairs(
            where_clause=where_clause,
            limit=None,
            host_contig_mode=host_contig_mode,
            phage_contig_mode=phage_contig_mode,
        )
        for start in range(0, len(df), batch_size):
            yield df.iloc[start:start + batch_size]

    def close(self):
        """Close the HTTP session."""
        self._session.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    # ── BLAST methods ──────────────────────────────────────────────────────

    def blast_search(
        self,
        sequence: str,
        program: str = "blastn",
        db: Optional[str] = None,
        max_hits: int = 10,
        evalue: float = 1e-5,
    ) -> pd.DataFrame:
        """Search a sequence against BLAST databases via the API.

        Parameters
        ----------
        sequence:
            Query sequence (nucleotide or protein depending on *program*).
        program:
            BLAST program: ``"blastn"``, ``"blastp"``, ``"blastx"``,
            ``"tblastn"``, ``"tblastx"``.
        db:
            Database to search: ``"phages"``, ``"proteins"``, ``"hosts"``,
            ``"private"``, ``"combined"``, or ``None`` for auto-selection.
        max_hits:
            Maximum number of hits to return (server caps at 100).
        evalue:
            E-value threshold.

        Returns
        -------
        pd.DataFrame
            BLAST results in tabular format.
        """
        payload = {
            "sequence": sequence,
            "program": program,
            "max_hits": max_hits,
            "evalue": evalue,
        }
        if db is not None:
            payload["db"] = db
        result = self._post("/blast/search", payload)
        return self._records_to_df(result.get("results", []))

    def list_blast_databases(self) -> Dict:
        """List available BLAST databases and their status.

        Returns
        -------
        dict
            Mapping of database name to info dict with keys:
            ``type`` (nucl/prot), ``exists`` (bool), ``path`` (str).
        """
        result = self._get("/blast/databases")
        return result.get("databases", {})

    def blast_status(self) -> Dict:
        """Get BLAST installation and database status.

        Returns
        -------
        dict
            Keys: ``blast_installed``, ``databases_built``,
            ``available_databases``.
        """
        return self._get("/blast/status")
