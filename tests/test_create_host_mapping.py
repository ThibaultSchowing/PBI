#!/usr/bin/env python3
"""Tests for create_host_mapping.create_host_mapping()."""

import json
import textwrap
from pathlib import Path

import pandas as pd
import pytest

from workflow.scripts.sequences.create_host_mapping import create_host_mapping


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_csv(path: Path, rows: list[dict]) -> None:
    """Write a list of dicts as a CSV."""
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)


def _touch(path: Path, size: int = 10) -> None:
    """Create a file with *size* bytes (non-empty)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("A" * size)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestCreateHostMapping:
    """Unit tests for the host-mapping creation logic."""

    def test_basic_mapping(self, tmp_path):
        """Hosts with existing FASTA files appear in the mapping."""
        input_dir = tmp_path / "fastas"
        input_dir.mkdir()

        csv_path = tmp_path / "hosts.csv"
        _write_csv(csv_path, [
            {"Host_ID": "GCF_001", "Download_Status": "success"},
            {"Host_ID": "GCF_002", "Download_Status": "success"},
        ])

        _touch(input_dir / "GCF_001.fna", 100)
        _touch(input_dir / "GCF_002.fna", 200)

        out = tmp_path / "mapping.json"
        create_host_mapping(
            metadata_csv=str(csv_path),
            private_mapping_json=str(tmp_path / "nonexistent.json"),
            input_dir=str(input_dir),
            output_mapping=str(out),
        )

        mapping = json.loads(out.read_text())
        assert set(mapping.keys()) == {"GCF_001", "GCF_002"}

    def test_failed_downloads_skipped(self, tmp_path):
        """Hosts with Download_Status='failed' are excluded from the mapping."""
        input_dir = tmp_path / "fastas"
        input_dir.mkdir()

        csv_path = tmp_path / "hosts.csv"
        _write_csv(csv_path, [
            {"Host_ID": "GCF_001", "Download_Status": "success"},
            {"Host_ID": "GCF_002", "Download_Status": "failed"},
            {"Host_ID": "GCF_003", "Download_Status": "success"},
        ])

        _touch(input_dir / "GCF_001.fna", 100)
        # GCF_002 has no file (download failed)
        _touch(input_dir / "GCF_003.fna", 100)

        out = tmp_path / "mapping.json"
        create_host_mapping(
            metadata_csv=str(csv_path),
            private_mapping_json=str(tmp_path / "nonexistent.json"),
            input_dir=str(input_dir),
            output_mapping=str(out),
        )

        mapping = json.loads(out.read_text())
        assert "GCF_001" in mapping
        assert "GCF_002" not in mapping
        assert "GCF_003" in mapping

    def test_missing_download_status_treated_as_present(self, tmp_path):
        """Hosts without Download_Status column are not skipped (backward compat)."""
        input_dir = tmp_path / "fastas"
        input_dir.mkdir()

        csv_path = tmp_path / "hosts.csv"
        _write_csv(csv_path, [
            {"Host_ID": "GCF_001"},
            {"Host_ID": "GCF_002"},
        ])

        _touch(input_dir / "GCF_001.fna", 100)
        _touch(input_dir / "GCF_002.fna", 200)

        out = tmp_path / "mapping.json"
        create_host_mapping(
            metadata_csv=str(csv_path),
            private_mapping_json=str(tmp_path / "nonexistent.json"),
            input_dir=str(input_dir),
            output_mapping=str(out),
        )

        mapping = json.loads(out.read_text())
        assert set(mapping.keys()) == {"GCF_001", "GCF_002"}

    def test_all_failed_raises(self, tmp_path):
        """When every host failed and no FASTA exists, a ValueError is raised."""
        input_dir = tmp_path / "fastas"
        input_dir.mkdir()

        csv_path = tmp_path / "hosts.csv"
        _write_csv(csv_path, [
            {"Host_ID": "GCF_001", "Download_Status": "failed"},
            {"Host_ID": "GCF_002", "Download_Status": "failed"},
        ])

        out = tmp_path / "mapping.json"
        with pytest.raises(ValueError, match="No valid host FASTA files found"):
            create_host_mapping(
                metadata_csv=str(csv_path),
                private_mapping_json=str(tmp_path / "nonexistent.json"),
                input_dir=str(input_dir),
                output_mapping=str(out),
            )

    def test_missing_file_still_included_when_status_success(self, tmp_path):
        """A successful download that produced no file is still skipped (file check)."""
        input_dir = tmp_path / "fastas"
        input_dir.mkdir()

        csv_path = tmp_path / "hosts.csv"
        _write_csv(csv_path, [
            {"Host_ID": "GCF_001", "Download_Status": "success"},
        ])

        # No .fna file created — file check should exclude it
        out = tmp_path / "mapping.json"
        with pytest.raises(ValueError, match="No valid host FASTA files found"):
            create_host_mapping(
                metadata_csv=str(csv_path),
                private_mapping_json=str(tmp_path / "nonexistent.json"),
                input_dir=str(input_dir),
                output_mapping=str(out),
            )

    def test_empty_csv_raises(self, tmp_path):
        """An empty CSV raises an error (pandas parse error or no valid files)."""
        input_dir = tmp_path / "fastas"
        input_dir.mkdir()

        csv_path = tmp_path / "hosts.csv"
        _write_csv(csv_path, [])

        out = tmp_path / "mapping.json"
        with pytest.raises(Exception):
            create_host_mapping(
                metadata_csv=str(csv_path),
                private_mapping_json=str(tmp_path / "nonexistent.json"),
                input_dir=str(input_dir),
                output_mapping=str(out),
            )

    def test_private_mapping_merged(self, tmp_path):
        """Private host entries are merged into the mapping."""
        input_dir = tmp_path / "fastas"
        input_dir.mkdir()

        csv_path = tmp_path / "hosts.csv"
        _write_csv(csv_path, [
            {"Host_ID": "GCF_001", "Download_Status": "success"},
        ])
        _touch(input_dir / "GCF_001.fna", 100)

        private_dir = tmp_path / "private"
        private_dir.mkdir()
        private_fasta = private_dir / "private_host.fna"
        _touch(private_fasta, 100)

        private_mapping = tmp_path / "private_mapping.json"
        private_mapping.write_text(json.dumps({"PRIV_001": str(private_fasta)}))

        out = tmp_path / "mapping.json"
        create_host_mapping(
            metadata_csv=str(csv_path),
            private_mapping_json=str(private_mapping),
            input_dir=str(input_dir),
            output_mapping=str(out),
        )

        mapping = json.loads(out.read_text())
        assert "GCF_001" in mapping
        assert "PRIV_001" in mapping
