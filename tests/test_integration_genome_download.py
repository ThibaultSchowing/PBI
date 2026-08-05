"""
Integration checks for genome download pipeline files.

These tests validate that required source files, configuration sections, and
environment specs are present. They do NOT make live network calls.
"""

import re
from pathlib import Path

BASE_DIR = Path(__file__).parent.parent


# ---------------------------------------------------------------------------
# GTDB identifier detection
# ---------------------------------------------------------------------------

_GTDB_PATTERN = re.compile(r'\bsp\d{9}\b', re.IGNORECASE)


def test_gtdb_valid_species_not_matched():
    assert not _GTDB_PATTERN.search("Escherichia coli")
    assert not _GTDB_PATTERN.search("Staphylococcus aureus")


def test_gtdb_identifiers_matched():
    assert _GTDB_PATTERN.search("Acidovorax sp000302535")
    assert _GTDB_PATTERN.search("sp001411535")
    assert _GTDB_PATTERN.search("Bacteria sp123456789")


# ---------------------------------------------------------------------------
# Configuration file
# ---------------------------------------------------------------------------

def test_genome_download_config_exists():
    assert (BASE_DIR / "workflow" / "config" / "genome_download_config.yaml").is_file()


def test_genome_download_config_required_sections():
    cfg_path = BASE_DIR / "workflow" / "config" / "genome_download_config.yaml"
    content = cfg_path.read_text()
    for section in ("download:", "cache:", "parsing:", "ncbi:", "validation:", "logging:"):
        assert section in content, f"Missing config section: {section}"


# ---------------------------------------------------------------------------
# Workflow scripts
# ---------------------------------------------------------------------------

def test_download_scripts_exist():
    for script in (
        "workflow/scripts/sequences/download_host_genomes_robust.py",
    ):
        assert (BASE_DIR / script).is_file(), f"Missing script: {script}"


# ---------------------------------------------------------------------------
# Environment/dependency files
# ---------------------------------------------------------------------------

def test_sequences_env_file_exists():
    assert (BASE_DIR / "workflow" / "envs" / "sequences.yaml").is_file()


def test_sequences_env_required_deps():
    content = (BASE_DIR / "workflow" / "envs" / "sequences.yaml").read_text()
    for dep in ("biopython", "pandas", "pyyaml", "aiohttp"):
        assert dep in content, f"Missing dependency: {dep}"
