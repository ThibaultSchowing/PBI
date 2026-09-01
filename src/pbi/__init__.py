"""
PBI - Phage Bacteria Interactions.

Main classes are exposed lazily so lightweight submodules (for example
``pbi.private_data``) can be imported without loading the full runtime stack.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .negative_examples import NegativeExampleGenerator
    from .sequence_retrieval import SequenceRetriever
    from .gff3_retrieval import GFF3Retriever
    from .blast_search import BlastSearcher
    from .api_client import APIClient
    from .streaming_dataset import (
        PhageHostIndexedDataset,
        PhageHostStreamingDataset,
        phage_host_collate_fn,
    )

# Package metadata
__version__ = "0.4.0"
__author__ = "Thibault Schowing, CI4CB"

# Public exports used by `from pbi import *`
__all__ = [
    'SequenceRetriever',
    'GFF3Retriever',
    'BlastSearcher',
    'APIClient',
    'NegativeExampleGenerator',
    'PhageHostStreamingDataset',
    'PhageHostIndexedDataset',
    'phage_host_collate_fn',
]


def __getattr__(name):
    # Lazy imports keep optional dependencies from loading at import time.
    if name == "SequenceRetriever":
        from .sequence_retrieval import SequenceRetriever
        return SequenceRetriever
    if name == "GFF3Retriever":
        from .gff3_retrieval import GFF3Retriever
        return GFF3Retriever
    if name == "APIClient":
        from .api_client import APIClient
        return APIClient
    if name == "NegativeExampleGenerator":
        from .negative_examples import NegativeExampleGenerator
        return NegativeExampleGenerator
    if name == "BlastSearcher":
        from .blast_search import BlastSearcher
        return BlastSearcher
    if name in {"PhageHostStreamingDataset", "PhageHostIndexedDataset", "phage_host_collate_fn"}:
        from .streaming_dataset import (
            PhageHostIndexedDataset,
            PhageHostStreamingDataset,
            phage_host_collate_fn,
        )
        exports = {
            "PhageHostStreamingDataset": PhageHostStreamingDataset,
            "PhageHostIndexedDataset": PhageHostIndexedDataset,
            "phage_host_collate_fn": phage_host_collate_fn,
        }
        return exports[name]
    raise AttributeError(
        f"module 'pbi' has no attribute {name!r}. "
        f"Available attributes: {', '.join(__all__)}"
    )


def get_default_paths():
    """
    Get default paths for database and FASTA files.

    Checks DATA_PATH first (container mode), then falls back to local project
    paths for development mode.

    Returns:
        dict: Paths to database and FASTA files.
    """
    import os
    from pathlib import Path

    data_path = os.environ.get('DATA_PATH')
    private_data_path = os.environ.get('PBI_PRIVATE_DATA_DIR')

    if data_path:
        # Container mode
        base_path = Path(data_path)
        private_base_path = Path(private_data_path) if private_data_path else Path('/private-data')
    else:
        # Local mode
        project_root = Path(__file__).parent.parent.parent
        base_path = project_root / 'data' / 'processed'
        private_base_path = Path(private_data_path) if private_data_path else project_root / 'private_data'

    return {
        'database': base_path / 'databases' / 'phage_database_optimized.duckdb',
        'phage_fasta': base_path / 'sequences' / 'all_phages.fasta',
        'protein_fasta': base_path / 'sequences' / 'all_proteins.fasta',
        'host_mapping': base_path / 'sequences' / 'host_fasta_mapping.json',
        'private_phage_mapping': private_base_path / 'private_phage_mapping.json',
        'gff3_dir': base_path / 'gff3',
        'gff3_index': base_path / 'gff3' / 'gff3_index.json',
        # Legacy path for backward compatibility
        'host_fasta': base_path / 'sequences' / 'all_hosts.fasta',
    }


def quick_connect():
    """
    Return a SequenceRetriever initialized with default data paths.

    Returns:
        SequenceRetriever: Connected retriever instance.
    """
    from .sequence_retrieval import SequenceRetriever
    paths = get_default_paths()

    # Prefer host mapping (current approach), fallback to legacy single-host FASTA.
    host_mapping = str(paths['host_mapping']) if paths['host_mapping'].exists() else None
    host_fasta = str(paths['host_fasta']) if paths['host_fasta'].exists() else None

    return SequenceRetriever(
        str(paths['database']),
        str(paths['phage_fasta']),
        str(paths['protein_fasta']),
        host_fasta_path=host_fasta,
        host_mapping_path=host_mapping,
        private_phage_mapping_path=(
            str(paths['private_phage_mapping'])
            if paths['private_phage_mapping'].exists()
            else None
        ),
    )
