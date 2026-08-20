#!/usr/bin/env python3
"""
Build the Host_ID → FASTA-file-path JSON mapping.

Reads the host assembly metadata CSV, matches each Host_ID to its downloaded
``.fna`` file under *params.input_dir*, merges any private-dataset entries,
and writes a JSON mapping file.

Used by the ``create_host_mapping`` Snakemake rule in ``hosts.smk``.
"""

import json
import logging
import os
import sys
import time
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Shared logging helper
# ---------------------------------------------------------------------------

_PRIVATE_CONFLICT_EXAMPLE_LIMIT = 20


def _setup_logging(log_file: str, also_stderr: bool = True) -> None:
    """Route root-logger output to *log_file* and optionally stderr."""
    try:
        _common = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'common')
        if _common not in sys.path:
            sys.path.insert(0, _common)
        from logging_utils import setup_logging  # noqa: PLC0415
        setup_logging(log_file, also_stderr=also_stderr)
    except Exception:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        root = logging.getLogger()
        root.handlers.clear()
        root.setLevel(logging.INFO)
        fh = logging.FileHandler(log_file)
        fh.setFormatter(fmt)
        root.addHandler(fh)
        if also_stderr:
            sh = logging.StreamHandler(sys.stderr)
            sh.setFormatter(fmt)
            root.addHandler(sh)


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def create_host_mapping(
    metadata_csv: str,
    private_mapping_json: str,
    input_dir: str,
    output_mapping: str,
) -> None:
    """Build Host_ID → FASTA path mapping and write it as JSON.

    Args:
        metadata_csv:         Path to host assembly metadata CSV.
        private_mapping_json: Path to private-dataset host mapping JSON.
        input_dir:            Directory that contains ``<Host_ID>.fna`` files.
        output_mapping:       Destination path for the output JSON.

    Raises:
        ValueError: if no valid host FASTA files are found.
    """
    t0 = time.time()
    logging.info("=" * 70)
    logging.info("🚀 Starting host mapping creation")
    logging.info(f"   Metadata CSV    : {metadata_csv}")
    logging.info(f"   Private mapping : {private_mapping_json}")
    logging.info(f"   Input dir       : {input_dir}")
    logging.info(f"   Output          : {output_mapping}")

    metadata_df = pd.read_csv(metadata_csv)
    logging.info(f"   Host rows in CSV: {len(metadata_df)}")

    output_path = Path(output_mapping)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    input_dir_path = Path(input_dir)

    host_mapping: dict = {}
    valid_count = 0
    missing_count = 0
    failed_count = 0
    missing_examples: list = []

    for _, row in metadata_df.iterrows():
        host_id = row['Host_ID']

        # Skip hosts whose genome download failed — they have no FASTA file
        # and would create orphan entries in the mapping.
        if str(row.get('Download_Status', '')).strip().lower() == 'failed':
            failed_count += 1
            if len(missing_examples) < 5:
                missing_examples.append(f"{host_id} (download failed)")
            continue

        fasta_file = input_dir_path / f"{host_id}.fna"

        if fasta_file.exists() and fasta_file.stat().st_size > 0:
            host_mapping[host_id] = str(fasta_file)
            valid_count += 1
        else:
            missing_count += 1
            if len(missing_examples) < 5:
                missing_examples.append(str(fasta_file))
            logging.warning(f"⚠️  Missing or empty file: {fasta_file}")

    if not host_mapping:
        raise ValueError("❌ No valid host FASTA files found!")

    if failed_count:
        logging.info(
            f"⏭️  Skipped {failed_count} host(s) with failed downloads"
        )
    if missing_count:
        logging.warning(
            f"⚠️  {missing_count} host file(s) were missing or empty"
            + (f" (e.g. {missing_examples[0]})" if missing_examples else "")
        )

    # Merge private host mapping entries (if any)
    private_mapping_path = Path(private_mapping_json)
    private_added = 0
    private_conflicts = 0
    conflict_ids: list = []

    if private_mapping_path.exists():
        with private_mapping_path.open("r") as f:
            private_mapping = json.load(f)

        for pid, fasta_path in private_mapping.items():
            if pid in host_mapping:
                private_conflicts += 1
                if len(conflict_ids) < _PRIVATE_CONFLICT_EXAMPLE_LIMIT:
                    conflict_ids.append(pid)
                continue
            fasta_file = Path(fasta_path)
            if fasta_file.exists() and fasta_file.stat().st_size > 0:
                host_mapping[pid] = str(fasta_file)
                private_added += 1

    # Write JSON mapping
    with open(output_mapping, 'w') as f:
        json.dump(host_mapping, f, indent=2)

    elapsed = time.time() - t0
    logging.info(f"✅ Created mapping for {valid_count} public host FASTA files")
    logging.info(f"✅ Added {private_added} private host FASTA files")
    if private_conflicts > 0:
        logging.warning(
            f"⚠️  Skipped {private_conflicts} conflicting private Host_ID entries"
        )
        logging.warning(
            f"   Conflicting examples: {', '.join(conflict_ids)}"
        )
    if missing_count > 0:
        logging.warning(f"⚠️  {missing_count} host file(s) were missing or empty")
    logging.info(f"   Total hosts in mapping : {len(host_mapping)}")
    logging.info(f"   Elapsed                : {elapsed:.1f}s")
    logging.info("=" * 70)


if __name__ == "__main__":
    # Snakemake integration
    _setup_logging(snakemake.log[0])  # noqa: F821

    create_host_mapping(
        metadata_csv=snakemake.input.metadata,         # noqa: F821
        private_mapping_json=snakemake.input.private_mapping,  # noqa: F821
        input_dir=snakemake.params.input_dir,          # noqa: F821
        output_mapping=snakemake.output.mapping,       # noqa: F821
    )
