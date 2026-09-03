#!/usr/bin/env python3
"""
Index individual host genome FASTA files using pyfaidx

This script reads the host mapping JSON file, runs FASTA quality-control
audits on each genome, and creates ``.fai`` index files for files that pass
the header-uniqueness check.

QC procedure (non-destructive — files are never modified):
1. **Header audit**: every sequence identifier must be unique.  If duplicate
   headers are detected the file is *rejected* (not indexed) and the event is
   recorded in the QC log.
2. **Sequence content audit**: identical sequences (same content, any header)
   trigger a *warning* that is recorded in the QC log, but the file is still
   indexed — no sequences are discarded.

The QC log (``host_fasta_qc_log``) is a CSV file that can be loaded as a
DataFrame for downstream analysis.
"""

import csv
import json
import logging
import os
import sys
from pathlib import Path

from pyfaidx import Fasta

from fasta_qc import run_qc

# ---------------------------------------------------------------------------
# Shared logging helper (delegates to workflow/scripts/common/logging_utils.py
# when available; falls back to an inline equivalent).
# ---------------------------------------------------------------------------

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


# Columns written to the CSV QC log — keep in sync with create_host_status_report.py
_QC_LOG_COLUMNS = [
    'host_id',
    'fasta_path',
    'index_existed',
    'total_sequences',
    'header_qc_status',
    'n_duplicate_headers',
    'duplicate_header_examples',
    'seq_qc_status',
    'n_identical_seq_groups',
    'identical_seq_group_examples',
    'index_status',
    'error_message',
]


def index_host_files(mapping_file: str, log_file: str, qc_log_file: str):
    """Index all host FASTA files listed in the mapping file.

    Args:
        mapping_file:  Path to JSON file mapping Host_ID -> FASTA file path.
        log_file:      Path to plain-text summary log (Snakemake convention).
        qc_log_file:   Path to CSV QC log (DataFrame-loadable).
    """
    with open(mapping_file, 'r') as f:
        host_mapping = json.load(f)

    total_files = len(host_mapping)
    logging.info(f"Indexing {total_files} host FASTA files...")

    indexed_count = 0
    skipped_count = 0
    rejected_count = 0
    error_count = 0

    qc_rows = []

    for host_id, fasta_path in host_mapping.items():
        fasta_file = Path(fasta_path)
        index_file = Path(str(fasta_file) + '.fai')

        row = {col: '' for col in _QC_LOG_COLUMNS}
        row['host_id'] = host_id
        row['fasta_path'] = fasta_path
        row['index_existed'] = index_file.exists()

        # --- Already indexed ---------------------------------------------------
        if index_file.exists():
            logging.debug(f"Index already exists for {host_id}")
            row['index_status'] = 'already_indexed'
            skipped_count += 1
            continue

        # --- Run QC audits -----------------------------------------------------
        try:
            qc = run_qc(fasta_path)
            row.update(qc)
        except Exception as exc:
            error_count += 1
            row['index_status'] = 'error'
            row['error_message'] = str(exc)
            logging.error(f"QC failed for {host_id}: {exc}")
            qc_rows.append(row)
            continue

        # --- Reject files with duplicate headers --------------------------------
        if row.get('header_qc_status', '').startswith('rejected'):
            rejected_count += 1
            row['index_status'] = 'rejected_duplicate_headers'
            logging.warning(
                f"Skipping index for {host_id} — duplicate FASTA headers detected"
            )
            qc_rows.append(row)
            continue

        # --- Index (headers are unique) ----------------------------------------
        try:
            with Fasta(str(fasta_file)):
                pass

            if index_file.exists():
                indexed_count += 1
                row['index_status'] = 'indexed'
                if indexed_count % 100 == 0:
                    logging.info(f"Indexed {indexed_count}/{total_files} files...")
            else:
                error_count += 1
                row['index_status'] = 'failed'
                row['error_message'] = 'index file not created after pyfaidx call'
                logging.warning(f"Failed to create index for {host_id}: {fasta_path}")

        except Exception as exc:
            error_count += 1
            row['index_status'] = 'error'
            row['error_message'] = str(exc)
            logging.error(f"Error indexing {host_id}: {exc}")

        qc_rows.append(row)

    # --- Write CSV QC log (DataFrame-loadable) --------------------------------
    Path(qc_log_file).parent.mkdir(parents=True, exist_ok=True)
    with open(qc_log_file, 'w', newline='') as csvf:
        writer = csv.DictWriter(csvf, fieldnames=_QC_LOG_COLUMNS)
        writer.writeheader()
        writer.writerows(qc_rows)

    # --- Write plain-text summary log (Snakemake convention) ------------------
    with open(log_file, 'w') as f:
        f.write(f"Total host FASTA files: {total_files}\n")
        f.write(f"Newly indexed: {indexed_count}\n")
        f.write(f"Already indexed (skipped): {skipped_count}\n")
        f.write(f"Rejected (duplicate headers): {rejected_count}\n")
        f.write(f"Errors: {error_count}\n")
        f.write(f"\nQC log: {qc_log_file}\n")

    logging.info(
        f"Indexing complete: {indexed_count} newly indexed, "
        f"{skipped_count} pre-existing, "
        f"{rejected_count} rejected, {error_count} errors"
    )
    if rejected_count > 0:
        logging.warning(
            f"{rejected_count} file(s) rejected due to duplicate FASTA headers — "
            f"see QC log: {qc_log_file}"
        )


if __name__ == "__main__":
    # Snakemake variables
    mapping_file = snakemake.input.mapping          # noqa: F821
    log_file = snakemake.log[0]                     # noqa: F821
    qc_log_file = snakemake.output.qc_log           # noqa: F821

    _setup_logging(log_file)
    index_host_files(mapping_file, log_file, qc_log_file)
