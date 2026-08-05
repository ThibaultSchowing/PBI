#!/usr/bin/env python

# merge_metadata.py — Generic metadata merger.
#
# Replaces the eight per-type merger scripts that were identical apart
# from three hard-coded values (schema contract path, dataset name,
# numerical columns list).  All three are now read from the schema
# contract YAML, so a single script handles every metadata type.
#
# Snakemake params:
#   schema_file  – filename of the schema YAML (e.g. "anti_crispr_metadata_merged.yaml")
#
# The script locates the schemas/ directory by walking up from its own
# file location until it finds the workflow/ parent, which works regardless
# of Docker mount points or temporary script copies.

import sys
import os
import logging
import csv

import numpy as np
import pandas as pd
from pathlib import Path

# Allow imports of shared helpers and schema contracts when run via Snakemake.
sys.path.append("scripts")
import utils
from schema_contracts import load_contract, normalize_df_schema

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def _find_schemas_dir():
    """Locate the schemas/ directory by walking up from this script's location.

    Works both locally and inside Docker containers where __file__ may point
    to a temporary copy.  Searches upward for a 'schemas' directory next to
    the workflow/ folder that contains this script.
    """
    # When run via Snakemake, __file__ points to a temp copy but the original
    # script is at workflow/scripts/preprocessing/mergers/merge_metadata.py.
    # Walk up from __file__ looking for the schemas/ directory.
    current = Path(__file__).resolve().parent
    for _ in range(10):  # safety limit
        # Check if schemas/ exists at this level
        schemas_dir = current / "schemas"
        if schemas_dir.is_dir():
            return schemas_dir
        # Also check inside workflow/ at this level (for Docker layouts)
        workflow_schemas = current / "workflow" / "schemas"
        if workflow_schemas.is_dir():
            return workflow_schemas
        parent = current.parent
        if parent == current:
            break
        current = parent

    raise FileNotFoundError(
        f"Cannot locate schemas/ directory starting from {Path(__file__).resolve()}"
    )


def _read_inputs(inputs):
    """Yield non-empty TSV input paths."""
    for infile in inputs:
        if utils.is_file_empty_or_invalid(infile):
            logging.warning(f"Skipping empty/invalid file: {infile}")
            continue
        yield infile


def _load_and_prepare(infile, contract, dataset_name):
    """Read a single TSV file, inject Source_DB if missing, normalize schema."""
    df = pd.read_csv(infile, sep="\t", quoting=csv.QUOTE_NONNUMERIC)

    if "Source_DB" not in df.columns:
        df["Source_DB"] = Path(infile).stem.split("_")[0]

    df, _ = normalize_df_schema(df, contract, dataset_name=dataset_name, logger=logging)
    return df


def main():
    """Merge all input TSV files into a single output TSV."""
    inputs = snakemake.input
    output = snakemake.output[0]
    schema_file = snakemake.params.schema_file

    # Resolve the full schema path relative to the discovered schemas/ directory.
    schemas_dir = _find_schemas_dir()
    schema_path = schemas_dir / schema_file

    # Load the schema contract — contains column spec AND numerical columns.
    contract = load_contract(schema_path)
    numerical_columns = contract.get("numerical_columns", [])

    # Derive a human-readable dataset name from the schema filename.
    dataset_name = schema_path.stem.replace("_merged", "")

    # Process each input file.
    dfs = []
    for infile in _read_inputs(inputs):
        df = _load_and_prepare(infile, contract, dataset_name)
        df = utils.convert_numerical_columns(df, numerical_columns)
        dfs.append(df)

    # Ensure a consistent column set before merging.
    if dfs:
        final = pd.concat([d.head(0) for d in dfs], ignore_index=True, sort=False)
        final, _ = normalize_df_schema(final, contract, dataset_name=f"{dataset_name}_merged", logger=logging)
        dfs = [d.reindex(columns=final.columns) for d in dfs]

    os.makedirs(os.path.dirname(output), exist_ok=True)
    total = utils.merge_dataframes_chunked(dfs, output)
    logging.info(f"Merged {len(inputs)} files into {output} with {total} rows")


if __name__ == "__main__":
    main()
else:
    # When Snakemake sources this file, execute immediately.
    main()
