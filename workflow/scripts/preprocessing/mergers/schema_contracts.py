#!/usr/bin/env python

from pathlib import Path
import logging
import pandas as pd
import yaml


def load_contract(path):
    """Load a schema contract YAML file.

    The contract defines the expected column structure for a dataset and
    includes the list of columns that should be coerced to numeric types.
    """
    contract_path = Path(path)
    with contract_path.open("r", encoding="utf-8") as handle:
        contract = yaml.safe_load(handle) or {}

    return {
        "required": list(contract.get("required", [])),
        "optional": list(contract.get("optional", [])),
        "aliases": dict(contract.get("aliases", {})),
        "defaults": dict(contract.get("defaults", {})),
        "numerical_columns": list(contract.get("numerical_columns", [])),
    }


def _log(logger, level, message):
    if logger is None:
        return
    getattr(logger, level)(message)


def normalize_df_schema(df, contract, *, dataset_name=None, logger=None):
    """Normalize a dataframe using a schema contract and return a report."""
    if logger is None:
        logger = logging.getLogger(__name__)

    normalized = df.copy()
    dataset = dataset_name or "dataset"

    required = list(contract.get("required", []))
    optional = list(contract.get("optional", []))
    aliases = dict(contract.get("aliases", {}))
    defaults = dict(contract.get("defaults", {}))

    report = {
        "missing_required": [],
        "missing_optional": [],
        "added_optional": [],
        "unknown_columns": [],
        "aliases_applied": [],
        "collisions": [],
    }

    # Strip whitespace in column names.
    normalized.columns = [str(col).strip() for col in normalized.columns]

    # Resolve duplicate columns after stripping by keeping first seen.
    seen_columns = set()
    duplicate_columns = []
    for col in normalized.columns:
        if col in seen_columns:
            duplicate_columns.append(col)
        seen_columns.add(col)
    if duplicate_columns:
        normalized = normalized.loc[:, ~normalized.columns.duplicated(keep="first")]
        report["collisions"].append(
            {
                "type": "duplicate_column_name",
                "columns": duplicate_columns,
                "resolution": "kept_first",
            }
        )
        _log(logger, "warning", f"[{dataset}] Duplicate columns removed: {duplicate_columns}")

    # Apply aliases with deterministic collision handling.
    for alias_name, canonical_name in aliases.items():
        if alias_name not in normalized.columns:
            continue

        if canonical_name in normalized.columns:
            normalized = normalized.drop(columns=[alias_name])
            collision = {
                "type": "alias_and_canonical_present",
                "alias": alias_name,
                "canonical": canonical_name,
                "resolution": "kept_canonical_dropped_alias",
            }
            report["collisions"].append(collision)
            _log(
                logger,
                "warning",
                f"[{dataset}] Found alias '{alias_name}' and canonical '{canonical_name}'. "
                "Keeping canonical, dropping alias.",
            )
            continue

        normalized = normalized.rename(columns={alias_name: canonical_name})
        report["aliases_applied"].append({"from": alias_name, "to": canonical_name})
        _log(logger, "info", f"[{dataset}] Applied alias '{alias_name}' -> '{canonical_name}'")

    existing_columns = list(normalized.columns)
    missing_required = [col for col in required if col not in existing_columns]
    report["missing_required"] = missing_required
    if missing_required:
        raise ValueError(
            f"[{dataset}] Missing required columns: {missing_required}. "
            f"Existing columns: {existing_columns}"
        )

    missing_optional = [col for col in optional if col not in existing_columns]
    report["missing_optional"] = missing_optional
    for col in missing_optional:
        normalized[col] = defaults.get(col, pd.NA)
        report["added_optional"].append(col)
        _log(logger, "info", f"[{dataset}] Added missing optional column '{col}'")

    contract_columns = set(required + optional)
    report["unknown_columns"] = [col for col in normalized.columns if col not in contract_columns]

    ordered_columns = required + optional + sorted(report["unknown_columns"])
    normalized = normalized.loc[:, ordered_columns]

    return normalized, report
