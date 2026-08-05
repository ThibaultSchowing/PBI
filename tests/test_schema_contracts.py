#!/usr/bin/env python3

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts" / "preprocessing" / "mergers"))
from schema_contracts import load_contract, normalize_df_schema


def _write_contract(tmp_path):
    contract_path = tmp_path / "contract.yaml"
    contract_path.write_text(
        """
required:
  - id
optional:
  - value
aliases:
  old_id: id
defaults:
  value: default_value
""".strip()
    )
    return load_contract(contract_path)


def test_alias_rename_works(tmp_path):
    contract = _write_contract(tmp_path)
    df = pd.DataFrame({" old_id ": ["A"]})

    normalized, report = normalize_df_schema(df, contract)

    assert "id" in normalized.columns
    assert normalized.loc[0, "id"] == "A"
    assert report["aliases_applied"] == [{"from": "old_id", "to": "id"}]


def test_missing_required_fails(tmp_path):
    contract = _write_contract(tmp_path)
    df = pd.DataFrame({"other": [1]})

    with pytest.raises(ValueError, match="Missing required columns"):
        normalize_df_schema(df, contract)


def test_missing_optional_added_as_default(tmp_path):
    contract = _write_contract(tmp_path)
    df = pd.DataFrame({"id": ["A"]})

    normalized, report = normalize_df_schema(df, contract)

    assert "value" in normalized.columns
    assert normalized.loc[0, "value"] == "default_value"
    assert report["added_optional"] == ["value"]


def test_unknown_columns_preserved(tmp_path):
    contract = _write_contract(tmp_path)
    df = pd.DataFrame({"id": ["A"], "new_field": ["x"]})

    normalized, report = normalize_df_schema(df, contract)

    assert "new_field" in normalized.columns
    assert report["unknown_columns"] == ["new_field"]


def test_deterministic_column_ordering(tmp_path):
    contract = _write_contract(tmp_path)
    df = pd.DataFrame({"zeta": [1], "id": ["A"], "alpha": [2]})

    normalized, _ = normalize_df_schema(df, contract)

    assert list(normalized.columns) == ["id", "value", "alpha", "zeta"]


def test_numerical_columns_loaded(tmp_path):
    """Verify that the numerical_columns key is read from the schema YAML."""
    contract_path = tmp_path / "contract.yaml"
    contract_path.write_text(
        """
required:
  - id
optional:
  - value
  - count
numerical_columns:
  - count
aliases: {}
defaults: {}
""".strip()
    )
    contract = load_contract(contract_path)
    assert contract["numerical_columns"] == ["count"]


def test_numerical_columns_defaults_to_empty(tmp_path):
    """Verify that a missing numerical_columns key defaults to an empty list."""
    contract = _write_contract(tmp_path)
    assert contract["numerical_columns"] == []
