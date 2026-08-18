#!/usr/bin/env python3
"""
Convert PBI_Perphect_Training.ipynb to a Python script.

Extracts only code cells from the notebook, producing a runnable .py file
that mirrors the notebook's logic. Useful for:
  - Verifying notebook and script produce identical results
  - Running the notebook workflow as a standalone script
  - Code review in a traditional diff-friendly format

Usage:
    python notebook_to_script.py
    python notebook_to_script.py --input Notebook.ipynb --output custom_output.py
"""

import argparse
import json
import re
import sys
from pathlib import Path


def convert_notebook(nb_path: str, py_path: str) -> int:
    """
    Convert a Jupyter notebook to a Python script.

    Args:
        nb_path: Path to the input .ipynb file.
        py_path: Path to the output .py file.

    Returns:
        Number of code cells extracted.
    """
    with open(nb_path) as f:
        nb = json.load(f)

    lines = [
        '#!/usr/bin/env python3',
        '"""',
        f'Auto-generated from {Path(nb_path).name}',
        '',
        'Run with: python {Path(py_path).name}',
        '"""',
        '',
    ]

    cell_count = 0
    for cell in nb["cells"]:
        if cell["cell_type"] == "code":
            source = "".join(cell["source"])
            if source.strip():
                lines.append(f"# --- Cell {cell_count + 1} ---")
                lines.append(source)
                lines.append("")
                cell_count += 1

    with open(py_path, "w", newline="\n") as f:
        f.write("\n".join(lines))

    return cell_count


def verify_scripts_match(nb_script: str, ref_script: str) -> list:
    """
    Compare two Python scripts for logical equivalence.

    Ignores comments, blank lines, and minor whitespace differences.
    Returns a list of differences found.
    """
    def normalize(path):
        with open(path) as f:
            lines = f.readlines()
        result = []
        for line in lines:
            stripped = line.strip()
            # Skip comments, blank lines, shebang, docstrings
            if not stripped or stripped.startswith("#") or stripped.startswith('"""'):
                continue
            # Normalize whitespace
            result.append(re.sub(r'\s+', ' ', stripped))
        return result

    nb_lines = normalize(nb_script)
    ref_lines = normalize(ref_script)

    differences = []
    if len(nb_lines) != len(ref_lines):
        differences.append(
            f"Line count differs: notebook script={len(nb_lines)}, "
            f"reference={len(ref_lines)}"
        )

    # Check first N lines (where N is min)
    min_len = min(len(nb_lines), len(ref_lines))
    for i in range(min_len):
        if nb_lines[i] != ref_lines[i]:
            differences.append(
                f"Line {i + 1} differs:\n"
                f"  notebook: {nb_lines[i][:80]}\n"
                f"  reference: {ref_lines[i][:80]}"
            )

    return differences


def main():
    parser = argparse.ArgumentParser(
        description="Convert Jupyter notebook to Python script"
    )
    parser.add_argument(
        "--input", "-i",
        default="PBI_Perphect_Training.ipynb",
        help="Input notebook (default: PBI_Perphect_Training.ipynb)",
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        help="Output script (default: train_notebook.py)",
    )
    parser.add_argument(
        "--verify",
        default=None,
        help="Verify against reference script (e.g., train.py)",
    )

    args = parser.parse_args()

    nb_path = Path(args.input)
    if not nb_path.exists():
        print(f"Error: {nb_path} not found")
        sys.exit(1)

    py_path = Path(args.output) if args.output else nb_path.with_suffix(".py").name.replace(
        "PBI_Perphect_Training", "train_notebook"
    )

    cell_count = convert_notebook(str(nb_path), str(py_path))
    print(f"Converted {nb_path} -> {py_path} ({cell_count} code cells)")

    if args.verify:
        ref_path = Path(args.verify)
        if not ref_path.exists():
            print(f"Warning: reference script {ref_path} not found, skipping verification")
            return

        differences = verify_scripts_match(str(py_path), str(ref_path))
        if differences:
            print(f"\nDifferences found ({len(differences)}):")
            for diff in differences[:10]:  # Show first 10
                print(f"  - {diff}")
            sys.exit(1)
        else:
            print("Verification passed: scripts are logically equivalent")


if __name__ == "__main__":
    main()
