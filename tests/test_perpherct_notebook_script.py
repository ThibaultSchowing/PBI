"""
Tests for PERPHECT notebook_to_script.py utility.
"""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "PERPHECT"))


def _create_test_notebook(tmp_path):
    """Create a minimal test notebook."""
    nb = {
        "cells": [
            {
                "cell_type": "markdown",
                "source": ["# Test Notebook\n"],
                "metadata": {},
            },
            {
                "cell_type": "code",
                "source": ["import os\n", "print('hello')"],
                "metadata": {},
                "outputs": [],
                "execution_count": None,
            },
            {
                "cell_type": "code",
                "source": ["x = 42\n", "print(x)"],
                "metadata": {},
                "outputs": [],
                "execution_count": None,
            },
            {
                "cell_type": "code",
                "source": [],
                "metadata": {},
                "outputs": [],
                "execution_count": None,
            },
        ],
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        },
        "nbformat": 4,
        "nbformat_minor": 4,
    }
    nb_path = tmp_path / "test_notebook.ipynb"
    with open(nb_path, "w") as f:
        json.dump(nb, f)
    return nb_path


class TestConvertNotebook:
    """Tests for convert_notebook function."""

    def test_extracts_code_cells(self, tmp_path):
        from notebook_to_script import convert_notebook
        nb_path = _create_test_notebook(tmp_path)
        py_path = tmp_path / "output.py"
        cell_count = convert_notebook(str(nb_path), str(py_path))
        assert cell_count == 2

    def test_output_file_exists(self, tmp_path):
        from notebook_to_script import convert_notebook
        nb_path = _create_test_notebook(tmp_path)
        py_path = tmp_path / "output.py"
        convert_notebook(str(nb_path), str(py_path))
        assert py_path.exists()

    def test_output_has_shebang(self, tmp_path):
        from notebook_to_script import convert_notebook
        nb_path = _create_test_notebook(tmp_path)
        py_path = tmp_path / "output.py"
        convert_notebook(str(nb_path), str(py_path))
        content = py_path.read_text()
        assert content.startswith("#!/usr/bin/env python3")

    def test_output_excludes_markdown(self, tmp_path):
        from notebook_to_script import convert_notebook
        nb_path = _create_test_notebook(tmp_path)
        py_path = tmp_path / "output.py"
        convert_notebook(str(nb_path), str(py_path))
        content = py_path.read_text()
        assert "# Test Notebook" not in content

    def test_output_includes_code(self, tmp_path):
        from notebook_to_script import convert_notebook
        nb_path = _create_test_notebook(tmp_path)
        py_path = tmp_path / "output.py"
        convert_notebook(str(nb_path), str(py_path))
        content = py_path.read_text()
        assert "import os" in content
        assert "x = 42" in content

    def test_with_real_notebook(self):
        """Test conversion with the actual PERPHECT notebook."""
        from notebook_to_script import convert_notebook
        nb_path = Path(__file__).parent.parent / "PERPHECT" / "PBI_Perphect_Training.ipynb"
        if not nb_path.exists():
            pytest.skip("Notebook not found")
        py_path = nb_path.parent / "train_notebook.py"
        cell_count = convert_notebook(str(nb_path), str(py_path))
        assert cell_count > 0
        content = py_path.read_text()
        assert "quick_connect" in content
        assert "PBIAdapter" in content


class TestVerifyScripts:
    """Tests for verify_scripts_match function."""

    def test_identical_scripts_pass(self, tmp_path):
        from notebook_to_script import verify_scripts_match
        script = "#!/usr/bin/env python3\nimport os\nx = 42\n"
        (tmp_path / "a.py").write_text(script)
        (tmp_path / "b.py").write_text(script)
        diffs = verify_scripts_match(str(tmp_path / "a.py"), str(tmp_path / "b.py"))
        assert diffs == []

    def test_different_scripts_fail(self, tmp_path):
        from notebook_to_script import verify_scripts_match
        (tmp_path / "a.py").write_text("x = 1\n")
        (tmp_path / "b.py").write_text("x = 2\n")
        diffs = verify_scripts_match(str(tmp_path / "a.py"), str(tmp_path / "b.py"))
        assert len(diffs) > 0

    def test_ignores_comments(self, tmp_path):
        from notebook_to_script import verify_scripts_match
        (tmp_path / "a.py").write_text("# comment\nx = 1\n")
        (tmp_path / "b.py").write_text("# different comment\nx = 1\n")
        diffs = verify_scripts_match(str(tmp_path / "a.py"), str(tmp_path / "b.py"))
        assert diffs == []
