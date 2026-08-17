"""
Tests for PERPHECT transforms module (translate_sequence_onehot).
"""

import sys
from pathlib import Path

import numpy as np
import pytest

# Add PERPHECT directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "PERPHECT"))


class TestTranslateSequenceOnehot:
    """Tests for the translate_sequence_onehot function."""

    def test_basic_atgc(self):
        from transforms import translate_sequence_onehot

        result = translate_sequence_onehot("ATGC")
        assert result.shape == (4, 4)
        assert result.dtype == np.uint8
        # A -> [1,0,0,0]
        np.testing.assert_array_equal(result[0], [1, 0, 0, 0])
        # T -> [0,1,0,0]
        np.testing.assert_array_equal(result[1], [0, 1, 0, 0])
        # G -> [0,0,1,0]
        np.testing.assert_array_equal(result[2], [0, 0, 1, 0])
        # C -> [0,0,0,1]
        np.testing.assert_array_equal(result[3], [0, 0, 0, 1])

    def test_lowercase_input(self):
        from transforms import translate_sequence_onehot

        result = translate_sequence_onehot("atgc")
        assert result.shape == (4, 4)
        np.testing.assert_array_equal(result[0], [1, 0, 0, 0])
        np.testing.assert_array_equal(result[3], [0, 0, 0, 1])

    def test_unknown_bases(self):
        from transforms import translate_sequence_onehot

        result = translate_sequence_onehot("ANGC")
        # N -> [0,0,0,0] (unknown)
        np.testing.assert_array_equal(result[0], [1, 0, 0, 0])
        np.testing.assert_array_equal(result[1], [0, 0, 0, 0])
        np.testing.assert_array_equal(result[2], [0, 0, 1, 0])
        np.testing.assert_array_equal(result[3], [0, 0, 0, 1])

    def test_empty_sequence(self):
        from transforms import translate_sequence_onehot

        result = translate_sequence_onehot("")
        assert result.shape == (0, 4)

    def test_single_base(self):
        from transforms import translate_sequence_onehot

        result = translate_sequence_onehot("G")
        assert result.shape == (1, 4)
        np.testing.assert_array_equal(result[0], [0, 0, 1, 0])

    def test_long_sequence(self):
        from transforms import translate_sequence_onehot

        seq = "ATCG" * 1000
        result = translate_sequence_onehot(seq)
        assert result.shape == (4000, 4)

    def test_all_ambiguous(self):
        from transforms import translate_sequence_onehot

        result = translate_sequence_onehot("NNNN")
        assert result.shape == (4, 4)
        assert np.sum(result) == 0  # all zeros
