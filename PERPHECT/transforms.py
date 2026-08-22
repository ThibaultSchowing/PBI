"""
Sequence transformation utilities for PERPHECT model.

Provides one-hot encoding of DNA sequences for use with the CNN architecture.
"""

import numpy as np


NUCLEOTIDE_MAP = {
    'A': 0,
    'T': 1,
    'G': 2,
    'C': 3,
}

# Lookup table: maps ASCII byte value (0-255) to nucleotide index.
# Non-ATGC bytes (e.g. N, ambiguous codes) map to 4 (treated as unknown).
_BASE_LUT = np.full(256, 4, dtype=np.int8)
_BASE_LUT[ord('A')] = 0
_BASE_LUT[ord('T')] = 1
_BASE_LUT[ord('G')] = 2
_BASE_LUT[ord('C')] = 3


def translate_sequence_onehot(sequence: str) -> np.ndarray:
    """
    Convert a DNA sequence string to a one-hot encoded numpy array.

    Uses vectorized numpy operations for performance on long sequences.

    Args:
        sequence: DNA sequence string containing characters A, T, G, C, N.
                  Unknown bases (N, ambiguous IUPAC codes, etc.) are encoded
                  as zero vectors.

    Returns:
        numpy array of shape (len(sequence), 4) with dtype uint8.
        Columns correspond to [A, T, G, C].

    Example:
        >>> translate_sequence_onehot("ATG")
        array([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 1, 0]], dtype=uint8)
    """
    seq_len = len(sequence)
    if seq_len == 0:
        return np.zeros((0, 4), dtype=np.uint8)

    # Convert string to uppercase byte array and look up nucleotide indices
    indices = _BASE_LUT[np.frombuffer(sequence.upper().encode('ascii'), dtype=np.uint8)]

    # One-hot encode: set (row, col) = 1 for valid bases (index < 4)
    onehot = np.zeros((seq_len, 4), dtype=np.uint8)
    valid = indices < 4
    onehot[np.arange(seq_len)[valid], indices[valid]] = 1

    return onehot
