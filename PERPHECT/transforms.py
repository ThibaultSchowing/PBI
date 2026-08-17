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


def translate_sequence_onehot(sequence: str) -> np.ndarray:
    """
    Convert a DNA sequence string to a one-hot encoded numpy array.

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
    onehot = np.zeros((seq_len, 4), dtype=np.uint8)

    for i, base in enumerate(sequence.upper()):
        idx = NUCLEOTIDE_MAP.get(base)
        if idx is not None:
            onehot[i, idx] = 1

    return onehot
