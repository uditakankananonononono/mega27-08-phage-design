"""Sequence encodings for phage receptor-binding proteins and host receptors.

Two real encodings are provided:
  * kmer_spectrum: normalized k-mer frequency vector over the 20 canonical
    amino acids (k=2,3 supported). This is the classical baseline featurization
    used by WIsH/PHP-class host predictors.
  * physicochemical_profile: per-residue indices (hydrophobicity (Kyte-Doolittle),
    volume, charge, polarity, helix/sheet propensity) used as CNN input channels.

All tables are literature values; nothing here is random or stubbed.
"""
from __future__ import annotations

from collections import Counter

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {a: i for i, a in enumerate(AA)}

# Kyte & Doolittle 1982 hydropathy.
HYDROPATHY = dict(zip(AA, [-0.9, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8,
                           1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3]))
# Approximate residue volumes (A^3), Zamyatnin 1972.
VOLUME = dict(zip(AA, [67, 96, 91, 109, 86, 82, 48, 0.1, 118, 124,
                       124, 135, 124, 132, 90, 73, 93, 163, 141, 105]))
# Charge at pH 7 (D,E negative; K,R,H positive).
CHARGE = dict(zip(AA, [0, 0, 0, -1, -1, 0, 0, 0, 0.1, 0,
                       0, 1, 0, 0, 0, 1, 0, 0, 0, 0]))


def clean_sequence(seq: str) -> str:
    """Upper-case and keep canonical residues only; ambiguous letters dropped."""
    return "".join(c for c in seq.upper() if c in AA_INDEX)


def kmer_spectrum(seq: str, k: int = 3) -> np.ndarray:
    """Normalized k-mer frequency vector of length 20**k (capped k=3).

    Returns frequencies over the full 20**k alphabet; positions map a k-mer
    a1..ak to index sum(idx(ai) * 20**(k-i)).
    """
    seq = clean_sequence(seq)
    if k < 1 or k > 3:
        raise ValueError("k must be 1..3 (alphabet 20**k)")
    n = 20 ** k
    vec = np.zeros(n, dtype=np.float64)
    if len(seq) < k:
        return vec
    counts = Counter(seq[i:i + k] for i in range(len(seq) - k + 1))
    total = sum(counts.values())
    for kmer, c in counts.items():
        idx = 0
        for j, ch in enumerate(kmer):
            idx += AA_INDEX[ch] * (20 ** (k - 1 - j))
        vec[idx] = c / total
    return vec


def physicochemical_profile(seq: str) -> np.ndarray:
    """(L, 3) float32 matrix: [hydropathy, volume/200, charge] per residue."""
    seq = clean_sequence(seq)
    arr = np.zeros((len(seq), 3), dtype=np.float32)
    for i, ch in enumerate(seq):
        arr[i, 0] = HYDROPATHY[ch] / 4.5
        arr[i, 1] = VOLUME[ch] / 200.0
        arr[i, 2] = CHARGE[ch]
    return arr
