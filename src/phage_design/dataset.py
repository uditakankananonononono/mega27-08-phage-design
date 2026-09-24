"""Supervised phage-host pair dataset.

Positive pair: phage x its GenBank-annotated host species receptor panel.
Negative pair: the same phage x the other species' panel (species-specificity
assumption; polyvalent exceptions are discussed in the paper, not hidden).

Each sample indexes into a shared protein table so batched training encodes
each distinct protein once per batch.
"""
from __future__ import annotations

import numpy as np
import torch
from torch.utils.data import Dataset

from .features.sequence import physicochemical_profile


class ProteinBank:
    """Padded physicochemical profiles for every protein in the task."""

    def __init__(self, sequences: list[str], max_len: int = 600):
        self.max_len = max_len
        self.seq_count = len(sequences)
        mats = []
        for s in sequences:
            m = physicochemical_profile(s[:max_len])
            pad = np.zeros((max_len, 3), dtype=np.float32)
            pad[:len(m)] = m
            mats.append(pad)
        self.tensor = torch.from_numpy(np.stack(mats))  # (N, max_len, 3)

    def __len__(self) -> int:
        return self.seq_count


class PairDataset(Dataset):
    """(phage RBP idx list, host receptor idx list, label) triples."""

    def __init__(self, samples: list[tuple[list[int], list[int], int]]):
        self.samples = samples

    def __len__(self) -> int:
        return len(self.samples)

    def __getitem__(self, i: int):
        return self.samples[i]


def collate_pairs(batch, bank: ProteinBank):
    """Gather unique proteins for the batch, one padded tensor per side."""
    p_idx = sorted({i for s in batch for i in s[0]})
    h_idx = sorted({i for s in batch for i in s[1]})
    p_map = {g: j for j, g in enumerate(p_idx)}
    h_map = {g: j for j, g in enumerate(h_idx)}
    p_groups = [torch.tensor([p_map[i] for i in s[0]]) for s in batch]
    h_groups = [torch.tensor([h_map[i] for i in s[1]]) for s in batch]
    labels = torch.tensor([s[2] for s in batch], dtype=torch.float32)
    return (bank.tensor[p_idx], bank.tensor[h_idx], p_groups, h_groups, labels)
