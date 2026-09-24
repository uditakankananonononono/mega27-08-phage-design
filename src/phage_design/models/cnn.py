"""1D CNN encoder over physicochemical residue profiles (torch, CPU).

Architecture follows the DeepHost-class design: stacked temporal convolutions
with increasing receptive fields capture local receptor-binding motifs
(e.g. the loop regions of tail fibers), global max pooling makes the encoder
length-invariant, and a linear head produces a fixed embedding.
"""
from __future__ import annotations

import torch
import torch.nn as nn


class ProteinCNNEncoder(nn.Module):
    def __init__(self, in_channels: int = 3, hidden: int = 64, emb_dim: int = 128):
        super().__init__()
        self.convs = nn.Sequential(
            nn.Conv1d(in_channels, hidden, kernel_size=5, padding=2),
            nn.ReLU(),
            nn.Conv1d(hidden, hidden, kernel_size=5, padding=2),
            nn.ReLU(),
            nn.Conv1d(hidden, hidden, kernel_size=7, padding=3),
            nn.ReLU(),
        )
        self.head = nn.Sequential(
            nn.Linear(2 * hidden, emb_dim),
            nn.ReLU(),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """x: (batch, L, C) -> embedding (batch, emb_dim) via masked max-pool."""
        mask = (x.abs().sum(dim=-1, keepdim=True) > 0).float()  # (B, L, 1)
        h = self.convs(x.transpose(1, 2))                       # (B, hidden, L)
        h = h.transpose(1, 2) * mask                            # zero padded rows
        lengths = mask.sum(dim=1).clamp(min=1.0)                # (B, 1)
        mean_pooled = h.sum(dim=1) / lengths                    # composition
        h_max = h.masked_fill(mask == 0, -1e9).max(dim=1).values
        max_pooled = torch.where(torch.isfinite(h_max), h_max, torch.zeros_like(h_max))
        return self.head(torch.cat([mean_pooled, max_pooled], dim=-1))
