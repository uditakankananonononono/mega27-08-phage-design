"""Batched grouped scoring: encode each protein once, mean-pool per sample."""
from __future__ import annotations

import torch
import torch.nn as nn


class GroupedInteractionModel(nn.Module):
    def __init__(self, phage_encoder: nn.Module, host_encoder: nn.Module, emb_dim: int = 64,
                 kmer_dim: int = 800):
        super().__init__()
        self.phage_encoder = phage_encoder
        self.host_encoder = host_encoder
        self.bilinear = nn.Bilinear(emb_dim, emb_dim, 1)
        # hybrid path: precomputed k-mer interaction features (Hadamard + diff)
        self.kmer_head = nn.Sequential(
            nn.Linear(kmer_dim, 128), nn.ReLU(), nn.Linear(128, 1))
        self.combine = nn.Linear(2, 1)

    @staticmethod
    def _group_mean(emb: torch.Tensor, groups: list[torch.Tensor]) -> torch.Tensor:
        out = torch.stack([emb[g].mean(dim=0) for g in groups])
        return out

    def forward(self, p_batch, h_batch, p_groups, h_groups,
                kmer_feats: torch.Tensor | None = None) -> torch.Tensor:
        ep = self._group_mean(self.phage_encoder(p_batch), p_groups)
        eh = self._group_mean(self.host_encoder(h_batch), h_groups)
        deep_logit = self.bilinear(ep, eh).squeeze(-1)
        if kmer_feats is None:
            return deep_logit
        kmer_logit = self.kmer_head(kmer_feats).squeeze(-1)
        return self.combine(torch.stack([deep_logit, kmer_logit], dim=-1)).squeeze(-1)
