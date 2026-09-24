"""Batched grouped scoring: encode each protein once, mean-pool per sample."""
from __future__ import annotations

import torch
import torch.nn as nn


class GroupedInteractionModel(nn.Module):
    def __init__(self, phage_encoder: nn.Module, host_encoder: nn.Module, emb_dim: int = 64):
        super().__init__()
        self.phage_encoder = phage_encoder
        self.host_encoder = host_encoder
        self.bilinear = nn.Bilinear(emb_dim, emb_dim, 1)

    @staticmethod
    def _group_mean(emb: torch.Tensor, groups: list[torch.Tensor]) -> torch.Tensor:
        out = torch.stack([emb[g].mean(dim=0) for g in groups])
        return out

    def forward(self, p_batch, h_batch, p_groups, h_groups) -> torch.Tensor:
        ep = self._group_mean(self.phage_encoder(p_batch), p_groups)
        eh = self._group_mean(self.host_encoder(h_batch), h_groups)
        return self.bilinear(ep, eh).squeeze(-1)
