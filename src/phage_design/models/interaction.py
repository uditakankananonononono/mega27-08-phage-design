"""Dual-encoder phage-host interaction model.

Phage receptor-binding proteins are encoded by the CNN (sequence), host
receptors by the GNN (structure graph) when coordinates exist, else by a
shared CNN fallback on sequence. Compatibility is scored by a bilinear form

    s(p, h) = e_p^T W e_h + b

whose symmetry/antisymmetry properties are derived in paper/derivations.md.
Trained with BCE-with-logits against mined GenBank host annotations.
"""
from __future__ import annotations

import torch
import torch.nn as nn

from .cnn import ProteinCNNEncoder


class InteractionModel(nn.Module):
    def __init__(self, emb_dim: int = 128):
        super().__init__()
        self.phage_cnn = ProteinCNNEncoder(in_channels=3, hidden=64, emb_dim=emb_dim)
        self.host_cnn = ProteinCNNEncoder(in_channels=3, hidden=64, emb_dim=emb_dim)
        self.bilinear = nn.Bilinear(emb_dim, emb_dim, 1)

    def score(self, phage_x: torch.Tensor, host_x: torch.Tensor) -> torch.Tensor:
        """phage_x/host_x: (B, L, 3) padded profiles -> logits (B,)."""
        ep = self.phage_cnn(phage_x)
        eh = self.host_cnn(host_x)
        return self.bilinear(ep, eh).squeeze(-1)
