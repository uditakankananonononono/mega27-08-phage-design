"""Graph neural network encoder over residue contact graphs (pure torch).

Graphs are built from experimentally determined or AlphaFold-DB C-alpha
coordinates: residues are nodes, edges join residues closer than a cutoff
(default 10 A), node features are the physicochemical profile. Message passing
follows the MPNN formulation (Gilmer et al. 2017):

    m_i^(t+1) = sum_{j in N(i)} M_t(h_i^t, h_j^t, e_ij)
    h_i^(t+1) = U_t(h_i^t, m_i^(t+1))

with M_t a linear map on concatenated endpoints and U_t a GRU-style gated
update (here a residual linear update, which keeps the 2-CPU budget honest).
A derivation of the permutation-invariance of this readout lives in
paper/derivations.md.
"""
from __future__ import annotations

import torch
import torch.nn as nn


def knn_edges(coords: torch.Tensor, cutoff: float = 10.0) -> torch.Tensor:
    """(L,3) coordinates -> (2, E) edge index with ||xi - xj|| < cutoff, i != j."""
    d = torch.cdist(coords, coords)
    mask = (d < cutoff) & (d > 0)
    return mask.nonzero(as_tuple=False).t().contiguous()


class ResidueGNNEncoder(nn.Module):
    def __init__(self, node_dim: int = 3, hidden: int = 64, emb_dim: int = 128,
                 n_layers: int = 3):
        super().__init__()
        self.in_proj = nn.Linear(node_dim, hidden)
        self.msg = nn.ModuleList(nn.Linear(2 * hidden, hidden) for _ in range(n_layers))
        self.upd = nn.ModuleList(nn.Linear(2 * hidden, hidden) for _ in range(n_layers))
        self.readout = nn.Sequential(nn.Linear(hidden, emb_dim), nn.ReLU())

    def forward(self, x: torch.Tensor, edge_index: torch.Tensor) -> torch.Tensor:
        """x: (L, node_dim); edge_index: (2, E). Returns (emb_dim,) embedding."""
        h = torch.relu(self.in_proj(x))
        for msg, upd in zip(self.msg, self.upd):
            src, dst = edge_index[0], edge_index[1]
            m = msg(torch.cat([h[src], h[dst]], dim=-1))       # (E, hidden)
            agg = torch.zeros_like(h).index_add_(0, dst, m)    # sum over N(i)
            deg = torch.zeros(h.size(0), 1, device=h.device)
            deg.index_add_(0, dst, torch.ones(dst.size(0), 1, device=h.device))
            agg = agg / deg.clamp(min=1)                       # mean-normalize
            h = torch.relu(h + upd(torch.cat([h, agg], dim=-1)))
        return self.readout(h.mean(dim=0, keepdim=True)).squeeze(0)
