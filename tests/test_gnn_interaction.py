import torch

from phage_design.models.gnn import ResidueGNNEncoder, knn_edges
from phage_design.models.interaction import InteractionModel


def test_knn_edges_respects_cutoff():
    coords = torch.tensor([[0., 0., 0.], [5., 0., 0.], [50., 0., 0.]])
    ei = knn_edges(coords, cutoff=10.0)
    pairs = set(map(tuple, ei.t().tolist()))
    assert (0, 1) in pairs and (1, 0) in pairs
    assert (0, 2) not in pairs


def test_gnn_permutation_invariant_readout():
    torch.manual_seed(0)
    gnn = ResidueGNNEncoder(node_dim=3, hidden=16, emb_dim=8, n_layers=2)
    x = torch.randn(12, 3)
    coords = torch.randn(12, 3) * 5
    ei = knn_edges(coords)
    perm = torch.randperm(12)
    e1 = gnn(x, ei)
    # remap edges under permutation
    inv = torch.empty_like(perm)
    inv[perm] = torch.arange(12)
    e2 = gnn(x[perm], inv[ei])
    assert torch.allclose(e1, e2, atol=1e-5)


def test_interaction_logits_shape():
    m = InteractionModel(emb_dim=32)
    px = torch.randn(4, 25, 23)
    hx = torch.randn(4, 60, 23)
    assert m.score(px, hx).shape == (4,)
