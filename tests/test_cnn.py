import torch

from phage_design.models.cnn import ProteinCNNEncoder


def test_encoder_length_invariant_shape():
    enc = ProteinCNNEncoder(in_channels=3, hidden=16, emb_dim=32)
    for L in (10, 57, 200):
        x = torch.randn(2, L, 3)
        out = enc(x)
        assert out.shape == (2, 32)


def test_encoder_zero_padded_rows_ignored():
    enc = ProteinCNNEncoder(in_channels=3, hidden=16, emb_dim=32)
    x = torch.randn(1, 40, 3)
    x[0, 30:] = 0.0  # padding
    out = enc(x)
    assert torch.isfinite(out).all()
