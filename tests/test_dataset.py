import torch

from phage_design.dataset import PairDataset, ProteinBank, collate_pairs
from phage_design.models.grouped import GroupedInteractionModel
from phage_design.models.cnn import ProteinCNNEncoder


def _bank():
    seqs = ["MKTAAVLLACDEFGHIK", "ACDEFGHIKLMNPQRST", "WYACDEFGHIKLMNPQR"]
    return ProteinBank(seqs, max_len=40)


def test_collate_dedups_proteins():
    bank = _bank()
    ds = PairDataset([([0, 1], [2], 1), ([1], [2], 0)])
    pb, hb, pg, hg, y = collate_pairs([ds[0], ds[1]], bank)
    assert pb.shape == (2, 40, 3)  # proteins 0,1 once
    assert hb.shape == (1, 40, 3)  # protein 2 once
    assert y.tolist() == [1.0, 0.0]


def test_grouped_model_logits():
    bank = _bank()
    ds = PairDataset([([0, 1], [2], 1), ([1], [2], 0)])
    model = GroupedInteractionModel(ProteinCNNEncoder(3, 16, 16),
                                    ProteinCNNEncoder(3, 16, 16), emb_dim=16)
    out = model(*collate_pairs([ds[0], ds[1]], bank)[:4])
    assert out.shape == (2,)
    assert torch.isfinite(out).all()
