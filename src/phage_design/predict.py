"""Inference library: host-range prediction and receptor screen for a new phage.

Loads the trained CNN-BILIN checkpoint (committed at results/interaction_cnn.pt)
and scores a phage genome's receptor-binding proteins against the receptor
panels. Offline: no network. The k-mer interaction features are z-scored with
corpus statistics saved at results/kmer_norm.json (written by
scripts/make_kmer_norm.py from the full model corpus).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import torch

from .data.genbank_parse import parse_genbank
from .dataset import ProteinBank, collate_pairs
from .features.sequence import kmer_spectrum
from .models.cnn import ProteinCNNEncoder
from .models.grouped import GroupedInteractionModel

EMB = 64
MAX_LEN = 600
MAX_RBP = 4


def load_panel(panel_csv: str | Path) -> dict[str, list[tuple[str, str]]]:
    """species -> [(receptor_name, sequence)] from the committed panel CSV."""
    import csv
    panel: dict[str, list[tuple[str, str]]] = {}
    with open(panel_csv) as fh:
        for row in csv.DictReader(fh):
            panel.setdefault(row["species"], []).append((row["receptor"], row["sequence"]))
    if not panel:
        raise ValueError(f"no receptors loaded from {panel_csv}")
    return panel


def load_model(checkpoint: str | Path) -> GroupedInteractionModel:
    model = GroupedInteractionModel(
        ProteinCNNEncoder(23, 48, EMB), ProteinCNNEncoder(23, 48, EMB), emb_dim=EMB)
    state = torch.load(checkpoint, map_location="cpu")
    model.load_state_dict(state)
    model.eval()
    return model


def load_norm(norm_json: str | Path) -> tuple[np.ndarray, np.ndarray]:
    d = json.loads(Path(norm_json).read_text())
    return np.asarray(d["mu"], dtype=np.float64), np.asarray(d["sd"], dtype=np.float64)


def _kmer_features(rbps: list[str], panel_seqs: list[str],
                   mu: np.ndarray, sd: np.ndarray) -> np.ndarray:
    phage_spec = np.mean([kmer_spectrum(x, k=2) for x in rbps], axis=0)
    host_spec = np.mean([kmer_spectrum(x, k=2) for x in panel_seqs], axis=0)
    feats = np.concatenate([phage_spec * host_spec, phage_spec - host_spec])
    return ((feats - mu) / sd).astype(np.float32)


def rbps_from_genbank(genbank_path: str | Path) -> list[tuple[str, str]]:
    rec = parse_genbank(Path(genbank_path).read_text(errors="replace"))
    if not rec.rbp_products:
        raise ValueError(f"no receptor-binding proteins annotated in {genbank_path}")
    return rec.rbp_products[:MAX_RBP]


def score_against(model, rbps: list[str], host_seqs: list[str],
                  mu: np.ndarray, sd: np.ndarray) -> float:
    """Logit of the (RBP set x host sequence set) pair."""
    proteins = list(rbps) + list(host_seqs)
    p_idx = list(range(len(rbps)))
    h_idx = list(range(len(rbps), len(proteins)))
    bank = ProteinBank(proteins, max_len=MAX_LEN)
    samples = [(p_idx, h_idx, 1)]
    kf = torch.from_numpy(_kmer_features(rbps, host_seqs, mu, sd)).unsqueeze(0)
    with torch.no_grad():
        (pb, hb, pg, hg, _) = collate_pairs(samples, bank)
        return float(model(pb, hb, pg, hg, kf).item())


def predict_host_range(model, panel: dict, rbps: list[str],
                       mu: np.ndarray, sd: np.ndarray) -> dict:
    """Score the phage against every species panel; return ranked probabilities."""
    logits = {sp: score_against(model, rbps, [s for _, s in seqs], mu, sd)
              for sp, seqs in panel.items()}
    exps = {sp: float(1.0 / (1.0 + np.exp(-l))) for sp, l in logits.items()}
    ranking = sorted(exps, key=exps.get, reverse=True)
    return {"logits": {k: round(v, 6) for k, v in logits.items()},
            "probabilities": {k: round(v, 6) for k, v in exps.items()},
            "ranking": ranking,
            "predicted_host": ranking[0]}


def receptor_screen(model, panel: dict, rbps: list[str],
                    mu: np.ndarray, sd: np.ndarray) -> list[dict]:
    """Score the RBP set against each individual receptor (discovery screen)."""
    rows = []
    for sp, seqs in panel.items():
        for name, seq in seqs:
            logit = score_against(model, rbps, [seq], mu, sd)
            rows.append({"species": sp, "receptor": name,
                         "logit": round(logit, 6),
                         "probability": round(float(1.0 / (1.0 + np.exp(-logit))), 6)})
    rows.sort(key=lambda r: r["logit"], reverse=True)
    return rows
