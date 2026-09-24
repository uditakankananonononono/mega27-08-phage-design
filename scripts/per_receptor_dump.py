"""Dump the full per-receptor score matrix (all screened phages x 14 receptors)
plus a per-receptor discrimination analysis, for the paper's per-receptor
section. Mirrors scripts/discovery_screen.py scoring exactly."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from phage_design.dataset import ProteinBank, collate_pairs
from phage_design.features.sequence import kmer_spectrum
from phage_design.models.cnn import ProteinCNNEncoder
from phage_design.models.grouped import GroupedInteractionModel

SPECIES = ["Escherichia coli", "Klebsiella pneumoniae"]
MAX_RBP, MAX_LEN, EMB = 4, 600, 64


def main():
    phages = pd.read_csv("data/processed/phage_rbp_table.csv")
    phages = phages[phages.host_species.isin(SPECIES)].reset_index(drop=True)
    panel = pd.read_csv("data/processed/host_receptor_panel.csv")

    model = GroupedInteractionModel(ProteinCNNEncoder(23, 48, EMB),
                                    ProteinCNNEncoder(23, 48, EMB), emb_dim=EMB)
    model.load_state_dict(torch.load("results/interaction_cnn.pt"))
    model.eval()

    receptors = panel.sequence.tolist()
    rec_meta = panel[["species", "receptor", "accession"]].to_dict("records")
    proteins, pidx = [], {}

    def idx_of(seq):
        if seq not in pidx:
            pidx[seq] = len(proteins)
            proteins.append(seq)
        return pidx[seq]

    rec_idx = [idx_of(s) for s in receptors]
    phage_rbps = {}
    for _, row in phages.iterrows():
        rbps = [s.strip() for s in row.rbp_sequences.split(" ||| ")][:MAX_RBP]
        phage_rbps[row.accession] = [idx_of(s) for s in rbps]
    bank = ProteinBank(proteins, max_len=MAX_LEN)

    def kmer_feat(rbps, host_seqs):
        ps = np.mean([kmer_spectrum(x, 2) for x in rbps], axis=0)
        hs = np.mean([kmer_spectrum(x, 2) for x in host_seqs], axis=0)
        return np.concatenate([ps * hs, ps - hs]).astype(np.float32)

    all_feats = np.stack([
        kmer_feat([s.strip() for s in r.rbp_sequences.split(" ||| ")][:MAX_RBP],
                  panel[panel.species == sp].sequence.tolist())
        for _, r in phages.iterrows() for sp in SPECIES
    ])
    mu, sd = all_feats.mean(0), all_feats.std(0) + 1e-8

    species_idx = {sp: [i for i, m in enumerate(rec_meta) if m["species"] == sp]
                   for sp in SPECIES}

    rows = []
    with torch.no_grad():
        for _, row in phages.iterrows():
            acc = row.accession
            p_idx = phage_rbps[acc]
            rbps = [s.strip() for s in row.rbp_sequences.split(" ||| ")][:MAX_RBP]
            out = {"accession": acc, "annotated_host": row.host_species}
            for sp in SPECIES:
                for i in species_idx[sp]:
                    feat_r = kmer_feat(rbps, [receptors[i]])
                    feat_r = torch.from_numpy(((feat_r - mu) / sd)).unsqueeze(0)
                    sample = [(p_idx, [rec_idx[i]], 0)]
                    pb, hb, pg, hg, _ = collate_pairs(sample, bank)
                    key = f"{sp.split()[0][:1]}-{rec_meta[i]['receptor']}"
                    out[key] = round(float(1 / (1 + np.exp(-model(pb, hb, pg, hg, feat_r).item()))), 4)
            rows.append(out)
    df = pd.DataFrame(rows)
    df.to_csv("results/per_receptor_scores.csv", index=False)

    # per-receptor discrimination: for each receptor of species S, compare
    # scores of phages annotated on S vs phages annotated on the other species
    analysis = {}
    for _, m in enumerate(rec_meta):
        name = f"{m['species'].split()[0][:1]}-{m['receptor']}"
        sp = m["species"]
        own = df.loc[df.annotated_host == sp, name].to_numpy()
        other = df.loc[df.annotated_host != sp, name].to_numpy()
        # Mann-Whitney AUC: P(own > other) + 0.5 P(==)
        gt = (own[:, None] > other[None, :]).mean()
        eq = (own[:, None] == other[None, :]).mean()
        analysis[name] = {
            "species": sp,
            "mean_own": round(float(own.mean()), 4),
            "mean_other": round(float(other.mean()), 4),
            "separation": round(float(own.mean() - other.mean()), 4),
            "auc": round(float(gt + 0.5 * eq), 4),
        }
    json.dump(analysis, open("results/per_receptor_analysis.json", "w"), indent=2)
    for k, v in analysis.items():
        print(f"{k:12s} {v['species']:22s} own={v['mean_own']:.3f} other={v['mean_other']:.3f} AUC={v['auc']:.3f}")


if __name__ == "__main__":
    main()
