"""Discovery screen: candidate host-switch and broad-host-range phages, with
per-receptor mediation calls.

For every phage with RBPs and a mapped host, score compatibility against:
  * each species panel (pooled), and
  * each individual receptor (per-receptor mediation).

Outputs (results/discovery_screen.csv):
  accession, annotated host, P(E. coli), P(Klebsiella), margin, call
  (confirmed-host / candidate-host-switch / candidate-broad-range), and the
  top predicted receptor per species. Every call is falsifiable by plaque
  assay on the named strain panel; candidates are stated as predictions,
  never as established host ranges.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from phage_design.dataset import PairDataset, ProteinBank, collate_pairs
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

    # protein bank: all phage RBPs + all panel receptors
    receptors = panel.sequence.tolist()
    rec_meta = panel[["species", "receptor", "accession"]].to_dict("records")
    rec_offset = 0
    proteins, pidx = [], {}

    def idx_of(seq):
        nonlocal proteins
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

    # k-mer features must match training standardization; recompute stats on
    # the full task (scoring-only use; no labels involved)
    panel_spec = {sp: np.mean([kmer_spectrum(x, 2) for x in
                               panel[panel.species == sp].sequence], axis=0)
                  for sp in SPECIES}

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
            species_prob, per_receptor = {}, {}
            for sp in SPECIES:
                # pooled panel score
                feat = kmer_feat(rbps, panel[panel.species == sp].sequence.tolist())
                feat = torch.from_numpy(((feat - mu) / sd)).unsqueeze(0)
                sample = [(p_idx, [rec_idx[i] for i in species_idx[sp]], 0)]
                pb, hb, pg, hg, _ = collate_pairs(sample, bank)
                logit = model(pb, hb, pg, hg, feat).item()
                species_prob[sp] = 1 / (1 + np.exp(-logit))
                # per-receptor mediation
                for i in species_idx[sp]:
                    feat_r = kmer_feat(rbps, [receptors[i]])
                    feat_r = torch.from_numpy(((feat_r - mu) / sd)).unsqueeze(0)
                    sample = [(p_idx, [rec_idx[i]], 0)]
                    pb, hb, pg, hg, _ = collate_pairs(sample, bank)
                    per_receptor[f"{sp}:{rec_meta[i]['receptor']}"] = \
                        1 / (1 + np.exp(-model(pb, hb, pg, hg, feat_r).item()))
            pe, pk = species_prob[SPECIES[0]], species_prob[SPECIES[1]]
            annotated = row.host_species
            own = pe if annotated == SPECIES[0] else pk
            other = pk if annotated == SPECIES[0] else pe
            if other > own + 0.2:
                call = "candidate-host-switch"
            elif own > 0.7 and other > 0.7:
                call = "candidate-broad-range"
            else:
                call = "confirmed-host"
            top_rec = {sp: max(((k, v) for k, v in per_receptor.items() if k.startswith(sp)),
                               key=lambda kv: kv[1]) for sp in SPECIES}
            rows.append({
                "accession": acc, "annotated_host": annotated,
                "p_ecoli": round(pe, 4), "p_klebsiella": round(pk, 4),
                "margin": round(own - other, 4), "call": call,
                "top_receptor_ecoli": f"{top_rec[SPECIES[0]][0].split(':')[1]}:{top_rec[SPECIES[0]][1]:.3f}",
                "top_receptor_klebsiella": f"{top_rec[SPECIES[1]][0].split(':')[1]}:{top_rec[SPECIES[1]][1]:.3f}",
            })
    df = pd.DataFrame(rows)
    df.to_csv("results/discovery_screen.csv", index=False)
    print(df.call.value_counts().to_string())
    print("\nTop host-switch candidates:")
    print(df[df.call == "candidate-host-switch"].sort_values("margin").head(10).to_string(index=False))


if __name__ == "__main__":
    main()
