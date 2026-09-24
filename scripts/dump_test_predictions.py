"""Dump per-phage held-out predictions of the seed-7 model for error analysis."""
import json
import sys
from pathlib import Path

import numpy as np
import torch
from sklearn.model_selection import GroupShuffleSplit

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from train_interaction import SPECIES, MAX_LEN, EMB, load_task
from phage_design.dataset import ProteinBank, collate_pairs
from phage_design.features.sequence import kmer_spectrum
from phage_design.models.cnn import ProteinCNNEncoder
from phage_design.models.grouped import GroupedInteractionModel

SEED = 7
torch.manual_seed(SEED)
np.random.seed(SEED)
phages, proteins, samples, groups, species_receptors, host_names, rbp_lists = load_task()
idx = np.arange(len(samples))
gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=SEED)
tr, te = next(gss.split(idx, groups=groups))
gss2 = GroupShuffleSplit(n_splits=1, test_size=0.125, random_state=SEED + 1)
tr_sub, va = next(gss2.split(tr, groups=[groups[i] for i in tr]))
tr_idx = tr[tr_sub]

bank = ProteinBank(proteins, max_len=MAX_LEN)
panel_spec = {sp: np.mean([kmer_spectrum(x, k=2) for x in seqs], axis=0)
              for sp, seqs in species_receptors.items()}
phage_spec = np.stack([np.mean([kmer_spectrum(x, k=2) for x in rbps], axis=0)
                       for rbps in rbp_lists])
host_spec = np.stack([panel_spec[sp] for sp in host_names])
kmer_feats = np.concatenate([phage_spec * host_spec, phage_spec - host_spec], axis=1)
mu, sd = kmer_feats[tr_idx].mean(0), kmer_feats[tr_idx].std(0) + 1e-8
kmer_feats = torch.from_numpy(((kmer_feats - mu) / sd).astype(np.float32))

model = GroupedInteractionModel(ProteinCNNEncoder(23, 48, EMB),
                                ProteinCNNEncoder(23, 48, EMB), emb_dim=EMB)
model.load_state_dict(torch.load("results/interaction_cnn_seed7.pt", map_location="cpu"))
model.eval()

# pair test rows into per-phage (pos, neg) records
by_phage = {}
for j in te:
    p_idx, _, label = samples[j]
    acc = groups[j]
    by_phage.setdefault(acc, {})[label] = j

rows = []
with torch.no_grad():
    for acc, pair in by_phage.items():
        logits = {}
        for label, j in pair.items():
            pb, hb, pg, hg, _ = collate_pairs([samples[j]], bank)
            logits[label] = float(model(pb, hb, pg, hg, kmer_feats[j:j+1]).item())
        truth = host_names[pair[1]]
        p_true = 1 / (1 + np.exp(-logits[1]))
        p_other = 1 / (1 + np.exp(-logits[0]))
        rows.append({"accession": acc, "truth": truth,
                     "p_true_species": round(p_true, 4), "p_other_species": round(p_other, 4),
                     "correct": bool(p_true > 0.5), "choice_correct": bool(p_true > p_other)})
json.dump(rows, open("results/own_test_predictions_seed7.json", "w"), indent=2)
wrong = [r for r in rows if not r["correct"]]
print(f"test phages={len(rows)} wrong={len(wrong)} choice_wrong={sum(1 for r in rows if not r['choice_correct'])}")
for r in sorted(wrong, key=lambda r: r["p_true_species"]):
    print(r["accession"], r["truth"], "p_true=%.3f p_other=%.3f" % (r["p_true_species"], r["p_other_species"]))
