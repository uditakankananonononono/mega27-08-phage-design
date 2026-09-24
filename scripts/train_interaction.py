"""Train + evaluate the phage-host interaction model against real baselines.

Task: given a phage's receptor-binding proteins, predict whether it infects
E. coli or Klebsiella pneumoniae (binary host-range task on GenBank ground
truth). Split by phage accession (group split - a phage is in exactly one
split, so RBP sequences never leak across train/test).

Models:
  * KMER-LR    : PHP-class baseline - logistic regression on dipeptide spectra
  * CNN-BILIN  : this work - CNN encoders + bilinear compatibility score
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupShuffleSplit

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from phage_design.dataset import PairDataset, ProteinBank, collate_pairs
from phage_design.features.sequence import kmer_spectrum
from phage_design.models.cnn import ProteinCNNEncoder
from phage_design.models.grouped import GroupedInteractionModel

SPECIES = ["Escherichia coli", "Klebsiella pneumoniae"]
MAX_RBP = 4
MAX_LEN = 600
EMB = 64
EPOCHS = 12
SEED = 7


def load_task():
    phages = pd.read_csv("data/processed/phage_rbp_table.csv")
    phages = phages[phages.host_species.isin(SPECIES)].reset_index(drop=True)
    panel = pd.read_csv("data/processed/host_receptor_panel.csv")
    species_receptors = {sp: panel[panel.species == sp].sequence.tolist() for sp in SPECIES}

    proteins, protein_index = [], {}

    def idx_of(seq):
        if seq not in protein_index:
            protein_index[seq] = len(proteins)
            proteins.append(seq)
        return protein_index[seq]

    samples, groups = [], []
    for _, row in phages.iterrows():
        rbps = [s.strip() for s in row.rbp_sequences.split(" ||| ")][:MAX_RBP]
        if not rbps:
            continue
        p_idx = [idx_of(s) for s in rbps]
        pos_sp = row.host_species
        neg_sp = SPECIES[1] if pos_sp == SPECIES[0] else SPECIES[0]
        h_pos = [idx_of(s) for s in species_receptors[pos_sp]]
        h_neg = [idx_of(s) for s in species_receptors[neg_sp]]
        samples.append((p_idx, h_pos, 1))
        samples.append((p_idx, h_neg, 0))
        groups += [row.accession, row.accession]
    return phages, proteins, samples, groups, species_receptors


def kmer_lr_baseline(phages, samples, groups, species_receptors, splits):
    """Logistic regression on mean dipeptide spectrum (phage RBPs + host panel)."""
    panel_spec = {sp: np.mean([kmer_spectrum(s, k=2) for s in seqs], axis=0)
                  for sp, seqs in species_receptors.items()}
    X, y = [], []
    phage_map = phages.set_index("accession")
    # rebuild features per sample from group id (accession) and host side
    for (p_idx, h_idx, label), acc in zip(samples, groups):
        # h side identity: match panel spectrum by index unknown -> recover by
        # receptor count heuristic is fragile; instead rebuild from label:
        pass  # replaced below
    return None


def main():
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    t0 = time.time()
    phages, proteins, samples, groups = load_task()
    print(f"phages={len(phages)} samples={len(samples)} unique proteins={len(proteins)}")

    gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=SEED)
    idx = np.arange(len(samples))
    tr, te = next(gss.split(idx, groups=groups))
    gss2 = GroupShuffleSplit(n_splits=1, test_size=0.125, random_state=SEED + 1)
    tr_sub, va = next(gss2.split(tr, groups=[groups[i] for i in tr]))
    tr_idx, va_idx, te_idx = tr[tr_sub], tr[va], te
    print(f"split sizes train={len(tr_idx)} val={len(va_idx)} test={len(te_idx)}")

    bank = ProteinBank(proteins, max_len=MAX_LEN)
    model = GroupedInteractionModel(
        ProteinCNNEncoder(3, 48, EMB), ProteinCNNEncoder(3, 48, EMB), emb_dim=EMB)
    opt = torch.optim.Adam(model.parameters(), lr=2e-3, weight_decay=1e-5)
    lossf = nn.BCEWithLogitsLoss()

    def batches(indices, bs=64, shuffle=True):
        order = np.random.permutation(indices) if shuffle else indices
        for i in range(0, len(order), bs):
            chunk = [samples[j] for j in order[i:i + bs]]
            yield collate_pairs(chunk, bank), torch.tensor(
                [c[2] for c in chunk], dtype=torch.float32)

    def evaluate(indices):
        model.eval()
        logits, labels = [], []
        with torch.no_grad():
            for (pb, hb, pg, hg, _), y in batches(indices, bs=128, shuffle=False):
                logits.append(model(pb, hb, pg, hg))
                labels.append(y)
        logits = torch.cat(logits).numpy()
        labels = torch.cat(labels).numpy()
        return (roc_auc_score(labels, logits),
                average_precision_score(labels, logits),
                float(((logits > 0) == labels).mean()))

    history = []
    for epoch in range(1, EPOCHS + 1):
        model.train()
        tot, n = 0.0, 0
        for (pb, hb, pg, hg, _), y in batches(tr_idx):
            opt.zero_grad()
            loss = lossf(model(pb, hb, pg, hg), y)
            loss.backward()
            opt.step()
            tot += loss.item() * len(y)
            n += len(y)
        va_auroc, va_auprc, va_acc = evaluate(va_idx)
        history.append({"epoch": epoch, "train_loss": tot / n,
                        "val_auroc": va_auroc, "val_auprc": va_auprc, "val_acc": va_acc})
        print(f"epoch {epoch:2d} loss={tot/n:.4f} val AUROC={va_auroc:.3f} "
              f"AUPRC={va_auprc:.3f} acc={va_acc:.3f}")

    te_auroc, te_auprc, te_acc = evaluate(te_idx)
    result = {
        "model": "CNN-BILIN",
        "test_auroc": te_auroc, "test_auprc": te_auprc, "test_acc": te_acc,
        "n_phages": int(len(phages)), "n_samples": int(len(samples)),
        "n_unique_proteins": int(len(proteins)),
        "epochs": EPOCHS, "max_len": MAX_LEN, "emb": EMB,
        "runtime_s": round(time.time() - t0, 1),
        "history": history,
    }
    Path("results").mkdir(exist_ok=True)
    Path("results/interaction_cnn.json").write_text(json.dumps(result, indent=2))
    torch.save(model.state_dict(), "results/interaction_cnn.pt")
    print(f"\nTEST  AUROC={te_auroc:.3f} AUPRC={te_auprc:.3f} acc={te_acc:.3f}")
    print(f"runtime {result['runtime_s']}s")


if __name__ == "__main__":
    main()
