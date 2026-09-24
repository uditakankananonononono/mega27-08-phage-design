"""Held-out-genera host-range benchmark: train on one genus, test on the other.

Both directions (Escherichia -> Klebsiella, Klebsiella -> Escherichia),
group-split by accession, same task/labels/model/baseline as
train_interaction.py. Seeds 7/42/123. Results -> results/genus_holdout.json.

This is the cross-genus generalization regime: test phages' RBPs may never
appear in training, and the model must score them against receptor panels
that were the *negative* class for every training phage.
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
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupShuffleSplit

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from phage_design.dataset import ProteinBank, collate_pairs
from phage_design.features.sequence import kmer_spectrum
from phage_design.models.cnn import ProteinCNNEncoder
from phage_design.models.grouped import GroupedInteractionModel
from phage_design.genus import genus_split, assert_no_leak, protein_novelty
from train_interaction import load_task, SPECIES, MAX_LEN, EMB, EPOCHS

SEEDS = [7, 42, 123]


def kmer_lr(samples, host_names, rbp_lists, species_receptors, tr_idx, va_idx, te_idx):
    panel_spec = {sp: np.mean([kmer_spectrum(x, k=2) for x in seqs], axis=0)
                  for sp, seqs in species_receptors.items()}
    phage_spec = np.stack([np.mean([kmer_spectrum(x, k=2) for x in rbps], axis=0)
                           for rbps in rbp_lists])
    host_spec = np.stack([panel_spec[sp] for sp in host_names])
    X = np.concatenate([phage_spec * host_spec, phage_spec - host_spec], axis=1)
    y = np.array([s[2] for s in samples])
    scaler = StandardScaler().fit(X[tr_idx])
    Xs = scaler.transform(X)
    clf = LogisticRegression(max_iter=5000, C=1.0)
    clf.fit(Xs[tr_idx], y[tr_idx])
    out = {}
    for name, idx in [("val", va_idx), ("test", te_idx)]:
        score = clf.decision_function(Xs[idx])
        out[f"{name}_auroc"] = float(roc_auc_score(y[idx], score))
        out[f"{name}_auprc"] = float(average_precision_score(y[idx], score))
        out[f"{name}_acc"] = float(((score > 0) == y[idx]).mean())
    return out


def run_seed(seed, samples, groups, proteins, species_receptors,
             host_names, rbp_lists, train_species):
    torch.manual_seed(seed)
    np.random.seed(seed)
    tr_all, te_idx = genus_split(groups, phage_species_map, train_species)
    # carve a validation slice from the training genus only (group-safe)
    gss = GroupShuffleSplit(n_splits=1, test_size=0.15, random_state=seed)
    tr_all = np.array(tr_all)
    tr_sub, va_sub = next(gss.split(tr_all, groups=[groups[i] for i in tr_all]))
    tr_idx, va_idx = tr_all[tr_sub], tr_all[va_sub]
    te_idx = np.array(te_idx)
    assert_no_leak(groups, tr_idx, te_idx)
    assert_no_leak(groups, tr_idx, va_idx)

    bank = ProteinBank(proteins, max_len=MAX_LEN)
    panel_spec = {sp: np.mean([kmer_spectrum(x, k=2) for x in seqs], axis=0)
                  for sp, seqs in species_receptors.items()}
    phage_spec = np.stack([np.mean([kmer_spectrum(x, k=2) for x in rbps], axis=0)
                           for rbps in rbp_lists])
    host_spec = np.stack([panel_spec[sp] for sp in host_names])
    kmer_feats = np.concatenate([phage_spec * host_spec, phage_spec - host_spec], axis=1)
    mu, sd = kmer_feats[tr_idx].mean(0), kmer_feats[tr_idx].std(0) + 1e-8
    kmer_feats = torch.from_numpy(((kmer_feats - mu) / sd).astype(np.float32))

    model = GroupedInteractionModel(
        ProteinCNNEncoder(23, 48, EMB), ProteinCNNEncoder(23, 48, EMB), emb_dim=EMB)
    opt = torch.optim.Adam(model.parameters(), lr=2e-3, weight_decay=1e-5)
    lossf = nn.BCEWithLogitsLoss()

    def batches(indices, bs=64, shuffle=True):
        order = np.random.permutation(indices) if shuffle else indices
        for i in range(0, len(order), bs):
            sel = order[i:i + bs]
            chunk = [samples[j] for j in sel]
            yield (collate_pairs(chunk, bank), kmer_feats[sel],
                   torch.tensor([c[2] for c in chunk], dtype=torch.float32))

    def evaluate(indices):
        model.eval()
        logits, labels = [], []
        with torch.no_grad():
            for (pb, hb, pg, hg, _), kf, y in batches(indices, bs=128, shuffle=False):
                logits.append(model(pb, hb, pg, hg, kf))
                labels.append(y)
        logits = torch.cat(logits).numpy()
        labels = torch.cat(labels).numpy()
        return (float(roc_auc_score(labels, logits)),
                float(average_precision_score(labels, logits)),
                float(((logits > 0) == labels).mean()))

    for epoch in range(EPOCHS):
        model.train()
        for (pb, hb, pg, hg, _), kf, y in batches(tr_idx):
            opt.zero_grad()
            loss = lossf(model(pb, hb, pg, hg, kf), y)
            loss.backward()
            opt.step()
        if epoch in (EPOCHS - 1,):
            va = evaluate(va_idx)

    te_auroc, te_auprc, te_acc = evaluate(te_idx)

    # per-phage choice accuracy on the held-out genus
    model.eval()
    phage_samples = {}
    for j, ((_, _, label), acc) in enumerate(zip(samples, groups)):
        phage_samples.setdefault(acc, {})[label] = j
    te_set = set(te_idx.tolist())
    ok, tot = 0, 0
    with torch.no_grad():
        for acc, d in phage_samples.items():
            if 1 not in d or 0 not in d or d[1] not in te_set or d[0] not in te_set:
                continue
            s_pos = model(*collate_pairs([samples[d[1]]], bank)[:4],
                          kmer_feats[d[1]].unsqueeze(0)).item()
            s_neg = model(*collate_pairs([samples[d[0]]], bank)[:4],
                          kmer_feats[d[0]].unsqueeze(0)).item()
            ok += int(s_pos > s_neg)
            tot += 1

    base = kmer_lr(samples, host_names, rbp_lists, species_receptors,
                   tr_idx, va_idx, te_idx)
    return {
        "test_auroc": te_auroc, "test_auprc": te_auprc, "test_acc": te_acc,
        "val_auroc": va[0],
        "test_choice_accuracy": ok / max(tot, 1), "test_choice_n": tot,
        "n_train_phages": len({groups[i] for i in tr_idx}),
        "n_test_phages": len({groups[i] for i in te_idx}),
        "baseline_kmer_lr": {k: v for k, v in base.items() if k.startswith("test_")},
        "protein_novelty": protein_novelty(samples, tr_idx, te_idx),
    }


if __name__ == "__main__":
    t0 = time.time()
    phages, proteins, samples, groups, species_receptors, host_names, rbp_lists = load_task()
    phage_species_map = dict(zip(phages.accession, phages.host_species))
    out = {"config": {"seeds": SEEDS, "epochs": EPOCHS, "max_len": MAX_LEN,
                      "emb": EMB, "task": "binary host-range (annotated host vs other species)",
                      "split": "genus holdout, group by accession, both directions"},
           "directions": {}}
    for train_sp in SPECIES:
        test_sp = SPECIES[1] if train_sp == SPECIES[0] else SPECIES[0]
        name = ("E2K" if train_sp == SPECIES[0] else "K2E")
        per_seed = {}
        for seed in SEEDS:
            r = run_seed(seed, samples, groups, proteins, species_receptors,
                         host_names, rbp_lists, train_sp)
            per_seed[str(seed)] = r
            print(f"{name} seed={seed} AUROC={r['test_auroc']:.3f} "
                  f"AUPRC={r['test_auprc']:.3f} choice={r['test_choice_accuracy']:.3f} "
                  f"({r['test_choice_n']} phages) baseAUROC={r['baseline_kmer_lr']['test_auroc']:.3f} "
                  f"novel={r['protein_novelty']['novel_fraction']:.3f}")
        aucs = [r["test_auroc"] for r in per_seed.values()]
        b_aucs = [r["baseline_kmer_lr"]["test_auroc"] for r in per_seed.values()]
        choices = [r["test_choice_accuracy"] for r in per_seed.values()]
        out["directions"][name] = {
            "train_species": train_sp, "test_species": test_sp,
            "seeds": per_seed,
            "cnn_auroc_mean": float(np.mean(aucs)), "cnn_auroc_sd": float(np.std(aucs)),
            "kmerlr_auroc_mean": float(np.mean(b_aucs)),
            "choice_acc_mean": float(np.mean(choices)),
        }
    out["runtime_s"] = round(time.time() - t0, 1)
    Path("results/genus_holdout.json").write_text(json.dumps(out, indent=2))
    print(json.dumps({d: {k: round(v, 4) for k, v in out["directions"][d].items()
                          if k.endswith("mean")} for d in out["directions"]}, indent=2))
    print(f"runtime {out['runtime_s']}s")
