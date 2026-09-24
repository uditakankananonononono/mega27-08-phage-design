"""3-genus holdout benchmark: train on Escherichia+Klebsiella phages,
test on Salmonella enterica phages - a true held-out genus.

Fixes the v1 degeneracy: with two genera in training, every receptor panel
appears as BOTH a positive and a negative, so no panel-identity prior can
form; the model must use RBP-host compatibility. Per Salmonella phage: 3
rows (Sal panel = positive; Ec and Kp panels = negatives). Metrics: AUROC/
AUPRC/acc over held-out rows, 3-way choice accuracy (Sal outscores both),
KMER-LR baseline, seeds 7/42/123. Results -> results/genus_holdout_v2.json.
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
from phage_design.predict import MAX_RBP
from train_interaction import load_task, SPECIES, MAX_LEN, EMB, EPOCHS

SEEDS = [7, 42, 123]
SAL_SP = "Salmonella enterica"


def load_salmonella():
    df = pd.read_csv("data/processed/salmonella_rbp_table.csv")
    panel = pd.read_csv("data/processed/salmonella_receptor_panel.csv")
    return df, panel.sequence.tolist()


def main():
    t0 = time.time()
    phages, proteins, samples, groups, species_receptors, host_names, rbp_lists = load_task()
    protein_index = {s: i for i, s in enumerate(proteins)}

    def idx_of(seq):
        if seq not in protein_index:
            protein_index[seq] = len(proteins)
            proteins.append(seq)
        return protein_index[seq]

    sal_df, sal_panel_seqs = load_salmonella()
    sal_h = [idx_of(s) for s in sal_panel_seqs]
    ec_h = [idx_of(s) for s in species_receptors[SPECIES[0]]]
    kp_h = [idx_of(s) for s in species_receptors[SPECIES[1]]]
    # test rows: per Salmonella phage, (Sal=1, Ec=0, Kp=0)
    sal_samples, sal_groups, sal_host_names, sal_rbp_lists = [], [], [], []
    for _, row in sal_df.iterrows():
        rbps = [s.strip() for s in row.rbp_sequences.split(" ||| ")][:MAX_RBP]
        if not rbps:
            continue
        p_idx = [idx_of(s) for s in rbps]
        for h_idx, hname, label in [(sal_h, SAL_SP, 1), (ec_h, SPECIES[0], 0),
                                    (kp_h, SPECIES[1], 0)]:
            sal_samples.append((p_idx, h_idx, label))
            sal_groups.append(row.accession)
            sal_host_names.append(hname)
            sal_rbp_lists.append(rbps)

    all_samples = samples + sal_samples
    all_host_names = host_names + sal_host_names
    all_rbp_lists = rbp_lists + sal_rbp_lists
    n_train_pool = len(samples)
    panel_spec = {sp: np.mean([kmer_spectrum(x, k=2) for x in seqs], axis=0)
                  for sp, seqs in {**species_receptors, SAL_SP: sal_panel_seqs}.items()}
    phage_spec = np.stack([np.mean([kmer_spectrum(x, k=2) for x in rbps], axis=0)
                           for rbps in all_rbp_lists])
    host_spec = np.stack([panel_spec[sp] for sp in all_host_names])
    kmer_all = np.concatenate([phage_spec * host_spec, phage_spec - host_spec], axis=1)
    y_all = np.array([s[2] for s in all_samples])
    te_idx = np.arange(n_train_pool, len(all_samples))
    print(f"train pool={n_train_pool} test rows={len(te_idx)} "
          f"({len(sal_df)} Salmonella phages) unique proteins={len(proteins)}")

    results = {}
    for seed in SEEDS:
        torch.manual_seed(seed)
        np.random.seed(seed)
        gss = GroupShuffleSplit(n_splits=1, test_size=0.15, random_state=seed)
        pool = np.arange(n_train_pool)
        tr_sub, va_sub = next(gss.split(pool, groups=groups))
        tr_idx, va_idx = pool[tr_sub], pool[va_sub]

        bank = ProteinBank(proteins, max_len=MAX_LEN)
        mu, sd = kmer_all[tr_idx].mean(0), kmer_all[tr_idx].std(0) + 1e-8
        kmer_feats = torch.from_numpy(((kmer_all - mu) / sd).astype(np.float32))
        model = GroupedInteractionModel(
            ProteinCNNEncoder(23, 48, EMB), ProteinCNNEncoder(23, 48, EMB), emb_dim=EMB)
        opt = torch.optim.Adam(model.parameters(), lr=2e-3, weight_decay=1e-5)
        lossf = nn.BCEWithLogitsLoss()

        def batches(indices, bs=64, shuffle=True):
            order = np.random.permutation(indices) if shuffle else indices
            for i in range(0, len(order), bs):
                sel = order[i:i + bs]
                chunk = [all_samples[j] for j in sel]
                yield (collate_pairs(chunk, bank), kmer_feats[sel], \
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

        te = evaluate(te_idx)
        # 3-way choice per Salmonella phage
        model.eval()
        ok, tot = 0, 0
        with torch.no_grad():
            for i in range(0, len(te_idx), 3):
                trip = te_idx[i:i + 3]
                logits = []
                for j in trip:
                    pb, hb, pg, hg, _ = collate_pairs([all_samples[j]], bank)
                    logits.append(model(pb, hb, pg, hg, kmer_feats[j].unsqueeze(0)).item())
                ok += int(np.argmax(logits) == 0)  # Sal row is first
                tot += 1

        # baseline
        scaler = StandardScaler().fit(kmer_all[tr_idx])
        Xs = scaler.transform(kmer_all)
        clf = LogisticRegression(max_iter=5000, C=1.0).fit(Xs[tr_idx], y_all[tr_idx])
        b_score = clf.decision_function(Xs[te_idx])
        b = {"test_auroc": float(roc_auc_score(y_all[te_idx], b_score)),
             "test_auprc": float(average_precision_score(y_all[te_idx], b_score))}

        results[str(seed)] = {
            "test_auroc": te[0], "test_auprc": te[1], "test_acc": te[2],
            "choice3_accuracy": ok / tot, "choice3_n": tot,
            "baseline_kmer_lr": b,
        }
        print(f"seed={seed} AUROC={te[0]:.3f} AUPRC={te[1]:.3f} "
              f"choice3={ok/tot:.3f} ({tot}) baseAUROC={b['test_auroc']:.3f}")

    aucs = [r["test_auroc"] for r in results.values()]
    ch = [r["choice3_accuracy"] for r in results.values()]
    b_aucs = [r["baseline_kmer_lr"]["test_auroc"] for r in results.values()]
    out = {"config": {"seeds": SEEDS, "epochs": EPOCHS, "max_len": MAX_LEN, "emb": EMB,
                      "task": "3-way host panel scoring; Sal=positive, Ec/Kp=negatives",
                      "split": "train Ec+Kp phages (group by accession); test = all Salmonella phages",
                      "n_salmonella_phages": int(len(sal_df)),
                      "n_test_rows": int(len(te_idx))},
           "seeds": results,
           "cnn_auroc_mean": float(np.mean(aucs)), "cnn_auroc_sd": float(np.std(aucs)),
           "choice3_acc_mean": float(np.mean(ch)),
           "kmerlr_auroc_mean": float(np.mean(b_aucs)),
           "runtime_s": round(time.time() - t0, 1)}
    Path("results/genus_holdout_v2.json").write_text(json.dumps(out, indent=2))
    print(json.dumps({k: round(v, 4) for k, v in out.items()
                      if k.endswith(("mean", "sd"))}, indent=2))
    print(f"runtime {out['runtime_s']}s")


if __name__ == "__main__":
    main()
