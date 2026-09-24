"""Figures 5-8 for the paper, generated only from committed result files."""
import csv
import json
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------- fig5: seed replication + baseline comparison ----------
seeds = (7, 42, 123)
runs = {s: json.load(open(f"results/interaction_cnn_seed{s}.json")) for s in seeds}
metrics = [("test_auroc", "AUROC"), ("test_acc", "Accuracy"), ("test_choice_accuracy", "Choice acc.")]
fig, ax = plt.subplots(figsize=(7.2, 4.0))
x = np.arange(len(metrics))
w = 0.22
for i, s in enumerate(seeds):
    vals = [runs[s][m] for m, _ in metrics]
    ax.bar(x + (i - 1) * w, vals, w, label=f"seed {s}")
base = [runs[7]["baseline_kmer_lr"]["test_auroc"], np.nan, np.nan]
ax.bar(x + 2 * w, [b if not np.isnan(b) else 0 for b in base], w,
       label="k-mer LR baseline (AUROC)", color="grey", alpha=0.7)
for i, b in enumerate(base):
    if not np.isnan(b):
        ax.text(i + 2 * w, b + 0.005, f"{b:.3f}", ha="center", fontsize=8)
for i, s in enumerate(seeds):
    for j, (m, _) in enumerate(metrics):
        ax.text(j + (i - 1) * w, runs[s][m] + 0.005, f"{runs[s][m]:.3f}", ha="center", fontsize=7)
ax.set_xticks(x)
ax.set_xticklabels([lbl for _, lbl in metrics])
ax.set_ylim(0.85, 1.0)
ax.set_ylabel("held-out test value")
ax.set_title("Seed replication on the same held-out 158-phage test set")
ax.legend(fontsize=8, loc="lower right")
fig.tight_layout()
fig.savefig("paper/figures/fig5_seed_replication.png", dpi=200)
plt.close(fig)

# ---------- fig6: PHP forced-binary score margins ----------
php = json.load(open("results/php_headtohead.json"))
pp = php["per_phage"]
margins_e, margins_k = [], []
for r in pp:
    m = r["php_maxscore_ecoli"] - r["php_maxscore_kleb"]
    (margins_e if r["truth"] == "E. coli" else margins_k).append(m)
fig, ax = plt.subplots(figsize=(7.2, 4.0))
bins = np.linspace(min(margins_e + margins_k), max(margins_e + margins_k), 40)
ax.hist(margins_e, bins=bins, alpha=0.65, label=f"true E. coli phages (n={len(margins_e)})")
ax.hist(margins_k, bins=bins, alpha=0.65, label=f"true Klebsiella phages (n={len(margins_k)})")
ax.axvline(0, color="k", lw=1, ls="--")
ax.set_xlabel("PHP forced-binary margin: max GMM score (E. coli refs) - (K. pneumoniae refs)")
ax.set_ylabel("phages")
ax.set_title("PHP forced-binary score margins on the 158 held-out test phages")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("paper/figures/fig6_php_margins.png", dpi=200)
plt.close(fig)

# ---------- fig7: discovery screen top-receptor usage ----------
rows = list(csv.DictReader(open("results/discovery_screen.csv")))
cands = [r for r in rows if r["call"] == "candidate-host-switch"]
receptors = sorted({r["top_receptor_ecoli"].split(":")[0] for r in cands}
                   | {r["top_receptor_klebsiella"].split(":")[0] for r in cands})
M = np.full((len(cands), len(receptors)), np.nan)
for i, r in enumerate(cands):
    for col, key in ((None, "top_receptor_ecoli"), (None, "top_receptor_klebsiella")):
        name, score = r[key].split(":")
        j = receptors.index(name)
        val = float(score)
        if np.isnan(M[i, j]) or val > M[i, j]:
            M[i, j] = val
fig, ax = plt.subplots(figsize=(7.2, 5.2))
masked = np.ma.masked_invalid(M)
cmap = plt.cm.viridis.copy()
cmap.set_bad("lightgrey")
im = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=0.0, vmax=1.0)
ax.set_xticks(range(len(receptors)))
ax.set_xticklabels(receptors, rotation=45, ha="right", fontsize=8)
ax.set_yticks(range(len(cands)))
ax.set_yticklabels([r["accession"] for r in cands], fontsize=8)
for i in range(len(cands)):
    for j in range(len(receptors)):
        if not np.isnan(M[i, j]):
            ax.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=6)
fig.colorbar(im, ax=ax, label="max pair score (top receptor)")
ax.set_title("Top-scoring receptor for each cross-host candidate (either direction)")
fig.tight_layout()
fig.savefig("paper/figures/fig7_receptor_heatmap.png", dpi=200)
plt.close(fig)

# ---------- fig8: PHP native top-1 predicted species distribution ----------
pred_species = Counter(r["php_top1_species"] for r in pp)
top = pred_species.most_common(12)
fig, ax = plt.subplots(figsize=(7.2, 4.0))
names = [t[0] for t in top][::-1]
counts = [t[1] for t in top][::-1]
ax.barh(names, counts, color="steelblue")
for i, c in enumerate(counts):
    ax.text(c + 0.3, i, str(c), va="center", fontsize=8)
ax.set_xlabel("times predicted as top-1 host (out of 158 test phages)")
ax.set_title("PHP native top-1: most predicted species")
fig.tight_layout()
fig.savefig("paper/figures/fig8_php_species.png", dpi=200)
plt.close(fig)

print("candidates:", len(cands), "receptors:", receptors)
print("unique predicted species:", len(pred_species), "top:", top[:5])
print("wrote fig5-fig8")
