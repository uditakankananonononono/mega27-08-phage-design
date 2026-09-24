"""Figure 11: cocktail member x targeted-receptor score heatmaps for
PD8-EC1 and PD8-KP1 (leg-A argmax design), from results/cocktail_design.json
and results/per_receptor_scores.csv. Discovery-candidate seats are starred."""
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

design = json.load(open("results/cocktail_design.json"))
scores = pd.read_csv("results/per_receptor_scores.csv").set_index("accession")

fig, axes = plt.subplots(2, 1, figsize=(9.5, 11),
                         gridspec_kw={"height_ratios": [18, 12]})
for ax, (sp, title) in zip(axes, [("ecoli", "PD8-EC1 (Escherichia coli panel)"),
                                  ("klebsiella", "PD8-KP1 (Klebsiella pneumoniae panel)")]):
    blk = design["species"][sp]
    a = blk["leg_a_argmax_design"]
    targeted = a["targeted_receptors"]
    prefix = "E-" if sp == "ecoli" else "K-"
    members = a["members"]
    M = np.array([[scores.loc[m["accession"], prefix + r] for r in targeted]
                  for m in members])
    im = ax.imshow(M, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    labels = [f"{m['member_id'].split('.')[-1]}. {m['phage_name']}"
              + (" *" if m["discovery_candidate"] else "") for m in members]
    ax.set_yticks(range(len(members)), labels, fontsize=7)
    ax.set_xticks(range(len(targeted)), targeted, fontsize=8, rotation=20,
                  ha="right")
    for i in range(len(members)):
        for j in range(len(targeted)):
            ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center",
                    fontsize=6, color="white" if M[i, j] < 0.6 else "black")
    ax.set_title(title, fontsize=10)
    # group separators every 3 members (k=3 per receptor)
    for g in range(3, len(members), 3):
        ax.axhline(g - 0.5, color="red", lw=0.8, alpha=0.6)
fig.suptitle("Cocktail member x targeted-receptor scores (argmax design; * = discovery-candidate seat)",
             fontsize=11)
fig.subplots_adjust(left=0.28, right=0.88, top=0.93, bottom=0.06, hspace=0.35)
cax = fig.add_axes([0.90, 0.35, 0.02, 0.30])
fig.colorbar(im, cax=cax, label="model receptor score")
fig.savefig("paper/figures/fig11_cocktail_heatmap.png", dpi=160)
print("wrote paper/figures/fig11_cocktail_heatmap.png")
