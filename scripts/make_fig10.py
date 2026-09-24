"""Figure 10: per-seed validation-AUROC trajectories, from committed seed JSONs only."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

seeds = (7, 42, 123)
runs = {s: json.load(open(f"results/interaction_cnn_seed{s}.json")) for s in seeds}
fig, ax = plt.subplots(figsize=(7.2, 4.0))
for s in seeds:
    hist = runs[s]["history"]
    epochs = [h["epoch"] for h in hist]
    val = [h["val_auroc"] for h in hist]
    ax.plot(epochs, val, marker="o", ms=3, lw=1.2,
            label=f"seed {s} (test {runs[s]['test_auroc']:.3f})")
    ax.annotate(f"{val[0]:.3f}", (epochs[0], val[0]), textcoords="offset points",
                xytext=(0, 6), fontsize=7, ha="center")
base = runs[7]["baseline_kmer_lr"]["test_auroc"]
ax.axhline(base, color="grey", ls="--", lw=1.0)
ax.text(1.1, base - 0.006, f"k-mer LR baseline test AUROC {base:.3f}", fontsize=7, color="grey")
ax.set_xlabel("epoch")
ax.set_ylabel("validation AUROC")
ax.set_title("Validation-AUROC trajectories across seeds (identical splits)")
ax.set_xticks([1, 5, 10, 15, 20])
ax.legend(fontsize=8, loc="lower right")
ax.set_ylim(0.88, 1.0)
fig.tight_layout()
fig.savefig("paper/figures/fig10_seed_curves.png", dpi=200)
print("wrote paper/figures/fig10_seed_curves.png")
for s in seeds:
    hist = runs[s]["history"]
    print(s, "e1", round(hist[0]["val_auroc"],3), "e10", round(hist[9]["val_auroc"],3), "e20", round(hist[-1]["val_auroc"],3))
