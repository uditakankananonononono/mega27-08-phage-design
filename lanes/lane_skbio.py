"""scikit-bio lane: RBP k-mer composition PERMANOVA vs host species."""
import csv, json
import numpy as np
import pandas as pd
from skbio import DistanceMatrix
from skbio.stats.distance import permanova
from scipy.spatial.distance import pdist, squareform

def kmer_freq(seq, k=3):
    from collections import Counter
    c = Counter(seq[i:i+k] for i in range(0, len(seq) - k + 1))
    tot = sum(c.values()) or 1
    return {m: v / tot for m, v in c.items()}

def main():
    df = pd.read_csv("data/processed/phage_rbp_table.csv")
    df = df[df.rbp_sequences.notna()].head(400)
    feats = [kmer_freq(s) for s in df.rbp_sequences]
    keys = sorted({k for f in feats for k in f})
    X = np.array([[f.get(k, 0.0) for k in keys] for f in feats])
    D = squareform(pdist(X, metric="braycurtis"))
    labels = df.host_species.fillna("unknown").values
    dm = DistanceMatrix(D, ids=[str(i) for i in range(len(D))])
    res = permanova(dm, labels, permutations=199)
    json.dump({"n_rbps": int(len(df)), "test": "PERMANOVA braycurtis k=3",
               "pseudo_F": float(res["test statistic"]), "p_value": float(res["p-value"]),
               "groups": sorted(set(labels))},
              open("results/skbio_rbp_permanova.json", "w"), indent=2)
    print(res["test statistic"], res["p-value"])

if __name__ == "__main__":
    main()
