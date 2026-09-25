"""MMseqs2 lane: cluster all corpus RBPs; compare clusters to receptor use."""
import json, subprocess, tempfile, os
import pandas as pd

MMSEQS = "/tmp/bin/mmseqs/bin/mmseqs"

def main():
    df = pd.read_csv("data/processed/phage_rbp_table.csv")
    rows = []
    for _, r in df.iterrows():
        if not isinstance(r.rbp_sequences, str):
            continue
        for i, s in enumerate(str(r.rbp_sequences).split(";")):
            s = s.strip()
            if len(s) >= 50:
                rows.append((f"{r.accession}_{i}", s, r.host_species))
    with tempfile.TemporaryDirectory() as td:
        faa = os.path.join(td, "rbps.faa")
        with open(faa, "w") as fh:
            for h, s, _ in rows:
                fh.write(f">{h}\n{s}\n")
        subprocess.run([MMSEQS, "easy-cluster", faa, os.path.join(td, "clu"),
                        os.path.join(td, "tmp"), "--min-seq-id", "0.4",
                        "-c", "0.5", "--cov-mode", "1"],
                       check=True, capture_output=True, timeout=900)
        clu = pd.read_csv(os.path.join(td, "clu_cluster.tsv"), sep="\t",
                          names=["rep", "member"])
    clu["cluster"] = clu.groupby("rep").ngroup()
    sizes = clu.groupby("rep").size().sort_values(ascending=False)
    host_map = dict(zip([h for h, _, _ in rows], [hs for _, _, hs in rows]))
    clu["host"] = clu.member.map(host_map)
    top = clu[clu.rep == sizes.index[0]]
    json.dump({"n_rbps": len(rows), "n_clusters": int(clu.cluster.nunique()),
               "largest_cluster_size": int(sizes.iloc[0]),
               "largest_cluster_hosts": sorted(set(top.host.dropna())),
               "n_singletons": int((sizes == 1).sum())},
              open("results/mmseqs_rbp_clusters.json", "w"), indent=2)
    clu.to_csv("results/mmseqs_rbp_clusters.tsv", sep="\t", index=False)
    top_members = list(top.member)
    with open("results/mmseqs_largest_cluster.faa", "w") as fh:
        seq = dict((h, s) for h, s, _ in rows)
        for m in top_members:
            fh.write(f">{m}\n{seq[m]}\n")
    print(len(rows), clu.cluster.nunique(), sizes.iloc[0])

if __name__ == "__main__":
    main()
