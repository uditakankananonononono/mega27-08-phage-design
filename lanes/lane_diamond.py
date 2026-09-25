"""DIAMOND lane: RBP->receptor blastp as an mmseqs-independent sensitivity check."""
import json, re, subprocess, tempfile, os
import pandas as pd

DIAMOND = "/tmp/bin/diamond"

def main():
    panel = pd.read_csv("data/processed/host_receptor_panel.csv")
    pairs = pd.read_csv("data/processed/rbp_host_pairs.csv")
    df = pd.read_csv("data/processed/phage_rbp_table.csv")
    with tempfile.TemporaryDirectory() as td:
        dbfaa = os.path.join(td, "receptors.faa")
        with open(dbfaa, "w") as fh:
            for _, r in panel.iterrows():
                seq=re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]","X",str(r.sequence)); fh.write(f">{r.receptor}|{r.species}\n{seq}\n")
        subprocess.run([DIAMOND, "makedb", "--in", dbfaa, "-d", os.path.join(td, "rec")],
                       check=True, capture_output=True)
        qfaa = os.path.join(td, "q.faa")
        n_q = 0
        with open(qfaa, "w") as fh:
            for _, r in df.iterrows():
                if not isinstance(r.rbp_sequences, str):
                    continue
                for i, s in enumerate(str(r.rbp_sequences).split(";")):
                    s = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]","X",s.strip())
                    if len(s) >= 50:
                        fh.write(f">{r.accession}_{i}\n{s}\n"); n_q += 1
        out = os.path.join(td, "hits.tsv")
        subprocess.run([DIAMOND, "blastp", "-d", os.path.join(td, "rec"), "-q", qfaa,
                        "-o", out, "--max-target-seqs", "1", "--evalue", "1e-3",
                        "--threads", "2"], check=True, capture_output=True, timeout=900)
        hits = pd.read_csv(out, sep="\t", names=["q", "s", "pident", "alen", "mism",
                                             "gapopen", "qs", "qe", "ss", "se", "evalue", "bits"])
    hits["receptor"] = hits.s.str.split("|").str[0]
    hits.to_csv("results/diamond_rbp_receptor_hits.tsv", sep="\t", index=False)
    json.dump({"n_queries": n_q, "n_hits": int(len(hits)),
               "top_receptors": hits.receptor.value_counts().head(5).to_dict(),
               "median_pident": float(hits.pident.median()) if len(hits) else None},
              open("results/diamond_rbp_receptor.json", "w"), indent=2)
    print(n_q, len(hits))

if __name__ == "__main__":
    main()
