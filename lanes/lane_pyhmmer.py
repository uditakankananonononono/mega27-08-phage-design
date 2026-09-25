"""pyhmmer lane: HMM from the muscle alignment scanned over all corpus RBPs."""
import json
import pandas as pd, re
import pyhmmer
from pyhmmer.easel import Alphabet, MSAFile, TextSequence
from pyhmmer.plan7 import Builder, Background, HMMFile

def main():
    abc = Alphabet.amino()
    seqs = []
    import io
    with MSAFile("results/muscle_largest_cluster_aln.faa", format="afa",
                 alphabet=abc) as msaf:
        msa = msaf.read().digitize(abc)
    msa.name = b"fiber_cluster"
    builder = Builder(alphabet=abc)
    hmm, _, _ = builder.build_msa(msa, Background(abc))
    df = pd.read_csv("data/processed/phage_rbp_table.csv")
    targets = []
    for _, r in df.iterrows():
        if not isinstance(r.rbp_sequences, str):
            continue
        for i, s in enumerate(str(r.rbp_sequences).split(";")):
            s = s.strip()
            if len(s) >= 50:
                targets.append((f"{r.accession}_{i}".encode(), s))
    hits = []
    ts = [TextSequence(name=n, sequence=re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]','X',s)).digitize(abc) for n, s in targets]
    from pyhmmer.hmmer import hmmsearch
    for top in hmmsearch(hmm, ts):
        for hit in top:
            if hit.included:
                hits.append({"target": hit.name, "evalue": hit.evalue,
                             "score": hit.score})
    hits.sort(key=lambda h: h["evalue"])
    pd.DataFrame(hits).to_csv("results/pyhmmer_fiber_scan.csv", index=False)
    json.dump({"n_targets": len(targets), "n_included_hits": len(hits),
               "best_evalue": hits[0]["evalue"] if hits else None},
              open("results/pyhmmer_fiber_scan.json", "w"), indent=2)
    print(len(targets), len(hits))

if __name__ == "__main__":
    main()
