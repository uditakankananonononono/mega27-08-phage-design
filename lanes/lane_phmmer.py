"""EBI HMMER phmmer lane: sequence homology search of two corpus RBPs vs PDB."""
import json, time, urllib.request, urllib.parse
import pandas as pd

UA = {"Content-Type": "application/json", "Accept": "application/json",
      "User-Agent": "mega27-research/1.0"}

def main():
    df = pd.read_csv("data/processed/phage_rbp_table.csv")
    seqs = [str(s).split(";")[0] for s in df.rbp_sequences.dropna().head(2)]
    outs = []
    for i, seq in enumerate(seqs):
        body = json.dumps({"database": "pdb", "seq": f">rbp{i}\n{seq}"}).encode()
        r = urllib.request.Request("https://www.ebi.ac.uk/Tools/hmmer/search/phmmer",
                                   data=body, headers=UA)
        try:
            with urllib.request.urlopen(r, timeout=60) as fh:
                res = json.load(fh)
        except Exception as e:
            outs.append({"query": f"rbp{i}", "error": str(e)})
            continue
        hits = res.get("results", {}).get("hits", [])[:10]
        outs.append({"query": f"rbp{i}", "n_hits": res.get("results", {}).get("nhits"),
                     "top": [{"target": h.get("acc"), "evalue": h.get("evalue"),
                              "name": h.get("name")} for h in hits]})
    json.dump({"note": "EBI HMMER phmmer vs PDB", "queries": outs},
              open("results/phmmer_rbp_pdb.json", "w"), indent=2)
    print(json.dumps(outs)[:300])

if __name__ == "__main__":
    main()
