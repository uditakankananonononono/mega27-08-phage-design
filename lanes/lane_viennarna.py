"""ViennaRNA lane: 5'UTR/RBS accessibility of RBP transcripts (expression risk
screen for cocktail candidates - low MFE hairpins over the RBS occlude
translation initiation, a published design consideration)."""
import json
import pandas as pd
import RNA
from Bio import SeqIO
import glob

def main():
    df = pd.read_csv("data/processed/phage_rbp_table.csv")
    accs = set(df.accession.head(12))
    rows = []
    for f in sorted(glob.glob("data/raw/genbank/*.gb*")):
        rec = SeqIO.read(f, "genbank")
        base = rec.id.split(".")[0]
        if base not in accs:
            continue
        for ft in rec.features:
            if ft.type != "CDS":
                continue
            start = max(int(ft.location.start) - 60, 0)
            window = str(rec.seq[start:int(ft.location.start) + 24])
            if len(window) < 40:
                continue
            ss, mfe = RNA.fold(window)
            rows.append({"accession": base, "gene": ft.qualifiers.get("gene", ["-"])[0],
                         "utr_mfe_kcal": round(mfe, 2),
                         "hairpin_over_rbs": ss[:60].count("(") > 8})
            break
    import csv
    with open("results/viennarna_rbs_accessibility.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    json.dump({"n_transcripts": len(rows),
               "median_utr_mfe": sorted(r["utr_mfe_kcal"] for r in rows)[len(rows)//2] if rows else None,
               "n_occluded_rbs": sum(r["hairpin_over_rbs"] for r in rows)},
              open("results/viennarna_rbs.json", "w"), indent=2)
    print(len(rows))

if __name__ == "__main__":
    main()
