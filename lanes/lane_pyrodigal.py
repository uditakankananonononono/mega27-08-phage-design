"""pyrodigal lane: independent gene-calling QC on corpus genomes."""
import csv, glob, json, random
from Bio import SeqIO
import pyrodigal

def main():
    random.seed(7)
    files = sorted(glob.glob("data/raw/genbank/*.gb*"))
    sample = random.sample(files, 8)
    rows = []
    for f in sample:
        rec = SeqIO.read(f, "genbank")
        n_gb = sum(1 for ft in rec.features if ft.type == "CDS")
        gf = pyrodigal.GeneFinder(meta=False)
        gf.train(bytes(rec.seq))
        genes = gf.find_genes(bytes(rec.seq))
        rows.append({"accession": rec.id, "genbank_cds": n_gb,
                     "pyrodigal_cds": len(genes),
                     "delta_pct": round(100 * (len(genes) - n_gb) / max(n_gb, 1), 2)})
        print(rec.id, n_gb, len(genes), flush=True)
    with open("results/pyrodigal_annotation_qc.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    json.dump({"n_genomes": len(rows),
               "median_abs_delta_pct": sorted(abs(r["delta_pct"]) for r in rows)[len(rows)//2],
               "note": "pyrodigal re-calls vs committed GenBank CDS counts"},
              open("results/pyrodigal_annotation_qc.json", "w"), indent=2)

if __name__ == "__main__":
    main()
