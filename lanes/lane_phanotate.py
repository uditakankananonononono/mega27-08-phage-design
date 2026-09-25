"""PHANOTATE lane: phage-specific gene-calling cross-check (third opinion)."""
import csv, glob, json, random, subprocess, tempfile, os
from Bio import SeqIO

def main():
    random.seed(7)
    files = sorted(glob.glob("data/raw/genbank/*.gb*"))
    sample = random.sample(files, 4)
    rows = []
    for f in sample:
        rec = SeqIO.read(f, "genbank")
        n_gb = sum(1 for ft in rec.features if ft.type == "CDS")
        with tempfile.NamedTemporaryFile("w", suffix=".fna", delete=False) as t:
            SeqIO.write(rec, t.name, "fasta")
            fna = t.name
        out = fna + ".phanotate"
        r = subprocess.run(["python3", os.path.expanduser("~/.local/bin/phanotate.py"), fna],
                           capture_output=True, text=True, timeout=300)
        n_ph = r.stdout.count("CDS") or r.stdout.count(">")
        os.unlink(fna)
        rows.append({"accession": rec.id, "genbank_cds": n_gb, "phanotate_cds": n_ph,
                     "delta_pct": round(100 * (n_ph - n_gb) / max(n_gb, 1), 2)})
        print(rec.id, n_gb, n_ph, flush=True)
    with open("results/phanotate_annotation_qc.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    json.dump({"n_genomes": len(rows),
               "median_abs_delta_pct": sorted(abs(r["delta_pct"]) for r in rows)[len(rows)//2]},
              open("results/phanotate_annotation_qc.json", "w"), indent=2)

if __name__ == "__main__":
    main()
