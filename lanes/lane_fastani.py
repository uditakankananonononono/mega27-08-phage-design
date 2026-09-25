"""fastANI lane: ANI species-demarcation check vs ICTV host labels."""
import glob, json, random, subprocess, tempfile, os
from Bio import SeqIO

FASTANI = "/tmp/bin/fastANI"

def main():
    random.seed(13)
    files = sorted(glob.glob("data/raw/genbank/*.gb*"))
    kleb = random.sample([f for f in files if "Klebsiella" in SeqIO.read(f, "genbank").annotations.get("organism", "")], 15)
    esco = random.sample([f for f in files if "Escherichia" in SeqIO.read(f, "genbank").annotations.get("organism", "")], 15)
    with tempfile.TemporaryDirectory() as td:
        fnas, meta = [], {}
        for grp, fs in (("Klebsiella", kleb), ("Escherichia", esco)):
            for f in fs:
                rec = SeqIO.read(f, "genbank")
                out = os.path.join(td, rec.id + ".fna")
                SeqIO.write(rec, out, "fasta")
                fnas.append(out); meta[rec.id] = grp
        lst = os.path.join(td, "list.txt")
        open(lst, "w").write("\n".join(fnas))
        out = os.path.join(td, "ani.out")
        subprocess.run([FASTANI, "--ql", lst, "--rl", lst, "-o", out, "-t", "2"],
                       check=True, capture_output=True, timeout=900)
        rows = [l.split("\t") for l in open(out)]
    same, diff = [], []
    for a, b, ani, *_ in rows:
        ga, gb = meta[os.path.basename(a)[:-4]], meta[os.path.basename(b)[:-4]]
        (same if ga == gb else diff).append(float(ani))
    import statistics as st
    json.dump({"n_pairs": len(rows),
               "same_host_genus_mean_ani": round(st.mean(same), 3) if same else None,
               "diff_host_genus_mean_ani": round(st.mean(diff), 3) if diff else None,
               "n_same": len(same), "n_diff": len(diff)},
              open("results/fastani_species_check.json", "w"), indent=2)
    print(len(same), len(diff))

if __name__ == "__main__":
    main()
