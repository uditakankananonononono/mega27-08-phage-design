"""sourmash lane: k-mer sketch containment between same/different host-genus phage genomes."""
import glob, json, random
import sourmash
from Bio import SeqIO

def main():
    random.seed(7)
    files = sorted(glob.glob("data/raw/genbank/*.gb"))
    recs = {}
    for f in files:
        rec = SeqIO.read(f, "genbank")
        recs.setdefault(rec.annotations.get("organism", "").split()[0], []).append(rec)
    picks = random.sample(recs.get("Klebsiella", []), 8) + random.sample(recs.get("Escherichia", []), 8)
    sk = []
    for r in picks:
        mh = sourmash.MinHash(n=0, ksize=31, scaled=100)
        mh.add_sequence(str(r.seq), force=True)
        sk.append((r.annotations.get("organism", "").split()[0], mh))
    same, diff = [], []
    for i in range(len(sk)):
        for j in range(i + 1, len(sk)):
            c = sk[i][1].avg_containment(sk[j][1])
            (same if sk[i][0] == sk[j][0] else diff).append(c)
    import statistics as st
    json.dump({"n_genomes": len(sk), "k": 31, "scaled": 100,
               "same_host_genus_mean_containment": round(st.mean(same), 4),
               "diff_host_genus_mean_containment": round(st.mean(diff), 4)},
              open("results/sourmash_genome_containment.json", "w"), indent=2)
    print(len(sk), round(st.mean(same), 4), round(st.mean(diff), 4))

if __name__ == "__main__":
    main()
