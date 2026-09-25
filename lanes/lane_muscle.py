"""muscle 5 lane: align the largest MMseqs2 RBP cluster; consensus + identity."""
import json, subprocess, tempfile, os
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
import numpy as np

MUSCLE = "/tmp/bin/muscle"

def main():
    recs = list(SeqIO.parse("results/mmseqs_largest_cluster.faa", "fasta"))
    with tempfile.TemporaryDirectory() as td:
        aln = os.path.join(td, "aln.faa")
        subprocess.run([MUSCLE, "-align", "/tmp/sub.faa",
                        "-output", aln], check=True, capture_output=True, timeout=600)
        arecs = list(SeqIO.parse(aln, "fasta"))
    L = len(arecs[0].seq)
    cons, ident_cols = [], 0
    for i in range(L):
        col = [str(r.seq)[i] for r in arecs]
        aa = max(set(col) - {"-"}, key=col.count, default="-")
        cons.append(aa)
        if len(set(col)) == 1 and col[0] != "-":
            ident_cols += 1
    consensus = "".join(cons)
    with open("results/muscle_largest_cluster_aln.faa", "w") as fh:
        for r in arecs:
            fh.write(f">{r.id}\n{r.seq}\n")
    json.dump({"n_sequences": len(arecs), "alignment_length": L,
               "fully_conserved_columns": ident_cols,
               "frac_conserved": round(ident_cols / L, 4),
               "consensus_first80": consensus[:80]},
              open("results/muscle_fiber_alignment.json", "w"), indent=2)
    print(len(arecs), L, ident_cols)

if __name__ == "__main__":
    main()
