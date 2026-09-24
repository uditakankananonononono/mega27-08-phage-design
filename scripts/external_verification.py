"""Local verification analyses using installed packages.
1. ProtParam cross-check of the 14 panel receptors (Biopython ProtParam).
2. Pairwise identity between E./K. homolog pairs (grounds conservation claims).
3. Mann-Whitney tests for the per-receptor probe (statsmodels/scipy)."""
import json
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import pairwise2
from scipy.stats import mannwhitneyu

panel = pd.read_csv("data/processed/host_receptor_panel.csv")
out = {"protparam": {}, "homolog_identity": {}, "probe_mannwhitney": {}}
for _, r in panel.iterrows():
    pa = ProteinAnalysis(r.sequence)
    out["protparam"][r.receptor] = {
        "species": r.species, "length": len(r.sequence),
        "gravy": round(pa.gravy(), 3),
        "instability": round(pa.instability_index(), 1),
        "aromaticity": round(pa.aromaticity(), 3)}

def ident(a, b):
    aln = pairwise2.align.globalxx(a, b, one_alignment_only=True)[0]
    matches = sum(1 for x, y in zip(aln.seqA, aln.seqB) if x == y and x != "-")
    return round(matches / max(len(aln.seqA.rstrip("-")), 1), 3)

seq = {(r.species, r.receptor): r.sequence for _, r in panel.iterrows()}
E, K = "Escherichia coli", "Klebsiella pneumoniae"
out["homolog_identity"]["LamB_vs_lamB"] = ident(seq[(E, "LamB")], seq[(K, "lamB")])
out["homolog_identity"]["OmpA_vs_OmpA"] = ident(seq[(E, "OmpA")], seq[(K, "OmpA")])
out["homolog_identity"]["OmpC_vs_ompK36"] = ident(seq[(E, "OmpC")], seq[(K, "ompK36")])
out["homolog_identity"]["OmpC_vs_ompK35"] = ident(seq[(E, "OmpC")], seq[(K, "ompK35")])

scores = pd.read_csv("results/per_receptor_scores.csv")
for col in scores.columns:
    if col in ("accession", "annotated_host"):
        continue
    sp = E if col.startswith("E-") else K
    own = scores.loc[scores.annotated_host == sp, col]
    other = scores.loc[scores.annotated_host != sp, col]
    U, p = mannwhitneyu(own, other, alternative="two-sided")
    out["probe_mannwhitney"][col] = {"U": float(U), "p_two_sided": float(p)}

json.dump(out, open("results/external_verification.json", "w"), indent=1)
print("homolog identity:", out["homolog_identity"])
sig = {k: v["p_two_sided"] for k, v in out["probe_mannwhitney"].items()}
print("max p across receptors:", max(sig.values()), "| all significant:", all(p < 1e-6 for p in sig.values()))
