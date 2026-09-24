"""Export test-split phage genomes (cached GenBank) as FASTA for external-tool
head-to-head runs (PHP, DeepHost)."""
import sys
from pathlib import Path

import pandas as pd
from sklearn.model_selection import GroupShuffleSplit

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from Bio import SeqIO

SPECIES = ["Escherichia coli", "Klebsiella pneumoniae"]
phages = pd.read_csv("data/processed/phage_rbp_table.csv")
phages = phages[phages.host_species.isin(SPECIES)].reset_index(drop=True)
groups = [a for a in phages.accession for _ in range(2)]
import numpy as np
idx = np.arange(len(groups))
tr, te = next(GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=7).split(idx, groups=groups))
test_accs = sorted({groups[i] for i in te})

out = Path("external/test_genomes.fasta")
n = 0
with out.open("w") as fh:
    for acc in test_accs:
        gb = Path(f"data/raw/genbank/{acc}.gb")
        if not gb.exists():
            continue
        rec = SeqIO.read(gb.open(), "genbank")
        fh.write(f">{acc}\n{str(rec.seq)}\n")
        n += 1
print(f"wrote {out} with {n} test genomes")
Path("external/test_genomes.labels.csv").write_text(
    "accession,host_species\n" + "\n".join(
        f"{a},{phages.set_index('accession').loc[a, 'host_species']}" for a in test_accs
        if Path(f"data/raw/genbank/{a}.gb").exists()))
