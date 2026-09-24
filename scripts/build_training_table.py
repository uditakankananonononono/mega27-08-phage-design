"""Rebuild a phage-level table with RBP AA sequences from cached GenBank files.

No network: reads data/raw/genbank/*.gb written by pull_phage_genomes.py.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from phage_design.data.genbank_parse import parse_genbank
from phage_design.data.hosts import normalize_host

rows = []
for gb in sorted(Path("data/raw/genbank").glob("*.gb")):
    try:
        rec = parse_genbank(gb.read_text(errors="replace"))
    except Exception:
        continue
    species, strain = normalize_host(rec.host)
    if species is None or not rec.rbp_products:
        continue
    rows.append({
        "accession": rec.accession,
        "organism": rec.organism,
        "host_species": species,
        "host_strain": strain or "",
        "genome_length": rec.genome_length,
        "n_rbps": len(rec.rbp_products),
        "rbp_sequences": " ||| ".join(aa for _, aa in rec.rbp_products),
        "rbp_products": " ||| ".join(p for p, _ in rec.rbp_products),
    })
df = pd.DataFrame(rows)
df.to_csv("data/processed/phage_rbp_table.csv", index=False)
print(f"phages with RBPs+mapped host: {len(df)}")
print(df.host_species.value_counts().to_string())
print("RBPs per phage:", df.n_rbps.describe()[["min", "50%", "max"]].to_dict())
