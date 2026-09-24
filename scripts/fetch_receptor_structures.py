"""Download AlphaFold DB models for every receptor in the panel (live)."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from phage_design.data.structures import ca_coords, fetch_alphafold_pdb

OUT = Path("data/raw/structures")
OUT.mkdir(parents=True, exist_ok=True)

panel = pd.read_csv("data/processed/host_receptor_panel.csv")
ok, fail = 0, []
for _, row in panel.iterrows():
    acc = row.accession
    pdb = fetch_alphafold_pdb(acc)
    if pdb is None:
        fail.append(acc)
        print(f"  ! no AlphaFold model for {acc} ({row.receptor})")
        continue
    (OUT / f"{acc}.pdb").write_text(pdb)
    n_ca = len(ca_coords(pdb))
    print(f"  {row.species} {row.receptor} ({acc}): {n_ca} C-alpha atoms")
    ok += 1
print(f"structures fetched: {ok}, missing: {fail}")
