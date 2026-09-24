"""Positive control: re-dock T5 pb5 onto E. coli FhuA (PDB 8A8C) and measure
interface RMSD of the top-scoring pose against the experimental complex."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from phage_design.data.structures import ca_coords, fetch_pdb
from phage_design.docking.rigid import dock, interface_rmsd

pdb_text = fetch_pdb("8A8C")
assert pdb_text, "8A8C fetch failed"
Path("data/raw/structures/8A8C.pdb").write_text(pdb_text)
rec = ca_coords(pdb_text, chain="A")   # FhuA
lig = ca_coords(pdb_text, chain="B")   # pb5
print(f"FhuA: {len(rec)} CA, pb5: {len(lig)} CA")

native_lig = np.array([c for _, c in lig])
rec_xyz = np.array([c for _, c in rec])

results = dock(rec, lig, n_rot=48, n_anchor=48, seed=0)
best = results[0]
# reconstruct top pose
lig_c = native_lig.mean(axis=0)
placed = (native_lig - lig_c) @ best[1].T
# translation: recover from stored pose info is messy; recompute via centroid offset
# dock() stores (score, R, t) where placed = (lig - lig_c) @ R.T + t_eff
# simpler: recompute using the same procedure - but dock() didn't return placed coords.
# We instead re-run the top pose placement explicitly here using best rotation and
# the anchor that produced it. To keep the module honest and simple, we patch: dock
# returns (score, R, t) with placed = (lig - lig_c) @ R.T + t
placed = (native_lig - lig_c) @ best[1].T + best[2]
irmsd = interface_rmsd(native_lig, placed, rec_xyz)
print(f"top pose score={best[0]:.1f}  iRMSD vs experimental complex = {irmsd:.2f} A")
out = {"pdb": "8A8C", "n_rec": len(rec), "n_lig": len(lig),
       "top_score": float(best[0]), "irmsd_A": irmsd,
       "n_poses": len(results), "top5_scores": [float(r[0]) for r in results[:5]]}
Path("results/docking_positive_control.json").write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
