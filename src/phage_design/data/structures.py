"""Fetch real 3D structures: AlphaFold DB models and PDB entries.

AlphaFold DB serves predicted C-alpha+ models by UniProt accession;
PDB serves experimental structures by PDB id (mmCIF/PDB format). Only C-alpha
coordinates are extracted for the residue graph.
"""
from __future__ import annotations

import time

import requests

AFDB = "https://alphafold.ebi.ac.uk/files"
PDB = "https://files.rcsb.org/download"


def fetch_alphafold_pdb(uniprot_acc: str, version: int = 4) -> str | None:
    """AlphaFold DB model (PDB format) for a UniProt accession, else None."""
    time.sleep(0.3)
    url = f"{AFDB}/AF-{uniprot_acc}-F1-model_v{version}.pdb"
    r = requests.get(url, timeout=60)
    return r.text if r.status_code == 200 and r.text.startswith("HEADER") else None


def fetch_pdb(pdb_id: str) -> str | None:
    time.sleep(0.3)
    r = requests.get(f"{PDB}/{pdb_id}.pdb", timeout=60)
    return r.text if r.status_code == 200 else None


def ca_coords(pdb_text: str, chain: str | None = None) -> list[tuple[str, tuple[float, float, float]]]:
    """Extract (resname, (x,y,z)) for C-alpha atoms, optionally one chain."""
    out = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            if chain and line[21] != chain:
                continue
            out.append((line[17:20].strip(),
                        (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
    return out
