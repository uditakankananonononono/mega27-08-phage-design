"""UniProt REST access for host receptor sequences (live; fixtures in tests)."""
from __future__ import annotations

import time

import requests

UNIPROT = "https://rest.uniprot.org/uniprotkb"


def fetch_fasta_by_accession(acc: str) -> tuple[str, str] | None:
    """Return (header, sequence) for one accession, or None if absent."""
    time.sleep(0.3)
    r = requests.get(f"{UNIPROT}/{acc}.fasta", timeout=30)
    if r.status_code != 200 or not r.text.strip():
        return None
    lines = r.text.strip().splitlines()
    return lines[0], "".join(lines[1:])


def search_receptor(gene: str, organism_id: int, max_results: int = 3) -> list[tuple[str, str, str]]:
    """Search UniProtKB (reviewed first) for gene+organism; return (acc, desc, seq)."""
    query = f"(gene:{gene}) AND (organism_id:{organism_id})"
    params = {"query": query, "format": "fasta", "size": max_results}
    time.sleep(0.3)
    r = requests.get(f"{UNIPROT}/search", params=params, timeout=30)
    if r.status_code != 200 or not r.text.strip():
        return []
    out, header, seq = [], None, []
    for line in r.text.splitlines():
        if line.startswith(">"):
            if header:
                out.append((header.split("|")[1], header, "".join(seq)))
            header, seq = line, []
        else:
            seq.append(line.strip())
    if header:
        out.append((header.split("|")[1], header, "".join(seq)))
    return out
