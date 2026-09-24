"""Live NCBI E-utilities access for phage genomes and host annotations.

All functions in this module hit the real NCBI endpoints. Tests never call
these; the test suite uses recorded fixtures instead (hermetic CI).
"""
from __future__ import annotations

import time
import urllib.parse
import xml.etree.ElementTree as ET

import requests

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def esearch(term: str, db: str = "nucleotide", retmax: int = 20) -> dict:
    """Search an NCBI database, returning count and id list."""
    url = f"{EUTILS}/esearch.fcgi"
    params = {"db": db, "term": term, "retmax": retmax, "retmode": "json"}
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    res = r.json()["esearchresult"]
    return {"count": int(res["count"]), "ids": res["idlist"]}


def esummary(ids: list[str], db: str = "nucleotide") -> list[dict]:
    """Fetch summary records for a list of ids."""
    if not ids:
        return []
    url = f"{EUTILS}/esummary.fcgi"
    params = {"db": db, "id": ",".join(ids), "retmode": "json"}
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    data = r.json()["result"]
    return [data[i] for i in data["uids"]]


def efetch_genbank(acc: str) -> str:
    """Fetch one GenBank record (flat text) by accession."""
    url = f"{EUTILS}/efetch.fcgi"
    params = {"db": "nucleotide", "id": acc, "rettype": "gb", "retmode": "text"}
    r = requests.get(url, params=params, timeout=60)
    r.raise_for_status()
    return r.text


def polite_pause(seconds: float = 0.34) -> None:
    """Respect NCBI's 3-requests/second limit without an API key."""
    time.sleep(seconds)
