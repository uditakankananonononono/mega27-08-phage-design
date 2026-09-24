"""Live NCBI E-utilities access for phage genomes and host annotations.

All functions in this module hit the real NCBI endpoints. Tests never call
these; the test suite uses recorded fixtures instead (hermetic CI). Rate
limits are respected: <=3 requests/second without an API key, with
exponential backoff on HTTP 429/5xx.
"""
from __future__ import annotations

import time

import requests

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _get(url: str, params: dict, timeout: int = 60, max_retries: int = 6) -> requests.Response:
    delay = 1.0
    for attempt in range(max_retries):
        r = requests.get(url, params=params, timeout=timeout)
        if r.status_code in (429, 500, 502, 503):
            time.sleep(delay)
            delay = min(delay * 2, 60)
            continue
        r.raise_for_status()
        return r
    r.raise_for_status()  # final attempt's error
    raise RuntimeError("unreachable")


def esearch(term: str, db: str = "nucleotide", retmax: int = 20, retstart: int = 0) -> dict:
    """Search an NCBI database, returning count and id list."""
    params = {"db": db, "term": term, "retmax": retmax, "retstart": retstart, "retmode": "json"}
    time.sleep(0.4)
    res = _get(f"{EUTILS}/esearch.fcgi", params).json()["esearchresult"]
    return {"count": int(res["count"]), "ids": res["idlist"]}


def efetch_genbank_multi(ids: list[str]) -> str:
    """Fetch multiple GenBank records in one call (id list -> concatenated flat text)."""
    params = {"db": "nucleotide", "id": ",".join(ids), "rettype": "gb", "retmode": "text"}
    time.sleep(0.4)
    return _get(f"{EUTILS}/efetch.fcgi", params, timeout=120).text


def efetch_genbank(acc: str) -> str:
    """Fetch one GenBank record (flat text) by accession/uid."""
    return efetch_genbank_multi([acc])
