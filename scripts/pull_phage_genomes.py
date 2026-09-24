"""Pull real phage genomes from NCBI and extract (phage, host, RBP) triples.

Live-NCBI step; results cached under data/raw/genbank so re-runs are cheap
and the test suite never touches the network.
"""
from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from phage_design.data.ncbi import efetch_genbank_multi, esearch
from phage_design.data.genbank_parse import parse_genbank_multi


def _single_record_text(text: str, accession: str) -> str:
    # split concatenated GenBank text and keep the record for this accession
    chunks = text.split("//\n")
    for c in chunks:
        if c.startswith("LOCUS") and accession in c.split("\n")[1][:200]:
            return c + "//\n"
    for c in chunks:
        if accession in c[:400]:
            return c + "//\n"
    return chunks[0] + "//\n"

RAW = Path("data/raw/genbank")
RAW.mkdir(parents=True, exist_ok=True)

QUERIES = {
    "klebsiella": 'klebsiella phage[Title] AND complete genome[Title]',
    "ecoli": '("Escherichia phage"[Title] OR coliphage[Title]) AND complete genome[Title]',
}


def main(retmax: int = 120) -> None:
    rows = []
    for label, term in QUERIES.items():
        res = esearch(term, retmax=retmax)
        ids = res["ids"]
        print(f"[{label}] {res['count']} hits; fetching {len(ids)}")
        for i in range(0, len(ids), 10):
            batch = ids[i:i + 10]
            try:
                text = efetch_genbank_multi(batch)
            except Exception as e:  # transient NCBI errors: skip, don't fabricate
                print(f"  ! efetch failed for batch {batch}: {e}")
                continue
            for rec in parse_genbank_multi(text):
                (RAW / f"{rec.accession}.gb").write_text(
                    text[text.index("LOCUS"):] if False else _single_record_text(text, rec.accession))
                for product, aa in rec.rbp_products:
                    rows.append({
                        "query": label, "accession": rec.accession,
                        "organism": rec.organism, "host": rec.host or "",
                        "genome_length": rec.genome_length,
                        "rbp_product": product, "rbp_aa_len": len(aa),
                    })
            print(f"  batch {i//10 + 1}: cumulative rows={len(rows)}")
    out = Path("data/processed/rbp_host_pairs.csv")
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else
                           ["query", "accession", "organism", "host", "genome_length",
                            "rbp_product", "rbp_aa_len"])
        w.writeheader()
        w.writerows(rows)
    print(f"WROTE {out} rows={len(rows)}")


if __name__ == "__main__":
    main(retmax=int(sys.argv[1]) if len(sys.argv) > 1 else 120)
