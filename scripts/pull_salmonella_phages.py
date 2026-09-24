"""Pull Salmonella phage genomes from NCBI for the 3-genus holdout benchmark.

Same machinery as pull_phage_genomes.py, separate cache dir so the two-species
model corpus is untouched. Writes data/processed/salmonella_rbp_table.csv
(phage-level, host-normalized, RBP sequences) for genus_holdout_v2.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from phage_design.data.ncbi import efetch_genbank_multi, esearch
from phage_design.data.genbank_parse import parse_genbank_multi
from phage_design.data.hosts import normalize_host
from pull_phage_genomes import _single_record_text

RAW = Path("data/raw/genbank_salmonella")
RAW.mkdir(parents=True, exist_ok=True)

QUERIES = {
    "salmonella": '("Salmonella phage"[Title] OR salmonella phage[Title]) AND complete genome[Title]',
}


def main(retmax: int = 120, retstart: int = 0) -> None:
    rows = []
    for label, term in QUERIES.items():
        res = esearch(term, retmax=retmax, retstart=retstart)
        ids = res["ids"]
        print(f"[{label}] {res['count']} hits; fetching {len(ids)}")
        for i in range(0, len(ids), 10):
            batch = ids[i:i + 10]
            try:
                text = efetch_genbank_multi(batch)
            except Exception as e:
                print(f"  ! efetch failed for batch {batch}: {e}")
                continue
            for rec in parse_genbank_multi(text):
                single = _single_record_text(text, rec.accession)
                if single:
                    (RAW / f"{rec.accession}.gb").write_text(single)
                species, strain = normalize_host(rec.host)
                for product, aa in rec.rbp_products:
                    rows.append({
                        "query": label, "accession": rec.accession,
                        "organism": rec.organism, "host": rec.host or "",
                        "host_species": species or "", "host_strain": strain or "",
                        "genome_length": rec.genome_length,
                        "rbp_product": product, "rbp_aa_len": len(aa),
                        "rbp_sequence": aa,
                    })
            print(f"  batch {i//10 + 1}: cumulative rows={len(rows)}")
    # collapse to phage level (like phage_rbp_table.csv)
    phages = {}
    for r in rows:
        if r["host_species"] != "Salmonella enterica":
            continue
        p = phages.setdefault(r["accession"], {
            "accession": r["accession"], "organism": r["organism"],
            "host_species": r["host_species"], "host_strain": r["host_strain"],
            "genome_length": r["genome_length"], "rbps": [], "prods": []})
        p["rbps"].append(r["rbp_sequence"])
        p["prods"].append(r["rbp_product"])
    out = Path("data/processed/salmonella_rbp_table.csv")
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "organism", "host_species",
                                           "host_strain", "genome_length",
                                           "n_rbps", "rbp_sequences", "rbp_products"])
        w.writeheader()
        for p in sorted(phages.values(), key=lambda x: x["accession"]):
            w.writerow({"accession": p["accession"], "organism": p["organism"],
                        "host_species": p["host_species"], "host_strain": p["host_strain"],
                        "genome_length": p["genome_length"], "n_rbps": len(p["rbps"]),
                        "rbp_sequences": " ||| ".join(p["rbps"]),
                        "rbp_products": " ||| ".join(p["prods"])})
    print(f"WROTE {out} phages={len(phages)} (of {len(rows)} RBP rows)")


if __name__ == "__main__":
    main(retmax=int(sys.argv[1]) if len(sys.argv) > 1 else 120,
         retstart=int(sys.argv[2]) if len(sys.argv) > 2 else 0)
