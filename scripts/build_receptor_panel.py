"""Build the host receptor panel: real receptor sequences for E. coli K-12/MG1655
and Klebsiella pneumoniae, from UniProtKB (live fetch, cached to CSV)."""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from phage_design.data.uniprot import fetch_fasta_by_accession, search_receptor

# E. coli K-12 (taxid 83333) canonical phage receptors (UniProt accessions).
ECOLI_PANEL = {
    "LamB": "P02943",  # maltoporin - lambda phage receptor
    "FhuA": "P06971",  # ferrichrome porin - T1/T5/phi80 receptor
    "OmpC": "P06996",  # outer membrane porin C - T4 receptor
    "OmpA": "P0A910",  # outer membrane protein A - K3/OX2-class
    "OmpF": "P02931",  # porin F - T2-class
    "Tsx":  "P0A927",  # nucleoside channel - T6 receptor
    "BtuB": "P06129",  # vitamin B12 transporter - BF23 receptor
    "FepA": "P05825",  # ferric enterobactin receptor
    "TolC": "P02930",  # outer membrane channel - TLS-class
    "YaeT_BamA": "P0A940",  # BamA beta-barrel assembly protein
}
# K. pneumoniae receptors resolved by gene search (accessions vary by strain).
KLEB_GENES = ["ompK36", "ompK35", "lamB", "fhuA"]
KLEB_DIRECT = {"OmpA": "P24017"}  # reviewed full-length K. pneumoniae OmpA
KLEB_TAXID = 573  # Klebsiella pneumoniae (species)


def main() -> None:
    rows = []
    for name, acc in ECOLI_PANEL.items():
        res = fetch_fasta_by_accession(acc)
        if not res:
            print(f"  ! missing {name} ({acc})")
            continue
        header, seq = res
        rows.append({"species": "Escherichia coli", "receptor": name,
                     "accession": acc, "length": len(seq), "sequence": seq})
        print(f"  E. coli {name}: {len(seq)} aa")
    for gene in KLEB_GENES:
        hits = search_receptor(gene, KLEB_TAXID, max_results=1)
        for acc, header, seq in hits:
            rows.append({"species": "Klebsiella pneumoniae", "receptor": gene,
                         "accession": acc, "length": len(seq), "sequence": seq})
            print(f"  K. pneumoniae {gene}: {len(seq)} aa ({acc})")
        if not hits:
            print(f"  ! no UniProt hit for {gene} in taxid {KLEB_TAXID}")
    for name, acc in KLEB_DIRECT.items():
        res = fetch_fasta_by_accession(acc)
        if res:
            header, seq = res
            rows.append({"species": "Klebsiella pneumoniae", "receptor": name,
                         "accession": acc, "length": len(seq), "sequence": seq})
            print(f"  K. pneumoniae {name}: {len(seq)} aa ({acc})")
    # prefer the longest hit per gene for K. pneumoniae (avoid fragments)
    best = {}
    for r in rows:
        if r["species"] == "Klebsiella pneumoniae":
            k = r["receptor"].lower()
            if k not in best or r["length"] > best[k]["length"]:
                best[k] = r
    rows = [r for r in rows if r["species"] != "Klebsiella pneumoniae"] + list(best.values())
    out = Path("data/processed/host_receptor_panel.csv")
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["species", "receptor", "accession", "length", "sequence"])
        w.writeheader()
        w.writerows(rows)
    print(f"WROTE {out} rows={len(rows)}")


if __name__ == "__main__":
    main()
