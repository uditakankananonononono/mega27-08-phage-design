"""Salmonella enterica receptor panel (live UniProt resolve, cached to CSV).

Genes resolved by search against S. enterica serovar Typhimurium LT2
(taxid 99287); accessions are resolved live, never memorized.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from phage_design.data.uniprot import search_receptor

# canonical Gram-negative phage receptor genes; reviewed LT2 entries preferred
GENES = ["lamB", "fhuA", "btuB", "ompA", "ompC", "ompF", "tsx", "tolC", "fepA", "bamA"]
TAXID = 99287  # Salmonella enterica serovar Typhimurium LT2


def main() -> None:
    rows = []
    for gene in GENES:
        hits = search_receptor(gene, TAXID, max_results=1)
        if not hits:
            print(f"  ! no UniProt hit for {gene}")
            continue
        acc, header, seq = hits[0]
        rows.append({"species": "Salmonella enterica", "receptor": gene,
                     "accession": acc, "length": len(seq), "sequence": seq})
        print(f"  {gene}: {acc} {len(seq)} aa")
    out = Path("data/processed/salmonella_receptor_panel.csv")
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["species", "receptor", "accession",
                                           "length", "sequence"])
        w.writeheader()
        w.writerows(rows)
    print(f"WROTE {out} receptors={len(rows)}")


if __name__ == "__main__":
    main()
