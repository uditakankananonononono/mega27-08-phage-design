from pathlib import Path

from phage_design.data.genbank_parse import parse_genbank

FIX = Path(__file__).parent / "fixtures" / "mini_phage.gb"


def test_parse_host_and_rbp():
    rec = parse_genbank(FIX.read_text())
    assert rec.host == "Escherichia coli"
    assert rec.genome_length == 50
    assert len(rec.rbp_products) == 1
    product, aa = rec.rbp_products[0]
    assert "tail fiber" in product.lower()
    assert aa.startswith("MKTAAV")


def test_hypothetical_protein_excluded():
    rec = parse_genbank(FIX.read_text())
    assert all("hypothetical" not in p.lower() for p, _ in rec.rbp_products)


def test_multi_record_parse():
    text = FIX.read_text()
    recs = list(__import__("phage_design.data.genbank_parse", fromlist=["parse_genbank_multi"]).parse_genbank_multi(text + "\n" + text))
    assert len(recs) == 2
    assert all(r.host == "Escherichia coli" for r in recs)
