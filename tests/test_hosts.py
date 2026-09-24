from phage_design.data.hosts import normalize_host


def test_species_only():
    assert normalize_host("Escherichia coli") == ("Escherichia coli", None)


def test_strain_extracted():
    assert normalize_host("Klebsiella pneumoniae INF359") == ("Klebsiella pneumoniae", "INF359")


def test_non_host_text_unmapped():
    assert normalize_host("hospital sewage") == (None, None)


def test_none_safe():
    assert normalize_host(None) == (None, None)
