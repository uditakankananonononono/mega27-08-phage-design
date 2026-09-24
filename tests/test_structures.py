from phage_design.data.structures import ca_coords

MINI_PDB = """HEADER    TEST
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00 20.00           C
ATOM      3  CA  GLY A   2      12.500  10.000  10.000  1.00 20.00           C
ATOM      4  CA  LYS B   1      99.000  99.000  99.000  1.00 20.00           C
END
"""


def test_ca_coords_all_chains():
    res = ca_coords(MINI_PDB)
    assert len(res) == 3
    assert res[0][0] == "ALA"


def test_ca_coords_chain_filter():
    res = ca_coords(MINI_PDB, chain="A")
    assert len(res) == 2
    assert all(r[0] in ("ALA", "GLY") for r in res)
