import numpy as np

from phage_design.features.sequence import (clean_sequence, kmer_spectrum,
                                            physicochemical_profile)


def test_clean_sequence_drops_ambiguous():
    assert clean_sequence("ACDXB Z") == "ACD"


def test_kmer_spectrum_normalized_and_indexed():
    v = kmer_spectrum("AAAA", k=1)  # only alanine
    assert np.isclose(v.sum(), 1.0)
    assert v[0] == 1.0  # 'A' is index 0


def test_kmer_spectrum_dipeptide_positions():
    v = kmer_spectrum("AC", k=2)
    # 'AC' -> idx(A)*20 + idx(C) = 0*20 + 1
    assert v[1] == 1.0 and v.sum() == 1.0


def test_physicochemical_shape_and_values():
    m = physicochemical_profile("DEKR")
    assert m.shape == (4, 3)
    assert m[0, 2] == -1  # Asp charged -1
    assert m[3, 2] == 1   # Arg charged +1


def test_literature_values_locked():
    from phage_design.features.sequence import AA, CHARGE, HYDROPATHY, VOLUME
    assert AA == "ACDEFGHIKLMNPQRSTVWY"
    assert HYDROPATHY["A"] == 1.8 and HYDROPATHY["I"] == 4.5 and HYDROPATHY["R"] == -4.5
    assert HYDROPATHY["V"] == 4.2 and HYDROPATHY["W"] == -0.9
    assert CHARGE["D"] == -1 and CHARGE["K"] == 1 and CHARGE["A"] == 0
    assert VOLUME["G"] == 48 and VOLUME["W"] == 163
