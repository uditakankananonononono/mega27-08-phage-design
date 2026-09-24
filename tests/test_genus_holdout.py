"""Hermetic tests for the genus-held-out split (no external data)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from phage_design.genus import genus_split, assert_no_leak, protein_novelty

# 3 Escherichia phages (EC1-3), 2 Klebsiella phages (KP1-2).
# Each phage contributes a positive and a negative sample row.
GROUPS = ["EC1", "EC1", "EC2", "EC2", "EC3", "EC3", "KP1", "KP1", "KP2", "KP2"]
SPECIES = {"EC1": "Escherichia coli", "EC2": "Escherichia coli",
           "EC3": "Escherichia coli", "KP1": "Klebsiella pneumoniae",
           "KP2": "Klebsiella pneumoniae"}
# samples: (phage protein idx tuple, host protein idx tuple, label)
SAMPLES = [((0,), (10,), 1), ((0,), (11,), 0),
           ((1,), (10,), 1), ((1,), (11,), 0),
           ((2,), (10,), 1), ((2,), (11,), 0),
           ((3,), (11,), 1), ((3,), (10,), 0),
           ((4,), (11,), 1), ((4,), (10,), 0)]


def test_genus_split_both_directions_partition():
    for train_sp in ("Escherichia coli", "Klebsiella pneumoniae"):
        tr, te = genus_split(GROUPS, SPECIES, train_sp)
        assert sorted(tr + te) == list(range(len(GROUPS)))
        assert_no_leak(GROUPS, tr, te)
        tr_acc = {GROUPS[i] for i in tr}
        expected = {"EC1", "EC2", "EC3"} if train_sp == "Escherichia coli" else {"KP1", "KP2"}
        assert tr_acc == expected


def test_genus_split_missing_accession_raises():
    with pytest.raises(KeyError):
        genus_split(["EC1", "XX"], SPECIES, "Escherichia coli")


def test_assert_no_leak_detects_overlap():
    bad_groups = ["EC1", "EC1", "KP1"]
    with pytest.raises(AssertionError):
        assert_no_leak(bad_groups, [0], [1, 2])


def test_protein_novelty_counts_test_only_proteins():
    tr, te = genus_split(GROUPS, SPECIES, "Escherichia coli")
    nov = protein_novelty(SAMPLES, tr, te)
    # train side proteins: phage {0,1,2}, host {10,11}; test: phage {3,4}, host {10,11}
    assert nov["n_test_proteins"] == 4  # {3,4,10,11}
    assert nov["n_novel_proteins"] == 2  # {3,4}
    assert nov["novel_fraction"] == pytest.approx(0.5)


def test_protein_novelty_empty_test_side():
    assert protein_novelty(SAMPLES, list(range(len(SAMPLES))), [])["novel_fraction"] == 0.0
