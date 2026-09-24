"""Genus-held-out split utilities for the cross-genus host-range benchmark.

The corpus carries exactly two host genera (Escherichia, Klebsiella). The
honest generalization test is train-on-one-genus / test-on-the-other, in
both directions. Group discipline is by phage accession: every sample
derived from one phage (positive and negative pair rows alike) stays on
one side, so no RBP sequence or host annotation can leak across the split.
"""


def genus_split(groups, phage_species, train_species):
    """Return (train_idx, test_idx) sample indices for one genus holdout.

    groups: per-sample phage accession (parallel to the samples list).
    phage_species: accession -> annotated host species (the positive host).
    train_species: the species whose phages form the training side.

    Every sample of a phage lands on exactly one side. The test side is the
    complement: phages annotated on any other species.
    """
    tr, te = [], []
    for i, acc in enumerate(groups):
        if acc not in phage_species:
            raise KeyError(f"accession {acc!r} missing from phage_species map")
        (tr if phage_species[acc] == train_species else te).append(i)
    return tr, te


def assert_no_leak(groups, tr_idx, te_idx):
    """Raise AssertionError if any accession appears on both sides."""
    tr_acc = {groups[i] for i in tr_idx}
    te_acc = {groups[i] for i in te_idx}
    overlap = tr_acc & te_acc
    assert not overlap, f"split leak: {len(overlap)} accessions on both sides"


def protein_novelty(samples, tr_idx, te_idx):
    """Fraction of (protein-set, side) memberships that are test-only.

    Reports, over unique protein indices: how many appear only in test
    samples (unseen during training). Characterizes the generalization
    regime: high novelty means the model must score proteins never
    encountered in training.
    """
    tr_prot, te_prot = set(), set()
    for i in tr_idx:
        tr_prot.update(samples[i][0])
        tr_prot.update(samples[i][1])
    for i in te_idx:
        te_prot.update(samples[i][0])
        te_prot.update(samples[i][1])
    novel = te_prot - tr_prot
    return {
        "n_test_proteins": len(te_prot),
        "n_novel_proteins": len(novel),
        "novel_fraction": (len(novel) / len(te_prot)) if te_prot else 0.0,
    }
