import numpy as np

from phage_design.docking.rigid import (dock, fibonacci_rotations,
                                        interface_rmsd, residue_charges,
                                        score_pose)


def _helix(n, offset=(0, 0, 0)):
    # straight C-alpha line along x (3.8 A rise/residue) - no self-curvature,
    # so inter-chain distances are exactly the offset plus axial differences.
    res = []
    for i in range(n):
        res.append(("ALA", (i * 3.8 + offset[0], offset[1], offset[2])))
    return res


def test_fibonacci_rotations_are_rotations():
    rots = fibonacci_rotations(8)
    for R in rots:
        assert np.allclose(R @ R.T, np.eye(3), atol=1e-6)
        assert np.isclose(np.linalg.det(R), 1.0, atol=1e-6)


def test_residue_charges():
    q = residue_charges(["ASP", "LYS", "ALA"])
    assert q.tolist() == [-1, 1, 0]


def test_score_prefers_touching_nonclashing():
    rec = _helix(10)
    lig_close = _helix(6, offset=(0, 7, 0))    # 7 A away: contact zone, no clashes
    lig_far = _helix(6, offset=(0, 60, 0))     # far away
    rx = np.array([c for _, c in rec]); lx1 = np.array([c for _, c in lig_close])
    lx2 = np.array([c for _, c in lig_far])
    rq = residue_charges([r for r, _ in rec]); lq = residue_charges([r for r, _ in lig_close])
    assert score_pose(rx, lx1, rq, lq) > score_pose(rx, lx2, rq, lq)


def test_dock_returns_sorted_poses():
    rec, lig = _helix(12), _helix(6)
    res = dock(rec, lig, n_rot=4, n_anchor=4)
    assert len(res) == 16
    scores = [r[0] for r in res]
    assert scores == sorted(scores, reverse=True)


def test_interface_rmsd_zero_for_native():
    rec, lig = _helix(12), _helix(6, offset=(0, 8, 0))
    rx = np.array([c for _, c in rec]); lx = np.array([c for _, c in lig])
    assert interface_rmsd(lx, lx, rx) == 0.0
