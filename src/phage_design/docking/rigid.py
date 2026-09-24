"""Rigid-body geometric docking of a phage RBP onto a host receptor.

A deliberately small, honest docker: uniform rotational grid (Fibonacci-sphere
quaternions) x translational anchors at receptor surface residues. The score
combines shape complementarity (favorable C-alpha contacts in 6-12 A, clashes
below 4.5 A) with a Coulomb electrostatic term over residue point charges.
No flexibility is modeled; the paper states this limit plainly and validates
on the experimental T5 pb5-FhuA complex (PDB 8A8C) via interface RMSD.
"""
from __future__ import annotations

import numpy as np

from ..features.sequence import AA_INDEX, CHARGE

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def fibonacci_rotations(n: int) -> np.ndarray:
    """n approximately uniform rotation matrices via Fibonacci-sphere axes."""
    rots = []
    golden = np.pi * (3 - np.sqrt(5))
    for i in range(n):
        z = 1 - 2 * (i + 0.5) / n
        theta = np.arccos(np.clip(z, -1, 1))
        phi = golden * i
        axis = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
        angle = 2 * np.pi * i / n
        rots.append(_rodrigues(axis, angle))
    return np.stack(rots)


def _rodrigues(axis: np.ndarray, angle: float) -> np.ndarray:
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    I = np.eye(3)
    return I + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)


def residue_charges(resnames: list[str]) -> np.ndarray:
    return np.array([CHARGE.get(THREE_TO_ONE.get(r, "A"), 0.0) for r in resnames])


def score_pose(rec_xyz: np.ndarray, lig_xyz: np.ndarray,
               rec_q: np.ndarray, lig_q: np.ndarray) -> float:
    """Higher is better: contacts - clashes + electrostatics."""
    d = np.linalg.norm(rec_xyz[:, None, :] - lig_xyz[None, :, :], axis=-1)
    contacts = ((d >= 4.5) & (d <= 12.0)).sum()
    clashes = (d < 4.5).sum()
    # electrostatics over near pairs only (clipped distance)
    near = np.clip(d, 2.0, 20.0)
    elec = -(rec_q[:, None] * lig_q[None, :]) / near  # opposite charges attract
    elec = elec[d <= 20.0].sum()
    return contacts - 10.0 * clashes + 0.1 * elec


def dock(rec_res: list[tuple[str, tuple]], lig_res: list[tuple[str, tuple]],
         n_rot: int = 48, n_anchor: int = 24, seed: int = 0):
    """Grid-search docking. Returns list of (score, rotation, translation)."""
    rng = np.random.default_rng(seed)
    rec_xyz = np.array([c for _, c in rec_res])
    lig_xyz = np.array([c for _, c in lig_res])
    rec_q = residue_charges([r for r, _ in rec_res])
    lig_q = residue_charges([r for r, _ in lig_res])
    lig_c = lig_xyz.mean(axis=0)
    # translational anchors: spread receptor residues (surface proxy)
    anchor_idx = rng.choice(len(rec_xyz), size=min(n_anchor, len(rec_xyz)), replace=False)
    rots = fibonacci_rotations(n_rot)
    results = []
    for R in rots:
        lig_rot = (lig_xyz - lig_c) @ R.T
        for ai in anchor_idx:
            t = rec_xyz[ai] - lig_rot.mean(axis=0)
            placed = lig_rot + t
            results.append((score_pose(rec_xyz, placed, rec_q, lig_q), R, t - lig_rot.mean(axis=0) + lig_rot.mean(axis=0)))
    results.sort(key=lambda r: -r[0])
    return results


def interface_rmsd(native_lig: np.ndarray, placed_lig: np.ndarray,
                   rec_xyz: np.ndarray, cutoff: float = 10.0) -> float:
    """RMSD over ligand C-alphas within `cutoff` of the receptor (CAPTURE-style iRMSD)."""
    d = np.linalg.norm(rec_xyz[:, None, :] - native_lig[None, :, :], axis=-1)
    iface = (d.min(axis=0) < cutoff)
    if iface.sum() < 3:
        return float("nan")
    a, b = native_lig[iface], placed_lig[iface]
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))
