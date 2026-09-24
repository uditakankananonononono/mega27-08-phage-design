"""Phage-cocktail design as k-redundant receptor multicover with exact
host-escape analysis.

Model. A phage p is predicted to use receptor set R_p (score >= theta).
Phage p infects a host iff at least one receptor in R_p is intact. The host
escapes phage p only by inactivating EVERY receptor in R_p, and escapes the
cocktail only by inactivating the union U = union of R_p over members
(Theorem: escape set K must satisfy K superset of U; min |K| = |U|).

Redundancy (k-cover depth) buys robustness to model error: if each positive
usage claim is independently wrong with probability eps, receptor r belongs
to the TRUE union with probability 1 - eps^{d_r}, so the expected true escape
depth is sum_r (1 - eps^{d_r}) -- increasing depth d_r raises it.

Escape dynamics. With per-receptor loss-of-function rate mu per generation,
a single lineage acquires all of U in G generations with probability
q ~ (G*mu)^|U| as mu -> 0; a population of N lineages escapes with
P_esc = 1 - (1 - q)^N <= N*q. We compute q exactly by DP over subsets.

All functions are pure and hermetic.
"""
from __future__ import annotations

import numpy as np


def usage_matrix(scores: np.ndarray, theta: float) -> np.ndarray:
    """Boolean phage x receptor usage matrix at score threshold theta."""
    scores = np.asarray(scores, dtype=float)
    if scores.ndim != 2:
        raise ValueError("scores must be 2-D (phage x receptor)")
    return scores >= theta


def coverage_depth(U: np.ndarray) -> np.ndarray:
    """Per-receptor number of phages using it (sum over rows)."""
    return np.asarray(U, dtype=int).sum(axis=0)


def greedy_k_cover(U: np.ndarray, k: int) -> list[int]:
    """Greedy minimum multicover: each receptor covered by >= k phages.

    f(S) = sum_r min(k, depth_S(r)) is monotone submodular; the greedy rule
    (add the row maximizing the marginal f-gain, tie-break by larger support
    then lower index) achieves |S| <= OPT * (1 + ln(max_p |R_p|)) by Wolsey's
    submodular-cover theorem.

    Returns selected row indices in selection order. If some receptor cannot
    reach depth k with the given pool, covers it as deeply as possible and
    stops when no row gives positive marginal gain (caller must check the
    achieved depth vector for feasibility).
    """
    U = np.asarray(U, dtype=bool)
    n, m = U.shape
    if k < 1:
        raise ValueError("k must be >= 1")
    depth = np.zeros(m, dtype=int)
    remaining = set(range(n))
    selected: list[int] = []
    while (depth < k).any():
        best, best_gain, best_breadth = None, 0, -1
        for i in sorted(remaining):
            deficit = np.maximum(k - depth, 0)
            gain = int((U[i] * deficit).sum())
            breadth = int(U[i].sum())
            if gain > best_gain or (gain == best_gain and gain > 0 and breadth > best_breadth):
                best, best_gain, best_breadth = i, gain, breadth
        if best is None or best_gain == 0:
            break
        selected.append(best)
        remaining.discard(best)
        depth += U[best].astype(int)
    return selected


def escape_union(U_sel: np.ndarray) -> np.ndarray:
    """Boolean vector of receptors in the union of member supports."""
    U_sel = np.atleast_2d(np.asarray(U_sel, dtype=bool))
    return U_sel.any(axis=0)


def exact_escape_probability(u: int, mu: float, generations: int,
                             population: float) -> dict:
    """Exact escape probability for a cocktail whose receptor union has
    size u, per-receptor LOF probability mu per generation, over the given
    number of generations and population size.

    Knockout acquisition is independent per receptor and monotone (no
    reversion), so receptor r is intact after G generations with probability
    (1-mu)^G, and the per-lineage full-escape probability is the closed form

        q = (1 - (1-mu)^G)^u.

    Population: P_esc = 1 - (1-q)^N (exact under lineage independence);
    union_bound = min(1, N*q) is the honest upper bound in the scaling
    theorem. The subset-DP used to cross-check this formula lives in
    tests/test_cocktail.py.
    """
    if u < 0 or not 0 <= mu <= 1:
        raise ValueError("bad u or mu")
    intact = (1.0 - mu) ** generations
    q = (1.0 - intact) ** u
    p_esc = 1.0 - (1.0 - q) ** population if q > 0 else 0.0
    return {"q_single_lineage": q, "p_escape_population": p_esc,
            "union_bound": min(1.0, population * q)}


def member_support_sizes(U_sel: np.ndarray) -> list[int]:
    """|R_p| per member: the minimum knockouts needed to escape that member."""
    return [int(x) for x in np.asarray(U_sel, dtype=int).sum(axis=1)]


def escape_ladder(U_sel: np.ndarray) -> list[dict]:
    """Greedy worst-case host escape order: repeatedly knock out the receptor
    whose loss escapes the most still-active members (ties: lower receptor
    index); if no single knockout escapes anyone, knock out the entire
    remaining support of the easiest remaining member. Member p is escaped
    once K superset of R_p. Returns ordered steps with cumulative escaped
    counts."""
    U_sel = np.atleast_2d(np.asarray(U_sel, dtype=bool))
    n_members, n_rec = U_sel.shape

    def still_active(knocked):
        return np.array([(U_sel[m] & ~knocked).any() for m in range(n_members)])

    knocked = np.zeros(n_rec, dtype=bool)
    steps = []
    while still_active(knocked).any():
        active = still_active(knocked)
        best_r, best_n = None, 0
        for r in range(n_rec):
            if knocked[r]:
                continue
            trial = knocked.copy()
            trial[r] = True
            newly = int((active & ~still_active(trial)).sum())
            if newly > best_n:
                best_r, best_n = r, newly
        if best_r is not None:
            knocked[best_r] = True
            steps.append({"receptors": [int(best_r)],
                          "escaped_after_step": int((~still_active(knocked)).sum())})
        else:
            rem = (U_sel & ~knocked).sum(axis=1)
            target = int(np.argmin(np.where(active, rem, n_rec + 1)))
            rs = [int(r) for r in range(n_rec) if U_sel[target, r] and not knocked[r]]
            for r in rs:
                knocked[r] = True
            steps.append({"receptors": rs,
                          "escaped_after_step": int((~still_active(knocked)).sum())})
    return steps


def leave_one_out(U_sel: np.ndarray, k: int) -> list[dict]:
    """Per member: coverage depth vector and escape-union size after removal."""
    U_sel = np.atleast_2d(np.asarray(U_sel, dtype=bool))
    out = []
    full_union = int(escape_union(U_sel).sum())
    for i in range(len(U_sel)):
        rest = np.delete(U_sel, i, axis=0)
        depth = coverage_depth(rest) if len(rest) else np.zeros(U_sel.shape[1], dtype=int)
        out.append({
            "removed_index": i,
            "depth": depth.tolist(),
            "receptors_below_k": int((depth < k).sum()),
            "escape_union_size": int(escape_union(rest).sum()) if len(rest) else 0,
            "escape_union_delta": full_union - (int(escape_union(rest).sum()) if len(rest) else 0),
        })
    return out


def expected_true_escape_depth(depth: np.ndarray, eps: float) -> float:
    """E[|U_true|] = sum_r (1 - eps^{d_r}) under independent claim-error eps."""
    depth = np.asarray(depth, dtype=float)
    return float((1.0 - np.power(eps, depth)).sum())
