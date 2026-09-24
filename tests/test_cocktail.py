"""Hermetic tests for phage_design.cocktail on toy matrices with known answers."""
import numpy as np
import pytest

from phage_design import cocktail as ck


def test_usage_matrix_threshold():
    S = np.array([[0.5, 0.49], [0.9, 0.1]])
    U = ck.usage_matrix(S, 0.5)
    assert U.tolist() == [[True, False], [True, False]]
    assert ck.usage_matrix(S, 0.9).tolist() == [[False, False], [True, False]]


def test_greedy_k_cover_minimal_and_feasible():
    # rows {0}, {1}, {0,1} with one copy each, k=2: any 2-row choice leaves
    # some receptor at depth 1, so OPT = 3; greedy must match.
    U = np.array([[1, 0], [0, 1], [1, 1]], dtype=bool)
    sel = ck.greedy_k_cover(U, 2)
    assert len(sel) == 3
    assert (ck.coverage_depth(U[sel]) >= 2).all()


def test_greedy_k_cover_respects_wolsey_bound():
    # |greedy| <= OPT * (1 + ln(max support)) must hold on this instance
    rng = np.random.RandomState(0)
    U = rng.rand(30, 6) > 0.5
    U[0] = True  # universal row keeps OPT small
    sel = ck.greedy_k_cover(U, 3)
    delta = U.sum(axis=1).max()
    # brute-force OPT by search over small subsets
    import itertools
    opt = None
    for n in range(1, 8):
        for combo in itertools.combinations(range(len(U)), n):
            if (ck.coverage_depth(U[list(combo)]) >= 3).all():
                opt = n
                break
        if opt:
            break
    assert opt is not None
    assert len(sel) <= opt * (1 + np.log(delta)) + 1e-9


def test_greedy_k_cover_infeasible_pool_stops():
    U = np.array([[1, 0], [1, 0]], dtype=bool)  # receptor 1 can never be covered
    sel = ck.greedy_k_cover(U, 2)
    assert (ck.coverage_depth(U[sel]) >= 2).sum() == 1  # only receptor 0 feasible


def _dp_escape_q(u, mu, G):
    """Brute-force subset DP cross-check of the closed form."""
    dist = np.zeros(1 << u)
    dist[0] = 1.0
    for _ in range(G):
        new = np.zeros(1 << u)
        for s in range(1 << u):
            if dist[s] == 0:
                continue
            missing = [r for r in range(u) if not (s >> r) & 1]
            for sub in range(1 << len(missing)):
                prob, t = dist[s], s
                for j, r in enumerate(missing):
                    if (sub >> j) & 1:
                        prob *= mu
                        t |= 1 << r
                    else:
                        prob *= (1 - mu)
                new[t] += prob
        dist = new
    return dist[-1]


def test_escape_probability_matches_dp():
    for u, mu, G in [(1, 0.1, 5), (3, 0.05, 8), (4, 0.2, 3)]:
        closed = ck.exact_escape_probability(u, mu, G, 1.0)["q_single_lineage"]
        assert np.isclose(closed, _dp_escape_q(u, mu, G), rtol=1e-10)


def test_escape_probability_monotone_in_u():
    shallow = ck.exact_escape_probability(2, 1e-6, 20, 1e8)["p_escape_population"]
    deep = ck.exact_escape_probability(6, 1e-6, 20, 1e8)["p_escape_population"]
    assert deep < shallow
    assert deep >= 0 and shallow <= 1


def test_union_bound_holds():
    r = ck.exact_escape_probability(3, 0.01, 10, 1e5)
    assert r["p_escape_population"] <= r["union_bound"] + 1e-12


def test_member_support_sizes_and_ladder():
    U = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=bool)  # 3 singletons
    assert ck.member_support_sizes(U) == [1, 1, 1]
    steps = ck.escape_ladder(U)
    assert len(steps) == 3  # one knockout escapes one member per step
    assert steps[-1]["escaped_after_step"] == 3
    knocked = sorted(r for s in steps for r in s["receptors"])
    assert knocked == [0, 1, 2]  # every receptor must fall: min escape = |U|


def test_ladder_multi_receptor_member_escapes_last():
    U = np.array([[1, 0], [1, 1]], dtype=bool)  # member 1 needs both knocked
    steps = ck.escape_ladder(U)
    assert steps[-1]["escaped_after_step"] == 2
    total_knocked = sorted(r for s in steps for r in s["receptors"])
    assert total_knocked == [0, 1]


def test_leave_one_out():
    U = np.array([[1, 0], [1, 0], [0, 1]], dtype=bool)
    loo = ck.leave_one_out(U, 2)
    # removing a {0}-member: depth (1,1) -> both receptors below k=2, union intact
    assert loo[0]["receptors_below_k"] == 2
    assert loo[0]["escape_union_delta"] == 0
    # removing the {1}-member: depth (2,0) -> receptor 1 leaves the escape union
    assert loo[2]["escape_union_delta"] == 1
    assert loo[2]["receptors_below_k"] == 1


def test_expected_true_escape_depth():
    depth = np.array([3, 3, 0])
    assert np.isclose(ck.expected_true_escape_depth(depth, 0.0), 2.0)
    assert np.isclose(ck.expected_true_escape_depth(depth, 1.0), 0.0)
    assert np.isclose(ck.expected_true_escape_depth(depth, 0.1),
                      2 * (1 - 0.1 ** 3))
