"""Tests of the Python front. The point defended is exactness, so the central test compares the
sampler with the conditional law obtained by complete enumeration."""
import math
from itertools import product
import exactcond as e


def enumerate_law(rs, cs, v):
    """Every table with these margins, weighted by the Fisher-Yates probability."""
    N = sum(rs)
    r, c = len(rs), len(cs)
    out = {}
    def rows(rem_cols, i, acc):
        if i == r - 1:
            if all(t >= 0 for t in rem_cols):
                yield acc + [list(rem_cols)]
            return
        for cell in product(*[range(t + 1) for t in rem_cols]):
            if sum(cell) != rs[i]:
                continue
            yield from rows([rem_cols[j] - cell[j] for j in range(c)], i + 1, acc + [list(cell)])
    for tab in rows(list(cs), 0, []):
        lw = (sum(math.lgamma(t + 1) for t in rs) + sum(math.lgamma(t + 1) for t in cs)
              - math.lgamma(N + 1) - sum(math.lgamma(x + 1) for row in tab for x in row))
        T = sum(v[i] * v[j] * tab[i][j] for i in range(r) for j in range(c))
        out[T] = out.get(T, 0.0) + math.exp(lw)
    tot = sum(out.values())
    return {k: w / tot for k, w in out.items()}


def test_sampler_matches_enumeration():
    rs, cs, v = [3, 2, 2], [3, 2, 2], [0, 1, 2]
    exact = enumerate_law(rs, cs, v)
    draws = e.ec_rc_null(rs, cs, "linlin", None, None, 200000, 21)
    for t, p in exact.items():
        got = sum(1 for d in draws if abs(d - t) < 1e-9) / len(draws)
        assert abs(got - p) < 0.005, (t, p, got)


def test_dispersion_components_are_the_contrasts():
    tab = [[3, 228, 72], [10, 548, 67], [7, 276, 131]]
    r = e.ec_rc_test(tab, "dispersion", None, [-1, 0, 1], 1000, 1)
    want = [(-row[0] + row[2]) / sum(row) for row in tab]
    for a, b in zip(r["row_means"], want):
        assert abs(a - b) < 1e-12
    gm = sum(-row[0] + row[2] for row in tab) / sum(sum(row) for row in tab)
    T = sum(sum(row) * (m - gm) ** 2 for row, m in zip(tab, r["row_means"]))
    assert abs(r["statistic"] - T) < 1e-9


def test_affine_rescoring_leaves_the_p_value_alone():
    tab = [[3, 228, 72], [10, 548, 67], [7, 276, 131]]
    a = e.ec_rc_test(tab, "dispersion", None, [-1, 0, 1], 20000, 41)
    b = e.ec_rc_test(tab, "dispersion", None, [4, 7, 10], 20000, 41)
    assert a["p_value"] == b["p_value"]


def test_identical_strong_tables_give_no_evidence():
    a = [[20, 8, 3], [9, 22, 7], [4, 8, 19]]
    z = e.ec_rc_two(a, a, None, None, 20000, 20000, 20, 0.05, 3)
    assert z["midp"] > 0.5
    assert z["S"] == 2 * z["T1"]


def test_opposed_association_is_detected():
    a = [[20, 8, 3], [9, 22, 7], [4, 8, 19]]
    z = e.ec_rc_two(a, [row[::-1] for row in a], None, None, 20000, 20000, 20, 0.05, 3)
    assert z["midp"] < 0.001


def test_sphere_draws_have_the_beta_law():
    n = 7
    M = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    A = [[1.0 if (i == 0 and j == 0) else 0.0 for j in range(n)] for i in range(n)]
    d = e.ec_sphere_maxquad(M, A, 40000, 51)
    mean = sum(d) / len(d)
    assert abs(mean - 1.0 / n) < 0.01          # E[u_1^2] = 1/n exactly
    assert all(0.0 <= t <= 1.0 for t in d)


def test_misuse_is_refused():
    import pytest
    with pytest.raises(ValueError):
        e.ec_rc_null([10, 10], [10, 11])
    with pytest.raises(ValueError):
        e.ec_rc_test([[1, 2], [3, 4]], "dispersion", None, [1])
    with pytest.raises(ValueError):
        e.ec_rc_two([[1, 2], [3, 4]], [[1, 2], [3, 4]], None, [0, 0.5])


def test_multi_one_binary_item_matches_the_table_test():
    import random
    random.seed(5)
    g = ["a"] * 40 + ["b"] * 50 + ["c"] * 45
    x = [1 if random.random() < (0.6 if v == "c" else 0.85) else 0 for v in g]
    om = e.ec_multi_test([x], g, None, "pearson", 40000, 9)
    tab = [[sum(1 for a, b in zip(g, x) if a == lev and b == c) for c in (0, 1)] for lev in ("a", "b", "c")]
    rc = e.ec_rc_test(tab, "pearson", None, None, 40000, 9)
    assert abs(om["sum"]["statistic"] - rc["statistic"]) < 1e-9
    assert abs(om["sum"]["p_value"] - rc["p_value"]) < 0.01


def test_multi_cluster_must_nest():
    g = ["a"] * 60 + ["b"] * 60
    x = [[i % 2 for i in range(120)], [(i // 3) % 2 for i in range(120)]]
    r = e.ec_multi_test(x, g, [i // 10 for i in range(120)], "pearson", 500, 3)
    assert r["clusters"] == 12 and r["exchangeable"] == "clusters"
    try:
        e.ec_multi_test(x, g, [i % 12 for i in range(120)], "pearson", 100, 3)
    except ValueError as exc:
        assert "inside one group" in str(exc)
    else:
        raise AssertionError("a straddling cluster was accepted")


def test_diffnet_finds_the_planted_edge():
    import random
    random.seed(11)
    n, p = 200, 5
    g = ["a"] * (n // 2) + ["b"] * (n // 2)
    cols = [[random.gauss(0, 1) for _ in range(n)] for _ in range(p)]
    for i in range(n // 2, n):                       # one edge differs in group b
        cols[1][i] = 0.8 * cols[0][i] + 0.6 * cols[1][i]
    r = e.ec_diffnet_test(cols, g, None, True, 999, 1)
    assert r["max"]["p_value"] < 0.01
    top = min(r["edges"], key=lambda d: d["p_value"])
    assert (top["from"], top["to"]) == (0, 1)


def test_diffnet_is_invariant_to_affine_rescoring():
    import random
    random.seed(12)
    n, p = 150, 4
    g = ["a"] * 75 + ["b"] * 75
    cols = [[random.gauss(0, 1) for _ in range(n)] for _ in range(p)]
    a = e.ec_diffnet_test(cols, g, None, True, 499, 5)
    cols2 = [list(c) for c in cols]
    cols2[1] = [3.5 * v - 7 for v in cols2[1]]
    b = e.ec_diffnet_test(cols2, g, None, True, 499, 5)
    assert a["max"]["p_value"] == b["max"]["p_value"]
    for x, y in zip(a["edges"], b["edges"]):
        assert abs(x["statistic"] - y["statistic"]) < 1e-9
