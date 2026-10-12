"""exactcond: Python (ctypes) binding to the shared C back-end. Pure standard library.

The same C sources back the R package 'exactcond'; the parity harness in the repository requires
every exported quantity to agree byte for byte between the two fronts, which holds because the
random stream lives in C and is seeded per replicate rather than drawn from either language's
generator.

LOADING has two paths, as in the other packages in this programme. An INSTALLED package carries a
library that pip built at install time, so no compiler is needed here. A SOURCE CHECKOUT has the
.c file but no library, so it is compiled on demand into a writable cache directory, which keeps
the edit-a-C-file-and-re-import loop working without writing into site-packages.
"""
import ctypes as _ct
import glob as _glob
import os as _os
import sys as _sys

_here = _os.path.dirname(_os.path.abspath(__file__))


def _find_built():
    for pat in ("_libexactcond*.so", "_libexactcond*.dylib", "_libexactcond*.pyd",
                "_libexactcond*.dll"):
        hits = sorted(_glob.glob(_os.path.join(_here, pat)))
        if hits:
            return hits[0]
    return None


def _cache_dir():
    if _os.access(_here, _os.W_OK):
        return _here
    base = _os.environ.get("XDG_CACHE_HOME") or _os.path.join(_os.path.expanduser("~"), ".cache")
    d = _os.path.join(base, "exactcond")
    _os.makedirs(d, exist_ok=True)
    return d


def _ext():
    return "dll" if _sys.platform == "win32" else ("dylib" if _sys.platform == "darwin" else "so")


def _compile_from_source():
    import subprocess as _sp
    srcs = sorted(f for f in _glob.glob(_os.path.join(_here, "*.c"))
                  if _os.path.basename(f) != "init.c")
    if not srcs:
        return None
    out = _os.path.join(_cache_dir(), "libexactcond." + _ext())
    cc = _os.environ.get("CC") or "cc"
    try:
        _sp.check_call([cc, "-O2", "-fPIC", "-shared", "-I", _here] + srcs + ["-lm", "-o", out])
    except (OSError, _sp.CalledProcessError) as exc:
        raise ImportError(
            "exactcond could not build its C back-end from the bundled sources (%s). In an "
            "installed copy this should never happen, because pip builds the library at install "
            "time. In a source checkout a C compiler is required; set CC if it is not 'cc'." % exc)
    return out


_libpath = _find_built()
if _libpath is None:
    _cached = _os.path.join(_cache_dir(), "libexactcond." + _ext())
    _libpath = _cached if _os.path.exists(_cached) else _compile_from_source()
if _libpath is None:
    raise ImportError("exactcond: no compiled back-end and no C sources to build one from.")

_lib = _ct.CDLL(_libpath)
_dp, _ip = _ct.POINTER(_ct.c_double), _ct.POINTER(_ct.c_int)

for _f, _sig in (
    ("C_ec_rc_null", [_ip, _ip, _ip, _ip, _ip, _ip, _dp, _dp, _ip, _dp]),
    ("C_ec_rc_test", [_ip, _ip, _ip, _ip, _ip, _dp, _dp, _ip, _dp]),
    ("C_ec_rc_two", [_ip, _ip, _ip, _ip, _ip, _ip, _ip, _dp, _dp, _dp, _ip, _dp]),
    ("C_ec_sphere_maxquad", [_dp, _dp, _ip, _ip, _ip, _ip, _dp]),
    ("C_ec_pvalue", [_dp, _dp, _ip, _dp]),
    ("C_ec_multi", [_ip, _ip, _ip, _ip, _ip, _ip, _ip, _ip, _ip, _ip, _ip, _dp, _dp, _dp]),
    ("C_ec_diffnet", [_dp, _ip, _ip, _ip, _ip, _ip, _ip, _ip, _ip, _dp, _dp, _dp]),
):
    getattr(_lib, _f).argtypes = _sig
    getattr(_lib, _f).restype = None

_STATS = {"linlin": 0, "dispersion": 1, "pearson": 2, "lr": 3}


def _d(xs):
    xs = list(xs)
    return (_ct.c_double * len(xs))(*[float(x) for x in xs])


def _i(xs):
    xs = list(xs)
    return (_ct.c_int * len(xs))(*[int(x) for x in xs])


def _one(x):
    return _ct.byref(_ct.c_int(int(x)))


def _oned(x):
    return _ct.byref(_ct.c_double(float(x)))


def _flat(x):
    """Accept a nested sequence or a flat one and return (rows, cols, flat row-major list)."""
    rows = list(x)
    if not rows:
        raise ValueError("empty table")
    if not hasattr(rows[0], "__len__"):
        raise ValueError("expected a two-dimensional table")
    c = len(rows[0])
    if any(len(r) != c for r in rows):
        raise ValueError("ragged table")
    return len(rows), c, [int(v) for r in rows for v in r]


def ec_scores(k, centre=False):
    """Integer scores for an ordinal margin.

    The dispersion and linear-by-linear statistics are unchanged in distribution by an affine
    transformation of the scores, so centring alters the reported statistic but never the p-value.
    """
    s = [float(i) for i in range(int(k))]
    if centre:
        m = sum(s) / len(s)
        s = [v - m for v in s]
    return s


def _scores(u, v, r, c):
    u = ec_scores(r) if u is None else [float(t) for t in u]
    v = ec_scores(c) if v is None else [float(t) for t in v]
    if len(u) != r:
        raise ValueError("u must have one score per row, so length %d" % r)
    if len(v) != c:
        raise ValueError("v must have one score per column, so length %d" % c)
    return u, v


def ec_rc_null(rows, cols, stat="dispersion", u=None, v=None, B=20000, seed=1):
    """Draws of a table statistic from the law conditional on BOTH margins."""
    rows = [int(t) for t in rows]
    cols = [int(t) for t in cols]
    if sum(rows) != sum(cols):
        raise ValueError("the row and column totals must agree")
    uu, vv = _scores(u, v, len(rows), len(cols))
    out = (_ct.c_double * int(B))()
    _lib.C_ec_rc_null(_i(rows), _i(cols), _one(len(rows)), _one(len(cols)), _one(B),
                      _one(_STATS[stat]), _d(uu), _d(vv), _one(seed), out)
    return list(out)


def ec_rc_test(x, stat="dispersion", u=None, v=None, B=20000, seed=1):
    """Exact conditional test in a two-way table, conditioning on both margins.

    Returns a dict with the observed statistic, the exact p-value and its mid-p companion, a
    Wilson interval reflecting the Monte Carlo error in B draws, the largest null draw, and the
    row mean scores, which for the dispersion statistic ARE the contrasts being reported.
    """
    r, c, flat = _flat(x)
    if any(t < 0 for t in flat):
        raise ValueError("x must hold non-negative counts")
    uu, vv = _scores(u, v, r, c)
    out = (_ct.c_double * 6)()
    _lib.C_ec_rc_test(_i(flat), _one(r), _one(c), _one(B), _one(_STATS[stat]),
                      _d(uu), _d(vv), _one(seed), out)
    n = [sum(flat[i * c:(i + 1) * c]) for i in range(r)]
    means = [sum(vv[j] * flat[i * c + j] for j in range(c)) / n[i] if n[i] else float("nan")
             for i in range(r)]
    return {"statistic": out[0], "p_value": out[1], "midp": out[2], "mc_lower": out[3],
            "mc_upper": out[4], "max_null": out[5], "row_means": means, "n": n,
            "B": int(B), "stat": stat}


def ec_rc_two(t1, t2, u=None, v=None, B=20000, burn=20000, thin=20, alpha=0.05, seed=1):
    """Exact conditional test that two tables share a linear-by-linear association parameter.

    Under the null of a common parameter the joint density carries it only through the total score,
    which the conditioning fixes, so it leaves the density altogether. The law is sampled by moves
    making a two by two swap in each table with equal and opposite score changes, starting at the
    observed pair -- which is what keeps it usable when both tables are strongly associated, the
    case in which convolving two separately estimated null distributions fails outright.
    """
    r, c, f1 = _flat(t1)
    r2, c2, f2 = _flat(t2)
    if (r, c) != (r2, c2):
        raise ValueError("the two tables must have the same dimension")
    if r < 2 or c < 2:
        raise ValueError("both tables must have at least two rows and two columns")
    uu, vv = _scores(u, v, r, c)
    if any(abs(t - round(t)) > 1e-9 for t in uu + vv):
        raise ValueError("ec_rc_two needs integer scores, so the total score lies on a lattice")
    out = (_ct.c_double * 5)()
    _lib.C_ec_rc_two(_i(f1), _i(f2), _one(r), _one(c), _one(B), _one(burn), _one(thin),
                     _d(uu), _d(vv), _oned(alpha), _one(seed), out)
    if out[0] < 0:
        raise RuntimeError("the sampler could not be initialised on these tables")
    return {"midp": out[0], "size": out[1], "T1": out[2], "S": out[3], "accept": out[4],
            "B": int(B), "burn": int(burn), "thin": int(thin), "alpha": float(alpha)}


def ec_sphere_maxquad(M, A, B=20000, seed=1):
    """Exact null of a quadratic form maximised over candidate structures.

    M is the residual maker; A is one symmetric matrix or a list of them. Each draw takes the
    maximum INSIDE the draw, which is what makes the reference valid after a structure has been
    chosen by fit.
    """
    rows = [list(r) for r in M]
    n = len(rows)
    if any(len(r) != n for r in rows):
        raise ValueError("M must be square")
    if not hasattr(A[0][0], "__len__"):
        A = [A]
    flat = []
    for a in A:
        ar = [list(r) for r in a]
        if len(ar) != n or any(len(r) != n for r in ar):
            raise ValueError("every candidate must be n by n")
        for i in range(n):
            for j in range(n):
                flat.append(float(ar[i][j]))
    Mf = [float(rows[i][j]) for i in range(n) for j in range(n)]
    out = (_ct.c_double * int(B))()
    _lib.C_ec_sphere_maxquad(_d(Mf), _d(flat), _one(n), _one(len(A)), _one(B), _one(seed), out)
    return list(out)


def ec_multi_test(x, group, cluster=None, stat="pearson", B=20000, seed=1):
    """One exact test for a battery of categorical items across groups.

    x is a sequence of items (each a sequence of integer codes from zero, or from one, over the
    respondents) or a sequence of respondent rows; group gives each respondent's group; cluster,
    if given, each respondent's cluster, every cluster inside one group, in which case the labels
    are permuted across clusters. The sum and the maximum over items of the per-item statistic
    (Pearson on the group-by-level table, or the dispersion of group mean scores) are each
    referred to the same permutations, giving one p-value for the battery, and the per-item
    statistics are referred to those permutations too, for a closed follow-up.
    """
    rows = [list(r) for r in x]
    n_g = len(list(group))
    if len(rows) == n_g and (not rows or hasattr(rows[0], "__len__")) and len(rows[0]) != n_g:
        items = [[int(r[j]) for r in rows] for j in range(len(rows[0]))]   # respondent rows
    else:
        items = [[int(v) for v in r] for r in rows]                          # item columns
    n = len(items[0]); m = len(items)
    if any(len(c) != n for c in items) or n_g != n:
        raise ValueError("every item and the group vector must have one entry per respondent")
    items = [[v - 1 for v in c] if min(c) == 1 else c for c in items]
    if any(min(c) < 0 for c in items):
        raise ValueError("items must be coded from zero or one")
    k = [max(c) + 1 for c in items]
    levels = sorted(set(group)); grp = [levels.index(v) for v in group]; g = len(levels)
    if g < 2:
        raise ValueError("group must have at least two levels")
    if cluster is None:
        cl = [0] * n; ncl = 0
    else:
        cluster = list(cluster)
        if len(cluster) != n:
            raise ValueError("cluster must have one entry per respondent")
        clev = sorted(set(cluster)); cl = [clev.index(v) for v in cluster]; ncl = len(clev)
        seen = {}
        for c, a in zip(cl, grp):
            if seen.setdefault(c, a) != a:
                raise ValueError("every cluster must lie inside one group")
    which = {"pearson": 0, "dispersion": 1}[stat]
    out = (_ct.c_double * 12)(); ist = (_ct.c_double * m)(); ip = (_ct.c_double * m)()
    _lib.C_ec_multi(_i([v for c in items for v in c]), _one(n), _one(m), _one(g), _i(grp), _i(k),
                    _i(cl), _one(ncl), _one(which), _one(B), _one(seed), out, ist, ip)
    def part(i, j, mx):   # obs at i; p, midp, lo, hi from j
        return {"statistic": out[i], "p_value": out[j], "midp": out[j + 1],
                "mc_lower": out[j + 2], "mc_upper": out[j + 3], "max_null": mx}
    return {"sum": part(0, 2, out[10]), "max": part(1, 6, out[11]),
            "items": [{"statistic": ist[j], "p_value": ip[j]} for j in range(m)],
            "n": [grp.count(a) for a in range(g)], "clusters": ncl, "B": int(B), "stat": stat,
            "exchangeable": "clusters" if ncl > 0 else "respondents"}


def _ranks(v):
    """Average ranks, ties shared, which is what makes the correlation Spearman's."""
    order = sorted(range(len(v)), key=lambda i: v[i])
    out = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        r = (i + j) / 2.0 + 1.0
        for t in range(i, j + 1):
            out[order[t]] = r
        i = j + 1
    return out


def ec_diffnet_test(x, group, cluster=None, rank=True, B=20000, seed=1):
    """A differential network without a Gaussian assumption.

    x is a sequence of variables (each a sequence over the respondents) or of respondent rows.
    The edge statistic is the heterogeneity of the Fisher-transformed correlations across groups,
    the same quantity the Gaussian theory uses, referred to the law obtained by permuting the group
    labels rather than to a chi-square; with rank=True the correlations are Spearman's and the
    procedure is distribution-free. Returns the sum and the maximum over edges, each with one
    p-value for the whole network, and the per-edge statistics and p-values from those same
    permutations for a closed follow-up.
    """
    rows = [list(r) for r in x]
    n_g = len(list(group))
    if len(rows) == n_g and rows and hasattr(rows[0], "__len__") and len(rows[0]) != n_g:
        cols = [[float(r[j]) for r in rows] for j in range(len(rows[0]))]
    else:
        cols = [[float(v) for v in r] for r in rows]
    p = len(cols); n = len(cols[0])
    if p < 2:
        raise ValueError("x must hold at least two variables")
    if any(len(c) != n for c in cols) or n_g != n:
        raise ValueError("every variable and the group vector must have one entry per respondent")
    if rank:
        cols = [_ranks(c) for c in cols]
    levels = sorted(set(group)); grp = [levels.index(v) for v in group]; g = len(levels)
    if g < 2:
        raise ValueError("group must have at least two levels")
    if min(grp.count(a) for a in range(g)) < 4:
        raise ValueError("every group needs at least four respondents")
    if cluster is None:
        cl = [0] * n; ncl = 0
    else:
        cluster = list(cluster)
        if len(cluster) != n:
            raise ValueError("cluster must have one entry per respondent")
        clev = sorted(set(cluster)); cl = [clev.index(v) for v in cluster]; ncl = len(clev)
        seen = {}
        for c, a in zip(cl, grp):
            if seen.setdefault(c, a) != a:
                raise ValueError("every cluster must lie inside one group")
    m = p * (p - 1) // 2
    out = (_ct.c_double * 12)(); est = (_ct.c_double * m)(); ep = (_ct.c_double * m)()
    _lib.C_ec_diffnet(_d([v for c in cols for v in c]), _one(n), _one(p), _one(g), _i(grp),
                      _i(cl), _one(ncl), _one(B), _one(seed), out, est, ep)
    def part(i, j, mx):
        return {"statistic": out[i], "p_value": out[j], "midp": out[j + 1],
                "mc_lower": out[j + 2], "mc_upper": out[j + 3], "max_null": mx}
    pairs = [(a, b) for a in range(p) for b in range(a + 1, p)]
    return {"sum": part(0, 2, out[10]), "max": part(1, 6, out[11]),
            "edges": [{"from": a, "to": b, "statistic": est[e], "p_value": ep[e]}
                      for e, (a, b) in enumerate(pairs)],
            "n": [grp.count(a) for a in range(g)], "clusters": ncl, "B": int(B),
            "ranked": bool(rank), "exchangeable": "clusters" if ncl > 0 else "respondents"}


def ec_pvalue(obs, null):
    """Tail probabilities from an observed value and a sample from the null."""
    null = [float(t) for t in null]
    out = (_ct.c_double * 4)()
    _lib.C_ec_pvalue(_oned(obs), _d(null), _one(len(null)), out)
    return {"p_value": out[0], "midp": out[1], "mc_lower": out[2], "mc_upper": out[3]}
