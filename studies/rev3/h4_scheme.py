"""A complete-by-construction sampling scheme for hybrid4.

The w-simplex (w >= 0, sum w <= 1) guarantees Kipping's conditions but covers
only 13% of the exact region P (h4_completeness.md). Two alternatives:

(I) The curvature simplex S*. Write g(xi) = dI/dxi, xi = mu^2. For hybrid4
        g   = 2 w1 xi + 2 w2 (1 - xi) + 2 a0 w3 (xi + e0)^-3 + 2 a1 w4 (xi + e1)^-3
        g'' = 24 a0 w3 (xi + e0)^-5 + 24 a1 w4 (xi + e1)^-5        (a_k = 1/N_k)
    g'' lies in a two-function span whose ratio is monotonic, so g'' >= 0 on
    [0, 1] iff g''(0) >= 0 and g''(1) >= 0. Then g' is increasing, so g' <= 0
    on [0, 1] iff g'(1) <= 0; then g is decreasing, so g >= 0 iff g(1) >= 0.
    Hence the hybrid4 profiles with
        I(0) >= 0,  I' >= 0,  I'' <= 0,  I''' >= 0   (derivatives in xi on [0,1])
    are exactly the 4-simplex of five linear inequalities
        sum w <= 1, g(1) >= 0, g'(1) <= 0, g''(0) >= 0, g''(1) >= 0,
    and a flat Dirichlet on its five vertices (stick-breaking from the unit
    hypercube, as for the w-simplex) samples it uniformly and completely.
    The two extra conditions are tested here against the tables, the raw
    model profiles and the fits that leave the w-simplex.

(II) The canonical map onto all of P. {xi, 1 - xi, (xi+e0)^-3, (xi+e1)^-3}
    is an extended Chebyshev system on [0, 1] (g'' has at most one zero, so
    g has at most three, counting multiplicity). Every g >= 0 in its span is
    then a sum of two extreme elements: E0(t), with zeros at 0 and a double
    zero at t, and E1(t), with a double zero at t and a zero at 1 (the
    analogue of Markov-Lukacs, p = x A^2 + (1 - x) B^2 for cubics). With the
    E normalised to sum w = 1,
        w = s [theta v0(t0) + (1 - theta) v1(t1)],   s, theta, t0, t1 in [0, 1]
    is onto P and physical by construction, but not uniform.

    .venv/bin/python -W ignore studies/rev3/h4_scheme.py   -> h4_scheme.md

Works in the monomial-pole basis (laws.hybrid4_monomial_law); ran at commit
ecfa090, when the tables were w-simplex fits in that basis. Since 2026-10-14
hybrid4 is tabulated in the curvature-simplex basis (ldhybrid.v1.HYBRID4_BASIS).
"""

import os
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from math import factorial
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).parent))
from ldhybrid import fit, laws, v1  # noqa: E402
from ldhybrid.interpolate import LDTables  # noqa: E402

HERE = Path(__file__).parent
OUT = HERE / "h4_scheme.md"
CDS = ROOT / "data" / "release" / "cds"
PROF = ROOT / "data" / "release" / "zenodo"
E0, E1 = v1.HYBRID4_EPS
A0, A1 = (1.0 / (e ** -2 - (1 + e) ** -2) for e in (E0, E1))
RNG = np.random.default_rng(20261014)


# ------------------------------------------------------------ g and derivatives
def u(x):
    """g(x) = u(x) . w, with g = dI/dxi."""
    x = np.asarray(x, float)
    return np.stack([2 * x, 2 * (1 - x), 2 * A0 * (x + E0) ** -3, 2 * A1 * (x + E1) ** -3], -1)


def u1(x):
    x = np.asarray(x, float)
    return np.stack([2 + 0 * x, -2 + 0 * x, -6 * A0 * (x + E0) ** -4, -6 * A1 * (x + E1) ** -4], -1)


def u2(x):
    x = np.asarray(x, float)
    return np.stack([0 * x, 0 * x, 24 * A0 * (x + E0) ** -5, 24 * A1 * (x + E1) ** -5], -1)


# ------------------------------------------------------------ (I) curvature simplex
def sstar():
    """S* as A w <= b, and its five vertices."""
    A = np.array([np.ones(4), -u(1.0), u1(1.0), -u2(0.0), -u2(1.0)])
    b = np.array([1.0, 0, 0, 0, 0])
    V = np.array([np.linalg.solve(A[list(S)], b[list(S)]) for S in combinations(range(5), 4)])
    return A, b, V


def q_to_w(q, V):
    """Unit hypercube -> S*: flat-Dirichlet stick-breaking onto the vertices."""
    q = np.atleast_2d(q)
    n = q.shape[1]
    lam = np.zeros((len(q), n + 1))
    rest = np.ones(len(q))
    for j in range(n):
        x = 1 - q[:, j] ** (1.0 / (n - j))
        lam[:, j] = rest * x
        rest = rest * (1 - x)
    lam[:, n] = rest
    return lam @ V


def w_to_q(w, V):
    w = np.atleast_2d(w)
    M = np.vstack([V.T, np.ones(len(V))])
    lam = np.linalg.solve(M, np.vstack([w.T, np.ones(len(w))])).T
    n = lam.shape[1] - 1
    q = np.zeros((len(w), n))
    rest = np.ones(len(w))
    for j in range(n):
        x = lam[:, j] / rest
        q[:, j] = (1 - x) ** (n - j)
        rest = rest - lam[:, j]
    return q


# ------------------------------------------------------------ (II) canonical map
def extreme(t, end):
    """Normalised extreme element: zeros at `end` (0 or 1) and a double zero at t."""
    M = np.array([u(float(end)), u(t), u1(t)])
    v = np.linalg.svd(M)[2][-1]
    return v / v.sum()


def canon(q):
    s, th, t0, t1 = q
    return s * (th * extreme(t0, 0) + (1 - th) * extreme(t1, 1))


def canon_inverse(w, n=161):
    """(s, theta, t0, t1) for w in P, by a grid search and a polish."""
    s = w.sum()
    b = w / s
    ts = np.linspace(1e-4, 1 - 1e-4, n)
    P0 = np.array([extreme(t, 0) for t in ts])
    P1 = np.array([extreme(t, 1) for t in ts])
    best = (np.inf, None)
    for i in range(n):
        D = P0[i][None, :] - P1                       # b - P1 = theta (P0 - P1)
        R = b[None, :] - P1
        th = np.einsum("ij,ij->i", D, R) / np.einsum("ij,ij->i", D, D)
        th = np.clip(th, 0, 1)
        res = np.linalg.norm(R - th[:, None] * D, axis=1)
        k = int(np.argmin(res))
        if res[k] < best[0]:
            best = (res[k], (th[k], ts[i], ts[k]))
    f = lambda p: p[0] * extreme(p[1], 0) + (1 - p[0]) * extreme(p[2], 1) - b
    r = least_squares(f, best[1], bounds=([0, 1e-6, 1e-6], [1, 1 - 1e-6, 1 - 1e-6]), xtol=1e-15, ftol=1e-15)
    return np.r_[s, r.x], float(np.max(np.abs(r.fun)))


# ------------------------------------------------------------ exact region P
def region_P():
    return fit._inequalities(laws.hybrid4_monomial_law(), "physical")


def sample_P(n, C, d):
    lo = np.array([-0.0922, -1.5798, -0.1355, -1.0902])
    hi = np.array([1.3612, 1.7916, 1.0806, 2.1109])
    out = []
    while sum(len(o) for o in out) < n:
        W = lo + (hi - lo) * RNG.random((50_000, 4))
        out.append(W[np.all(W @ C.T <= d, axis=1)])
    return np.vstack(out)[:n]


# ------------------------------------------------------------ fits with S*
def _ineq(law, constraint, _orig=fit._inequalities):
    if constraint != "sstar":
        return _orig(law, constraint)
    A, b, _ = sstar()
    n = np.max(np.abs(A), axis=1)
    return A / n[:, None], b / n


_G = {}


def task(job):
    fit._inequalities = _ineq
    from ldhybrid.profile import Profile
    mu, y, measure, do_lc = job
    L = laws.hybrid4_monomial_law()
    pr = Profile(mu, y)
    out = {}
    for con in ("simplex", "physical", "sstar"):
        f = fit.fit_linear(L, pr, measure, v1.RESAMPLE, con, v1.FCM)
        out[con] = dict(w=np.asarray(f.w), rms=f.rms)
    if do_lc:
        if "g" not in _G:
            from ldhybrid import lightcurve as lc
            _G["g"] = lc.LCGrid()
            _G["cols"] = _G["g"].linear_columns(L)
        g, (one, cols) = _G["g"], _G["cols"]
        truth = g.columns(pr)
        for con in out:
            out[con]["lc"] = g.errors(*g.combine(one, cols, out[con]["w"]), *truth)["max_ppm"]
    return out


def bound_jobs():
    import h5py
    d = LDTables(CDS, cache=False).table("phoenix", "hybrid4")
    with h5py.File(PROF / "ldhybrid_v1.0_band_profiles_phoenix.h5", "r") as h:
        P = {k: h[k][()] for k in ("teff", "logg", "mh", "mu", "y")}
        bands = [b.decode() for b in h["bands"][()]]
    key = {(round(a, 3), round(b, 3), round(c, 3)): i
           for i, (a, b, c) in enumerate(zip(P["teff"], P["logg"], P["mh"]))}
    jobs, meta = [], []
    for method in ("r", "mu"):
        for r in np.flatnonzero((d["method"] == method) & (d["Bound"] == 1)):
            i = key[(round(d["teff"][r], 3), round(d["logg"][r], 3), round(d["mh"][r], 3))]
            j = bands.index(d["band"][r])
            ok = np.isfinite(P["mu"][i, j])
            jobs.append((P["mu"][i, j][ok].astype(float), P["y"][i, j][ok].astype(float),
                         method, method == "r"))
            meta.append(method)
    return jobs, np.array(meta)


# ------------------------------------------------------------ raw-profile curvature
def profile_curvature(grid, cuts=(0.0, 0.05, 0.1)):
    """Fraction of (model, band) profiles whose node finite differences in xi
    have I'' > 1e-3 or I''' < -1e-2 anywhere at mu >= cut."""
    import h5py
    with h5py.File(PROF / f"ldhybrid_v1.0_band_profiles_{grid}.h5", "r") as h:
        mu, y = h["mu"][()], h["y"][()]
    bad = {c: [0, 0] for c in cuts}
    tot = 0
    for i in range(mu.shape[0]):
        for j in range(mu.shape[1]):
            m = np.isfinite(mu[i, j])
            if m.sum() < 6:
                continue
            tot += 1
            x0, Y0 = mu[i, j][m].astype(float), y[i, j][m].astype(float)
            o = np.argsort(x0)
            x0, Y0 = x0[o], Y0[o]
            for c in cuts:
                k = x0 >= c
                x, Y = x0[k] ** 2, Y0[k]
                d1 = np.diff(Y) / np.diff(x)
                xm = (x[1:] + x[:-1]) / 2
                d2 = np.diff(d1) / np.diff(xm)
                xm2 = (xm[1:] + xm[:-1]) / 2
                d3 = np.diff(d2) / np.diff(xm2)
                bad[c][0] += bool(np.any(d2 > 1e-3))
                bad[c][1] += bool(np.any(d3 < -1e-2))
    return tot, bad


# ------------------------------------------------------------ report
def main():
    A, b, V = sstar()
    C, d = region_P()
    L = ["# A complete sampling scheme for hybrid4", "",
         "Generated by `studies/rev3/h4_scheme.py`; see its docstring for the derivations.", ""]

    # (I) geometry
    vol = abs(np.linalg.det(V[1:] - V[0])) / factorial(4)
    L += ["## I. The curvature simplex S*", "",
          "Conditions (derivatives in xi = mu^2 on [0, 1]): I(0) >= 0, I' >= 0, I'' <= 0, I''' >= 0.",
          "They are equivalent to five linear inequalities, so S* is a 4-simplex with vertices", "",
          "| vertex | w1 | w2 | w3 | w4 | profile |", "|---|---|---|---|---|---|"]
    names = {}
    for v in V:
        Ix = lambda x, v=v: 1 - (laws.hybrid4_monomial_law().T(np.sqrt(x)) @ v)
        tag = ("I = 1 (no darkening)" if np.allclose(v, 0, atol=1e-12) else
               "I = 1 - (1 - mu^2)^2" if np.allclose(v, [0, 1, 0, 0], atol=1e-12) else
               "I = mu^2" if np.allclose(v, [.5, .5, 0, 0], atol=1e-12) else
               f"I(0) = 0, pole-dominated (I(mu=0.1) = {Ix(0.01):.3f}, I(0.5) = {Ix(0.25):.3f})")
        L.append("| | " + " | ".join(f"{x:+.5f}" for x in v) + f" | {tag} |")
    W = sample_P(400_000, C, d)
    inS = np.all(W @ A.T <= b + 1e-12, axis=1)
    inW = np.all(W >= 0, axis=1)
    Wd = RNG.dirichlet(np.ones(5), 400_000)
    old_in_S = np.mean(np.all(Wd[:, :4] @ A.T <= b + 1e-12, axis=1))
    Q = RNG.random((200_000, 4))
    Ws = q_to_w(Q, V)
    phys = np.all(Ws @ C.T <= d + 1e-10, axis=1).mean()
    rt = np.max(np.abs(w_to_q(Ws, V) - Q))
    L += ["", f"Volume {vol:.4f} (w-simplex 1/24 = {1 / 24:.4f}; exact region P = 0.322). "
          f"Of P, S* holds {100 * inS.mean():.1f}% and the w-simplex {100 * inW.mean():.1f}%; "
          f"{100 * np.mean(inS & ~inW) / max(inS.mean(), 1e-30):.1f}% of S* lies outside the w-simplex, and "
          f"{100 * (1 - old_in_S):.1f}% of the w-simplex lies outside S* (profiles convex in mu^2, "
          "e.g. the vertex I = mu^4).", "",
          f"q-map: {len(Q)} uniform q -> S*: {100 * phys:.4f}% satisfy the exact conditions I-II; "
          f"round trip q -> w -> q max error {rt:.1e}.", ""]

    # (I) data support
    t = LDTables(CDS, cache=False)
    L += ["### Support in the tables", "",
          "| grid | method | rows | in S* | min margin -I'' | min I''' |", "|---|---|---|---|---|---|"]
    xi = np.linspace(0, 1, 2001)
    for grid in ("phoenix", "mps1", "mps2"):
        D = t.table(grid, "hybrid4")
        for m in ("r", "mu"):
            Wt = laws.hybrid4_monomial(D["coeffs"][D["method"] == m])
            ins = np.all(Wt @ A.T <= b + 2e-6, axis=1)       # 2e-6: 6-decimal rounding
            L.append(f"| {grid} | {m} | {len(Wt)} | {100 * ins.mean():.4f}% | "
                     f"{-(Wt @ u1(xi).T).max():.4f} | {(Wt @ u2(xi).T).min():.4f} |")
    L += ["", "On the w-simplex I''' >= 0 holds automatically (it depends only on w3, w4 >= 0); "
          "the MPS-ATLAS fits never meet the simplex constraint, so for them both conditions are "
          "properties of the unconstrained fits.", ""]
    L += ["### Support in the raw model profiles (node finite differences)", "",
          "| grid | profiles | I'' > 0 somewhere: all mu / mu >= 0.05 / mu >= 0.1 | "
          "I''' < 0 somewhere: all mu / mu >= 0.05 / mu >= 0.1 |", "|---|---|---|---|"]
    for grid in ("phoenix", "mps2"):
        tot, bad = profile_curvature(grid)
        L.append(f"| {grid} | {tot} | " + " / ".join(f"{100 * bad[c][0] / tot:.2f}%" for c in bad)
                 + " | " + " / ".join(f"{100 * bad[c][1] / tot:.2f}%" for c in bad) + " |")
    L += ["", "Violations sit at the extreme limb (mu < 0.1: the PHOENIX spherical limb, and "
          "a few MPS-ATLAS blue bands at mu = 0.01-0.03), a scale hybrid4 does not resolve.", ""]

    # (I) bound PHOENIX fits
    jobs, meta = bound_jobs()
    with ProcessPoolExecutor(min(4, os.cpu_count() or 1)) as ex:
        res = list(ex.map(task, jobs, chunksize=20))
    L += ["### Fits where the w-simplex binds (PHOENIX)", "",
          "| method | rows | exact-region fit in S* | S* fit outside w-simplex | "
          "rms I: simplex / S* / exact (median) |", "|---|---|---|---|---|"]
    for m in ("r", "mu"):
        R = [r for r, mm in zip(res, meta) if mm == m]
        Wp = np.array([r["physical"]["w"] for r in R])
        Wst = np.array([r["sstar"]["w"] for r in R])
        L.append(f"| {m} | {len(R)} | {100 * np.mean(np.all(Wp @ A.T <= b + 1e-7, axis=1)):.1f}% | "
                 f"{100 * np.mean(np.min(Wst, axis=1) < -1e-8):.1f}% | "
                 + " / ".join(f"{np.median([r[c]['rms'] for r in R]):.2e}" for c in ("simplex", "sstar", "physical"))
                 + " |")
    R = [r for r, mm in zip(res, meta) if mm == "r"]
    lc = {c: np.array([r[c]["lc"] for r in R]) for c in ("simplex", "sstar", "physical")}
    L += ["", "Light-curve error of the r-weighted fits (max ppm, production grid): "
          + "; ".join(f"{c} median {np.median(v):.1f}, max {v.max():.1f}" for c, v in lc.items()) + ".", ""]

    # (II) canonical map
    Qc = RNG.random((20_000, 4))
    Wc = np.array([canon(q) for q in Qc])
    physc = np.all(Wc @ C.T <= d + 1e-9, axis=1).mean()
    Wt = sample_P(1000, C, d)
    errs = np.array([canon_inverse(w)[1] for w in Wt])
    frac_simplex = np.mean(np.all(Wc >= 0, axis=1))
    frac_S = np.mean(np.all(Wc @ A.T <= b + 1e-12, axis=1))
    L += ["## II. The canonical map onto the whole exact region", "",
          f"Forward: {len(Qc)} uniform q -> {100 * physc:.3f}% physical. "
          f"Completeness: of {len(Wt)} points drawn uniformly in P, {np.sum(errs < 1e-6)} are "
          f"reproduced to 1e-6 (max residual {errs.max():.1e}).",
          f"Non-uniformity: uniform q puts {100 * frac_simplex:.1f}% of draws in the w-simplex "
          f"(volume share {100 * inW.mean():.1f}%) and {100 * frac_S:.1f}% in S* "
          f"(volume share {100 * inS.mean():.1f}%).", ""]
    OUT.write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
