"""Double check: does any hybrid4 profile on the simplex (w >= 0, sum w <= 1)
violate the Kipping (2013) physicality conditions,
  (I)  I(mu) >= 0 for 0 <= mu <= 1, and
  (II) dI/dmu >= 0 (monotonically brighter towards the centre)?

Checks: the shapes; dense random and extreme points of the simplex; the q-map;
every tabulated hybrid4 row (all grids, both methods); every corrected set
(w + Δw); interpolated coefficients; the standalone files.

    .venv/bin/python studies/rev3/h4_kipping_check.py
"""
import sys
import warnings
from itertools import product
from pathlib import Path

import numpy as np
from squishierplanet.laws import simplex_from_q

from ldhybrid import laws, v1
from ldhybrid.interpolate import LDTables
from ldhybrid.paths import DATA

ROOT = Path(__file__).resolve().parents[2]
L = laws.get("hybrid4")
mu = np.concatenate([np.geomspace(1e-9, 1e-2, 1000), np.linspace(1e-2, 1, 3001)])
mu = np.unique(np.r_[0.0, mu])
xi = mu * mu
T = L.T(mu)                     # (n_mu, 4)
D = L.dT(xi)                    # dT/dxi; dI/dmu = -2 mu (D @ w)
out = []


def check(W, label, tol=1e-12, chunk=500):
    W = np.atleast_2d(np.asarray(W, float))
    neg, nonmono, minI, minD = [], [], np.inf, np.inf
    for i in range(0, len(W), chunk):
        w = W[i:i + chunk]
        I = 1.0 - w @ T.T
        dIdxi = -(w @ D.T)
        neg.append((I < -tol).any(axis=1))
        nonmono.append((dIdxi < -tol).any(axis=1))
        minI, minD = min(minI, I.min()), min(minD, dIdxi.min())
    neg, nonmono = np.concatenate(neg), np.concatenate(nonmono)
    out.append(f"{label}: {len(W)} sets; violate (I) {neg.sum()}, violate (II) {nonmono.sum()}; "
               f"min I = {minI:.3g}, min dI/dxi = {minD:.3g}")
    return neg | nonmono


# 1. the shapes themselves
out.append(f"poles: {v1.HYBRID4_EPS}")
out.append(f"T_j(0) = {T[0].round(12)}, T_j(1) = {T[-1].round(12)}; "
           f"max dT_j/dxi over [0, 1] = {D.max():.3g} (must be <= 0)")
# 2. simplex: vertices, edge midpoints, random interior and boundary points
V = np.vstack([np.zeros(4), np.eye(4)])
check(V, "simplex vertices")
mids = np.array([(a + b) / 2 for a, b in product(V, V)])
check(mids, "edge midpoints")
rng = np.random.default_rng(0)
check(rng.dirichlet(np.ones(5), 100000)[:, :4], "random interior (flat Dirichlet)")
B = rng.dirichlet(np.ones(4), 100000)                     # sum w = 1 face (I(0) = 0)
check(B, "random points on the face sum w = 1")
Z = rng.dirichlet(np.ones(5), 100000)[:, :4]
Z[np.arange(len(Z)), rng.integers(0, 4, len(Z))] = 0.0    # faces w_j = 0
check(Z, "random points on the faces w_j = 0")
# 3. q-map: q in [0,1]^4 -> simplex
Q = np.vstack([rng.random((50000, 4)), np.array(list(product([0.0, 1.0], repeat=4)))])
Wq = np.array([simplex_from_q(q) for q in Q])
out.append(f"q-map: min w = {Wq.min():.3g}, max sum w = {Wq.sum(1).max():.15f}")
check(Wq, "q-map images of q in [0,1]^4 (incl. corners)")
# 4. every tabulated row
t = LDTables(DATA / "release" / "cds")
for grid in ("phoenix", "mps1", "mps2"):
    d = t.table(grid, "hybrid4")
    W = d["coeffs"]
    out.append(f"{grid} table: {len(W)} rows; min w = {W.min():.3g}, max sum w = {W.sum(1).max():.6f}, "
               f"rows with sum w = 1 (zero limb intensity): {(W.sum(1) > 1 - 1e-6).sum()}")
    bad = check(W, f"{grid} tabulated hybrid4 (r and mu methods)", tol=1e-9)
# 5. corrected sets: w(mps2) + Δw for every MURaM state, every grid point and band
d = t.table("mps2", "hybrid4")
m = d["method"] == "r"
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for s in ("ssd", "B100G", "B200G", "B300G"):
        Wc = []
        for band in np.unique(d["band"][m]):
            k = m & (d["band"] == band)
            try:
                dw = t.delta_w("hybrid4", band, d["teff"][k], s)
            except ValueError:
                continue
            Wc.append(d["coeffs"][k] + dw)
        Wc = np.vstack(Wc)
        out.append(f"mps2 + {s}: outside simplex {((Wc < -1e-9).any(1) | (Wc.sum(1) > 1 + 1e-9)).sum()} of {len(Wc)}")
        check(Wc, f"mps2 + {s} corrected sets", tol=1e-9)
# 6. interpolated coefficients (random points, all three grids, TESS and Kepler)
for grid, (lo, hi) in {"phoenix": ((2300, 3.5, -2), (7400, 5.5, 1)),
                       "mps2": ((3500, 3.5, -2), (7500, 5.0, 1))}.items():
    Wi = []
    for _ in range(3000):
        x = [rng.uniform(a, b) for a, b in zip(lo, hi)]
        band = ("TESS", "Kepler", "JWST_NIRSpec_PRISM")[rng.integers(3)]
        try:
            Wi.append(list(t.coefficients("hybrid4", *x, band, grid=grid)["coeffs"].values()))
        except ValueError:
            pass
    check(Wi, f"{grid} interpolated (random points)", tol=1e-9)
# 7. standalone files
sys.path.insert(0, str(ROOT / "standalone"))
import ldhybrid_hybrid4_mps2_tess as S1, ldhybrid_hybrid4_phoenix_kepler as S2  # noqa: E402
Ws = []
for _ in range(2000):
    Ws.append(S1.coefficients(rng.uniform(3500, 7500), rng.uniform(3.5, 5.0), rng.uniform(-2, 1),
                              correction=("B100G" if rng.random() < 0.5 else None))["w"])
    try:
        Ws.append(S2.coefficients(rng.uniform(2300, 7400), rng.uniform(3.5, 5.5), rng.uniform(-2, 1))["w"])
    except ValueError:
        pass
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    check(Ws, "standalone files (random points)", tol=1e-9)

txt = "# hybrid4 vs the Kipping (2013) conditions\n\n" + "\n".join(f"* {o}" for o in out) + "\n"
(Path(__file__).parent / "h4_kipping_check.md").write_text(txt)
print(txt)
