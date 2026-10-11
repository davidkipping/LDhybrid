"""Shared set-up for the revision studies (studies/rev).

P2S: the P2 subset (34 PHOENIX v3.0 + 16 MPS-ATLAS Set 2 models) × 8 bands.
T3:  the 15 cases of the paper's Table 3 (studies/lc/lc_section.py).
Fits use the production (v1) settings; light-curve errors use the production
grid products.PROD_LC (p = 0.1, b ∈ {0, 0.3, 0.6, 0.9}) unless stated.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
import importlib.util as _iu  # noqa: E402
_spec = _iu.spec_from_file_location("p2common", ROOT / "studies" / "p2" / "common.py")
p2common = _iu.module_from_spec(_spec)
_spec.loader.exec_module(p2common)

from ldhybrid import bands, fit, laws, v1  # noqa: E402
from ldhybrid import lightcurve as lc  # noqa: E402
from ldhybrid.products import PROD_LC, constraint_for  # noqa: E402
from ldhybrid.profile import band_profile  # noqa: E402

HERE = Path(__file__).parent
LAWS6 = ["hybrid5", "hybrid4", "claret4", "power2", "hybrid2", "quadratic"]
P2S_BANDS = p2common.BANDS
T3_STARS = [("late M", (2800, 5.0)), ("mid M", (3300, 5.0)), ("early M", (3800, 5.0)),
            ("K dwarf", (4600, 4.5)), ("Sun analogue", (5800, 4.5))]
T3_BANDS = [("TESS", "TESS"), ("HST_WFC3_G141", "WFC3 G141"),
            ("JWST_NIRSpec_PRISM", "NIRSpec PRISM")]


def p2s_models():
    return p2common.models()


def load(kind, key):
    return p2common.load(kind, key)


def load_t3(teff, logg):
    from ldhybrid.grids import phoenix
    return phoenix.load(teff, logg, "-0.0")


def profile(spec, band, limb=v1.LIMB):
    lam, S = bands.response(band)
    if spec.coverage(lam, S) < 0.999:
        return None
    if spec.meta.get("geometry") == "spherical":
        return band_profile(spec, lam, S, limb=limb)
    return band_profile(spec, lam, S)


def poly_law(n):
    """Even polynomial in μ² without poles: T_k = 1 − μ^{2k}, k = 1..n."""
    from squishierplanet import laws as sl
    shapes = [sl.Law([1.0] + [0.0] * (k - 1) + [-1.0]) for k in range(1, n + 1)]
    return laws._from_shapes(f"poly{n}", [f"w{k}:1-mu{2 * k}" for k in range(1, n + 1)], shapes)


def get_law(name):
    if name.startswith("poly"):
        return poly_law(int(name[4:]))
    return laws.get(name)


def fit_law(name, pr, measure=v1.MEASURE, constraint=None):
    """(I_of_mu, coefficients, Fit) with production settings."""
    if name == "power2":
        f = fit.fit_power2(pr, measure, v1.RESAMPLE)
        c, a = f.w
        return (lambda mu, c=c, a=a: laws.power2(c, a, mu)), f.w, f
    L = get_law(name)
    con = constraint or (constraint_for(name) if not name.startswith("poly") else "physical")
    f = fit.fit_linear(L, pr, measure, v1.RESAMPLE, con, v1.FCM)
    return (lambda mu, L=L, w=f.w: L.intensity(w, mu)), f.w, f


_GRIDS = {}


def lcgrid(**kw):
    key = tuple(sorted((k, tuple(v) if isinstance(v, (list, tuple)) else v) for k, v in kw.items()))
    if key not in _GRIDS:
        _GRIDS[key] = lc.LCGrid(**(kw or PROD_LC))
    return _GRIDS[key]


def score(model, truth_cols, g=None):
    """Worst |ΔF| (ppm) of a profile I(μ) against precomputed truth columns."""
    g = g or lcgrid()
    return g.errors(*g.columns(model), *truth_cols)["max_ppm"]


def pct(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return np.median(x), np.percentile(x, 99), np.max(x)


def fmt3(x):
    m, p, mx = pct(x)
    f = lambda v: f"{v:.2g}" if v < 10 else f"{v:.0f}" if v >= 100 else f"{v:.1f}" if v < 10 else f"{v:.0f}"
    return f"{f(m)} / {f(p)} / {f(mx)}"


def run_parallel(task, jobs, workers=5):
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(workers) as ex:
        out = []
        for i, r in enumerate(ex.map(task, jobs, chunksize=1), 1):
            out.append(r)
            if i % 20 == 0:
                print(f"  {i}/{len(jobs)}", flush=True)
    return out
