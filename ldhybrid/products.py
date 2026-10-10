"""P3 production: coefficients and diagnostics for one model atmosphere over
all bands, laws and methods, as arrays (assembled into tables elsewhere).

Light-curve diagnostic in production: the worst |ΔF| (ppm) against the exact
light curve of the tabulated profile on a reduced grid, p = 0.1,
b ∈ {0, 0.3, 0.6, 0.9}, 60 points per half-chord (PROD_LC). The P2 studies used
the full grid; this one is ~12× cheaper and tracks it closely.
"""

from __future__ import annotations

import numpy as np

from . import bands as B
from . import fit, laws, v1
from . import lightcurve as lc
from .profile import band_profile

LAWS = ("hybrid2", "hybrid4", "hybrid5", "quadratic", "claret4", "power2")
METHODS = ("r", "mu")                     # headline (v1.MEASURE), secondary
PROD_LC = dict(p_grid=(0.1,), b_grid=(0.0, 0.3, 0.6, 0.9), n=60)
NPAR = {"hybrid2": 2, "hybrid4": 4, "hybrid5": 5, "quadratic": 2, "claret4": 4, "power2": 2}

_GRID = None
_COLS = {}


def _grid():
    global _GRID
    if _GRID is None:
        _GRID = lc.LCGrid(**PROD_LC)
    return _GRID


def _cols(name):
    if name not in _COLS:
        _COLS[name] = _grid().linear_columns(laws.get(name))
    return _COLS[name]


def constraint_for(name):
    return v1.CONSTRAINT.get(name, "physical")


def model_products(spec, band_names, limb: str = v1.LIMB):
    """Arrays over bands for one model. Bands the model does not cover get NaN."""
    nb = len(band_names)
    out = {"r_c": np.full(nb, np.nan), "mu_min": np.full(nb, np.nan),
           "covered": np.zeros(nb, bool)}
    for name in LAWS:
        for m in METHODS:
            k = f"{name}|{m}"
            out[k + "|w"] = np.full((nb, NPAR[name]), np.nan)
            for f in ("rms_I", "max_abs_I", "flux_err", "lc_ppm"):
                out[k + "|" + f] = np.full(nb, np.nan, np.float32)
            out[k + "|physical"] = np.zeros(nb, bool)
            out[k + "|bound"] = np.zeros(nb, bool)
    g = _grid()
    spherical = spec.meta.get("geometry") == "spherical"
    for i, b in enumerate(band_names):
        lam, S = B.response(b)
        if spec.coverage(lam, S) < 0.999:
            continue
        pr = band_profile(spec, lam, S, limb=limb) if spherical else band_profile(spec, lam, S)
        out["covered"][i] = True
        out["r_c"][i] = pr.meta.get("r_c", 1.0)
        out["mu_min"][i] = pr.mu[0]
        truth = g.columns(pr)
        for name in LAWS:
            for m in METHODS:
                k = f"{name}|{m}"
                if name == "power2":
                    f = fit.fit_power2(pr, m, v1.RESAMPLE)
                    c, a = f.w
                    model = g.columns(lambda mu: laws.power2(c, a, mu))
                    phys = bool(0 <= c <= 1 and a > 0)
                    bound = False
                else:
                    f = fit.fit_linear(laws.get(name), pr, m, v1.RESAMPLE, constraint_for(name), v1.FCM)
                    model = g.combine(*_cols(name), f.w)
                    phys = bool(f.extra["physical"])
                    bound = bool(f.active)
                out[k + "|w"][i] = f.w
                out[k + "|rms_I"][i] = f.rms
                out[k + "|max_abs_I"][i] = f.max_abs
                out[k + "|flux_err"][i] = f.flux_err
                out[k + "|lc_ppm"][i] = g.errors(*model, *truth)["max_ppm"]
                out[k + "|physical"][i] = phys
                out[k + "|bound"][i] = bound
    return out
