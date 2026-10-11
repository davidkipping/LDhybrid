"""User-facing API.

    import ldhybrid as ld
    r = ld.coefficients("hybrid5", teff=5600, logg=4.4, mh=0.1, band="TESS")
    r["w"], r["q"], r["physical"]
    # recommended for FGK priors (P4): MPS-ATLAS Set 2 + MURaM 100 G
    r = ld.coefficients("hybrid5", 5600, 4.4, 0.1, "TESS", grid="mps2", correction="B100G")
    # any throughput or spectroscopic bin (PHOENIX coefficient spectra)
    r = ld.coefficients_custom("hybrid4", 4500, 4.6, 0.0, lam_angstrom, throughput)
    law = ld.to_law("hybrid5", r["w"])      # squishierplanet.Law, for exact light curves
    # interpolate the band profiles and refit, instead of interpolating coefficients
    r = ld.coefficients("hybrid5", 5600, 4.4, 0.1, "TESS", refit=True)
    # prior for uncertain stellar parameters (+ grid-to-grid systematic)
    p = ld.prior("hybrid5", 5600, 100, 4.4, 0.1, 0.1, 0.1, "TESS", systematic=True)
    # apply an empirical offset (observed − model) in h'1 = I(2/3), with its scatter
    m, C = ld.condition(p["mean"], p["cov"], "hybrid5", "h1", +0.008, 0.006)

Coefficients are interpolated linearly in (Teff, log g, [M/H]) within one grid
(never across grids), which keeps them inside the convex physical region.
"""

from __future__ import annotations

import itertools
import json
from functools import lru_cache

import numpy as np

from . import fit, laws, predict, qspace, spectra
from .paths import DATA, REDUCED


def observables(law, w=None):
    """Exact linear functionals h1, h2, ell, delta of a linear law's weights
    (interpolate.observables)."""
    from . import interpolate
    return interpolate.observables(law, w)


def condition(mean, cov, law, observable, delta, sigma):
    """Gaussian update of a weight prior by an empirical offset in one
    observable (interpolate.condition)."""
    from . import interpolate
    return interpolate.condition(mean, cov, law, observable, delta, sigma)

GRIDS = ("phoenix", "mps1", "mps2", "kurucz")


def coefficients(law: str, teff: float, logg: float, mh: float, band: str,
                 grid: str = "mps2", method: str = "r", correction: str | None = None,
                 refit: bool = False) -> dict:
    """Tabulated coefficients of `law` for a star in a registered passband.

    correction: None, or a MURaM magnetic state ("ssd", "B100G", "B200G",
    "B300G"; hybrid laws only), added as Δw interpolated in Teff.
    refit: interpolate the band-integrated profiles (data/products/profiles)
    and refit with the production method, rather than interpolating the
    tabulated coefficients (r-weighted method only).
    Returns dict(w, q (hybrids), physical, grid, method, correction).
    """
    if grid not in GRIDS:
        raise ValueError(f"grid must be one of {GRIDS}")
    if refit:
        if method != "r":
            raise ValueError("refit reproduces the headline (r-weighted) method only")
        r = _refitter().refit(law, teff, logg, mh, band, grid)
        w = np.array(list(r["coeffs"].values()))
    else:
        w = predict.weights(grid, law, band, teff, logg, mh, method)
    if np.any(np.isnan(w)):
        raise ValueError(f"({teff}, {logg}, {mh}) is outside the {grid} grid for {band}")
    if correction:
        if not law.startswith("hybrid"):
            raise ValueError("corrections are tabulated for the hybrid laws only")
        dw = predict.magnetic(law, band, teff, correction)
        if dw is None:
            raise ValueError(f"no {correction} correction for {band}")
        w = w + dw
    out = dict(law=law, w=w, grid=grid, method=method, correction=correction)
    if law == "power2":
        out["physical"] = bool(0 <= w[0] <= 1 and w[1] > 0)
    else:
        out["physical"] = fit.is_physical(laws.get(law), w)
    if law.startswith("hybrid"):
        out["q"] = qspace.to_q(law, w)
    return out


@lru_cache(maxsize=1)
def _refitter():
    from . import interpolate
    return interpolate.LDTables(DATA / "products", cache=False,
                                profiles=DATA / "products" / "profiles")


def prior(law: str, teff: float, sig_teff: float, logg: float, sig_logg: float, mh: float,
          sig_mh: float, band: str, grid: str = "mps2", method: str = "r", n: int = 2000,
          systematic: bool = False, seed: int = 0) -> dict:
    """Monte-Carlo prior on a law's weights for Gaussian uncertainties in
    (Teff, log g, [M/H]): mean and covariance of the interpolated weights over
    n draws (draws outside the grid are dropped). systematic=True adds the
    covariance across the grids (phoenix, mps1, mps2) that cover the star at
    the nominal parameters. Returns dict(mean, cov, n_used, grids)."""
    rng = np.random.default_rng(seed)
    f, _, _ = predict._interp(grid, law, method, band)
    X = np.column_stack([rng.normal(teff, sig_teff, n), rng.normal(logg, sig_logg, n),
                         rng.normal(mh, sig_mh, n)]) / predict.SCALE
    Y = np.asarray(f(X))
    Y = Y[np.all(np.isfinite(Y), axis=1)]
    if len(Y) == 0:
        raise ValueError("no draws inside the grid")
    cov = np.atleast_2d(np.cov(Y.T, ddof=1))
    grids = [grid]
    if systematic:
        W = []
        for g in ("phoenix", "mps1", "mps2"):
            try:
                w = predict.weights(g, law, band, teff, logg, mh, method)
            except ValueError:
                continue
            if np.all(np.isfinite(w)):
                W.append(w)
                if g != grid:
                    grids.append(g)
        if len(W) > 1:
            cov = cov + np.atleast_2d(np.cov(np.array(W).T, ddof=1))
    return dict(law=law, mean=Y.mean(0), cov=cov, n_used=len(Y), grids=grids, band=band,
                method=method)


@lru_cache(maxsize=64)
def _spectrum(teff, logg, mhs):
    from .grids.phoenix import file_name
    z = np.load(DATA / "products" / "spectra" / "phoenix" /
                file_name(teff, logg, mhs).replace("SpecInt.fits", "coeffspec.npz"))
    return spectra.CoeffSpectrum(z["lam_edges"], z["I1"].astype(float),
                                 {k: z[f"a_{k}"].astype(float) for k in spectra.SPECTRA_LAWS},
                                 float(z["r_c"]), {})


@lru_cache(maxsize=1)
def _phoenix_catalogue():
    return [(e["teff"], e["logg"], e["mhs"]) for e in
            json.load(open(REDUCED / "phoenix" / "catalogue.json"))]


def _corners(values, x):
    v = np.array(sorted(set(values)))
    if x <= v[0]:
        return [(v[0], 1.0)]
    if x >= v[-1]:
        return [(v[-1], 1.0)]
    i = np.searchsorted(v, x) - 1
    t = (x - v[i]) / (v[i + 1] - v[i])
    return [(v[i], 1 - t), (v[i + 1], t)]


def coefficients_custom(law: str, teff: float, logg: float, mh: float, lam, throughput) -> dict:
    """Hybrid coefficients for any photon-counting throughput (λ in Å), from
    the PHOENIX coefficient spectra (R = 3000; one stellar edge per model),
    interpolated trilinearly between grid models (missing corners dropped and
    the weights renormalised)."""
    if law not in spectra.SPECTRA_LAWS:
        raise ValueError(f"coefficient spectra exist for {spectra.SPECTRA_LAWS}")
    cat = _phoenix_catalogue()
    have = {(t, g, float(m)): m for t, g, m in cat}
    terms = []
    for (t, a), (g, b), (m, c) in itertools.product(
            _corners([k[0] for k in cat], teff), _corners([k[1] for k in cat], logg),
            _corners([float(k[2]) for k in cat], mh)):
        key = (int(t), float(g), float(m))
        if key in have and a * b * c > 0:
            w, _ = _spectrum(int(t), float(g), have[key]).band_coefficients(lam, throughput, law)
            terms.append((a * b * c, w))
    if not terms:
        raise ValueError("no PHOENIX models bracket these parameters")
    s = sum(f for f, _ in terms)
    w = sum(f * x for f, x in terms) / s
    return dict(law=law, w=w, grid="phoenix-spectra", physical=fit.is_physical(laws.get(law), w),
                q=qspace.to_q(law, w))


def to_law(law: str, w):
    """squishierplanet.Law for closed-form light curves (hybrid laws)."""
    return laws.get(law).to_law(np.asarray(w, float))
