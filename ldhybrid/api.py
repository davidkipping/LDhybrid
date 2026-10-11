"""User-facing API.

    import ldhybrid as ld
    r = ld.coefficients("hybrid5", teff=5600, logg=4.4, mh=0.1, band="TESS")
    r["w"], r["q"], r["physical"]
    # recommended for FGK priors (P4): MPS-ATLAS Set 2 + MURaM 100 G
    r = ld.coefficients("hybrid5", 5600, 4.4, 0.1, "TESS", grid="mps2", correction="B100G")
    # any throughput or spectroscopic bin: PHOENIX coefficient spectra (hybrid
    # laws), or MPS-ATLAS models band-integrated on the fly (any law; needs the
    # MPS-ATLAS archive, ~9.5 GB per set, in data/raw/mpsatlas)
    r = ld.coefficients_custom("hybrid4", 4500, 4.6, 0.0, lam_angstrom, throughput)
    r = ld.coefficients_custom("claret4", 5600, 4.4, 0.1, lam_angstrom, throughput, grid="mps2")
    law = ld.to_law("hybrid5", r["w"])      # squishierplanet.Law, for exact light curves
    F = ld.light_curve("hybrid5", r["w"], z, p=0.1)   # numerical, any law
    # interpolate the band profiles and refit, instead of interpolating coefficients
    r = ld.coefficients("hybrid5", 5600, 4.4, 0.1, "TESS", refit=True)
    # prior for uncertain stellar parameters (+ grid-to-grid systematic), in
    # the recommended correction state, then the empirical Kepler offsets
    # measured in the SAME state (h'1 and h'2 jointly)
    p = ld.prior("hybrid5", 5600, 100, 4.4, 0.1, 0.1, 0.1, "Kepler", correction="B100G",
                 systematic=True)
    m, C = ld.apply_offsets(p, "Kepler")

Coefficients are interpolated linearly in (Teff, log g, [M/H]) within one grid
(never across grids), which keeps them inside the convex physical region.
Priors, offsets and light curves are those of the standalone
ldhybrid.interpolate, run on the CDS release (data/release/cds).
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


def condition(mean, cov, law, observable, delta, sigma, rho=None):
    """Gaussian update of a weight prior by empirical offsets in one or
    several observables (interpolate.condition)."""
    from . import interpolate
    return interpolate.condition(mean, cov, law, observable, delta, sigma, rho)


def light_curve(law, coeffs, z, p):
    """Circular-planet transit light curve of any tabulated law, by the
    numerical quadrature of the reference light curves (interpolate.light_curve)."""
    from . import interpolate
    return interpolate.light_curve(law, coeffs, z, p)


@lru_cache(maxsize=1)
def _release():
    from . import interpolate
    return interpolate.LDTables(DATA / "release" / "cds", profiles=DATA / "products" / "profiles")


def apply_offsets(prior, regime, statistics=None, route="hybrid5"):
    """Condition a prior from prior() on the tabulated empirical offsets of
    the same grid and correction state (interpolate.LDTables.apply_offsets)."""
    return _release().apply_offsets(prior, regime, statistics, route)


GRIDS = ("phoenix", "mps1", "mps2")


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
          systematic: bool = False, seed: int = 0, correction: str | None = None,
          m_dwarf_floor: bool = True, cov_floor=None) -> dict:
    """Monte-Carlo prior on a law's weights for Gaussian uncertainties in
    (Teff, log g, [M/H]), optionally in a MURaM correction state, with the
    grid-to-grid systematic (and its M-dwarf floor) and a variance floor;
    see interpolate.LDTables.prior. Returns dict(mean, cov, n_used, dropped,
    grids, correction, ...)."""
    return _release().prior(law, teff, sig_teff, logg, sig_logg, mh, sig_mh, band, grid,
                            method, n, systematic, seed, correction, m_dwarf_floor, cov_floor)


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


def coefficients_custom(law: str, teff: float, logg: float, mh: float, lam, throughput,
                        grid: str = "phoenix") -> dict:
    """Coefficients for any photon-counting throughput (λ in Å).

    grid="phoenix": hybrid laws from the PHOENIX coefficient spectra (R = 3000;
    one stellar edge per model; unconstrained per-bin fits summed), interpolated
    trilinearly between grid models (missing corners dropped and the weights
    renormalised).
    grid="mps1"/"mps2": any law, with the production method. The MPS-ATLAS
    models at the corners of the Delaunay simplex containing the star are read
    from the MPS-ATLAS archive (Kostogryz et al. 2022, Edmond doi
    10.17617/3.NJ56TR; set1.zip and set2.zip, ~9.5 GB each, expected in
    data/raw/mpsatlas; scripts/fetch_archives.py mpsatlas), band-integrated over
    the throughput, combined with the barycentric weights (all MPS-ATLAS
    profiles share the same μ nodes) and fitted (uniform in r, exact
    quadrature, the law's production constraint). At a grid model this
    reproduces the tabulated coefficients of a registered band."""
    if grid in ("mps1", "mps2"):
        return _custom_mps(law, teff, logg, mh, np.asarray(lam, float),
                           np.asarray(throughput, float), grid)
    if grid != "phoenix":
        raise ValueError("grid must be phoenix, mps1 or mps2")
    if law not in spectra.SPECTRA_LAWS:
        raise ValueError(f"PHOENIX coefficient spectra exist for {spectra.SPECTRA_LAWS}; "
                         "use grid='mps1' or 'mps2' for the other laws")
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


@lru_cache(maxsize=2)
def _mps_triangulation(set_name):
    from scipy.spatial import Delaunay
    from .grids import mpsatlas
    cat = [(t, g, m) for m, t, g in mpsatlas.catalogue(set_name)
           if 3500 <= t <= 7500 and 3.5 <= g <= 5.0 and -2.0 - 1e-9 <= m <= 1.0 + 1e-9]
    X = np.array(cat, float) / predict.SCALE
    return Delaunay(X), cat


def _custom_mps(law, teff, logg, mh, lam, S, grid):
    from . import v1
    from .grids import mpsatlas
    from .products import constraint_for
    from .profile import Profile
    set_name = "set" + grid[-1]
    tri, cat = _mps_triangulation(set_name)
    x = np.array([teff, logg, mh]) / predict.SCALE
    k = int(tri.find_simplex(x))
    if k < 0:
        raise ValueError(f"({teff}, {logg}, {mh}) is outside the {grid} grid")
    T = tri.transform[k]
    b = T[:3] @ (x - T[3])
    lam_w = np.r_[b, 1 - b.sum()]
    mu, y, corners = None, 0.0, []
    for wgt, c in zip(lam_w, tri.simplices[k]):
        if wgt <= 1e-12:
            continue
        t_, g_, m_ = cat[c]
        spec = mpsatlas.load(set_name, m_, t_, g_)
        if spec.coverage(lam, S) < 0.999:
            raise ValueError("the throughput extends beyond the MPS-ATLAS wavelength range")
        I = spec.band(lam, S)
        mu = spec.mu                                  # the same 24 nodes for every model
        y = y + wgt * I / I[np.argmax(mu)]            # normalised at the disc centre
        corners.append(((t_, g_, m_), float(wgt)))
    pr = Profile(mu, y)
    if law == "power2":
        w = fit.fit_power2(pr, v1.MEASURE, v1.RESAMPLE).w
        phys = bool(0 <= w[0] <= 1 and w[1] > 0)
    else:
        w = fit.fit_linear(laws.get(law), pr, v1.MEASURE, v1.RESAMPLE, constraint_for(law),
                           v1.FCM).w
        phys = fit.is_physical(laws.get(law), w)
    out = dict(law=law, w=w, grid=grid, physical=phys, corners=corners)
    if law.startswith("hybrid"):
        out["q"] = qspace.to_q(law, w)
    return out


def to_law(law: str, w):
    """squishierplanet.Law for closed-form light curves (hybrid laws)."""
    return laws.get(law).to_law(np.asarray(w, float))
