"""User-facing API.

    import ldhybrid as ld
    r = ld.coefficients("hybrid5", teff=5600, logg=4.4, mh=0.1, band="TESS")
    r["w"], r["q"], r["physical"]
    # recommended for FGK priors (P4): MPS-ATLAS Set 2 + MURaM 100 G
    r = ld.coefficients("hybrid5", 5600, 4.4, 0.1, "TESS", grid="mps2", correction="B100G")
    # any throughput or spectroscopic bin (PHOENIX coefficient spectra)
    r = ld.coefficients_custom("hybrid4", 4500, 4.6, 0.0, lam_angstrom, throughput)
    law = ld.to_law("hybrid5", r["w"])      # squishierplanet.Law, for exact light curves

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

GRIDS = ("phoenix", "mps1", "mps2", "kurucz")


def coefficients(law: str, teff: float, logg: float, mh: float, band: str,
                 grid: str = "mps2", method: str = "r", correction: str | None = None) -> dict:
    """Tabulated coefficients of `law` for a star in a registered passband.

    correction: None, or a MURaM magnetic state ("ssd", "B100G", "B200G",
    "B300G"; hybrid laws only), added as Δw interpolated in Teff.
    Returns dict(w, q (hybrids), physical, grid, method, correction).
    """
    if grid not in GRIDS:
        raise ValueError(f"grid must be one of {GRIDS}")
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
