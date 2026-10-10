"""Grids in the ExoTiC-LD data format (Grant & Wakeford 2024; data bundle
v3.1.2): 1D Kurucz ATLAS9 (17 μ) and 3D Stagger (Magic et al. 2015; 10 μ).

Each `<model>_spectra.dat` holds a title line, a line of μ values, then rows of
λ [Å] followed by I(λ, μ) in photon units [photons s⁻¹ cm⁻² Å⁻¹ sr⁻¹].
Both grids are plane-parallel: no limb cut is applied.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..paths import RAW
from ..specint import SpecInt, bin_photons

ROOT = RAW / "exoticld" / "unzipped"


def path(grid: str, mh: float, teff: int, logg: float) -> Path:
    mh = 0.0 if mh == 0 else mh
    return ROOT / grid / f"MH{round(mh, 2)}" / f"teff{int(teff)}" / f"logg{round(logg, 1)}" / f"{grid}_spectra.dat"


def load(grid: str, mh: float, teff: int, logg: float, refine: int = 8) -> SpecInt:
    """Native-resolution SpecInt. Each interval between consecutive samples is
    split into `refine` bins, so that the bin-mean throughput approximation in
    SpecInt.band stays accurate on coarse (Kurucz ODF) sampling; the photon
    integrals are exact for the piecewise-linear spectrum."""
    p = path(grid, mh, teff, logg)
    mu = np.loadtxt(p, skiprows=1, max_rows=1)
    d = np.loadtxt(p, skiprows=2)
    lam, I = d[:, 0], d[:, 1:].T                     # (nmu, nlam)
    o = np.argsort(mu)
    mu, I = mu[o], I[o]
    mu[-1] = 1.0                                     # stored as float32 1.0
    edges = (lam[:-1, None] + np.diff(lam)[:, None] * np.arange(refine) / refine).ravel()
    edges = np.append(edges, lam[-1])
    P = bin_photons(lam, I, edges, photon_units=True)
    return SpecInt(edges, mu, P, dict(grid=grid, geometry="plane-parallel", mh=mh,
                                    teff=float(teff), logg=float(logg), source_file=str(p)))
