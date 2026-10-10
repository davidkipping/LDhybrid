"""MPS-ATLAS library (Witzke et al. 2021; Kostogryz et al. 2022), Set 1 and
Set 2, read straight from the Edmond zip archives.

Plane-parallel, 24 μ (1 … 0.01), 1221 ODF wavelengths from 9 nm to 160 µm.
Intensities are I_ν [erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹]. The photon density per unit
wavelength is I_λ λ / (hc) = I_ν / (h λ), so we bin I_ν / λ in photon units
(constant factors drop out of normalised profiles).

  Set 1: Grevesse & Sauval (1998) abundances, mixing length 1.25.
  Set 2: Asplund et al. (2009) abundances, mixing length from Viani et al. (2018).
"""

from __future__ import annotations

import zipfile
from functools import lru_cache

import numpy as np

from ..paths import RAW
from ..specint import SpecInt, bin_photons


@lru_cache(maxsize=2)
def _zip(set_name: str) -> zipfile.ZipFile:
    return zipfile.ZipFile(RAW / "mpsatlas" / f"{set_name}.zip")


@lru_cache(maxsize=None)
def _members(set_name: str) -> dict:
    """(mh, teff, logg) -> archive member name (directory spellings vary,
    e.g. MH-1 and MH-0.15)."""
    out = {}
    for n in _zip(set_name).namelist():
        if n.endswith("mpsa_intensity_spectra.dat"):
            _, mh, teff, logg, _ = n.split("/")
            key = (round(float(mh[2:]), 2) + 0.0, int(teff[4:]), round(float(logg[4:]), 1))
            out[key] = n
    return out


def member(set_name: str, mh: float, teff: int, logg: float) -> str:
    return _members(set_name)[(round(mh, 2) + 0.0, int(teff), round(logg, 1))]


def catalogue(set_name: str) -> list[tuple[float, int, float]]:
    """(mh, teff, logg) of every model in the archive."""
    return sorted(_members(set_name))


def load(set_name: str, mh: float, teff: int, logg: float, refine: int = 8) -> SpecInt:
    name = member(set_name, mh, teff, logg)
    lines = _zip(set_name).read(name).decode().splitlines()
    head = lines[0].split()
    mu = np.array(lines[1].split()[2:], float)
    d = np.array([l.split() for l in lines[2:] if l.strip()], float)
    lam = d[:, 0] * 10.0                                 # nm -> Å
    I = (d[:, 1:] / lam[:, None]).T                      # photon density ∝ I_ν/λ
    o = np.argsort(mu)
    mu, I = mu[o], I[o]
    edges = (lam[:-1, None] + np.diff(lam)[:, None] * np.arange(refine) / refine).ravel()
    edges = np.append(edges, lam[-1])
    P = bin_photons(lam, I, edges, photon_units=True)
    return SpecInt(edges, mu, P, dict(grid=f"MPS-ATLAS-{set_name}", geometry="plane-parallel",
                                     mh=mh, teff=float(teff), logg=float(logg),
                                     header=" ".join(head), source=f"{set_name}.zip:{name}"))
