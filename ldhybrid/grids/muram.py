"""MURaM 3D-MHD magnetic limb-darkening library (Kostogryz et al. 2026,
arXiv:2606.21912; spectra synthesised with MPS-ATLAS).

Nine main-sequence stars F3–M4 at [M/H] = 0 (G2 also at −1.0 and +0.5); five
magnetisation states (hydro, ssd, B100G, B200G, B300G); 10 μ from 0.1 to 1;
978 wavelengths 200 nm–10 µm (ODF, R ≈ 250–400); I_λ in erg s⁻¹ cm⁻² sr⁻¹ Å⁻¹,
vacuum wavelengths. There are no intensities below μ = 0.1.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np

from ..paths import RAW
from ..specint import SpecInt, bin_photons

ROOT = RAW / "muram"
STATES = ("hydro", "ssd", "B100G", "B200G", "B300G")


def stars() -> list[str]:
    return sorted(p.stem for p in ROOT.glob("*.h5"))


def params(star: str) -> dict:
    with h5py.File(ROOT / f"{star}.h5") as h:
        g = h[list(h.keys())[0]]
        states = [s.decode() for s in g["magnetizations"][()]]
        return dict(star=star, logg=float(g.attrs["logg"]), mh=float(g.attrs["MH"]),
                    teff={s: float(t) for s, t in zip(states, g["teff"][()])})


def load(star: str, state: str, refine: int = 8) -> SpecInt:
    with h5py.File(ROOT / f"{star}.h5") as h:
        g = h[list(h.keys())[0]]
        states = [s.decode() for s in g["magnetizations"][()]]
        k = states.index(state)
        lam = g["wavelengths"][()] * 10.0                     # nm -> Å (vacuum)
        mu = g["mu"][()]
        I = g["spectra"][k].T                                 # (nmu, nlam), I_λ
        teff = float(g["teff"][k])
        logg, mh = float(g.attrs["logg"]), float(g.attrs["MH"])
    o = np.argsort(mu)
    mu, I = mu[o], I[o]
    edges = (lam[:-1, None] + np.diff(lam)[:, None] * np.arange(refine) / refine).ravel()
    edges = np.append(edges, lam[-1])
    P = bin_photons(lam, I, edges)                            # energy units -> × λ
    return SpecInt(edges, mu, P, dict(grid="MURaM", geometry="plane-parallel", star=star,
                                     state=state, teff=teff, logg=logg, mh=mh))
