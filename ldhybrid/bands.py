"""Passband registry: name -> (λ in Å, photon-counting response S).

Every curve is returned as a *photon* response, ready for
SpecInt.band(lam, S), which forms ∫ I_λ λ S dλ:
  * HST and JWST mode throughputs come from the ExoTiC-LD v3.1.2 bundle
    (instrument throughputs, i.e. photon fractions);
  * photometric bands come from the SVO Filter Profile Service. Each entry
    states whether its curve is a photon-counting response or an energy
    response; an energy response R_E becomes the photon response R_E / λ,
    since ∫ F_λ R_E dλ = ∫ (F_λ λ)(R_E / λ) dλ.

    We set this per band rather than trusting SVO's DetectorType flag: SVO
    flags TESS and Spitzer/IRAC as energy counters (0), but their curves are
    the missions' photon-counting responses (the SVO TESS curve is identical to
    ExoTiC-LD's, and treating it as energy-type shifts TESS profiles by ~1%
    against Kostogryz et al. 2022). Only the Bessell Johnson–Cousins curves
    are energy-type;
  * Spitzer/IRAC uses the official SSC response files (see IRAC below).

Fetch the curves first with scripts/fetch_bands.py.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from functools import lru_cache

import numpy as np

from .paths import RAW

SVO_URL = "http://svo2.cab.inta-csic.es/theory/fps/fps.php?ID={}"

# name: SVO filter ID
SVO = {
    "Kepler": "Kepler/Kepler.K",
    "TESS": "TESS/TESS.Red",
    "CHEOPS": "CHEOPS/CHEOPS.band",
    "NGTS": "Paranal/NGTS.NGTS",
    "Johnson_B": "Generic/Johnson.B",
    "Johnson_V": "Generic/Johnson.V",
    "Cousins_R": "Generic/Cousins.R",
    "Cousins_I": "Generic/Cousins.I",
    "Sloan_g": "SLOAN/SDSS.gprime_filter",
    "Sloan_r": "SLOAN/SDSS.rprime_filter",
    "Sloan_i": "SLOAN/SDSS.iprime_filter",
    "Sloan_z": "SLOAN/SDSS.zprime_filter",
    "2MASS_J": "2MASS/2MASS.J",
    "2MASS_H": "2MASS/2MASS.H",
    "2MASS_Ks": "2MASS/2MASS.Ks",
}

# SVO curves that are energy-type responses (Bessell 1990); all others photon
ENERGY_TYPE = {"Johnson_B", "Johnson_V", "Cousins_R", "Cousins_I"}

# Spitzer/IRAC: official SSC spectral responses (electrons per photon).
# Most transit/eclipse photometry used the subarray, whose response differs
# markedly from the full array's (ch1 by factors 0.4–2.4 across the band);
# SVO's IRAC curves follow the full array.
IRAC_URL = "https://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/calibrationfiles/spectralresponse/{}"
IRAC = {
    "Spitzer_IRAC1": "080924ch1trans_sub.txt",
    "Spitzer_IRAC2": "080924ch2trans_sub.txt",
    "Spitzer_IRAC1_full": "201125ch1trans_full.txt",
    "Spitzer_IRAC2_full": "201125ch2trans_full.txt",
}

# name: ExoTiC-LD sensitivity-file stem
EXOTICLD = {
    "HST_WFC3_G102": "HST_WFC3_G102",
    "HST_WFC3_G141": "HST_WFC3_G141",
    "HST_WFC3_G280p1": "HST_WFC3_G280p1",
    "HST_WFC3_G280n1": "HST_WFC3_G280n1",
    "HST_STIS_G430L": "HST_STIS_G430L",
    "HST_STIS_G750L": "HST_STIS_G750L",
    "JWST_NIRISS_SOSS_o1": "JWST_NIRISS_SOSSo1",
    "JWST_NIRISS_SOSS_o2": "JWST_NIRISS_SOSSo2",
    "JWST_NIRSpec_PRISM": "JWST_NIRSpec_prism",
    "JWST_NIRSpec_G140H_F070LP": "JWST_NIRSpec_G140H-f070",
    "JWST_NIRSpec_G140H_F100LP": "JWST_NIRSpec_G140H-f100",
    "JWST_NIRSpec_G140M_F070LP": "JWST_NIRSpec_G140M-f070",
    "JWST_NIRSpec_G140M_F100LP": "JWST_NIRSpec_G140M-f100",
    "JWST_NIRSpec_G235H": "JWST_NIRSpec_G235H",
    "JWST_NIRSpec_G235M": "JWST_NIRSpec_G235M",
    "JWST_NIRSpec_G395H": "JWST_NIRSpec_G395H",
    "JWST_NIRSpec_G395M": "JWST_NIRSpec_G395M",
    "JWST_NIRCam_F322W2": "JWST_NIRCam_F322W2",
    "JWST_NIRCam_F444W": "JWST_NIRCam_F444",
    "JWST_MIRI_LRS": "JWST_MIRI_LRS",
}

BANDS = sorted(SVO) + sorted(IRAC) + sorted(EXOTICLD)


def svo_path(name):
    return RAW / "bands" / "svo" / (SVO[name].replace("/", "__") + ".xml")


def _parse_svo(path):
    root = ET.parse(path).getroot()
    ns = {"v": root.tag.split("}")[0].strip("{")} if root.tag.startswith("{") else {}
    pre = "v:" if ns else ""
    params = {p.get("name"): p.get("value") for p in root.iter(f"{{{ns['v']}}}PARAM" if ns else "PARAM")}
    rows = [[float(td.text) for td in tr] for tr in root.iter(f"{{{ns['v']}}}TR" if ns else "TR")]
    a = np.array(rows)
    return a[:, 0], a[:, 1], params


@lru_cache(maxsize=None)
def response(name: str) -> tuple[np.ndarray, np.ndarray]:
    """(λ [Å], photon response) for a registered band, sorted in λ, S >= 0."""
    if name in SVO:
        lam, S, params = _parse_svo(svo_path(name))
        if params.get("WavelengthUnit", "Angstrom") != "Angstrom":
            raise ValueError(f"{name}: unexpected unit {params.get('WavelengthUnit')}")
        if name in ENERGY_TYPE:
            S = S / lam
    elif name in IRAC:
        a = np.loadtxt(RAW / "bands" / "irac" / IRAC[name], comments="#")
        lam, S = a[:, 0] * 1e4, a[:, 1]                  # µm -> Å
    elif name in EXOTICLD:
        path = RAW / "bands" / "exoticld" / f"{EXOTICLD[name]}_throughput.csv"
        a = np.loadtxt(path, delimiter=",", comments="#", skiprows=0 if _has_hash(path) else 1)
        lam, S = a[:, 0], a[:, 1]
    else:
        raise KeyError(f"unknown band {name!r}; known: {BANDS}")
    o = np.argsort(lam)
    lam, S = lam[o], np.clip(S[o], 0.0, None)
    return lam, S


def _has_hash(path):
    with open(path) as f:
        return f.readline().lstrip().startswith("#")
