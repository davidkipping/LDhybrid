"""PHOENIX-ACES-AGSS-COND SpecInt v3.0 (Husser et al. 2013; v3.0 extends the
2011 grid to 5.5 µm).

Raw files: 78 μ nodes (spherical models, μ measured from the outermost model
layer), 1 Å sampling over 500 Å – 5.5 µm, float64, 34 MB each. We reduce each
file on download to photon integrals in log bins (default R = 3000 over
0.1–5.5 µm, float32, ~3.5 MB) and keep the raw file only for the P2 subset.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
from astropy.io import fits

from ..download import fetch, list_ftp, record, sha256
from ..paths import RAW, REDUCED
from ..specint import SpecInt, bin_photons, log_edges

BASE = "ftp://ftp.gwdg.de/pub/misc/phoenix/v3.0/SpecIntFITS"
METALLICITIES = ["-2.0", "-1.5", "-1.0", "-0.5", "-0.0", "+0.5", "+1.0"]
TEFF_RANGE = (2300, 7500)
LOGG_RANGE = (3.5, 5.5)

LAM_MIN, LAM_MAX, R_BIN = 1000.0, 54990.0, 3000.0

_NAME = re.compile(r"lte(\d{5})([+-]\d\.\d\d)([+-]\d\.\d)(?:\.Alpha=([+-]\d\.\d\d))?"
                   r"\.PHOENIX-ACES-AGSS-COND-2011-SpecInt\.fits$")


def parse_name(name: str) -> dict | None:
    m = _NAME.match(name)
    if not m:
        return None
    return dict(teff=int(m.group(1)), logg=-float(m.group(2)), mh=float(m.group(3)),
                alpha=float(m.group(4) or 0.0))


def file_name(teff: int, logg: float, mh: str | float) -> str:
    mh = mh if isinstance(mh, str) else ("-0.0" if mh == 0 else f"{mh:+.1f}")
    g = "-0.00" if logg == 0 else f"{-logg:+.2f}"
    return f"lte{teff:05d}{g}{mh}.PHOENIX-ACES-AGSS-COND-2011-SpecInt.fits"


def zdir(mh: str) -> str:
    return f"Z{mh}"


# The P2 method-study subset (raw files kept): an M-to-F dwarf sequence at
# solar metallicity, plus metallicity and gravity offsets.
P2_SUBSET = (
    [(t, 5.0 if t <= 3800 else 4.5, "-0.0")
     for t in (2600, 3000, 3400, 3800, 4200, 4600, 5000, 5400, 5800, 6200, 6600, 7000, 7400)]
    + [(t, 5.0 if t <= 3800 else 4.5, mh) for t in (3400, 4600, 5800, 6600) for mh in ("-1.0", "+0.5")]
    + [(t, g, "-0.0") for t in (4600, 5800, 6600) for g in (3.5, 4.0, 5.0)]
    + [(t, g, "-0.0") for t in (3000, 3400) for g in (4.5, 5.5)]
)


def in_scope(p: dict) -> bool:
    return (TEFF_RANGE[0] <= p["teff"] <= TEFF_RANGE[1]
            and LOGG_RANGE[0] <= p["logg"] <= LOGG_RANGE[1] and p["alpha"] == 0.0)


def remote_catalogue(metallicities=METALLICITIES) -> list[dict]:
    """In-scope remote files with parameters, URL and size."""
    out = []
    for mh in metallicities:
        for name, size in list_ftp(f"{BASE}/{zdir(mh)}/"):
            p = parse_name(name)
            if p and in_scope(p):
                out.append(dict(p, mhs=mh, name=name, size=size,
                                url=f"{BASE}/{zdir(mh)}/{name}"))
    return out


def reduce_file(path: Path, R: float = R_BIN, lam_min: float = LAM_MIN,
                lam_max: float = LAM_MAX) -> SpecInt:
    with fits.open(path, memmap=False) as h:
        H = h[0].header
        I = np.asarray(h[0].data, dtype=np.float64)          # (nmu, nlam)
        mu = np.asarray(h["MU"].data, dtype=np.float64)
    lam = H["CRVAL1"] + H["CDELT1"] * (np.arange(H["NAXIS1"]) + 1 - H["CRPIX1"])
    if not (np.all(np.diff(mu) > 0) and mu[-1] == 1.0):
        raise ValueError(f"{path.name}: unexpected μ ordering")
    edges = log_edges(lam_min, lam_max, R)
    P = bin_photons(lam, I, edges).astype(np.float32)
    meta = dict(grid="PHOENIX-ACES-AGSS-COND-SpecInt-v3.0", geometry="spherical",
                teff=float(H["PHXTEFF"]), logg=float(H["PHXLOGG"]),
                mh=float(H["PHXM_H"]), alpha=0.0, vturb=float(H["PHXXI_L"]),
                reff_cm=float(H["PHXREFF"]), mixing_length=float(H["PHXMXLEN"]),
                phoenix_version=str(H["PHXVER"]), model_date=str(H["DATE"]),
                source_file=path.name, source_sha256=sha256(path),
                binning=dict(R=R, lam_min=lam_min, lam_max=lam_max,
                             units="P = ∫ I_λ λ dλ over bin; I_λ in erg/s/cm²/cm/sr, λ in Å"))
    return SpecInt(edges, mu, P, meta)


def reduced_path(entry: dict) -> Path:
    return REDUCED / "phoenix" / zdir(entry["mhs"]) / entry["name"].replace(".fits", ".npz")


def raw_path(entry: dict) -> Path:
    return RAW / "phoenix" / zdir(entry["mhs"]) / entry["name"]


def keep_raw(entry: dict) -> bool:
    return (entry["teff"], entry["logg"], entry["mhs"]) in set(P2_SUBSET)


def fetch_one(entry: dict) -> str:
    """Download, reduce and (outside the P2 subset) delete the raw file.
    Idempotent: skips models whose reduced file exists."""
    out = reduced_path(entry)
    raw = raw_path(entry)
    if out.exists() and (raw.exists() or not keep_raw(entry)) and _current(out):
        return "skip"
    fetch(entry["url"], raw, source="phoenix-v3.0", expected_size=entry["size"],
          log=keep_raw(entry))
    s = reduce_file(raw)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.stem + ".tmp.npz")
    s.save(tmp)
    tmp.replace(out)
    record("phoenix-v3.0-reduced", entry["url"], out,
           raw_sha256=s.meta["source_sha256"], raw_size=entry["size"],
           raw_kept=keep_raw(entry))
    if not keep_raw(entry):
        raw.unlink()
    return "done"


def _current(path: Path) -> bool:
    """True if a reduced file was made with the current binning."""
    try:
        b = SpecInt.load(path).meta["binning"]
    except Exception:
        return False
    return (b["R"], b["lam_min"], b["lam_max"]) == (R_BIN, LAM_MIN, LAM_MAX)


def load(teff: int, logg: float, mh: str | float = "-0.0") -> SpecInt:
    mhs = mh if isinstance(mh, str) else ("-0.0" if mh == 0 else f"{mh:+.1f}")
    p = REDUCED / "phoenix" / zdir(mhs) / file_name(teff, logg, mhs).replace(".fits", ".npz")
    return SpecInt.load(p)
