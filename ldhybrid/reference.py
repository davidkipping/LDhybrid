"""Readers for published reference tables (mirrored by scripts/fetch_tables.py)."""

from __future__ import annotations

import gzip

import numpy as np

from .paths import RAW

# Kostogryz et al. 2022 (J/A+A/666/A60): CLVs at 23 μ (plus I(1)), 4 bands
KOSTOGRYZ_MU = np.array([0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.35, 0.3, 0.25, 0.22, 0.2, 0.17,
                         0.15, 0.12, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.03, 0.02, 0.01])


def _rows(path):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as f:
        return [line.split() for line in f if line.strip()]


def kostogryz_clv(table: str = "tablee4"):
    """Structured array: mh, teff, logg, band, I1, y[23] (I/I(1) at
    KOSTOGRYZ_MU). tablee3 = Set 1, tablee4 = Set 2."""
    rows = _rows(RAW / "vizier" / "666_A60" / f"{table}.dat.gz")
    mh = np.array([float(r[0]) for r in rows])
    teff = np.array([int(r[1]) for r in rows])
    logg = np.array([float(r[2]) for r in rows])
    band = np.array([r[3] for r in rows])
    I1 = np.array([float(r[4]) for r in rows])
    y = np.array([[float(v) for v in r[5:28]] for r in rows])
    return dict(mh=mh, teff=teff, logg=logg, band=band, I1=I1, y=y)


def kostogryz_coeffs(table: str = "tablee5"):
    """mh, teff, logg, band, power-2 (c, alpha), 4-term (a1..a4)."""
    rows = _rows(RAW / "vizier" / "666_A60" / f"{table}.dat.gz")
    a = np.array([[float(v) for v in r[4:10]] for r in rows])
    return dict(mh=np.array([float(r[0]) for r in rows]),
                teff=np.array([int(r[1]) for r in rows]),
                logg=np.array([float(r[2]) for r in rows]),
                band=np.array([r[3] for r in rows]),
                p2=a[:, :2], claret4=a[:, 2:])
