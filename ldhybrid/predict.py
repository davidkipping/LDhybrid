"""Predictions from the v1.0 tables for observed stars (P4).

* weights(): linear interpolation of a law's tabulated weights at (Teff,
  log g, [M/H]) in one band, on scaled coordinates (Teff/100 K, log g/0.5,
  [M/H]/0.5) with a Delaunay triangulation of the available models (so grid
  holes are bridged, never filled from another grid). Linear interpolation
  keeps hybrid weights inside their convex physical region.
* magnetic(): MURaM correction Δw, interpolated in Teff between its stars.
* Law-agnostic statistics of a profile (Maxted 2023: h'1 = I(2/3)/I(1),
  h'2 = h'1 − I(1/3)/I(1)), and SPAM-style refits (Howarth 2011): the law a
  light-curve fit would return for this profile and geometry.
"""

from __future__ import annotations

import csv
from functools import lru_cache

import h5py
import numpy as np
from scipy.interpolate import LinearNDInterpolator

from . import laws
from . import lightcurve as lc
from .paths import DATA

P = DATA / "products"
SCALE = np.array([100.0, 0.5, 0.5])


@lru_cache(maxsize=None)
def _h5(grid):
    return h5py.File(P / f"{grid}.h5")


@lru_cache(maxsize=None)
def _interp(grid, law, method, band):
    h = _h5(grid)
    bands = [b.decode() for b in h["bands"][()]]
    j = bands.index(band)
    X = np.column_stack([h["teff"][()], h["logg"][()], h["mh"][()]]) / SCALE
    ok = h["covered"][:, j]
    W = h[f"{law}|{method}|w"][:, j, :]
    return LinearNDInterpolator(X[ok], W[ok]), X[ok].min(0), X[ok].max(0)


def weights(grid, law, band, teff, logg, mh, method="r"):
    """Interpolated weights, or NaN outside the grid's convex hull."""
    f, lo, hi = _interp(grid, law, method, band)
    return np.asarray(f(np.array([teff, logg, mh]) / SCALE))[0]


@lru_cache(maxsize=None)
def _corrections():
    return list(csv.DictReader(open(P / "corrections.csv")))


def magnetic(law, band, teff, state="B100G", mh=0.0):
    """MURaM Δw for a field state, linear in Teff between the [M/H] = 0 stars
    (clamped at the ends)."""
    rows = [r for r in _corrections() if r["source"] == "MURaM" and r["kind"] == "magnetic"
            and r["state"] == state and r["law"] == law and r["band"] == band
            and float(r["mh"]) == mh]
    if not rows:
        return None
    rows.sort(key=lambda r: float(r["teff"]))
    T = np.array([float(r["teff"]) for r in rows])
    m = laws.get(law).m
    D = np.array([[float(r[f"dw{k + 1}"]) for k in range(m)] for r in rows])
    return np.array([np.interp(teff, T, D[:, k]) for k in range(m)])


def h_prime(I_of_mu):
    """Maxted (2023): h'1 = I(2/3)/I(1), h'2 = h'1 − I(1/3)/I(1)."""
    i1 = I_of_mu(np.array([1.0]))[0]
    a, b = I_of_mu(np.array([2 / 3, 1 / 3])) / i1
    return a, a - b


def spam(I_of_mu, p, b, law="quadratic", n=200):
    """Coefficients of a *linear* law fitted to the light curve of profile
    I(μ) for one transit chord (p, b): what a light-curve fit would return."""
    L = laws.get(law)
    z = lc.chord_z(p, b, n)
    F = lc.light_curve(I_of_mu, z, p)
    # F = 1 − dF/F★ with dF = dF1 − Σ w dF_j and F★ = F1 − Σ w F_j: solve the
    # linear equations (1 − F) F★ = dF for w by least squares
    I_r = lambda f: (lambda r: f(np.sqrt(np.clip(1 - r * r, 0, None))))
    one = lambda mu: np.ones_like(mu)
    dF1, F1 = lc.blocked(I_r(one), z, p), lc.star_flux(one)
    cols = [(lc.blocked(I_r(lambda mu, j=j: L.T(mu)[..., j]), z, p),
             lc.star_flux(lambda mu, j=j: L.T(mu)[..., j])) for j in range(L.m)]
    # (1 − F)(F1 − Σ w Fj) = dF1 − Σ w dFj  →  Σ w [dFj − (1 − F) Fj] = dF1 − (1 − F) F1
    A = np.column_stack([c[0] - (1 - F) * c[1] for c in cols])
    y = dF1 - (1 - F) * F1
    return np.linalg.lstsq(A, y, rcond=None)[0]
