"""Interpolate the LDhybrid v1.0 coefficient tables (the CDS release).

Self-contained: needs only numpy and scipy, so this file can be copied and
used on its own with the downloaded tables (`<grid>_<law>.dat[.gz]`).

    from ldhybrid.interpolate import LDTables
    t = LDTables("path/to/cds")                      # directory with the tables
    r = t.coefficients("hybrid5", teff=5600, logg=4.4, mh=0.1, band="TESS",
                       grid="mps2")
    r["coeffs"]          # {'w1': ..., ..., 'w5': ...}
    r["physical"]        # positive and monotonic?
    t.intensity("hybrid5", r["coeffs"], mu)          # I(mu)/I(1)

Command line:

    python -m ldhybrid.interpolate hybrid5 5600 4.4 0.1 TESS --grid mps2 --tables path/to/cds

Interpolation is linear within one grid, on a Delaunay triangulation of the
models in scaled coordinates (Teff/100 K, log g/0.5, [M/H]/0.5), so grid holes
are bridged and the result is exact at the models. For the hybrid laws,
linear interpolation keeps the weights inside their convex physical region.
Points outside the grid raise ValueError unless extrapolate="nearest".
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

GRIDS = ("phoenix", "mps1", "mps2", "kurucz")
COEFFS = {"hybrid2": ("w1", "w2"), "hybrid4": ("w1", "w2", "w3", "w4"),
          "hybrid5": ("w1", "w2", "w3", "w4", "w5"), "quadratic": ("u1", "u2"),
          "claret4": ("a1", "a2", "a3", "a4"), "power2": ("c", "alpha")}
LAWS = tuple(COEFFS)
#: v1.0 pole positions (frozen; ldhybrid/v1.py)
POLES = {"hybrid2": (0.36,), "hybrid4": (0.01176198, 0.30407926),
         "hybrid5": (0.00159655, 0.04127504, 0.38743472)}
SCALE = np.array([100.0, 0.5, 0.5])


# ----------------------------------------------------------------- the laws

def _pole_shape(mu, eps):
    """Π_ε(μ) = [(μ²+ε)⁻² − (1+ε)⁻²] / [ε⁻² − (1+ε)⁻²]: 1 at the limb, 0 at the centre."""
    return ((mu * mu + eps) ** -2 - (1 + eps) ** -2) / (eps ** -2 - (1 + eps) ** -2)


def shapes(law, mu):
    """The terms T_j(μ) of I(μ)/I(1) = 1 − Σ_j c_j T_j(μ), shape (len(mu), n)."""
    mu = np.asarray(mu, float)
    if law == "hybrid2":
        T = [1 - mu ** 2] + [_pole_shape(mu, e) for e in POLES[law]]
    elif law in ("hybrid4", "hybrid5"):
        T = [1 - mu ** 4, (1 - mu ** 2) ** 2] + [_pole_shape(mu, e) for e in POLES[law]]
    elif law == "quadratic":
        T = [1 - mu, (1 - mu) ** 2]
    elif law == "claret4":
        T = [1 - mu ** (k / 2) for k in (1, 2, 3, 4)]
    else:
        raise ValueError(f"no linear shapes for {law!r}")
    return np.stack(T, axis=-1)


def intensity(law, coeffs, mu):
    """I(μ)/I(1) for a law and its coefficients (dict or sequence)."""
    c = np.array([coeffs[k] for k in COEFFS[law]] if isinstance(coeffs, dict) else coeffs, float)
    mu = np.asarray(mu, float)
    if law == "power2":
        return 1 - c[0] * (1 - mu ** c[1])
    return 1 - shapes(law, mu) @ c


def is_physical(law, coeffs, n=2001):
    """Positive and non-increasing towards the limb on a dense μ grid."""
    mu = np.linspace(0.0, 1.0, n)
    with np.errstate(all="ignore"):
        I = intensity(law, coeffs, mu)
    return bool(np.all(np.isfinite(I)) and I[0] >= -1e-10 and np.all(np.diff(I) >= -1e-10))


# ----------------------------------------------------------------- the tables

class LDTables:
    """Reader and interpolator for the CDS tables in one directory."""

    def __init__(self, directory, cache=True):
        """cache: keep a parsed copy (<table>.npz) next to each table, if the
        directory is writable, so later sessions load in a fraction of a second."""
        self.dir = Path(directory)
        self.cache = cache
        self._data = {}
        self._interp = {}

    def _file(self, grid, law):
        for name in (f"{grid}_{law}.dat", f"{grid}_{law}.dat.gz"):
            if (self.dir / name).exists():
                return self.dir / name
        raise FileNotFoundError(f"{grid}_{law}.dat[.gz] not found in {self.dir}")

    def table(self, grid, law):
        """Parsed table as a dict of numpy arrays (cached)."""
        key = (grid, law)
        if key in self._data:
            return self._data[key]
        path = self._file(grid, law)
        npz = path.with_name(path.name.split(".")[0] + ".npz")
        if self.cache and npz.exists() and npz.stat().st_mtime >= path.stat().st_mtime:
            with np.load(npz) as z:
                self._data[key] = {k: z[k] for k in z.files}
            return self._data[key]
        op = gzip.open if path.suffix == ".gz" else open
        with op(path, "rt") as f:
            rows = [ln.split() for ln in f if ln.strip()]
        n = len(COEFFS[law])
        a = np.array(rows, dtype=object)
        d = dict(teff=a[:, 0].astype(float), logg=a[:, 1].astype(float),
                 mh=a[:, 2].astype(float), band=a[:, 3].astype(str),
                 method=a[:, 4].astype(str), mucri=a[:, 5].astype(float),
                 coeffs=a[:, 6:6 + n].astype(float))
        for k, name in enumerate(("sigma", "dImax", "dFlux", "LCppm", "Phys", "Bound")):
            d[name] = a[:, 6 + n + k].astype(float)
        self._data[key] = d
        if self.cache:
            try:
                np.savez(npz, **d)
            except OSError:
                pass
        return self._data[key]

    def bands(self, grid, law="hybrid5"):
        return sorted(set(self.table(grid, law)["band"]))

    def _interpolator(self, grid, law, band, method, extrapolate):
        key = (grid, law, band, method, extrapolate)
        if key not in self._interp:
            d = self.table(grid, law)
            m = (d["band"] == band) & (d["method"] == method)
            if not m.any():
                raise ValueError(f"no {grid} {law} rows for band {band!r}, method {method!r}")
            X = np.column_stack([d["teff"][m], d["logg"][m], d["mh"][m]]) / SCALE
            Y = np.column_stack([d["coeffs"][m], d["mucri"][m], d["LCppm"][m]])
            f = LinearNDInterpolator(X, Y)
            g = NearestNDInterpolator(X, Y) if extrapolate == "nearest" else None
            self._interp[key] = (f, g)
        return self._interp[key]

    def coefficients(self, law, teff, logg, mh, band, grid="mps2", method="r",
                     extrapolate=None):
        """Interpolated coefficients of `law` at (Teff, log g, [M/H]) in `band`.

        Returns dict(coeffs, physical, mucri, lc_ppm, inside, grid, method, band).
        lc_ppm is the interpolated light-curve error of the tabulated fits."""
        if law not in COEFFS:
            raise ValueError(f"law must be one of {LAWS}")
        if grid not in GRIDS:
            raise ValueError(f"grid must be one of {GRIDS}")
        f, g = self._interpolator(grid, law, band, method, extrapolate)
        x = np.array([[teff, logg, mh]]) / SCALE
        y = f(x)[0]
        inside = bool(np.all(np.isfinite(y)))
        if not inside:
            if g is None:
                raise ValueError(f"({teff}, {logg}, {mh}) is outside the {grid} grid for "
                                 f"{band}; pass extrapolate='nearest' to use the nearest model")
            y = g(x)[0]
        n = len(COEFFS[law])
        coeffs = dict(zip(COEFFS[law], (float(v) for v in y[:n])))
        return dict(law=law, coeffs=coeffs, physical=is_physical(law, coeffs),
                    mucri=float(y[n]), lc_ppm=float(y[n + 1]), inside=inside,
                    grid=grid, method=method, band=band)

    intensity = staticmethod(intensity)


def main(argv=None):
    p = argparse.ArgumentParser(description="Interpolate LDhybrid v1.0 coefficient tables.")
    p.add_argument("law", choices=LAWS)
    p.add_argument("teff", type=float)
    p.add_argument("logg", type=float)
    p.add_argument("mh", type=float)
    p.add_argument("band")
    p.add_argument("--grid", default="mps2", choices=GRIDS)
    p.add_argument("--method", default="r", choices=("r", "mu"))
    p.add_argument("--tables", default=".", help="directory with the CDS tables")
    p.add_argument("--nearest", action="store_true", help="nearest model outside the grid")
    a = p.parse_args(argv)
    r = LDTables(a.tables).coefficients(a.law, a.teff, a.logg, a.mh, a.band, a.grid,
                                        a.method, "nearest" if a.nearest else None)
    print(" ".join(f"{k} = {v:.6f}" for k, v in r["coeffs"].items()))
    print(f"physical = {r['physical']}, mu_cri = {r['mucri']:.5f}, "
          f"fit LC error ≈ {r['lc_ppm']:.2f} ppm, inside grid = {r['inside']}")


if __name__ == "__main__":
    main()
