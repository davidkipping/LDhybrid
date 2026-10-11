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

Priors and empirical offsets:

    p = t.prior("hybrid5", 5600, 100, 4.4, 0.1, 0.1, 0.1, "TESS", systematic=True)
    observables("hybrid5")["h1"]                     # (const, grad): h'1 = const + grad·w
    m, C = condition(p["mean"], p["cov"], "hybrid5", "h1", delta, sigma)

Refitting interpolated intensity profiles (needs the profile files and h5py):

    t = LDTables("path/to/cds", profiles="path/to/profiles")
    r = t.coefficients("hybrid5", 5600, 4.4, 0.1, "TESS", refit=True)

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


def shape_slopes(law, x):
    """dT_j/dx on the law's natural variable x (ξ = μ² for the hybrids, μ for
    the others), shape (len(x), n); monotonicity is −Σ_j c_j dT_j/dx ≥ 0, as in
    the production fits (ldhybrid.laws)."""
    x = np.asarray(x, float)
    if law.startswith("hybrid"):
        def dpole(e):
            norm = e ** -2 - (1 + e) ** -2
            return -2.0 * (x + e) ** -3 / norm
        if law == "hybrid2":
            D = [-np.ones_like(x)] + [dpole(e) for e in POLES[law]]
        else:
            D = [-2 * x, -2 * (1 - x)] + [dpole(e) for e in POLES[law]]
    elif law == "quadratic":
        D = [-np.ones_like(x), -2 * (1 - x)]
    elif law == "claret4":
        D = [-(k / 2) * x ** (k / 2 - 1) for k in (1, 2, 3, 4)]
    else:
        raise ValueError(f"no linear shapes for {law!r}")
    return np.stack(D, axis=-1)


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


# ----------------------------------------------------------------- observables

#: μ grid of the quadratic fit that defines ℓ and δ (Sing et al. 2026: μ ≥ 0.2,
#: least squares uniform in μ; 161 points as in the comparison of the paper)
MU_ELL = np.linspace(0.2, 1.0, 161)


def observables(law, coeffs=None):
    """Exact linear functionals of a linear law's coefficients c:
    obs = const + grad · c.

      h1    = I(2/3)                     (Maxted 2023, h'1)
      h2    = I(2/3) − I(1/3)            (Maxted 2023, h'2)
      ell   = 1 − (c1 + c2)              (Sing et al. 2026, ℓ)
      delta = c2 / 4                     (Sing et al. 2026, δ)

    where c1, c2 are the quadratic-law coefficients of a least-squares fit
    uniform in μ over 0.2 ≤ μ ≤ 1 (MU_ELL). Every term is linear in the law's
    coefficients, because I = 1 − Σ c_j T_j(μ) and the quadratic fit is a fixed
    linear map. Returns {name: (const, grad)}, or {name: value} if coeffs are
    given. Not defined for power-2 (nonlinear)."""
    if law == "power2":
        raise ValueError("observables are linear functionals only for the linear laws")
    t1, t2 = shapes(law, [2 / 3, 1 / 3])
    A = np.column_stack([1 - MU_ELL, (1 - MU_ELL) ** 2])
    M = np.linalg.pinv(A) @ shapes(law, MU_ELL)            # (c1, c2) = M · c
    out = {"h1": (1.0, -t1), "h2": (0.0, t2 - t1),
           "ell": (1.0, -(M[0] + M[1])), "delta": (0.0, M[1] / 4)}
    if coeffs is None:
        return out
    c = np.array([coeffs[k] for k in COEFFS[law]] if isinstance(coeffs, dict) else coeffs, float)
    return {k: float(a + g @ c) for k, (a, g) in out.items()}


def condition(mean, cov, law, observable, delta, sigma):
    """Gaussian update of a coefficient prior N(mean, cov) by an empirical
    offset in one observable (observed − model = delta, with scatter sigma,
    e.g. a row of the empirical-offsets table). The observable is linear in
    the coefficients (observables()), so the update is exact:

        g = ∇obs,  K = C g / (gᵀ C g + σ²),  mean' = mean + K δ,  C' = C − K gᵀ C.

    Apply once per observable (treating them as independent). Returns
    (mean', cov'). The result is Gaussian: for the hybrid laws it can put
    probability outside the physical region; truncate there if needed."""
    m = np.asarray(mean, float)
    C = np.asarray(cov, float)
    g = observables(law)[observable][1]
    K = C @ g / (g @ C @ g + sigma ** 2)
    return m + K * delta, C - np.outer(K, g @ C)


# ----------------------------------------------------------------- refitting

#: production fitting set-up (LDhybrid v1.0): least squares uniform in r,
#: integrated exactly over the PCHIP profile by N_GL-point Gauss–Legendre on
#: every interval between nodes; constraints per law
N_GL = 12
CONSTRAINT = {"hybrid2": "physical", "hybrid4": "simplex", "hybrid5": "simplex",
              "quadratic": "physical", "claret4": "physical", "power2": None}


def _constraint_rows(law, kind, n=2001):
    """C c ≤ d: Σc ≤ 1 (I(0) ≥ 0) plus c ≥ 0 (simplex) or dI/dx ≥ 0 on a dense
    grid of x (ξ = μ² for the hybrids, μ ≥ 10⁻⁶ for the others), exactly as the
    production fits."""
    m = len(COEFFS[law])
    rows, rhs = [np.ones(m)], [1.0]
    if kind == "simplex":
        rows += list(-np.eye(m))
        rhs += [0.0] * m
    else:
        x = np.linspace(0.0 if law.startswith("hybrid") else 1e-6, 1.0, n)
        rows += list(shape_slopes(law, x))                 # dI/dx = −D c ≥ 0
        rhs += [0.0] * n
    C = np.array(rows)
    s = np.abs(C).max(1)
    s[s == 0] = 1
    return C / s[:, None], np.array(rhs) / s


def fit_profile(law, mu_q, y_q, q=None):
    """Fit `law` to y(μ) at quadrature points μ_q with weights q (from
    quadrature_r; default equal weights), with the production constraint.
    Returns the coefficients."""
    from scipy.optimize import least_squares, lsq_linear, minimize
    q = np.ones_like(mu_q) if q is None else np.asarray(q, float)
    sq = np.sqrt(q)
    if law == "power2":
        f = lambda p: sq * ((1 - p[0] * (1 - mu_q ** p[1])) - y_q)
        return least_squares(f, (0.7, 0.6), method="lm", xtol=1e-14, ftol=1e-14).x
    A = shapes(law, mu_q) * sq[:, None]
    b = (1 - y_q) * sq
    s = np.linalg.norm(A, axis=0)
    s[s == 0] = 1
    As = A / s
    C, d = _constraint_rows(law, CONSTRAINT[law])
    Cs = C / s
    v = np.linalg.lstsq(As, b, rcond=None)[0]
    if np.all(Cs @ v <= d + 1e-12):
        return v / s
    if CONSTRAINT[law] == "simplex":
        v = lsq_linear(As, b, bounds=(0, np.inf), method="bvls", tol=1e-14).x
        if np.sum(v / s) <= 1 + 1e-12:
            return v / s
    cons = [{"type": "ineq", "fun": lambda v: d - Cs @ v, "jac": lambda v: -Cs}]
    obj = lambda v: (np.sum((As @ v - b) ** 2), 2 * As.T @ (As @ v - b))
    r = minimize(obj, np.clip(v, 0, None) * 0.99 if CONSTRAINT[law] == "simplex" else v * 0.99,
                 jac=True, method="SLSQP", constraints=cons,
                 options=dict(ftol=1e-15, maxiter=500))
    return r.x / s


def quadrature_r(breaks, n=N_GL):
    """Points and weights for ∫ (·) dr, r = √(1 − μ²), over [μ_min, 1], exact for
    integrands that are smooth on each interval between the μ breakpoints
    (Gauss–Legendre in μ with dr = μ/r dμ where μ ≤ 0.5, in r elsewhere), as
    in the production fits."""
    t, wt = np.polynomial.legendre.leggauss(n)
    b = np.unique(np.asarray(breaks, float))
    mus, qs = [], []
    for a, c in zip(b[:-1], b[1:]):
        if a >= 0.5:
            ra, rc = np.sqrt(1 - c * c), np.sqrt(1 - a * a)
            r = 0.5 * (rc - ra) * t + 0.5 * (ra + rc)
            mus.append(np.sqrt(np.clip(1 - r * r, 0.0, None)))
            qs.append(0.5 * (rc - ra) * wt)
        else:
            mu = 0.5 * (c - a) * t + 0.5 * (a + c)
            mus.append(mu)
            qs.append(0.5 * (c - a) * wt * mu / np.sqrt(1 - mu * mu))
    return np.concatenate(mus), np.concatenate(qs)


# ----------------------------------------------------------------- the tables

class LDTables:
    """Reader and interpolator for the CDS tables in one directory."""

    def __init__(self, directory, cache=True, profiles=None):
        """cache: keep a parsed copy (<table>.npz) next to each table, if the
        directory is writable, so later sessions load in a fraction of a second.
        profiles: directory (or file) of band profiles, for refit()."""
        self.dir = Path(directory)
        self.cache = cache
        self.profiles = profiles
        self._data = {}
        self._interp = {}
        self._prof = {}
        self._tri = {}

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
                     extrapolate=None, refit=False):
        """Interpolated coefficients of `law` at (Teff, log g, [M/H]) in `band`.

        Returns dict(coeffs, physical, mucri, lc_ppm, inside, grid, method, band).
        lc_ppm is the interpolated light-curve error of the tabulated fits.
        refit=True interpolates the band profiles instead and refits (refit();
        needs the profile files and h5py; r-weighted method only)."""
        if refit:
            if method != "r":
                raise ValueError("refit reproduces the headline (r-weighted) method only")
            return self.refit(law, teff, logg, mh, band, grid)
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
    observables = staticmethod(observables)
    condition = staticmethod(condition)

    # ---------------------------------------------------------- priors

    def prior(self, law, teff, sig_teff, logg, sig_logg, mh, sig_mh, band, grid="mps2",
              method="r", n=2000, systematic=False, seed=0):
        """Monte-Carlo prior on the coefficients for uncertain stellar
        parameters: draws (Teff, log g, [M/H]) from independent Gaussians,
        interpolates the table at each draw (draws outside the grid are
        dropped) and returns dict(mean, cov, n_used, grids).

        systematic=True adds the model systematic: the covariance, across the
        grids phoenix, mps1 and mps2 that cover the star in this band, of their
        interpolated coefficients at the nominal parameters. With σ → 0 and no
        systematic, mean equals coefficients() and cov → 0."""
        rng = np.random.default_rng(seed)
        f, _ = self._interpolator(grid, law, band, method, None)
        X = np.column_stack([rng.normal(teff, sig_teff, n), rng.normal(logg, sig_logg, n),
                             rng.normal(mh, sig_mh, n)]) / SCALE
        Y = f(X)[:, :len(COEFFS[law])]
        Y = Y[np.all(np.isfinite(Y), axis=1)]
        if len(Y) == 0:
            raise ValueError("no draws inside the grid")
        mean = Y.mean(0)
        cov = np.cov(Y.T, ddof=1) if len(Y) > 1 else np.zeros((Y.shape[1],) * 2)
        grids = [grid]
        if systematic:
            W = []
            for gname in ("phoenix", "mps1", "mps2"):
                try:
                    r = self.coefficients(law, teff, logg, mh, band, gname, method)
                except (ValueError, FileNotFoundError):
                    continue
                W.append([r["coeffs"][k] for k in COEFFS[law]])
                if gname != grid:
                    grids.append(gname)
            if len(W) > 1:
                cov = cov + np.cov(np.array(W).T, ddof=1)
        return dict(law=law, mean=mean, cov=np.atleast_2d(cov), n_used=len(Y), grids=grids,
                    band=band, method=method)

    # ---------------------------------------------------------- refitting

    def _profiles(self, grid):
        """Band profiles after the limb cut (<grid>.h5 with datasets teff, logg,
        mh, bands, mu, y; NaN-padded), looked for in `profiles` (constructor
        argument), then in <tables>/profiles/. Needs h5py."""
        if grid in self._prof:
            return self._prof[grid]
        import h5py
        cands = ([Path(self.profiles)] if self.profiles else []) + [self.dir / "profiles"]
        for c in cands:
            p = c if c.suffix == ".h5" else c / f"{grid}.h5"
            if p.exists():
                with h5py.File(p, "r") as h:
                    d = {k: h[k][()] for k in ("teff", "logg", "mh", "mu", "y")}
                    d["bands"] = [b.decode() for b in h["bands"][()]]
                self._prof[grid] = d
                return d
        raise FileNotFoundError(f"no profile file for {grid} in {cands}")

    def refit(self, law, teff, logg, mh, band, grid="mps2"):
        """Coefficients from interpolated *intensity profiles* rather than
        interpolated coefficients: the profiles of the Delaunay corner models
        are combined (PCHIP in μ, as in production) on the union of their nodes,
        combined with the barycentric weights, and the law is refitted with
        the production weighting and constraint. Exact at the models (to the
        float32 storage of the profiles); between them it removes the part of
        the interpolation error that comes from interpolating coefficients."""
        from scipy.spatial import Delaunay
        P = self._profiles(grid)
        j = P["bands"].index(band)
        ok = np.isfinite(P["y"][:, j]).any(axis=1)
        idx = np.flatnonzero(ok)
        key = (grid, band)
        if key not in self._tri:
            X = np.column_stack([P["teff"][idx], P["logg"][idx], P["mh"][idx]]) / SCALE
            self._tri[key] = Delaunay(X)
        tri = self._tri[key]
        x = np.array([teff, logg, mh]) / SCALE
        s = int(tri.find_simplex(x))
        if s < 0:
            raise ValueError(f"({teff}, {logg}, {mh}) is outside the {grid} grid for {band}")
        T = tri.transform[s]
        b = T[:3] @ (x - T[3])
        lam = np.r_[b, 1 - b.sum()]
        corners = idx[tri.simplices[s]]
        nodes = []
        for k in corners:
            m = np.isfinite(P["mu"][k, j])
            nodes.append((P["mu"][k, j][m].astype(float), P["y"][k, j][m].astype(float)))
        from scipy.interpolate import PchipInterpolator
        used = [(l, mu, y) for l, (mu, y) in zip(lam, nodes) if l > 1e-14]
        mu_min = max(mu[0] for _, mu, _ in used)
        breaks = np.concatenate([mu[mu >= mu_min] for _, mu, _ in used] + [[mu_min]])
        mu_q, q = quadrature_r(breaks)
        y_q = sum(l * PchipInterpolator(mu, y)(mu_q) for l, mu, y in used)
        c = fit_profile(law, mu_q, y_q, q)
        coeffs = dict(zip(COEFFS[law], (float(v) for v in c)))
        return dict(law=law, coeffs=coeffs, physical=is_physical(law, coeffs), grid=grid,
                    method="r", band=band, refit=True,
                    corners=[(float(P["teff"][k]), float(P["logg"][k]), float(P["mh"][k]))
                             for k in corners], weights=lam.tolist())


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
