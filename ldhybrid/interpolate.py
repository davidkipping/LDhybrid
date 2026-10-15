"""Interpolate the LDhybrid v1.0 coefficient tables (the CDS release).

Self-contained: needs only numpy and scipy, so this file can be copied and
used on its own with the downloaded tables (`<grid>_<law>.dat[.gz]`, and for
corrections and offsets `corrections.dat`, `offsets.dat`).

    from ldhybrid.interpolate import LDTables
    t = LDTables("path/to/cds")                      # directory with the tables
    r = t.coefficients("hybrid5", teff=5600, logg=4.4, mh=0.1, band="TESS",
                       grid="mps2")
    r["coeffs"]          # {'w1': ..., ..., 'w5': ...}
    r["physical"]        # positive and monotonic?
    t.intensity("hybrid5", r["coeffs"], mu)          # I(mu)/I(1)

The recommended FGK prior (MPS-ATLAS Set 2 + MURaM 100 G), and an empirical
offset measured in the SAME correction state (never apply an offset measured
against uncorrected predictions to a corrected prior: that counts the
correction twice):

    p = t.prior("hybrid5", 5600, 100, 4.4, 0.1, 0.1, 0.1, "Kepler",
                correction="B100G", systematic=True)
    m, C = t.apply_offsets(p, regime="Kepler")       # joint (h'1, h'2) update
    # or by hand, one or several observables (linear in the weights):
    m, C = condition(p["mean"], p["cov"], "hybrid5", ["h1", "h2"],
                     [-0.0023, -0.0039], [0.0052, 0.0126], rho=0.04)

Transit light curve of any tabulated law (circular planet; the same
numerical quadrature as the paper's reference light curves):

    F = light_curve("hybrid5", r["coeffs"], z, p=0.1)

Refitting interpolated intensity profiles (needs the profile files and h5py):

    t = LDTables("path/to/cds", profiles="path/to/profiles")
    r = t.coefficients("claret4", 5600, 4.4, 0.1, "TESS", refit=True, mu_min=0.02)

Command line:

    python -m ldhybrid.interpolate hybrid5 5600 4.4 0.1 TESS --grid mps2 --tables path/to/cds

Interpolation is linear within one grid, on a Delaunay triangulation of the
models in scaled coordinates (Teff/100 K, log g/0.5, [M/H]/0.5), so grid holes
are bridged and the result is exact at the models. For the hybrid laws,
linear interpolation keeps the weights inside their convex physical region.
Points outside the grid raise ValueError unless extrapolate="nearest".

hybrid4 is tabulated on its curvature-simplex basis (H4_BASIS): w ≥ 0,
Σw ≤ 1 is exactly the set of hybrid4 profiles that are positive, monotonic,
concave in μ² and of convex slope in μ². from_q()/to_q() map the hybrid
weights to and from sampling coordinates q in the unit hypercube (for
hybrid4, q1 = (1 − I(0))⁴, the analogue of Kipping 2013's q1); a uniform q
is a uniform prior over the region. Weights on SquishierPlanet's monomial-pole
hybrid4 shapes are w @ H4_BASIS.
"""

from __future__ import annotations

import argparse
import gzip
import warnings
from pathlib import Path

import numpy as np
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

GRIDS = ("phoenix", "mps1", "mps2")
COEFFS = {"hybrid2": ("w1", "w2"), "hybrid4": ("w1", "w2", "w3", "w4"),
          "hybrid5": ("w1", "w2", "w3", "w4", "w5"), "quadratic": ("u1", "u2"),
          "claret4": ("a1", "a2", "a3", "a4"), "power2": ("c", "alpha")}
LAWS = tuple(COEFFS)
#: v1.0 pole positions (frozen; ldhybrid/v1.py)
POLES = {"hybrid2": (0.36,), "hybrid4": (0.028601207249486824, 0.38582417584735224),
         "hybrid5": (0.0015965451792942786, 0.04127504477274302, 0.38743471584956696)}
SCALE = np.array([100.0, 0.5, 0.5])
#: layout of the parsed-table cache (<table>.npz); bump when the parser changes
TABLE_FORMAT = 2
#: MURaM magnetic states in corrections.dat
STATES = ("ssd", "B100G", "B200G", "B300G")
#: fixed-width layouts of the auxiliary tables (field widths, single spaces between)
_CORR_W = [("Source", 7), ("Kind", 8), ("State", 5), ("Star", 9), ("Teff", 5), ("logg", 5),
           ("mh", 5), ("Band", 25), ("Law", 7), ("dw1", 9), ("dw2", 9), ("dw3", 9), ("dw4", 9),
           ("dw5", 9), ("dLC", 7), ("Ext", 1), ("Base", 22)]
_OFF_W = [("Regime", 15), ("Bands", 26), ("Grid", 7), ("Corr", 11), ("Route", 7), ("Stat", 5),
          ("Offset", 7), ("e_Offset", 6), ("Scatter", 6), ("rho", 6), ("N", 2), ("Ref", 30)]
#: offsets-table statistic -> observables() name
_STAT = {"h1p": "h1", "h2p": "h2", "ell": "ell", "delta": "delta"}
#: lowest Teff shared by PHOENIX and MPS-ATLAS (for the M-dwarf systematic floor)
TEFF_SHARED = 3500.0


# ----------------------------------------------------------------- the laws

def _pole_shape(mu, eps):
    """Π_ε(μ) = [(μ²+ε)⁻² − (1+ε)⁻²] / [ε⁻² − (1+ε)⁻²]: 1 at the limb, 0 at the centre."""
    return ((mu * mu + eps) ** -2 - (1 + eps) ** -2) / (eps ** -2 - (1 + eps) ** -2)


def _hybrid4_basis(eps):
    """hybrid4's shapes as rows of weights on (1 − μ⁴, (1 − μ²)², Π_ε0, Π_ε1): the
    four non-trivial vertices of the curvature simplex (ldhybrid/v1.py). With
    g = dI/dξ, ξ = μ²: S1 = 1 − μ², S2 = (1 − μ²)², S3 with I(0) = 0 and
    g(1) = g'(1) = g''(1) = 0, S4 with I(0) = 0, g(1) = g'(1) = 0, g''(0) = 0."""
    e0, e1 = eps
    a0, a1 = (1.0 / (e ** -2 - (1 + e) ** -2) for e in eps)
    g = lambda x: [2 * x, 2 * (1 - x), 2 * a0 * (x + e0) ** -3, 2 * a1 * (x + e1) ** -3]
    g1 = lambda x: [2.0, -2.0, -6 * a0 * (x + e0) ** -4, -6 * a1 * (x + e1) ** -4]
    g2 = lambda x: [0.0, 0.0, 24 * a0 * (x + e0) ** -5, 24 * a1 * (x + e1) ** -5]
    rhs = [1.0, 0.0, 0.0, 0.0]
    S3 = np.linalg.solve(np.array([[1.0] * 4, g(1.0), g1(1.0), g2(1.0)]), rhs)
    S4 = np.linalg.solve(np.array([[1.0] * 4, g(1.0), g1(1.0), g2(0.0)]), rhs)
    return np.array([[0.5, 0.5, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], S3, S4])


#: hybrid4 basis (since 2026-10-14): S_k = Σ_j H4_BASIS[k, j] T_j. w ≥ 0, Σw ≤ 1 is
#: then exactly the set of hybrid4 profiles that are positive, monotonic, concave
#: in μ² and of convex slope in μ². Weights on the monomial-pole shapes
#: (SquishierPlanet's hybrid4) are w @ H4_BASIS.
H4_BASIS = _hybrid4_basis(POLES["hybrid4"])
#: hybrid2's triangle vertices (exact physical region; O is the third)
H2_VERTICES = ((-0.124615, 1.124615), (1.201039, -0.201039))


def from_q(law, q):
    """Hybrid weights from sampling coordinates q in the unit hypercube (a
    uniform q gives a uniform prior over the law's region). hybrid4 (the
    analogue of Kipping 2013's q1, q2): s = Σw = 1 − I(0) = q1^(1/4),
    w1 = s(1 − q2^(1/3)), w2 = s q2^(1/3)(1 − q3^(1/2)),
    w3 = s q2^(1/3) q3^(1/2)(1 − q4), w4 = s q2^(1/3) q3^(1/2) q4.
    hybrid5: stick-breaking w_j = (1 − q_j^(1/(n+1−j))) Π_{i<j} q_i^(1/(n+1−i)).
    hybrid2: w = √q1 [(1 − q2) V_c + q2 V_ℓ]."""
    q = np.asarray(q, float)
    if law == "hybrid4":
        s = q[..., 0] ** 0.25
        a, b = q[..., 1] ** (1 / 3), q[..., 2] ** 0.5
        return np.stack([s * (1 - a), s * a * (1 - b), s * a * b * (1 - q[..., 3]),
                         s * a * b * q[..., 3]], axis=-1)
    if law == "hybrid5":
        n = q.shape[-1]
        out, stick = [], 1.0
        for j in range(n):
            r = q[..., j] ** (1.0 / (n - j))
            out.append(stick * (1 - r))
            stick = stick * r
        return np.stack(out, axis=-1)
    if law == "hybrid2":
        Vc, Vl = (np.array(v) for v in H2_VERTICES)
        s = np.sqrt(q[..., 0])[..., None]
        return s * ((1 - q[..., 1])[..., None] * Vc + q[..., 1][..., None] * Vl)
    raise ValueError("sampling coordinates are defined for the hybrid laws")


def to_q(law, coeffs):
    """Inverse of from_q for one set of weights (dict or sequence). Weights
    outside the region give q outside [0, 1]; rounding noise within 1e-9 of 0
    or 1 is put on the bound."""
    q = _to_q(law, coeffs)
    q = np.where((q < 0) & (q > -1e-9), 0.0, q)
    return np.where((q > 1) & (q < 1 + 1e-9), 1.0, q)


def _to_q(law, coeffs):
    """Inverse of from_q for one set of weights (dict or sequence). Weights
    outside the region give q outside [0, 1]."""
    w = np.array([coeffs[k] for k in COEFFS[law]] if isinstance(coeffs, dict) else coeffs, float)
    if law == "hybrid4":
        s = w.sum()
        if s == 0:
            return np.zeros(4)
        f = w / s
        r1 = 1 - f[0]
        r2 = 1 - f[1] / r1 if r1 != 0 else 0.0
        f34 = f[2] + f[3]
        return np.array([np.sign(s) * abs(s) ** 4, np.sign(r1) * abs(r1) ** 3,
                         r2 * r2 if r2 >= 0 else r2, f[3] / f34 if f34 != 0 else 0.0])
    if law == "hybrid5":
        n, q, stick = w.size, np.empty(w.size), 1.0
        for j in range(n):
            f = 1.0 - w[j] / stick if stick > 0 else 1.0
            q[j] = abs(f) ** (n - j) if f >= 0 else f
            stick = stick * f if f > 0 else 0.0
        return q
    if law == "hybrid2":
        Vc, Vl = (np.array(v) for v in H2_VERTICES)
        s, t = np.linalg.solve(np.column_stack([Vc, Vl - Vc]), w)
        # s|s|, not s²: a negative depth (limb brighter than centre) gives q1 < 0
        return np.array([s * abs(s), t / s if s != 0 else 0.0])
    raise ValueError("sampling coordinates are defined for the hybrid laws")


def shapes(law, mu):
    """The terms T_j(μ) of I(μ)/I(1) = 1 − Σ_j c_j T_j(μ), shape (len(mu), n)."""
    mu = np.asarray(mu, float)
    if law == "hybrid2":
        T = [1 - mu ** 2] + [_pole_shape(mu, e) for e in POLES[law]]
    elif law == "hybrid4":
        T0 = np.stack([1 - mu ** 4, (1 - mu ** 2) ** 2] + [_pole_shape(mu, e) for e in POLES[law]], -1)
        return T0 @ H4_BASIS.T
    elif law == "hybrid5":
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
            if law == "hybrid4":
                return np.stack(D, axis=-1) @ H4_BASIS.T
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


def condition(mean, cov, law, observable, delta, sigma, rho=None):
    """Gaussian update of a coefficient prior N(mean, cov) by empirical offsets
    (observed − model) in one or several observables. Each observable is linear
    in the coefficients (observables()), so the update is exact:

        G = ∇obs (k × n),  S = G C Gᵀ + R,  K = C Gᵀ S⁻¹,
        mean' = mean + K δ,  C' = C − K G C.

    observable: a name ("h1", "h2", "ell", "delta") or a list of names.
    delta: the offset(s). sigma: their 1σ scatter (scalar or one per
    observable), or a full k × k covariance R. rho: correlation between two
    observables (e.g. the h'1–h'2 correlation in the offsets table), used when
    sigma gives the standard deviations. Use offsets measured in the same
    correction state as the prior. Offsets are relative to the prior mean, so
    apply several observables jointly (a list); a sequential second update
    would have to subtract the change of its observable from its delta.
    Returns (mean', cov'). The result is
    Gaussian: for the hybrid laws it can put probability outside the physical
    region; truncate there if needed."""
    m = np.asarray(mean, float)
    C = np.asarray(cov, float)
    names = [observable] if isinstance(observable, str) else list(observable)
    obs = observables(law)
    G = np.array([obs[k][1] for k in names])
    d = np.atleast_1d(np.asarray(delta, float))
    sig = np.asarray(sigma, float)
    if sig.ndim == 2:
        R = sig
    else:
        sd = np.broadcast_to(np.atleast_1d(sig), (len(names),)).astype(float)
        R = np.diag(sd ** 2)
        if rho is not None and len(names) == 2:
            R[0, 1] = R[1, 0] = rho * sd[0] * sd[1]
    S = G @ C @ G.T + R
    K = np.linalg.solve(S.T, (C @ G.T).T).T
    return m + K @ d, C - K @ G @ C


# ----------------------------------------------------------------- light curves

def _gl(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * np.pi * (x + 1.0), 0.5 * np.pi * w


def _star_flux(I_of_mu, n=4000):
    th, wt = _gl(n)
    mu = 0.5 * (1.0 - np.cos(th))
    return float(2.0 * np.pi * np.sum(wt * 0.5 * np.sin(th) * I_of_mu(mu) * mu))


def _blocked(I_of_r, z, p, n_sub=4, n_gl=48):
    th, wt = _gl(n_gl)
    z = np.atleast_1d(np.asarray(z, float))
    lo = np.clip(z - p, 0.0, None)
    hi = np.minimum(1.0, z + p)
    bp = np.stack([lo, np.clip(p - z, 0.0, 1.0), np.clip(z - p, 0.0, 1.0),
                   np.minimum(1.0, z + p), hi], axis=1)
    bp = np.sort(np.clip(bp, lo[:, None], hi[:, None]), axis=1)
    edges = [bp[:, :1]]
    for j in range(bp.shape[1] - 1):
        a, b = bp[:, j:j + 1], bp[:, j + 1:j + 2]
        edges.append(a + (b - a) * np.arange(1, n_sub + 1) / n_sub)
    E = np.concatenate(edges, axis=1)
    a, b = E[:, :-1, None], E[:, 1:, None]
    r = a + (b - a) * 0.5 * (1.0 - np.cos(th))
    dr = (b - a) * 0.5 * np.sin(th) * wt
    zz = z[:, None, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        cosk = (r * r + zz * zz - p * p) / (2.0 * r * zz)
    kap = np.where(zz == 0, np.where(r < p, np.pi, 0.0),
                   np.arccos(np.clip(np.nan_to_num(cosk, nan=-1.0), -1.0, 1.0)))
    return np.sum(I_of_r(r) * 2.0 * r * kap * dr, axis=(1, 2))


def light_curve(law, coeffs, z, p):
    """Normalised flux F(z) = 1 − ΔF/F★ of a circular planet of radius p
    (stellar radii) at sky separations z, for any tabulated law. Numerical
    quadrature, identical to the reference light curves of the paper
    (ldhybrid.lightcurve): Gauss–Legendre with endpoint substitutions on each
    segment between the chord breakpoints (agrees with the closed-form hybrid5
    light curve to ~1e-8 ppm). For speed and differentiability in a fitting
    code, use a dedicated transit code with the same law."""
    I_mu = lambda mu: intensity(law, coeffs, mu)
    I_r = lambda r: I_mu(np.sqrt(np.clip(1.0 - r * r, 0.0, None)))
    return 1.0 - _blocked(I_r, z, p) / _star_flux(I_mu)


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
                if "_format" in z.files and int(z["_format"]) == TABLE_FORMAT:
                    self._data[key] = {k: z[k] for k in z.files if k != "_format"}
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
        extra = a.shape[1] - (6 + n + 6)          # hybrid tables: q and expanded form
        if extra > 0:
            d["q"] = a[:, 6 + n:6 + 2 * n].astype(float)
            d["expanded"] = a[:, 6 + 2 * n:6 + n + extra].astype(float)
        for k, name in enumerate(("sigma", "dImax", "dFlux", "LCppm", "Phys", "Bound")):
            d[name] = a[:, 6 + n + extra + k].astype(float)
        self._data[key] = d
        if self.cache:
            try:
                np.savez(npz, _format=TABLE_FORMAT, **d)
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
                     extrapolate=None, refit=False, correction=None, mu_min=0.0):
        """Interpolated coefficients of `law` at (Teff, log g, [M/H]) in `band`.

        Returns dict(coeffs, physical, mucri, lc_ppm, inside, grid, method, band,
        correction, and for the hybrid laws q: the sampling coordinates, to_q()). lc_ppm is the interpolated light-curve error of the
        tabulated fits. refit=True interpolates the band profiles instead and
        refits (refit(); needs the profile files and h5py; r-weighted method
        only; mu_min as in refit()). correction: None or a MURaM state
        ("ssd", "B100G", "B200G", "B300G"; hybrid laws only), whose Δw from
        corrections.dat is added (delta_w())."""
        if refit:
            if method != "r":
                raise ValueError("refit reproduces the headline (r-weighted) method only")
            r = self.refit(law, teff, logg, mh, band, grid, mu_min=mu_min)
        else:
            if mu_min:
                raise ValueError("mu_min applies to refit=True only")
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
            r = dict(law=law, coeffs=coeffs, mucri=float(y[n]), lc_ppm=float(y[n + 1]),
                     inside=inside, grid=grid, method=method, band=band)
        if correction:
            dw = self.delta_w(law, band, teff, correction)
            r["coeffs"] = dict(zip(COEFFS[law], (v + d for v, d in
                                                  zip(r["coeffs"].values(), dw))))
        r["correction"] = correction
        r["physical"] = is_physical(law, r["coeffs"])
        if law.startswith("hybrid"):
            r["q"] = tuple(float(v) for v in to_q(law, r["coeffs"]))
        return r

    # ---------------------------------------------------------- corrections, offsets

    def _fixed(self, name, widths):
        """Rows of a fixed-width auxiliary table as dicts of strings ('' = null)."""
        path = self.dir / name
        if not path.exists() and (self.dir / (name + ".gz")).exists():
            path = self.dir / (name + ".gz")
        op = gzip.open if path.suffix == ".gz" else open
        out = []
        with op(path, "rt") as f:
            for ln in f:
                if not ln.strip():
                    continue
                d, pos = {}, 0
                for lab, w in widths:
                    d[lab] = ln[pos:pos + w].strip()
                    pos += w + 1
                out.append(d)
        return out

    def corrections(self):
        """corrections.dat as a list of dicts (cached)."""
        if "_corr" not in self._data:
            self._data["_corr"] = self._fixed("corrections.dat", _CORR_W)
        return self._data["_corr"]

    def delta_w(self, law, band, teff, state="B100G"):
        """MURaM magnetic Δw for `state`, interpolated linearly in Teff between
        the [M/H] = 0 MURaM stars and clamped outside them (as ldhybrid.predict).
        MURaM exists only at these Teff, at one log g each and [M/H] = 0, and
        has no intensities below μ = 0.1 (the ratio is held constant there):
        applying Δw at other log g and [M/H] is an extrapolation. teff may be
        an array."""
        if not law.startswith("hybrid"):
            raise ValueError("corrections are tabulated for the hybrid laws only")
        if state not in STATES:
            raise ValueError(f"state must be one of {STATES}")
        rows = [r for r in self.corrections() if r["Source"] == "MURaM" and r["Kind"] == "magnetic"
                and r["State"] == state and r["Law"] == law and r["Band"] == band
                and float(r["mh"]) == 0.0]
        if not rows:
            raise ValueError(f"no {state} correction for {law} in {band}")
        rows.sort(key=lambda r: float(r["Teff"]))
        T = np.array([float(r["Teff"]) for r in rows])
        D = np.array([[float(r[f"dw{k + 1}"]) for k in range(len(COEFFS[law]))] for r in rows])
        t = np.asarray(teff, float)
        if np.any((t < T[0]) | (t > T[-1])):
            warnings.warn(f"Teff outside the MURaM range {T[0]:.0f}-{T[-1]:.0f} K: the {state} "
                          "correction is clamped to the nearest MURaM star", stacklevel=2)
        out = np.stack([np.interp(t, T, D[:, k]) for k in range(D.shape[1])], axis=-1)
        return out

    def offsets(self):
        """offsets.dat as a list of dicts (cached)."""
        if "_off" not in self._data:
            self._data["_off"] = self._fixed("offsets.dat", _OFF_W)
        return self._data["_off"]

    def apply_offsets(self, prior, regime, statistics=None, route="hybrid5"):
        """Condition a prior (from prior()) on the tabulated empirical offsets
        for `regime` ("Kepler", "TESS", "JWST_0.86-5.5um"), measured against
        predictions of the same grid and correction state as the prior (an
        offset row in another state is never used: that would count the
        magnetic correction twice). Kepler/TESS h'1 and h'2 are applied jointly
        with their tabulated correlation; the star-to-star scatter is the
        width. route="hybrid5" uses the offsets of the hybrid5 fits' statistics
        (exact for conditioning weights). Returns (mean', cov')."""
        corr = {None: "none"}.get(prior.get("correction"), "MURaM_" + str(prior.get("correction")))
        stats = statistics or (["ell"] if regime.startswith("JWST") else ["h1p", "h2p"])
        rows = {r["Stat"]: r for r in self.offsets() if r["Regime"] == regime
                and r["Grid"] == prior["grid"] and r["Corr"] == corr and r["Route"] == route}
        miss = [k for k in stats if k not in rows]
        if miss:
            raise ValueError(f"no {route} offsets {miss} for {regime}, grid {prior['grid']}, "
                             f"correction state {corr}")
        delta = [float(rows[k]["Offset"]) for k in stats]
        sd = [float(rows[k]["Scatter"]) for k in stats]
        rho = rows[stats[0]]["rho"]
        return condition(prior["mean"], prior["cov"], prior["law"], [_STAT[k] for k in stats],
                         delta, sd, rho=float(rho) if (rho and len(stats) == 2) else None)

    intensity = staticmethod(intensity)
    observables = staticmethod(observables)
    condition = staticmethod(condition)

    # ---------------------------------------------------------- priors

    def prior(self, law, teff, sig_teff, logg, sig_logg, mh, sig_mh, band, grid="mps2",
              method="r", n=2000, systematic=False, seed=0, correction=None,
              m_dwarf_floor=True, cov_floor=None):
        """Monte-Carlo prior on the coefficients for uncertain stellar
        parameters: draws (Teff, log g, [M/H]) from independent Gaussians,
        interpolates the table at each draw and returns dict(mean, cov, n_used,
        dropped, grids, ...). Draws outside the grid are dropped (with a
        warning: near a grid edge this biases the mean).

        correction: MURaM state added as Δw(Teff) to each draw (hybrid laws).
        systematic=True adds the model systematic: the covariance, across the
        grids phoenix, mps1 and mps2 that cover the star in this band, of their
        interpolated coefficients at the nominal parameters. Where only one
        grid covers the star (Teff < 3500 K: PHOENIX only) and m_dwarf_floor is
        True, the spread at the nearest point shared by all three grids
        (Teff = 3500 K, log g and [M/H] clipped to the MPS-ATLAS range) is used
        instead, with a warning.
        cov_floor: minimum variance (scalar, or one per coefficient) on the
        diagonal, for the uncertainties the tables do not contain (limb
        convention, μ sampling, model-observation offsets); the interpolated
        covariance alone is nearly singular along directions the stellar
        parameters do not move. With σ → 0, no systematic and no floor, mean
        equals coefficients() and cov → 0."""
        rng = np.random.default_rng(seed)
        f, _ = self._interpolator(grid, law, band, method, None)
        T = rng.normal(teff, sig_teff, n)
        X = np.column_stack([T, rng.normal(logg, sig_logg, n), rng.normal(mh, sig_mh, n)]) / SCALE
        Y = f(X)[:, :len(COEFFS[law])]
        ok = np.all(np.isfinite(Y), axis=1)
        if not ok.any():
            raise ValueError("no draws inside the grid")
        dropped = 1.0 - ok.mean()
        if dropped > 0:
            warnings.warn(f"{100 * dropped:.1f}% of the draws fall outside the {grid} grid and "
                          "were dropped; the prior mean is biased towards the grid interior",
                          stacklevel=2)
        Y = Y[ok]
        if correction:
            Y = Y + self.delta_w(law, band, T[ok], correction)
        mean = Y.mean(0)
        cov = np.atleast_2d(np.cov(Y.T, ddof=1)) if len(Y) > 1 else np.zeros((Y.shape[1],) * 2)
        grids, floor_point = [grid], None
        if systematic:
            def spread(t, g, m, extrap=None):
                W, names = [], []
                for gname in ("phoenix", "mps1", "mps2"):
                    try:
                        r = self.coefficients(law, t, g, m, band, gname, method, extrap)
                    except (ValueError, FileNotFoundError):
                        continue
                    W.append([r["coeffs"][k] for k in COEFFS[law]])
                    names.append(gname)
                return np.array(W), names
            W, names = spread(teff, logg, mh)
            grids = [grid] + [g_ for g_ in names if g_ != grid]
            if len(W) < 2:
                msg = f"only {names} cover this star in {band}: no grid-to-grid systematic"
                if m_dwarf_floor:
                    floor_point = (max(teff, TEFF_SHARED), float(np.clip(logg, 3.5, 5.0)),
                                   float(np.clip(mh, -2.0, 1.0)))
                    W, names = spread(*floor_point, extrap="nearest")
                    msg += (f"; using the spread at the nearest shared point {floor_point} "
                            f"({names})")
                warnings.warn(msg, stacklevel=2)
            if len(W) > 1:
                cov = cov + np.atleast_2d(np.cov(W.T, ddof=1))
        if cov_floor is not None:
            fl = np.broadcast_to(np.asarray(cov_floor, float), (cov.shape[0],))
            cov = cov + np.diag(np.clip(fl - np.diag(cov), 0.0, None))
        return dict(law=law, mean=mean, cov=np.atleast_2d(cov), n_used=int(ok.sum()),
                    dropped=float(dropped), grids=grids, grid=grid, band=band, method=method,
                    correction=correction, floor_point=floor_point)

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

    def refit(self, law, teff, logg, mh, band, grid="mps2", mu_min=0.0):
        """Coefficients from interpolated *intensity profiles* rather than
        interpolated coefficients: the profiles of the Delaunay corner models
        are combined (PCHIP in μ, as in production) on the union of their nodes,
        combined with the barycentric weights, and the law is refitted with
        the production weighting and constraint. Exact at the models (to the
        float32 storage of the profiles); between them it removes the part of
        the interpolation error that comes from interpolating coefficients.
        mu_min > 0 leaves the profile at μ < mu_min out of the fit (e.g. 0.02,
        roughly the outer pressure scale height, which the four-term and
        power-2 laws cannot follow and which costs them accuracy elsewhere)."""
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
        cut = max([mu[0] for _, mu, _ in used] + [float(mu_min)])
        breaks = np.concatenate([mu[mu >= cut] for _, mu, _ in used] + [[cut]])
        mu_q, q = quadrature_r(breaks)
        y_q = sum(l * PchipInterpolator(mu, y)(mu_q) for l, mu, y in used)
        c = fit_profile(law, mu_q, y_q, q)
        coeffs = dict(zip(COEFFS[law], (float(v) for v in c)))
        return dict(law=law, coeffs=coeffs, physical=is_physical(law, coeffs), grid=grid,
                    method="r", band=band, refit=True, mu_min=float(mu_min),
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
    p.add_argument("--correction", default=None, choices=STATES,
                   help="add a MURaM magnetic correction (hybrid laws)")
    a = p.parse_args(argv)
    r = LDTables(a.tables).coefficients(a.law, a.teff, a.logg, a.mh, a.band, a.grid,
                                        a.method, "nearest" if a.nearest else None,
                                        correction=a.correction)
    print(" ".join(f"{k} = {v:.6f}" for k, v in r["coeffs"].items()))
    print(f"physical = {r['physical']}, mu_cri = {r['mucri']:.5f}, "
          f"fit LC error ≈ {r['lc_ppm']:.2f} ppm, inside grid = {r['inside']}, "
          f"correction = {r['correction']}")


if __name__ == "__main__":
    main()
