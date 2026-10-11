"""Band intensity profiles, limb treatment and fitting quadrature.

A Profile is a normalised centre-to-limb curve y(μ) = I(μ)/I(1) on nodes.

Spherical models (PHOENIX) carry an extended, optically thin outer layer:
their μ is measured from the outermost model shell, and the intensity
collapses well inside it. limb_cut() puts the stellar radius at the steepest
intensity drop (Wittkowski et al. 2004; Espinoza & Jordán 2015's r-method;
Claret & Southworth 2023's M1/M2) and rescales r' = r/r_c,
μ' = √(1 − r'²):

  * "M1"   r_c at the node of max |dI/dr|            (Claret & Southworth M1)
  * "M2"   r_c midway to the next node inward        (Claret & Southworth M2;
           their tabulated μ_cri(M2) > μ_cri(M1), i.e. the cut moves inward)
  * "peak" r_c at the parabolic refinement of the max (squishierplanet)

The cut is evaluated on the continuous (PCHIP-in-r) profile, so the rescaled
profile always reaches μ' = 0, with I(r_c) interpolated.

Plane-parallel models need no cut; their profile stops at the smallest
tabulated μ (no extrapolation is used in the fits).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import PchipInterpolator


@dataclass
class Profile:
    mu: np.ndarray                # increasing, last = 1
    y: np.ndarray                 # I/I(1)
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        o = np.argsort(self.mu)
        self.mu, self.y = np.asarray(self.mu, float)[o], np.asarray(self.y, float)[o]
        if not np.isclose(self.mu[-1], 1.0):
            raise ValueError("profile must include μ = 1")
        self.y = self.y / self.y[-1]

    def interpolant(self):
        """PCHIP in μ (monotone, C¹) on the tabulated range."""
        return PchipInterpolator(self.mu, self.y, extrapolate=False)

    def __call__(self, mu):
        """y(μ); below the smallest node, linear continuation to μ = 0,
        clipped at 0 (used only for light-curve truths, never for fits)."""
        mu = np.asarray(mu, float)
        f = self.interpolant()
        out = f(np.clip(mu, self.mu[0], 1.0))
        lo = mu < self.mu[0]
        if np.any(lo):
            s = (self.y[1] - self.y[0]) / (self.mu[1] - self.mu[0])
            out = np.where(lo, np.clip(self.y[0] + s * (mu - self.mu[0]), 0.0, None), out)
        return out


def limb_cut(mu, I, method: str = "M1", r_c: float | None = None):
    """Rescale a spherical-model profile at the limb. Returns (Profile, r_c).
    Pass r_c to impose a radius found elsewhere (e.g. one per model)."""
    mu = np.asarray(mu, float)
    I = np.asarray(I, float)
    o = np.argsort(mu)[::-1]                 # increasing r
    r = np.sqrt(np.clip(1.0 - mu[o] ** 2, 0.0, None))
    Ir = I[o]
    if r_c is None:
        dIdr = np.gradient(Ir, r)
        k = int(np.argmax(np.abs(dIdr)))
        if method == "M1":
            r_c = r[k]
        elif method == "M2":
            r_c = 0.5 * (r[k] + r[max(k - 1, 0)])
        elif method == "peak":
            r_c = r[k]
            if 0 < k < len(r) - 1:
                y0, y1, y2 = np.abs(dIdr[k - 1:k + 2])
                x0, x1, x2 = r[k - 1:k + 2]
                den = (x0 - x1) * (x0 - x2) * (x1 - x2)
                A = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / den
                B = (x2 ** 2 * (y0 - y1) + x1 ** 2 * (y2 - y0) + x0 ** 2 * (y1 - y2)) / den
                if A < 0:
                    r_c = float(np.clip(-B / (2 * A), r[k - 1], r[k + 1]))
        else:
            raise ValueError(f"unknown limb method {method!r}")
    fr = PchipInterpolator(r, Ir)
    keep = r < r_c
    rk = np.append(r[keep], r_c)
    Ik = np.append(Ir[keep], fr(r_c))
    mup = np.sqrt(np.clip(1.0 - (rk / r_c) ** 2, 0.0, None))
    return Profile(mup, Ik, dict(limb=method, r_c=float(r_c))), float(r_c)


def band_profile(spec, lam, S, geometry: str | None = None, limb: str = "M1",
                 r_c: float | None = None, meta: dict | None = None) -> Profile:
    """Photon-weighted band profile from a SpecInt, with the limb treatment
    appropriate to its geometry."""
    I = spec.band(lam, S)
    geometry = geometry or spec.meta.get("geometry", "plane-parallel")
    m = dict(spec.meta, **(meta or {}))
    if geometry == "spherical":
        p, rc = limb_cut(spec.mu, I, limb, r_c)
        p.meta = dict(m, **p.meta)
        return p
    return Profile(spec.mu, I, dict(m, limb="none", r_c=1.0))


# ---------------------------------------------------------------- quadrature

MEASURES = ("mu", "r", "area", "nodes")


def measure_coordinate(mu, measure):
    mu = np.asarray(mu, float)
    if measure == "mu":
        return mu
    if measure == "r":
        return np.sqrt(np.clip(1.0 - mu * mu, 0.0, None))
    if measure == "area":
        return 1.0 - mu * mu                  # r²
    raise ValueError(f"unknown measure {measure!r}")


N_EXACT = 12          # Gauss–Legendre points per interval between profile nodes


def _exact_quadrature(profile: Profile, measure: str):
    """Gauss–Legendre on every interval between the profile's nodes, so that
    ∫ (·)² dx over the PCHIP profile is integrated to rounding however steep
    the limb (the PCHIP is cubic in μ on each interval). For x = r the
    variable is μ (dr = μ/r dμ) on intervals with μ ≤ 0.5, where r → 1, and r
    on the others, where μ → 1; for x = μ and x = r² it is μ throughout."""
    t, wt = np.polynomial.legendre.leggauss(N_EXACT)
    nodes = profile.mu
    mus, qs = [], []
    for a, b in zip(nodes[:-1], nodes[1:]):
        if measure == "r" and a >= 0.5:
            ra, rb = np.sqrt(1 - b * b), np.sqrt(1 - a * a)        # r increasing
            r = 0.5 * (rb - ra) * t + 0.5 * (ra + rb)
            mus.append(np.sqrt(np.clip(1 - r * r, 0.0, None)))
            qs.append(0.5 * (rb - ra) * wt)
            continue
        mu = 0.5 * (b - a) * t + 0.5 * (a + b)
        jac = {"mu": np.ones_like(mu), "area": 2.0 * mu,
               "r": mu / np.sqrt(np.clip(1 - mu * mu, 1e-300, None))}[measure]
        mus.append(mu)
        qs.append(0.5 * (b - a) * wt * jac)
    mu = np.concatenate(mus)
    return mu, profile.interpolant()(mu), np.concatenate(qs)


def quadrature(profile: Profile, measure: str = "r", resample: int | str | None = 2000,
               mu_min: float = 0.0):
    """Nodes μ_i, targets y_i and weights q_i for ∫ (·)² dx, x ∈ {μ, r, r²},
    over the tabulated μ range.

    resample=None: trapezoid weights on the native nodes (ExoTETHyS-style).
    resample=n:    n nodes uniform in x, profile by PCHIP in μ.
    measure="nodes": equal weights on the native nodes (Claret 2000,
    Kostogryz et al. 2022), whatever resample is.
    mu_min: fit only μ ≥ mu_min (ExoTiC-LD uses 0.1; Sing et al. 2026 0.2)."""
    if mu_min > profile.mu[0]:
        keep = profile.mu >= mu_min
        profile = Profile(profile.mu[keep], profile.y[keep], profile.meta)
    if measure == "nodes":
        return profile.mu, profile.y, np.ones_like(profile.mu)
    if resample == "exact":
        return _exact_quadrature(profile, measure)
    if resample is None:
        mu, y = profile.mu, profile.y
        x = measure_coordinate(mu, measure)
        o = np.argsort(x)
        xs = x[o]
        q = np.zeros_like(xs)
        dx = np.diff(xs)
        q[:-1] += 0.5 * dx
        q[1:] += 0.5 * dx
        w = np.empty_like(q)
        w[o] = q
        return mu, y, w
    x0, x1 = sorted(measure_coordinate(np.array([profile.mu[0], 1.0]), measure))
    x = np.linspace(x0, x1, resample)
    if measure == "mu":
        mu = x
    elif measure == "r":
        mu = np.sqrt(np.clip(1.0 - x * x, 0.0, None))
    else:
        mu = np.sqrt(np.clip(1.0 - x, 0.0, None))
    mu = np.clip(mu, profile.mu[0], 1.0)
    q = np.full(resample, (x1 - x0) / (resample - 1))
    q[[0, -1]] *= 0.5
    return mu, profile.interpolant()(mu), q
