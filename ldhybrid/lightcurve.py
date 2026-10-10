"""Circular-planet transit light curves of an arbitrary radial profile, for
scoring fits in light-curve space.

The blocked flux of a planet of radius p at separation z is

    ΔF(z) = ∫ I(r) · 2r·κ(r) dr,   κ(r) = arccos[(r² + z² − p²) / (2rz)] ∈ [0, π],

over r ∈ [max(0, z − p), min(1, z + p)]. The integrand has square-root
endpoints at r = |z − p|, z + p and at the limb, so each segment between
breakpoints is integrated by Gauss–Legendre in θ with r = a + (b − a)(1 − cos θ)/2,
which removes them. The model and the truth use the same quadrature.

    F(z) = 1 − ΔF(z) / F★,   F★ = 2π ∫₀¹ I(μ) μ dμ.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

N_GL = 48


@lru_cache(maxsize=None)
def _gl(n):
    x, w = np.polynomial.legendre.leggauss(n)
    th = 0.5 * np.pi * (x + 1.0)                 # θ ∈ (0, π)
    return th, 0.5 * np.pi * w


_TH, _WT = _gl(N_GL)


def star_flux(I_of_mu, n: int = 4000) -> float:
    """F★ = 2π ∫₀¹ I(μ) μ dμ, with cosine-substituted Gauss–Legendre in μ
    (dense, because tabulated profiles are only C¹)."""
    th, wt = _gl(n)
    mu = 0.5 * (1.0 - np.cos(th))
    return float(2.0 * np.pi * np.sum(wt * 0.5 * np.sin(th) * I_of_mu(mu) * mu))


def blocked(I_of_r, z, p, n_sub: int = 4) -> np.ndarray:
    """ΔF(z) for a vector of separations z (planet radius p)."""
    z = np.atleast_1d(np.asarray(z, float))
    lo = np.clip(z - p, 0.0, None)
    hi = np.minimum(1.0, z + p)
    bp = np.stack([lo, np.clip(p - z, 0.0, 1.0), np.clip(z - p, 0.0, 1.0),
                   np.minimum(1.0, z + p), hi], axis=1)
    bp = np.sort(np.clip(bp, lo[:, None], hi[:, None]), axis=1)
    # sub-divide each segment uniformly so tabulated-profile kinks are resolved
    edges = [bp[:, :1]]
    for j in range(bp.shape[1] - 1):
        a, b = bp[:, j:j + 1], bp[:, j + 1:j + 2]
        edges.append(a + (b - a) * np.arange(1, n_sub + 1) / n_sub)
    E = np.concatenate(edges, axis=1)                     # (nz, nseg + 1)
    a, b = E[:, :-1, None], E[:, 1:, None]
    r = a + (b - a) * 0.5 * (1.0 - np.cos(_TH))           # (nz, nseg, n)
    dr = (b - a) * 0.5 * np.sin(_TH) * _WT
    zz = z[:, None, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        cosk = (r * r + zz * zz - p * p) / (2.0 * r * zz)
    kap = np.where(zz == 0, np.where(r < p, np.pi, 0.0),
                   np.arccos(np.clip(np.nan_to_num(cosk, nan=-1.0), -1.0, 1.0)))
    return np.sum(I_of_r(r) * 2.0 * r * kap * dr, axis=(1, 2))


def light_curve(I_of_mu, z, p, F_star: float | None = None) -> np.ndarray:
    """Normalised flux F(z) = 1 − ΔF/F★ for a profile given as I(μ)."""
    F_star = star_flux(I_of_mu) if F_star is None else F_star
    I_of_r = lambda r: I_of_mu(np.sqrt(np.clip(1.0 - r * r, 0.0, None)))
    return 1.0 - blocked(I_of_r, z, p) / F_star


# ------------------------------------------------------------ scoring grid

P_GRID = (0.03, 0.1, 0.15)
B_GRID = tuple(np.round(np.arange(0.0, 0.91, 0.1), 2))
N_CHORD = 120


def chord_z(p: float, b: float, n: int = N_CHORD) -> np.ndarray:
    """Separations along half a transit chord (symmetric), uniform in time
    from mid-transit to last contact."""
    xmax = np.sqrt(max((1.0 + p) ** 2 - b * b, 0.0))
    x = np.linspace(0.0, xmax, n)
    return np.sqrt(x * x + b * b)


def lc_errors(model_of_mu, truth_of_mu, p_grid=P_GRID, b_grid=B_GRID,
              n: int = N_CHORD) -> dict:
    """Light-curve error of a model profile against a truth profile (ppm)
    over the (p, b) grid: worst |ΔF| per (p, b), and the rms along the chord."""
    Fm, Ft = star_flux(model_of_mu), star_flux(truth_of_mu)
    worst = np.zeros((len(p_grid), len(b_grid)))
    rms = np.zeros_like(worst)
    for i, p in enumerate(p_grid):
        for j, b in enumerate(b_grid):
            z = chord_z(p, b, n)
            d = (light_curve(model_of_mu, z, p, Fm) - light_curve(truth_of_mu, z, p, Ft)) * 1e6
            worst[i, j] = np.max(np.abs(d))
            rms[i, j] = np.sqrt(np.mean(d * d))
    return dict(worst=worst, rms=rms, max_ppm=float(worst.max()),
                mean_rms_ppm=float(rms.mean()), p_grid=p_grid, b_grid=b_grid)


class LCGrid:
    """The (p, b) scoring grid flattened once, so that the light curve of any
    profile is one vectorised quadrature and a linear law's light curve is a
    combination of precomputed shape columns.

        g = LCGrid()
        dF_t, Fs_t = g.columns(truth)                 # truth, once per profile
        dF_j, Fs_j = g.columns(shape_j)               # per shape, once per law
        g.errors(dF_model, Fs_model, dF_t, Fs_t)      # same dict as lc_errors
    """

    def __init__(self, p_grid=P_GRID, b_grid=B_GRID, n: int = N_CHORD):
        self.p_grid, self.b_grid, self.n = tuple(p_grid), tuple(b_grid), n
        self.z = [[chord_z(p, b, n) for b in b_grid] for p in p_grid]

    def columns(self, I_of_mu):
        """Blocked flux on the whole grid, shape (np, nb, n), and F★."""
        I_of_r = lambda r: I_of_mu(np.sqrt(np.clip(1.0 - r * r, 0.0, None)))
        dF = np.array([[blocked(I_of_r, z, p) for z in row]
                       for p, row in zip(self.p_grid, self.z)])
        return dF, star_flux(I_of_mu)

    def errors(self, dF_m, Fs_m, dF_t, Fs_t) -> dict:
        d = (dF_t / Fs_t - dF_m / Fs_m) * 1e6           # F_m − F_t in ppm
        worst = np.max(np.abs(d), axis=-1)
        rms = np.sqrt(np.mean(d * d, axis=-1))
        return dict(worst=worst, rms=rms, max_ppm=float(worst.max()),
                    mean_rms_ppm=float(rms.mean()), p_grid=self.p_grid, b_grid=self.b_grid)

    def linear_columns(self, law):
        """Columns for I = 1 − Σ w_j T_j: (dF of 1, F★ of 1, dF_j, F★_j)."""
        one = self.columns(lambda mu: np.ones_like(mu))
        cols = [self.columns(lambda mu, j=j: law.T(mu)[..., j]) for j in range(law.m)]
        return one, cols

    @staticmethod
    def combine(one, cols, w):
        dF = one[0] - sum(wj * c[0] for wj, c in zip(w, cols))
        Fs = one[1] - sum(wj * c[1] for wj, c in zip(w, cols))
        return dF, Fs
