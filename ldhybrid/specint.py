"""Binned specific intensities and photon-weighted band integration.

A SpecInt holds, for one model atmosphere, the *photon* integral of the
specific intensity over wavelength bins,

    P[k, j] = ∫_{bin k} I_λ(λ, μ_j) λ dλ      (∝ photons s⁻¹ cm⁻² sr⁻¹),

on the model's native μ nodes. Bins are contiguous, so any throughput can be
applied bin by bin:

    I_band(μ) ∝ Σ_k P[k, μ] S̄_k ,

where S̄_k is the throughput averaged over bin k. Because the bins carry
integrals (not samples), coarser binning loses only the within-bin variation
of S, never flux.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class SpecInt:
    lam_edges: np.ndarray          # (nbin + 1,) Angstrom, increasing
    mu: np.ndarray                 # (nmu,) native nodes, increasing
    P: np.ndarray                  # (nbin, nmu) photon integrals per bin
    meta: dict = field(default_factory=dict)

    @property
    def lam_mid(self) -> np.ndarray:
        return 0.5 * (self.lam_edges[1:] + self.lam_edges[:-1])

    def bin_throughput(self, lam, S) -> np.ndarray:
        """Throughput (lam in Angstrom, S >= 0) averaged over each bin, by
        exact integration of its piecewise-linear interpolant; zero outside
        the curve."""
        lam = np.asarray(lam, float)
        S = np.asarray(S, float)
        o = np.argsort(lam)
        lam, S = lam[o], S[o]
        # cumulative integral of the linear interpolant, evaluated at the edges
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (S[1:] + S[:-1]) * np.diff(lam))])

        def C(x):
            x = np.clip(x, lam[0], lam[-1])
            i = np.clip(np.searchsorted(lam, x, side="right") - 1, 0, len(lam) - 2)
            t = x - lam[i]
            dl = lam[i + 1] - lam[i]
            slope = np.where(dl > 0, (S[i + 1] - S[i]) / np.where(dl > 0, dl, 1), 0.0)
            return cum[i] + S[i] * t + 0.5 * slope * t * t

        return np.diff(C(self.lam_edges)) / np.diff(self.lam_edges)

    def coverage(self, lam, S) -> float:
        """Fraction of ∫S dλ that falls inside the model's wavelength range."""
        lam = np.asarray(lam, float)
        S = np.asarray(S, float)
        inside = np.diff(self.lam_edges) @ self.bin_throughput(lam, S)
        total = np.trapezoid(S, lam)
        return float(inside / total) if total > 0 else 0.0

    def band(self, lam, S, min_coverage: float = 0.999) -> np.ndarray:
        """Photon-weighted band intensity on the native μ nodes (arbitrary
        units; normalise by the μ = 1 value for a limb-darkening profile).
        Raises if less than min_coverage of the throughput lies inside the
        model's wavelength range."""
        cov = self.coverage(lam, S)
        if cov < min_coverage:
            raise ValueError(f"only {cov:.4f} of the throughput lies inside the model's "
                             f"wavelength range {self.lam_edges[0]:.0f}–{self.lam_edges[-1]:.0f} Å")
        return self.bin_throughput(lam, S) @ self.P

    def save(self, path) -> None:
        np.savez_compressed(path, lam_edges=self.lam_edges, mu=self.mu, P=self.P,
                            meta=np.array(self.meta, dtype=object))

    @classmethod
    def load(cls, path) -> "SpecInt":
        with np.load(path, allow_pickle=True) as z:
            return cls(z["lam_edges"], z["mu"], z["P"], z["meta"].item())


def log_edges(lam_min: float, lam_max: float, R: float) -> np.ndarray:
    """Bin edges with constant λ/Δλ = R from lam_min, never beyond lam_max."""
    n = int(np.floor(np.log(lam_max / lam_min) * R))
    return lam_min * np.exp(np.arange(n + 1) / R)


def bin_photons(lam: np.ndarray, I: np.ndarray, edges: np.ndarray,
                photon_units: bool = False) -> np.ndarray:
    """Exact bin integrals of the piecewise-linear interpolant of I(λ)·λ
    (energy units), or of I(λ) itself when photon_units is True.

    lam: (n,) increasing sample wavelengths; I: (nmu, n); edges inside
    [lam[0], lam[-1]]. Returns (nbin, nmu)."""
    if edges[0] < lam[0] or edges[-1] > lam[-1]:
        raise ValueError("bin edges extend beyond the sampled wavelengths")
    f = I if photon_units else I * lam            # photon weighting
    dl = np.diff(lam)
    cum = np.concatenate([np.zeros((f.shape[0], 1)),
                          np.cumsum(0.5 * (f[:, 1:] + f[:, :-1]) * dl, axis=1)], axis=1)
    i = np.clip(np.searchsorted(lam, edges, side="right") - 1, 0, len(lam) - 2)
    t = (edges - lam[i]) / dl[i]
    f0, f1 = f[:, i], f[:, i + 1]
    # integral of the linear segment from lam[i] to the edge
    Ce = cum[:, i] + dl[i] * (f0 * t + 0.5 * (f1 - f0) * t * t)
    return np.diff(Ce, axis=1).T
