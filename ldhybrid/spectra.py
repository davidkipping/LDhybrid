"""Coefficient spectra: hybrid-law coefficients per wavelength bin that
integrate *exactly* to the coefficients of any passband.

For each PHOENIX model (v1: one stellar edge r_c per model, from the
photon-bolometric profile, M1), every reduced bin k (R = 3000) has photon
intensities I_k(μ_i) on the native nodes μ_i inside the edge. We fit the
*unnormalised* shape form

    I_k(μ) ≈ I_k(1) − Σ_j a_kj T_j(μ)

by r-weighted least squares with trapezoid weights on the native nodes. The
fit is a fixed linear map of the data, so for a passband with bin-averaged
throughput S̄_k,

    w_j(band) = Σ_k S̄_k a_kj / Σ_k S̄_k I_k(1)

reproduces the band fit to rounding. The per-bin fits are unconstrained (single
bins can be unphysical inside strong lines); physicality is a property of the
band result, which band_coefficients() reports.

Differences from the band tables (PCHIP resample, r_c per band) are a few ppm
in light-curve terms (P2d: one r_c per model 2–6 ppm; native vs resample ≲ 2).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import fit, laws
from .profile import limb_cut, quadrature
from .specint import SpecInt

SPECTRA_LAWS = ("hybrid2", "hybrid4", "hybrid5")


@dataclass
class CoeffSpectrum:
    lam_edges: np.ndarray            # (nbin + 1,) Å
    I1: np.ndarray                   # (nbin,) photon intensity at μ = 1 per bin
    a: dict                          # law -> (nbin, m)
    r_c: float
    meta: dict

    def band_coefficients(self, lam, S, law: str):
        """Weights of `law` for a passband (λ in Å, photon response S)."""
        w_bin = SpecInt(self.lam_edges, np.array([0.0, 1.0]),
                        np.zeros((len(self.I1), 2))).bin_throughput(lam, S)
        w = (w_bin @ self.a[law]) / (w_bin @ self.I1)
        return w, fit.is_physical(laws.get(law), w)

    def rebin(self, factor: int) -> "CoeffSpectrum":
        """Merge `factor` adjacent bins (exact: the a_kj and I_k(1) are bin
        integrals)."""
        n = (len(self.I1) // factor) * factor
        s = lambda x: x[:n].reshape(n // factor, factor, *x.shape[1:]).sum(1)
        return CoeffSpectrum(self.lam_edges[: n + 1 : factor], s(self.I1),
                             {k: s(v) for k, v in self.a.items()}, self.r_c, dict(self.meta))


def model_spectrum(spec: SpecInt, measure: str = "r") -> CoeffSpectrum:
    lam_all = np.array([spec.lam_edges[0], spec.lam_edges[-1]])
    bol = spec.band(lam_all, np.ones(2))
    prof, r_c = limb_cut(spec.mu, bol, "M1")
    # nodes kept inside the edge, in the order limb_cut returns (μ' increasing)
    o = np.argsort(spec.mu)[::-1]                          # increasing r
    r = np.sqrt(np.clip(1.0 - spec.mu[o] ** 2, 0.0, None))
    keep = o[r <= r_c + 1e-15]                             # M1: r_c is a node
    mup = np.sqrt(np.clip(1.0 - (r[r <= r_c + 1e-15] / r_c) ** 2, 0.0, None))
    order = np.argsort(mup)
    nodes, mup = keep[order], mup[order]
    _, _, q = quadrature(type(prof)(mup, np.ones_like(mup)), measure, None)
    P = spec.P.astype(np.float64)[:, nodes]                # (nbin, nnode)
    I1 = P[:, -1]                                          # μ' = 1 is the last node
    D = I1[:, None] - P                                    # data: I(1) − I(μ)
    a = {}
    for name in SPECTRA_LAWS:
        A = laws.get(name).T(mup)                          # (nnode, m)
        Aw = A * q[:, None]
        Mmap = np.linalg.solve(A.T @ Aw, Aw.T)             # (m, nnode)
        a[name] = D @ Mmap.T                               # (nbin, m)
    return CoeffSpectrum(spec.lam_edges, I1, a, float(r_c),
                         dict(spec.meta, nodes=nodes.tolist(), measure=measure))


def band_fit_same_settings(spec: SpecInt, lam, S, name: str, r_c: float, measure: str = "r"):
    """Reference: the direct band fit with the spectra's settings (one r_c,
    native-node trapezoid, unconstrained) — for the closure test."""
    prof, _ = limb_cut(spec.mu, spec.band(lam, S), "M1", r_c=r_c)
    return fit.fit_linear(laws.get(name), prof, measure, None, "none").w
