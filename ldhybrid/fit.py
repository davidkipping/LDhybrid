"""Fitting engine: weighted least squares in the shape basis, with linear
constraints, as a convex QP.

    minimise  Σ_i q_i (y_i − 1 + Σ_j w_j T_j(μ_i))²

Constraints (all linear in w):
  "none"      unconstrained;
  "simplex"   w ≥ 0, Σ w ≤ 1 (hybrid4/5's sufficient region);
  "physical"  the exact region of conditions I–II: Σ w ≤ 1 and
              −Σ_j w_j T_j'(x) ≥ 0 on a dense grid of the law's variable
              (ξ = μ² for hybrids, μ for controls); for hybrid2 this is the
              exact triangle;
  fcm=True    adds flux conservation, ∫ I μ dμ (law) = ∫ I μ dμ (profile),
              over the profile's tabulated range.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import cvxpy as cp
import numpy as np
from scipy.optimize import least_squares

from .laws import LinearLaw, power2
from .profile import Profile, quadrature

N_CONSTRAINT = 2001          # dense grid for the monotonicity constraint
TOL_ACTIVE = 1e-9


@dataclass
class Fit:
    law: str
    w: np.ndarray
    measure: str
    constraint: str
    fcm: bool
    rms: float                      # weighted rms intensity residual
    max_abs: float                  # max |Δy| on the tabulated range
    flux_err: float                 # relative error in ∫ I μ dμ
    active: list = field(default_factory=list)
    extra: dict = field(default_factory=dict)


def _flux_terms(law: LinearLaw | None, profile: Profile, n: int = 4001):
    """∫ T_j μ dμ (per shape) and ∫ (1 − y) μ dμ over [μ_min, 1], by
    Simpson on a dense grid."""
    from scipy.integrate import simpson
    mu = np.linspace(profile.mu[0], 1.0, n)
    y = profile.interpolant()(mu)
    data = simpson((1.0 - y) * mu, x=mu)
    if law is None:
        return None, data
    return simpson(law.T(mu) * mu[:, None], x=mu, axis=0), data


def _diagnostics(law_eval, profile, mu_q, y_q, q):
    r = law_eval(mu_q) - y_q
    rms = float(np.sqrt(np.sum(q * r * r) / np.sum(q)))
    mu = np.linspace(profile.mu[0], 1.0, 4001)
    max_abs = float(np.max(np.abs(law_eval(mu) - profile.interpolant()(mu))))
    from scipy.integrate import simpson
    F_law = simpson(law_eval(mu) * mu, x=mu)
    F_dat = simpson(profile.interpolant()(mu) * mu, x=mu)
    return rms, max_abs, float(F_law / F_dat - 1.0)


def fit_linear(law: LinearLaw, profile: Profile, measure: str = "r",
               resample: int | str | None = "exact", constraint: str = "none",
               fcm: bool = False, mu_min: float = 0.0) -> Fit:
    mu_q, y_q, q = quadrature(profile, measure, resample, mu_min)
    A = law.T(mu_q)
    b = 1.0 - y_q
    sq = np.sqrt(q)
    As, bs = A * sq[:, None], b * sq
    scale = np.linalg.norm(As, axis=0)
    scale[scale == 0] = 1.0

    if constraint == "none" and not fcm:
        w = np.linalg.lstsq(As / scale, bs, rcond=None)[0] / scale
        active = []
    else:
        C, d = _inequalities(law, constraint)
        E = e = None
        if fcm:
            Tm, data = _flux_terms(law, profile)
            E, e = Tm[None, :], np.array([data])
        # exact shortcut: if the (equality-constrained) least-squares optimum
        # already satisfies the inequalities, it is the constrained optimum
        w = _eq_lstsq(As / scale, bs, (E / scale) if E is not None else np.zeros((0, law.m)),
                      e if e is not None else np.zeros(0)) / scale
        if len(d) and np.any(C @ w > d + 1e-12):
            w = _qp(As, bs, scale, C, d, E, e)
        active = _active(law, w, constraint)
        if fcm:
            active.append("fcm")

    rms, max_abs, ferr = _diagnostics(lambda m: law.intensity(w, m), profile, mu_q, y_q, q)
    return Fit(law.name, w, measure, constraint, fcm, rms, max_abs, ferr, active,
               extra=dict(physical=is_physical(law, w)))


def _inequalities(law: LinearLaw, constraint: str):
    """Rows of C w <= d, each row normalised to unit max-abs."""
    if constraint == "none":
        return np.zeros((0, law.m)), np.zeros(0)
    rows, rhs = [np.ones(law.m)], [1.0]                     # limb intensity >= 0
    if constraint == "simplex":
        rows += list(-np.eye(law.m))
        rhs += [0.0] * law.m
    elif constraint == "physical":
        lo = 1e-6 if law.variable == "mu" else 0.0
        G = law.dT(np.linspace(lo, 1.0, N_CONSTRAINT))      # dI/dx = −G w ≥ 0
        rows += list(G)
        rhs += [0.0] * len(G)
    else:
        raise ValueError(f"unknown constraint {constraint!r}")
    C, d = np.array(rows), np.array(rhs)
    n = np.max(np.abs(C), axis=1)
    n[n == 0] = 1.0
    return C / n[:, None], d / n


def _eq_lstsq(R, g, M, m):
    """min ||R x − g|| subject to M x = m (M may be rank-deficient)."""
    if M.shape[0] == 0:
        return np.linalg.lstsq(R, g, rcond=None)[0]
    U, S, Vt = np.linalg.svd(M)
    rank = int(np.sum(S > 1e-10 * S[0]))
    x0 = Vt[:rank].T @ ((U[:, :rank].T @ m) / S[:rank])
    N = Vt[rank:].T
    if N.shape[1] == 0:
        return x0
    z = np.linalg.lstsq(R @ N, g - R @ x0, rcond=None)[0]
    return x0 + N @ z


def _qp(As, bs, scale, C, d, E, e, tol=1e-8):
    """min ||As w − bs||² s.t. C w ≤ d, E w = e. Clarabel on the scaled,
    QR-reduced problem, then an exact polish: an equality-constrained least
    squares on the active set, kept if it is feasible and no worse."""
    Q, R = np.linalg.qr(As / scale)
    g = Q.T @ bs
    Cs = C / scale
    v = cp.Variable(As.shape[1])
    cons = [Cs @ v <= d] if len(d) else []
    if E is not None:
        cons.append((E / scale) @ v == e)
    prob = cp.Problem(cp.Minimize(cp.sum_squares(R @ v - g)), cons)
    prob.solve(solver=cp.CLARABEL, tol_gap_abs=1e-12, tol_gap_rel=1e-12,
               tol_feas=1e-12, max_iter=500)
    if prob.status not in ("optimal", "optimal_inaccurate"):
        raise RuntimeError(f"QP failed: {prob.status}")
    v0 = np.asarray(v.value)
    obj0 = np.sum((R @ v0 - g) ** 2)
    act = (d - Cs @ v0) <= tol if len(d) else np.zeros(0, bool)
    M = np.vstack([Cs[act]] + ([E / scale] if E is not None else []))
    m = np.concatenate([d[act]] + ([e] if E is not None else []))
    v1 = _eq_lstsq(R, g, M, m)
    feasible = (not len(d) or np.all(Cs @ v1 <= d + 1e-12))
    if feasible and np.sum((R @ v1 - g) ** 2) <= obj0 * (1 + 1e-9) + 1e-30:
        v0 = v1
    return v0 / scale


def is_physical(law: LinearLaw, w, tol: float = 1e-10) -> bool:
    """Conditions I–II on a dense grid of the law's variable."""
    lo = 1e-6 if law.variable == "mu" else 0.0
    x = np.linspace(lo, 1.0, N_CONSTRAINT)
    return bool(1.0 - np.sum(w) >= -tol and np.all(-law.dT(x) @ np.asarray(w) >= -tol))


def _active(law, w, constraint):
    act = []
    if constraint in ("simplex", "physical") and 1.0 - np.sum(w) < 1e-7:
        act.append("limb=0")
    if constraint == "simplex":
        act += [f"w{j + 1}=0" for j in np.flatnonzero(w < 1e-8)]
    if constraint == "physical":
        lo = 1e-6 if law.variable == "mu" else 0.0
        x = np.linspace(lo, 1.0, N_CONSTRAINT)
        s = -law.dT(x) @ w
        if np.min(s) < 1e-7 * max(1.0, np.max(np.abs(s))):
            act.append(f"slope=0@{law.variable}={x[np.argmin(s)]:.3g}")
    return act


def fit_power2(profile: Profile, measure: str = "r", resample: int | str | None = "exact",
               x0=(0.7, 0.6), mu_min: float = 0.0) -> Fit:
    """Power-2 law by weighted Levenberg–Marquardt (as Claret & Southworth
    2022/2023), with the same quadrature as the linear laws."""
    mu_q, y_q, q = quadrature(profile, measure, resample, mu_min)
    sq = np.sqrt(q)
    res = least_squares(lambda p: sq * (power2(p[0], p[1], mu_q) - y_q), x0,
                        method="lm", xtol=1e-14, ftol=1e-14)
    c, a = res.x
    rms, max_abs, ferr = _diagnostics(lambda m: power2(c, a, m), profile, mu_q, y_q, q)
    return Fit("power2", np.array([c, a]), measure, "none", False, rms, max_abs, ferr,
               extra=dict(success=bool(res.success)))
