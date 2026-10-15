"""Sampler coordinates q ∈ [0, 1]^n for the hybrid weights; a uniform q gives
a uniform prior over the law's region.

* hybrid2: w = √q1 [(1 − q2) V_c + q2 V_ℓ] on the exact triangle O–V_c–V_ℓ
  (inverse of squishierplanet.laws.hybrid2_from_q).
* hybrid4 (curvature-simplex basis, ldhybrid.v1.HYBRID4_BASIS): the analogue
  of Kipping's (2013) q1, q2. The depth s = Σw = 1 − I(0) = q1^{1/4}, and
  the direction on the face Σw = s by stick-breaking of a flat Dirichlet:
      w1 = s (1 − q2^{1/3}),       w2 = s q2^{1/3} (1 − q3^{1/2}),
      w3 = s q2^{1/3} q3^{1/2} (1 − q4),   w4 = s q2^{1/3} q3^{1/2} q4.
  Every q in the unit hypercube gives a profile that is positive, monotonic,
  concave in μ² and of convex slope in μ², and every such profile has a q.
* hybrid5: stick-breaking of a flat Dirichlet on the w-simplex,
  w_j = (1 − q_j^{1/(n+1−j)}) Π_{i<j} q_i^{1/(n+1−i)}
  (inverse of squishierplanet.laws.simplex_from_q).

Points outside the region return q outside [0, 1] (no clipping), so a table
user can see it.
"""

from __future__ import annotations

import numpy as np
from squishierplanet.laws import hybrid2_vertices

from . import v1


def hybrid2_to_q(w, eps: float = v1.HYBRID2_EPS):
    Vc, Vl = hybrid2_vertices(eps)
    M = np.column_stack([Vc, Vl - Vc])                 # w = s Vc + t (Vl − Vc)
    s, t = np.linalg.solve(M, np.asarray(w, float))
    # q1 = s|s|, not s², so that a negative depth (limb brighter than centre)
    # gives q1 < 0 rather than mirroring into the cube
    return _tidy([_spow(s, 2), t / s if s != 0 else 0.0])


def _tidy(q, tol=1e-9):
    """Rounding noise within tol of 0 or 1 is put on the bound (a point on a face
    of the region has q exactly 0 or 1); larger excursions are kept, so points
    outside the region still show it."""
    q = np.asarray(q, float)
    q = np.where((q < 0) & (q > -tol), 0.0, q)
    return np.where((q > 1) & (q < 1 + tol), 1.0, q)


def simplex_to_q(w):
    w = np.asarray(w, float)
    n = w.size
    q = np.empty(n)
    stick = 1.0
    for j in range(n):
        f = 1.0 - w[j] / stick if stick > 0 else 1.0
        q[j] = np.sign(f) * abs(f) ** (n - j) if f >= 0 else f
        stick = stick * f if f > 0 else 0.0
    return _tidy(q)


def _spow(x, p):
    return np.sign(x) * np.abs(x) ** p


def hybrid4_from_q(q):
    """Weights (curvature-simplex basis) from q ∈ [0, 1]^4; vectorised over rows."""
    q = np.asarray(q, float)
    s = q[..., 0] ** 0.25
    a, b = q[..., 1] ** (1 / 3), q[..., 2] ** 0.5
    return np.stack([s * (1 - a), s * a * (1 - b), s * a * b * (1 - q[..., 3]),
                     s * a * b * q[..., 3]], axis=-1)


def hybrid4_to_q(w):
    """Inverse of hybrid4_from_q; vectorised over rows. At s = 0 the direction
    is undefined and q2..q4 are set to 0."""
    w = np.asarray(w, float)
    s = w.sum(-1)
    safe = np.where(s != 0, s, 1.0)
    f = w / safe[..., None]
    r1 = 1 - f[..., 0]                                  # = q2^{1/3}
    safe1 = np.where(r1 != 0, r1, 1.0)
    r2 = 1 - f[..., 1] / safe1                          # = q3^{1/2}
    f34 = f[..., 2] + f[..., 3]
    q = np.stack([_spow(s, 4), _spow(r1, 3), np.where(r2 >= 0, r2 * r2, r2),
                  np.where(f34 != 0, f[..., 3] / np.where(f34 != 0, f34, 1.0), 0.0)], axis=-1)
    q[..., 2] = np.where(r1 != 0, q[..., 2], 0.0)
    q[..., 1:] = np.where((s != 0)[..., None], q[..., 1:], 0.0)
    return _tidy(q)


def to_q(law: str, w):
    if law == "hybrid2":
        return hybrid2_to_q(w)
    if law == "hybrid4":
        return hybrid4_to_q(w)
    if law == "hybrid5":
        return simplex_to_q(w)
    raise ValueError(law)
