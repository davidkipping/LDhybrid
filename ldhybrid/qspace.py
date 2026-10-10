"""Sampler coordinates q ∈ [0, 1]^n for the hybrid weights (inverse maps of
squishierplanet.laws.hybrid2_from_q and simplex_from_q).

* hybrid2: w = √q1 [(1 − q2) V_c + q2 V_ℓ] on the exact triangle O–V_c–V_ℓ.
* hybrid4/5: stick-breaking of a flat Dirichlet on the simplex
  w_j = (1 − q_j^{1/(n+1−j)}) Π_{i<j} q_i^{1/(n+1−i)}.

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
    return np.array([s * s, t / s if s != 0 else 0.0])


def simplex_to_q(w):
    w = np.asarray(w, float)
    n = w.size
    q = np.empty(n)
    stick = 1.0
    for j in range(n):
        f = 1.0 - w[j] / stick if stick > 0 else 1.0
        q[j] = np.sign(f) * abs(f) ** (n - j) if f >= 0 else f
        stick = stick * f if f > 0 else 0.0
    return q


def to_q(law: str, w):
    if law == "hybrid2":
        return hybrid2_to_q(w)
    if law in ("hybrid4", "hybrid5"):
        return simplex_to_q(w)
    raise ValueError(law)
