"""Limb-darkening laws in shape form, I(μ)/I(1) = 1 − Σ_j w_j T_j(μ).

Each shape T_j equals 1 at the limb and 0 at the centre, so w_j is the share
of the centre-to-limb drop its term carries and I(0)/I(1) = 1 − Σ w.

* The hybrids take their shape functions from squishierplanet.laws, with the
  pole positions frozen in ldhybrid.v1 (hybrid2: ε = 0.36, re-tuned in P2c;
  squishierplanet still carries the J-band 0.208 until it is updated).
* The control laws (linear, quadratic, Claret 4-term) are linear too and use
  the same engine. Power-2 is nonlinear and handled in fit.py.

Physicality (conditions I and II) is linear in w. With x the law's natural
variable (ξ = μ² for the hybrids, μ for the controls), monotonicity is
dI/dx = −Σ_j w_j T_j'(x) ≥ 0, and positivity then reduces to 1 − Σ w ≥ 0.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from squishierplanet import laws as sl

from . import v1


def _law_value(L: sl.Law, xi):
    return L.intensity(np.sqrt(np.clip(xi, 0.0, None)))


def _law_dxi(L: sl.Law, xi):
    xi = np.asarray(xi, float)
    n = np.arange(1, L.c.size)
    out = (np.polynomial.polynomial.polyval(xi, n * L.c[1:]) if L.c.size > 1
           else np.zeros_like(xi))
    for dk, ek in zip(L.d, L.eps):
        out = out - 2.0 * dk / (xi + ek) ** 3
    return out


@dataclass(frozen=True)
class LinearLaw:
    name: str
    labels: tuple[str, ...]
    variable: str                                  # "xi" (= μ²) or "mu"
    _T: Callable                                    # μ -> (..., m)
    _dT: Callable                                   # x -> (..., m), d/dx
    shapes: tuple = ()                              # squishierplanet Laws (hybrids)

    @property
    def m(self) -> int:
        return len(self.labels)

    def T(self, mu) -> np.ndarray:
        return self._T(np.asarray(mu, float))

    def dT(self, x) -> np.ndarray:
        return self._dT(np.asarray(x, float))

    def intensity(self, w, mu) -> np.ndarray:
        return 1.0 - self.T(mu) @ np.asarray(w, float)

    def to_law(self, w) -> sl.Law:
        """squishierplanet Law (for exact light curves); hybrids only."""
        if not self.shapes:
            raise TypeError(f"{self.name} is not in the elementary class")
        law = sl.Law([1.0])
        for wj, Tj in zip(np.asarray(w, float), self.shapes):
            law = law + (-wj) * Tj
        return law


def _from_shapes(name, labels, shapes):
    shapes = tuple(shapes)
    return LinearLaw(
        name, tuple(labels), "xi",
        lambda mu: np.stack([_law_value(S, mu * mu) for S in shapes], axis=-1),
        lambda xi: np.stack([_law_dxi(S, xi) for S in shapes], axis=-1),
        shapes)


def hybrid2(eps: float = v1.HYBRID2_EPS) -> LinearLaw:
    return _from_shapes(f"hybrid2(eps={eps:g})" if eps != v1.HYBRID2_EPS else "hybrid2",
                        ("w1:1-mu2", f"w2:pole{eps:g}"),
                        [sl.Law([1.0, -1.0]), sl.pole_shape(eps)])


def _monomial_pole(eps):
    return [sl.Law([1.0, 0.0, -1.0]), sl.Law([1.0, -2.0, 1.0])] + [sl.pole_shape(e) for e in eps]


def hybrid(name: str) -> LinearLaw:
    if name == "hybrid2":
        return hybrid2()
    if name == "hybrid4":
        # curvature-simplex basis: S_k = Σ_j B[k, j] T_j (ldhybrid.v1.HYBRID4_BASIS)
        T = _monomial_pole(v1.HYBRID4_EPS)
        shapes = []
        for row in v1.HYBRID4_BASIS:
            S = sl.Law([0.0])
            for b, Tj in zip(row, T):
                if b != 0.0:
                    S = S + float(b) * Tj
            shapes.append(S)
        labels = ["w1:1-mu2", "w2:(1-mu2)^2", "w3:vertex-inner", "w4:vertex-outer"]
        return _from_shapes(name, labels, shapes)
    shapes = _monomial_pole(v1.HYBRID5_EPS)
    labels = ["w1:1-mu4", "w2:(1-mu2)^2"] + [f"w{3 + k}:pole{e:.4g}" for k, e in enumerate(v1.HYBRID5_EPS)]
    return _from_shapes(name, labels, shapes)


def hybrid4_monomial_law() -> LinearLaw:
    """hybrid4 on the monomial-pole shapes (SquishierPlanet's basis; the LDhybrid
    basis before 2026-10-14), for studies and conversions."""
    eps = v1.HYBRID4_EPS
    return _from_shapes("hybrid4_monomial",
                        ["w1:1-mu4", "w2:(1-mu2)^2"] + [f"w{3 + k}:pole{e:.4g}" for k, e in enumerate(eps)],
                        _monomial_pole(eps))


def hybrid4_monomial(w):
    """hybrid4 weights on SquishierPlanet's monomial-pole shapes
    (1 − μ⁴, (1 − μ²)², Π_ε0, Π_ε1), from weights on the curvature-simplex basis."""
    return np.asarray(w, float) @ v1.HYBRID4_BASIS


def hybrid4_from_monomial(w_mp):
    return np.linalg.solve(v1.HYBRID4_BASIS.T, np.asarray(w_mp, float).T).T


def linear() -> LinearLaw:
    return LinearLaw("linear", ("u",), "mu",
                     lambda mu: (1.0 - mu)[..., None],
                     lambda mu: -np.ones_like(mu)[..., None])


def quadratic() -> LinearLaw:
    return LinearLaw("quadratic", ("u1", "u2"), "mu",
                     lambda mu: np.stack([1.0 - mu, (1.0 - mu) ** 2], axis=-1),
                     lambda mu: np.stack([-np.ones_like(mu), -2.0 * (1.0 - mu)], axis=-1))


def claret4() -> LinearLaw:
    k = np.arange(1, 5) / 2.0
    return LinearLaw("claret4", ("a1", "a2", "a3", "a4"), "mu",
                     lambda mu: 1.0 - np.asarray(mu)[..., None] ** k,
                     lambda mu: -k * np.asarray(mu)[..., None] ** (k - 1.0))


LAWS = {"hybrid2": lambda: hybrid("hybrid2"), "hybrid4": lambda: hybrid("hybrid4"),
        "hybrid5": lambda: hybrid("hybrid5"), "linear": linear,
        "quadratic": quadratic, "claret4": claret4}


def get(name: str) -> LinearLaw:
    return LAWS[name]()


def power2(c, alpha, mu):
    """Power-2 law I/I(1) = 1 − c(1 − μ^α) (Hestroffer 1997)."""
    return 1.0 - c * (1.0 - np.asarray(mu, float) ** alpha)
