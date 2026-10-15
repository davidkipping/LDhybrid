"""LDhybrid v1.0 definitions, frozen at Gate G2 (2026-10-10). Never change
these after tables are published; a change means a new version.

Evidence: studies/p2/results.md.
"""

import numpy as np
from squishierplanet.laws import ladder

# --- laws -----------------------------------------------------------------
#: hybrid2's pole, re-tuned over 7,524 model-bands (P2c; was 0.208, J band)
HYBRID2_EPS = 0.36
#: hybrid4 / hybrid5 poles: the analytic ladder, alpha_eff = 1 (unchanged)
#: hybrid4: the analytic ladder with alpha_eff = 25/16 (ln eps scaled by 0.8),
#: adopted 2026-10-11 after studies/rev3 (was alpha_eff = 1: 0.011762, 0.30408)
HYBRID4_ALPHA = 25.0 / 16.0
HYBRID4_EPS = tuple(float(e) for e in ladder(2, HYBRID4_ALPHA))  # 0.028601, 0.38582
HYBRID5_EPS = tuple(float(e) for e in ladder(3))      # 0.0015966, 0.041275, 0.38743


def _hybrid4_vertices(eps=HYBRID4_EPS):
    """hybrid4's basis: the four non-trivial vertices of the curvature simplex,
    as rows of weights on the monomial-pole shapes (1 − μ⁴, (1 − μ²)², Π_ε0, Π_ε1).

    With g(ξ) = dI/dξ (ξ = μ²), the profiles with I(0) ≥ 0, g ≥ 0, g' ≤ 0 and
    g'' ≥ 0 on [0, 1] are exactly the 4-simplex Σw ≤ 1, g(1) ≥ 0, g'(1) ≤ 0,
    g''(0) ≥ 0, g''(1) ≥ 0 (g'' spans two poles with a monotonic ratio; see
    studies/rev3/h4_scheme.py). Its vertices are I = 1 and the rows below:
      S1 = 1 − μ²             (I = μ²: g' = 0, g'' = 0)
      S2 = (1 − μ²)²          (I = 1 − (1 − μ²)²: g'' = 0, g(1) = 0)
      S3 = inner-pole vertex  (I(0) = 0, g(1) = g'(1) = g''(1) = 0)
      S4 = outer-pole vertex  (I(0) = 0, g(1) = g'(1) = g''(0) = 0)
    """
    e0, e1 = eps
    a0, a1 = (1.0 / (e ** -2 - (1 + e) ** -2) for e in eps)
    g = lambda x: [2 * x, 2 * (1 - x), 2 * a0 * (x + e0) ** -3, 2 * a1 * (x + e1) ** -3]
    g1 = lambda x: [2.0, -2.0, -6 * a0 * (x + e0) ** -4, -6 * a1 * (x + e1) ** -4]
    g2 = lambda x: [0.0, 0.0, 24 * a0 * (x + e0) ** -5, 24 * a1 * (x + e1) ** -5]
    rhs = [1.0, 0.0, 0.0, 0.0]
    S3 = np.linalg.solve(np.array([[1.0] * 4, g(1.0), g1(1.0), g2(1.0)]), rhs)
    S4 = np.linalg.solve(np.array([[1.0] * 4, g(1.0), g1(1.0), g2(0.0)]), rhs)
    return np.array([[0.5, 0.5, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], S3, S4])


#: hybrid4 shapes S_k = Σ_j HYBRID4_BASIS[k, j] T_j (adopted 2026-10-14): w ≥ 0,
#: Σw ≤ 1 is then exactly the set of hybrid4 profiles that are positive,
#: monotonic, concave in μ² and of convex slope in μ². Monomial-pole weights
#: (SquishierPlanet's hybrid4 basis) are w @ HYBRID4_BASIS.
HYBRID4_BASIS = _hybrid4_vertices()

# --- headline fitting method (P2a, P2b) -------------------------------------
MEASURE = "r"              # least squares uniform in projected radius
#: quadrature of the weighted least squares: Gauss–Legendre on every interval
#: between profile nodes (exact for the PCHIP profile). Was 2000 points uniform
#: in r until 2026-10-10, which under-resolved the steep limb of spherical
#: models (studies/rev/val_resample_check.md).
RESAMPLE = "exact"
CONSTRAINT = {"hybrid2": "physical",     # exact triangle
              "hybrid4": "simplex", "hybrid5": "simplex"}
FCM = False                # flux conservation not imposed (flux error tabulated)
#: secondary table, comparable with Claret: uniform in μ
SECONDARY_MEASURE = "mu"

# --- limb treatment of spherical models (P2d) -------------------------------
LIMB = "M1"                # stellar edge at max |dI/dr| (Claret & Southworth M1)
LIMB_BAND_TABLES = "per band"     # r_c found per passband (Claret convention)
LIMB_SPECTRA = "per model"        # one r_c per model (photon-bolometric profile),
                                  # so coefficient spectra integrate exactly

# --- scope ----------------------------------------------------------------
#: MIRI LRS: FGK only (Kurucz 3500–6500 K); no M dwarfs in v1.0
MIRI_TEFF_MIN = 3500

# guard against any change to squishierplanet's ladder
assert np.allclose(HYBRID4_EPS, (0.02860121, 0.38582418), rtol=1e-6)
assert np.allclose(HYBRID5_EPS, (0.00159655, 0.04127504, 0.38743472), rtol=1e-6)
assert np.allclose(HYBRID4_BASIS[2:], [[0.00063, 0.00140, 1.02096, -0.02299],
                                       [-0.08116, -0.25685, -0.00059, 1.33860]], atol=6e-6)
