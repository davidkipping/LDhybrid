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
HYBRID4_EPS = tuple(float(e) for e in ladder(2))      # 0.011762, 0.30408
HYBRID5_EPS = tuple(float(e) for e in ladder(3))      # 0.0015966, 0.041275, 0.38743

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
assert np.allclose(HYBRID4_EPS, (0.01176198, 0.30407926), rtol=1e-6)
assert np.allclose(HYBRID5_EPS, (0.00159655, 0.04127504, 0.38743472), rtol=1e-6)
