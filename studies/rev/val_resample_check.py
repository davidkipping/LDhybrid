"""Side check found during D8: convergence of the production resample
(2000 points uniform in r) for steep PHOENIX limbs. LC error (production grid)
of r-weighted fits against the model profile vs the number of resample points.

    .venv/bin/python studies/rev/val_resample_check.py
"""
import numpy as np
from ldhybrid import bands, fit, laws
from ldhybrid import lightcurve as lc
from ldhybrid.grids import phoenix
from ldhybrid.products import PROD_LC, constraint_for
from ldhybrid.profile import band_profile

rng = np.random.default_rng(3)
cases = [(2600, 4.0, "JWST_NIRCam_F444W"), (3300, 5.0, "TESS"), (3300, 5.0, "JWST_NIRSpec_PRISM"),
         (4600, 4.5, "HST_WFC3_G141"), (5800, 4.5, "TESS"), (5800, 4.5, "JWST_NIRSpec_PRISM"),
         (6500, 4.0, "Kepler"), (2800, 5.0, "HST_WFC3_G141")]
NS = (2000, 8000, 30000, 100000)
print("case | law | " + " | ".join(f"n={n}" for n in NS))
for t, g, b in cases:
    pr = band_profile(phoenix.load(t, g, "-0.0"), *bands.response(b), limb="M1")
    for law in ("hybrid5", "claret4", "power2"):
        e = []
        for n in NS:
            if law == "power2":
                c, a = fit.fit_power2(pr, "r", n).w
                f = lambda mu, c=c, a=a: laws.power2(c, a, mu)
            else:
                L = laws.get(law)
                w = fit.fit_linear(L, pr, "r", n, constraint_for(law)).w
                f = lambda mu, L=L, w=w: L.intensity(w, mu)
            e.append(lc.lc_errors(f, pr, **PROD_LC)["max_ppm"])
        print(f"{t} {g} {b} | {law} | " + " | ".join(f"{x:.1f}" for x in e), flush=True)
