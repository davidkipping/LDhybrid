"""D13: photon- vs energy-weighted band integration. For the five Table 3 stars
in TESS, Kepler, WFC3 G141, NIRSpec PRISM and 2MASS J, fit each law to the
photon-weighted profile ∫ I λ S dλ (as tabulated) and to the energy-weighted
one ∫ I S dλ, and compare the coefficients and the light curves (p = 0.1,
production grid).

    .venv/bin/python studies/rev/d13_photon_energy.py   # writes d13.json
"""
import json

import numpy as np

import rc
from ldhybrid import bands
from ldhybrid import lightcurve as lc
from ldhybrid.profile import band_profile

BANDS = ["TESS", "Kepler", "HST_WFC3_G141", "JWST_NIRSpec_PRISM", "2MASS_J"]
LAWS = ["hybrid5", "claret4", "quadratic", "power2"]


def main():
    out = {}
    for sname, (t, g) in rc.T3_STARS:
        spec = rc.load_t3(t, g)
        for b in BANDS:
            lam, S = bands.response(b)
            prP = band_profile(spec, lam, S, limb="M1")
            prE = band_profile(spec, lam, S / lam, limb="M1")
            e = dict(r_c=[prP.meta["r_c"], prE.meta["r_c"]],
                     profile_lc=lc.lc_errors(prE, prP, p_grid=(0.1,), b_grid=(0, .3, .6, .9), n=60)["max_ppm"],
                     dI_max=float(np.max(np.abs(prE(np.linspace(0, 1, 2001)) - prP(np.linspace(0, 1, 2001))))))
            for n in LAWS:
                fP, wP, _ = rc.fit_law(n, prP)
                fE, wE, _ = rc.fit_law(n, prE)
                e[n] = dict(dw=(np.asarray(wE) - np.asarray(wP)).tolist(),
                            lc=lc.lc_errors(fE, fP, p_grid=(0.1,), b_grid=(0, .3, .6, .9), n=60)["max_ppm"])
            out[f"{sname}|{b}"] = e
            print(sname, b, f"profile {e['profile_lc']:.1f} ppm", flush=True)
    (rc.HERE / "d13.json").write_text(json.dumps(out))


if __name__ == "__main__":
    main()
