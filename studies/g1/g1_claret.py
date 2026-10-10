"""Gate G1(c): compare with Claret & Southworth (2023, J/A+A/674/A63),
power-2 from spherical PHOENIX with the M1/M2 limb cuts.

Not a reproduction: their PHOENIX-COND intensities (vturb 2 km/s) are not
public, ours are PHOENIX-ACES SpecInt v3.0. We apply their recipe (cut at
max |dI/dr| (M1) or midway outward (M2), rescale μ, 900 uniform μ points,
unweighted Levenberg–Marquardt) and report the differences, which measure
model-version systematics rather than pipeline error.

    .venv/bin/python studies/g1/g1_claret.py
"""

import json
from pathlib import Path

import numpy as np

from ldhybrid import bands, fit
from ldhybrid.grids import phoenix
from ldhybrid.paths import RAW
from ldhybrid.profile import band_profile

HERE = Path(__file__).parent
BANDS = {"Kepler": 3, "TESS": 4, "CHEOPS": 5}         # column offsets in their table
STARS = [(t, 5.0 if t <= 3800 else 4.5) for t in
         (2600, 3000, 3400, 3800, 4200, 4600, 5000, 5400, 5800, 6200, 6600, 7000, 7400)]


def claret_table(method):
    t = {"M1": "table1.dat", "M2": "table2.dat"}[method]
    rows = {}
    for line in open(RAW / "vizier" / "674_A63" / t):
        v = line.split()
        logg, teff = float(v[0]), float(v[1])
        a = np.array(v[4:], float)
        rows[(int(teff), logg)] = dict(g=a[0:6], h=a[6:12], mucri=a[12:18], chi2=a[18:24])
    return rows


def main():
    out = []
    for method in ("M1", "M2"):
        cl = claret_table(method)
        for teff, logg in STARS:
            if (teff, logg) not in cl:
                continue
            s = phoenix.load(teff, logg, "-0.0")
            for b, j in BANDS.items():
                lam, S = bands.response(b)
                pr = band_profile(s, lam, S, limb=method)
                f = fit.fit_power2(pr, measure="mu", resample=900)
                c = cl[(teff, logg)]
                mucri = float(np.sqrt(1 - pr.meta["r_c"] ** 2))
                out.append(dict(method=method, teff=teff, logg=logg, band=b,
                                ours=f.w.tolist(), claret=[c["g"][j], c["h"][j]],
                                mucri_ours=mucri, mucri_claret=float(c["mucri"][j])))
                r = out[-1]
                print(f"{method} {teff} {logg} {b:6s} g {r['ours'][0]:.4f}/{r['claret'][0]:.4f} "
                      f"h {r['ours'][1]:.4f}/{r['claret'][1]:.4f}  "
                      f"mu_cri {mucri:.4f}/{r['mucri_claret']:.4f}", flush=True)
    (HERE / "g1_claret.json").write_text(json.dumps(out, indent=1))
    for m in ("M1", "M2"):
        d = np.array([np.subtract(r["ours"], r["claret"]) for r in out if r["method"] == m])
        print(m, "median |Δg|, |Δh|:", np.median(np.abs(d), 0), " max:", np.abs(d).max(0))


if __name__ == "__main__":
    main()
