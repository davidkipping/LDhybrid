"""Paper section: transit light-curve error of each law against the exact light
curve of the model's own intensity profile.

Stars (PHOENIX v3.0, [M/H] = 0): late M (2800 K, log g 5.0), mid M (3300 K,
5.0), early M (3800 K, 5.0), K dwarf (4600 K, 4.5), Sun analogue (5800 K, 4.5).
Bands: TESS, HST WFC3 G141, JWST NIRSpec PRISM. Central transit (b = 0) of a
planet with p = 0.1; coefficients as tabulated (r-weighted, constrained fits).
Truth: the exact light curve of the tabulated profile (PCHIP between nodes).
Also: the worst error for b = 0.3, 0.5, 0.7, 0.85, to show it barely depends on b.

    .venv/bin/python studies/lc/lc_section.py
"""

import json
from pathlib import Path

import numpy as np

from ldhybrid import bands, fit, laws, v1
from ldhybrid import lightcurve as lc
from ldhybrid.grids import phoenix
from ldhybrid.products import constraint_for
from ldhybrid.profile import band_profile

HERE = Path(__file__).parent
STARS = [("late M", (2800, 5.0)), ("mid M", (3300, 5.0)), ("early M", (3800, 5.0)),
         ("K dwarf", (4600, 4.5)), ("Sun analogue", (5800, 4.5))]
BANDS = [("TESS", "TESS"), ("HST_WFC3_G141", "WFC3 G141"), ("JWST_NIRSpec_PRISM", "NIRSpec PRISM")]
LAWS = ["quadratic", "power2", "claret4", "hybrid2", "hybrid4", "hybrid5"]
P = 0.1
B_EXTRA = (0.3, 0.5, 0.7, 0.85)
N = 1201


def law_profile(name, pr):
    if name == "power2":
        c, a = fit.fit_power2(pr, "r", v1.RESAMPLE).w
        return lambda mu: laws.power2(c, a, mu), [c, a]
    L = laws.get(name)
    w = fit.fit_linear(L, pr, "r", v1.RESAMPLE, constraint_for(name)).w
    return (lambda mu, L=L, w=w: L.intensity(w, mu)), w.tolist()


def main():
    x = np.linspace(-(1 + P) * 1.08, (1 + P) * 1.08, N)       # chord position, b = 0
    z = np.abs(x)
    out = dict(p=P, b=0.0, x=x.tolist(), stars=[s for s, _ in STARS],
               bands=[b for _, b in BANDS], laws=LAWS, resid={}, truth={}, worst={},
               worst_b={}, coeffs={})
    for sname, (t, g) in STARS:
        spec = phoenix.load(t, g, "-0.0")
        for bkey, blab in BANDS:
            pr = band_profile(spec, *bands.response(bkey), limb="M1")
            Ft = lc.light_curve(pr, z, P)
            key = f"{sname}|{blab}"
            out["truth"][key] = Ft.tolist()
            for n in LAWS:
                f, coef = law_profile(n, pr)
                d = (lc.light_curve(f, z, P) - Ft) * 1e6
                out["resid"][f"{key}|{n}"] = d.tolist()
                out["worst"][f"{key}|{n}"] = float(np.max(np.abs(d)))
                out["coeffs"][f"{key}|{n}"] = coef
                e = lc.lc_errors(f, pr, p_grid=(P,), b_grid=(0.0,) + B_EXTRA)["worst"][0]
                out["worst_b"][f"{key}|{n}"] = e.tolist()
            print(sname, blab, " ".join(f"{n}={out['worst'][f'{key}|{n}']:.1f}" for n in LAWS),
                  flush=True)
    (HERE / "lc_section.json").write_text(json.dumps(out))

    lines = ["# Light-curve error, central transit (b = 0), p = 0.1, PHOENIX [M/H] = 0", "",
             "Worst |ΔF| (ppm) of each law, with its tabulated coefficients, against the "
             "exact light curve of the model profile.", "",
             "| star | band | " + " | ".join(LAWS) + " |", "|---|---|" + "---|" * len(LAWS)]
    for sname, _ in STARS:
        for _, blab in BANDS:
            k = f"{sname}|{blab}"
            lines.append(f"| {sname} | {blab} | " +
                         " | ".join(f"{out['worst'][f'{k}|{n}']:.1f}" for n in LAWS) + " |")
    ratio = [max(v) / min(v) for v in out["worst_b"].values() if min(v) > 0]
    lines += ["", f"Dependence on b (0–0.85): max/min of the worst error per case, median "
              f"{np.median(ratio):.2f}, max {np.max(ratio):.2f}."]
    (HERE / "lc_section.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
