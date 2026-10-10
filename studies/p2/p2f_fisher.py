"""P2f — how many limb-darkening combinations does a transit constrain?

Fisher analysis of the hybrid4/5 weights for realistic transits, with the
geometry fixed (optimistic) and free (p, b, a/R★, t0; realistic). The
fiducial profile is the law fitted to a PHOENIX model, so the derivatives are
taken at a realistic point. Prior scale: the simplex, i.e. ~0.1–0.3 per weight.

Cases: JWST white-light-like (σ = 100 ppm per 1-min point) and a TESS-like
2-min cadence (σ = 1000 ppm), one transit, p = 0.1 and 0.15, b = 0.3 and 0.7,
P = 3 d, a/R★ = 10.

Reports, per case, the posterior standard deviations along the Fisher
eigen-directions of w, the number of directions constrained below 0.05 and
0.02, and the leading eigenvectors.

    .venv/bin/python studies/p2/p2f_fisher.py
"""

import json
from pathlib import Path

import numpy as np

from ldhybrid import bands, fit, laws
from ldhybrid import lightcurve as lc
from ldhybrid.grids import phoenix
from ldhybrid.profile import band_profile

HERE = Path(__file__).parent
STARS = [((5800, 4.5, "-0.0"), "HST_WFC3_G141"), ((5800, 4.5, "-0.0"), "TESS"),
         ((3400, 5.0, "-0.0"), "JWST_NIRSpec_PRISM"), ((4600, 4.5, "-0.0"), "JWST_NIRISS_SOSS_o1")]
CASES = [dict(name="JWST-like", sigma=100e-6, cadence_min=1.0),
         dict(name="TESS-like", sigma=1000e-6, cadence_min=2.0)]
GEOM = [(0.1, 0.3), (0.1, 0.7), (0.15, 0.3), (0.15, 0.7)]
PERIOD_D, AR = 3.0, 10.0


def flux(law, w, p, b, ar, t0, t):
    """Circular orbit, small-angle: z(t) = sqrt(b² + (2π ar (t − t0)/P)²)."""
    x = 2 * np.pi * ar * (t - t0) / PERIOD_D
    z = np.sqrt(b * b + x * x)
    out = np.ones_like(t)
    m = z < 1 + p
    out[m] = lc.light_curve(lambda mu: law.intensity(w, mu), z[m], p)
    return out


def fisher(law, w, p, b, case):
    T14 = PERIOD_D / np.pi * np.arcsin(np.sqrt((1 + p) ** 2 - b * b) / AR)
    t = np.arange(-0.75 * T14, 0.75 * T14, case["cadence_min"] / 1440.0)
    theta = dict(p=p, b=b, ar=AR, t0=0.0)
    f0 = lambda **kw: flux(law, kw.get("w", w), kw.get("p", p), kw.get("b", b),
                           kw.get("ar", AR), kw.get("t0", 0.0), t)
    J = []
    for j in range(law.m):                                  # weights: central difference
        h = 1e-4
        wp, wm = w.copy(), w.copy()
        wp[j] += h
        wm[j] -= h
        J.append((f0(w=wp) - f0(w=wm)) / (2 * h))
    for k, h in (("p", 1e-5), ("b", 1e-4), ("ar", 1e-3), ("t0", 1e-6)):
        J.append((f0(**{k: theta[k] + h}) - f0(**{k: theta[k] - h})) / (2 * h))
    J = np.array(J).T / case["sigma"]
    return J.T @ J


def analyse(Fm, m):
    """Posterior covariance of the weights: geometry fixed (top-left block
    inverted) and marginalised (full inverse); add a weak prior so that the
    matrices are invertible (σ_prior = 1 per weight, far wider than the simplex)."""
    P = np.zeros_like(Fm)
    P[:m, :m] = np.eye(m)
    out = {}
    for label, Cfull in (("fixed", np.linalg.inv(Fm[:m, :m] + np.eye(m))),
                         ("free", np.linalg.inv(Fm + P)[:m, :m])):
        ev, V = np.linalg.eigh(Cfull)
        sd = np.sqrt(ev)
        out[label] = dict(sd=sd.tolist(), n_lt_005=int(np.sum(sd < 0.05)),
                          n_lt_002=int(np.sum(sd < 0.02)),
                          best_dir=V[:, 0].tolist(), sd_marginal=np.sqrt(np.diag(Cfull)).tolist())
    return out


def main():
    res = []
    lines = ["# P2f: how many weight combinations does a transit constrain?", "",
             "Posterior s.d. along Fisher eigen-directions of the hybrid weights "
             "(smallest first), one transit, a/R★ = 10, P = 3 d. Geometry fixed / free "
             "(p, b, a/R★, t0). Prior scale for comparison: the simplex, ~0.1–0.3 "
             "per weight.", ""]
    for (key, band) in STARS:
        pr = band_profile(phoenix.load(*key), *bands.response(band), limb="M1")
        for name in ("hybrid4", "hybrid5"):
            L = laws.get(name)
            w = fit.fit_linear(L, pr, "r", 2000, "simplex").w
            for case in CASES:
                for p, b in GEOM:
                    a = analyse(fisher(L, w, p, b, case), L.m)
                    res.append(dict(star=key, band=band, law=name, case=case["name"],
                                    p=p, b=b, w=w.tolist(), **a))
                    lines.append(
                        f"* {key[0]} K {band}, {name}, {case['name']}, p={p}, b={b}: "
                        f"fixed sd {np.round(a['fixed']['sd'], 3).tolist()} "
                        f"(<0.05: {a['fixed']['n_lt_005']}); free sd "
                        f"{np.round(a['free']['sd'], 3).tolist()} (<0.05: {a['free']['n_lt_005']})")
                    print(lines[-1], flush=True)
    (HERE / "p2f_fisher.json").write_text(json.dumps(res, indent=1))
    (HERE / "p2f_fisher.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
