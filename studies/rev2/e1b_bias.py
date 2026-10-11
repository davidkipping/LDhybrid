"""E1b: parameter bias from the limb-darkening law (noiseless), on an equal
footing, and with an illustrative prior.

Set-up. Truth: the light curve of the model band profile (PHOENIX: the 15
Table-3 cases; MPS-ATLAS Set 2 at [M/H] = 0 for 3800 K/log g 5.0, 4600 K and
5800 K/log g 4.5), planet p = 0.1 at b = 0.3 on a straight chord
z(t) = sqrt(b² + (v t)²), v = 1, 600 points over |t| ≤ 1.2(1 + p), no noise,
no exposure-time integration, unweighted least squares. Mid-time and
out-of-transit baseline are fixed. Free geometry: p, b, v.
Modes:
  fixed — coefficients at their tabulated (r-weighted) values;
  free  — coefficients free within the law's EXACT physical region for every
          linear law (I(0) ≥ 0 and dI/dx ≥ 0 on the 2001-point grid, as fit.py
          "physical"; for hybrid4/5 this is larger than the tabulated simplex),
          solved by SLSQP; power-2 within 0 ≤ c ≤ 1, α > 0;
  prior — as free, plus an independent Gaussian prior of width 0.05 on each
          coefficient, centred on the tabulated values, with 180 ppm per point
          (≈100 ppm per minute at this sampling for a ~3 h transit).
Reported: Δ(p²) (ppm), Δb, Δv/v, rms residual (ppm), whether a constraint binds.

    .venv/bin/python studies/rev2/e1b_bias.py
"""
import json

import numpy as np
from scipy.optimize import least_squares, minimize

import rc
from e1_bias_v1 import B0, P0, V0, Model, z_of
from ldhybrid import fit, laws
from ldhybrid import lightcurve as lc
from ldhybrid.grids import mpsatlas

SIG_PT = 180.0
SIG_W = 0.05
GB = ([0.01, 0.0, 0.2], [0.5, 1.5, 5.0])


def ineq(name):
    C, d = fit._inequalities(laws.get(name), "physical")
    return C, d


def run_fixed(M, w_tab, F_true):
    res = least_squares(lambda x: (M.flux(*x, w_tab) - F_true) * 1e6, [P0, B0, V0],
                        bounds=GB, x_scale=[0.01, 0.05, 0.05], diff_step=1e-7,
                        xtol=1e-12, ftol=1e-12, gtol=1e-12)
    return res.x, np.array(w_tab, float), res.fun, False


def run_linear_free(M, name, w_tab, F_true, prior):
    C, d = ineq(name)
    n = len(w_tab)
    x_fixed, _, _, _ = run_fixed(M, w_tab, F_true)
    x0 = np.r_[x_fixed, w_tab]
    scale = np.r_[[0.01, 0.05, 0.05], np.full(n, 0.05)]

    def resid(y):
        x = y * scale
        r = (M.flux(*x[:3], x[3:]) - F_true) * 1e6
        if prior:
            r = np.r_[r / SIG_PT, (x[3:] - w_tab) / SIG_W]
        return r

    obj = lambda y: float(np.sum(resid(y) ** 2))
    cons = [{"type": "ineq", "fun": lambda y: d - C @ (y[3:] * scale[3:])}]
    bnds = [(GB[0][i] / scale[i], GB[1][i] / scale[i]) for i in range(3)] + [(None, None)] * n
    r = minimize(obj, x0 / scale, method="SLSQP", constraints=cons, bounds=bnds,
                 options=dict(ftol=1e-16, maxiter=2000, eps=1e-9))
    y = r.x
    # polish with LM if no constraint binds
    x = y * scale
    w = x[3:]
    slack = d - C @ w
    binds = bool(np.min(slack) < 1e-7)
    if not binds:
        lm = least_squares(lambda yy: resid(yy), y, method="lm", xtol=1e-13, ftol=1e-13, gtol=1e-13)
        if np.all(d - C @ (lm.x[3:] * scale[3:]) >= -1e-10) and np.sum(lm.fun ** 2) <= obj(y):
            y = lm.x
            x = y * scale
            w = x[3:]
    fun = (M.flux(*x[:3], w) - F_true) * 1e6
    return x[:3], w, fun, binds


def run_power2_free(M, w_tab, F_true, prior):
    c0, a0 = w_tab
    x_fixed, _, _, _ = run_fixed(M, w_tab, F_true)
    x0 = np.r_[x_fixed, min(max(c0, 1e-6), 1 - 1e-6), a0]

    def resid(x):
        r = (M.flux(*x[:3], x[3:]) - F_true) * 1e6
        if prior:
            r = np.r_[r / SIG_PT, (x[3:] - np.array(w_tab)) / SIG_W]
        return r

    res = least_squares(resid, x0, bounds=(np.r_[GB[0], 0.0, 1e-4], np.r_[GB[1], 1.0, np.inf]),
                        diff_step=1e-7, xtol=1e-12, ftol=1e-12, gtol=1e-12)
    x = res.x
    binds = bool(x[3] > 1 - 1e-6 or x[3] < 1e-6)
    return x[:3], x[3:], (M.flux(*x[:3], x[3:]) - F_true) * 1e6, binds


def summarise(x, w, fun, binds):
    p, b, v = x
    return dict(dp2_ppm=float((p * p - P0 * P0) * 1e6), db=float(b - B0), dv=float(v / V0 - 1),
                rms_ppm=float(np.sqrt(np.mean(fun ** 2))), binds=binds, w=[float(t) for t in w])


def task(job):
    grid, sname, spec_key = job
    spec = rc.load_t3(*spec_key) if grid == "phoenix" else mpsatlas.load("set2", 0.0, *spec_key)
    out = {}
    for bkey, blab in rc.T3_BANDS:
        pr = rc.profile(spec, bkey)
        if pr is None:
            continue
        F_true = lc.light_curve(pr, z_of(B0, V0), P0)
        for name in rc.LAWS6:
            _, w_tab, _ = rc.fit_law(name, pr, "r")
            M = Model(name)
            res = {"fixed": summarise(*run_fixed(M, w_tab, F_true))}
            for mode in ("free", "prior"):
                if name == "power2":
                    res[mode] = summarise(*run_power2_free(M, w_tab, F_true, mode == "prior"))
                else:
                    res[mode] = summarise(*run_linear_free(M, name, np.array(w_tab, float), F_true, mode == "prior"))
            out[f"{grid}|{sname}|{blab}|{name}"] = res
        print(grid, sname, blab, flush=True)
    return out


def main():
    jobs = [("phoenix", n, k) for n, k in rc.T3_STARS]
    jobs += [("mps2", "early M", (3800, 5.0)), ("mps2", "K dwarf", (4600, 4.5)),
             ("mps2", "Sun analogue", (5800, 4.5))]
    res = {k: v for r in rc.run_parallel(task, jobs) for k, v in r.items()}
    (rc.HERE / "e1b_bias.json").write_text(json.dumps(res))
    L = ["## E1b — parameter bias, equal footing (noiseless; see script docstring for the set-up)", ""]
    for grid in ("phoenix", "mps2"):
        keys = [k for k in res if k.startswith(grid + "|")]
        cases = sorted({k.split("|")[1] + "|" + k.split("|")[2] for k in keys})
        L += [f"### {grid} ({len(cases)} cases)", "",
              "| law | mode | |Δ(p²)| med / max (ppm) | |Δb| med / max | max |Δv/v| | rms med (ppm) | chromatic spread max (ppm) | binds |",
              "|---|---|---|---|---|---|---|---|"]
        for name in rc.LAWS6:
            for mode in ("fixed", "free", "prior"):
                rr = [res[f"{grid}|{c}|{name}"][mode] for c in cases]
                dp = np.abs([r["dp2_ppm"] for r in rr])
                db = np.abs([r["db"] for r in rr])
                dv = np.abs([r["dv"] for r in rr])
                rms = [r["rms_ppm"] for r in rr]
                stars = sorted({c.split("|")[0] for c in cases})
                spread = max(np.ptp([res[f"{grid}|{c}|{name}"][mode]["dp2_ppm"] for c in cases if c.startswith(s + "|")])
                             for s in stars)
                nb = sum(r["binds"] for r in rr)
                L.append(f"| {name} | {mode} | {np.median(dp):.2f} / {dp.max():.1f} | {np.median(db):.4f} / {db.max():.4f} | "
                         f"{dv.max():.4f} | {np.median(rms):.2f} | {spread:.1f} | {nb}/{len(rr)} |")
        L.append("")
    (rc.HERE / "e1b_bias.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
