"""C7: convergence of the reference light curve and of the light-curve error.

Variants of lightcurve's quadrature (Gauss–Legendre nodes per segment N_GL,
sub-segments n_sub of each integration interval) and of the chord sampling
(points per half chord, production 60). Reported: the largest change, over
laws and cases, of the light-curve error of the tabulated fits (production
geometry p = 0.1, b ∈ {0, 0.3, 0.6, 0.9}), and of the model light curve itself.
T3 cases plus 30 PHOENIX P2S model-bands.

    .venv/bin/python studies/rev2/c7_convergence.py
"""
import json

import numpy as np

import rc
from ldhybrid import lightcurve as lc

VARIANTS = [("base", 48, 4, 60), ("GL96", 96, 4, 60), ("nsub16", 48, 16, 60),
            ("GL96+nsub16", 96, 16, 60), ("chord240", 48, 4, 240), ("chord960", 48, 4, 960),
            ("all", 96, 16, 960)]
_orig_blocked = lc.blocked


def set_quad(n_gl, n_sub):
    lc._TH, lc._WT = lc._gl(n_gl)
    lc.blocked = lambda I_of_r, z, p, n_sub_=n_sub: _orig_blocked(I_of_r, z, p, n_sub=n_sub_)


def case(pr):
    fits = {n: rc.fit_law(n, pr)[0] for n in rc.LAWS6}
    out = {}
    ref_curve = None
    for name, ngl, nsub, nch in VARIANTS:
        set_quad(ngl, nsub)
        g = lc.LCGrid(p_grid=(0.1,), b_grid=(0.0, 0.3, 0.6, 0.9), n=nch)
        tr = g.columns(pr)
        errs = {n: g.errors(*g.columns(f), *tr)["max_ppm"] for n, f in fits.items()}
        # model light curve on a common fine chord (b = 0) to compare quadratures
        z = lc.chord_z(0.1, 0.0, 400)
        curve = lc.light_curve(pr, z, 0.1)
        if ref_curve is None:
            ref_curve = curve
        out[name] = dict(errs=errs, dcurve_ppm=float(np.max(np.abs(curve - ref_curve)) * 1e6))
    set_quad(48, 4)
    return out


def task(job):
    kind, key, bandlist, label = job
    spec = rc.load_t3(*key) if kind == "t3" else rc.load(kind, key)
    rows = []
    for b in bandlist:
        pr = rc.profile(spec, b)
        if pr is not None:
            rows.append(dict(case=label, band=b, res=case(pr)))
    return rows


def main():
    rng = np.random.default_rng(3)
    jobs = [("t3", key, [b for b, _ in rc.T3_BANDS], name) for name, key in rc.T3_STARS]
    ph = [key for kind, key in rc.p2s_models() if kind == "phoenix"]
    for i in rng.choice(len(ph), 15, replace=False):
        jobs.append(("phoenix", ph[i], list(rng.choice(rc.P2S_BANDS, 2, replace=False)), str(ph[i])))
    rows = [r for res in rc.run_parallel(task, jobs) for r in res]
    (rc.HERE / "c7_convergence.json").write_text(json.dumps(rows))
    L = ["## C7 — convergence of the light-curve error (ppm)", "",
         f"{len(rows)} model-bands (T3 + 30 PHOENIX P2S); change relative to the production "
         "quadrature (N_GL = 48, n_sub = 4, 60 points per half chord).", "",
         "| variant | max |Δ lc error| over laws & cases | median |Δ| | max |Δ model curve| (b = 0) |",
         "|---|---|---|---|"]
    for name, *_ in VARIANTS[1:]:
        d = [abs(r["res"][name]["errs"][n] - r["res"]["base"]["errs"][n]) for r in rows for n in rc.LAWS6]
        dc = [r["res"][name]["dcurve_ppm"] for r in rows]
        L.append(f"| {name} | {max(d):.3f} | {np.median(d):.4f} | {max(dc):.2e} |")
    worst = max(((abs(r["res"]["chord960"]["errs"][n] - r["res"]["base"]["errs"][n]), r["case"], r["band"], n,
                  r["res"]["base"]["errs"][n]) for r in rows for n in rc.LAWS6))
    L.append("")
    L.append(f"Largest chord-sampling change: {worst[0]:.2f} ppm ({worst[1]} {worst[2]} {worst[3]}, "
             f"production {worst[4]:.2f} ppm).")
    rel = [abs(r["res"]["all"]["errs"][n] / r["res"]["base"]["errs"][n] - 1) for r in rows for n in rc.LAWS6]
    L.append(f"Relative change (all refinements): median {100 * np.median(rel):.2f}%, max {100 * max(rel):.1f}%.")
    (rc.HERE / "c7_convergence.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
