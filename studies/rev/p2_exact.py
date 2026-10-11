"""Re-run of the P2 method studies with the exact quadrature (v1.RESAMPLE =
"exact"), on the P2 subset (34 PHOENIX + 16 MPS-ATLAS Set 2 models × 8 bands)
and the P2 light-curve metric (LCGrid default: p ∈ {0.03, 0.1, 0.15},
b ∈ {0, …, 0.9}, worst |ΔF|):

  weighting  (P2a / Table 1): measures mu, r, area with the exact quadrature,
             and equal weights on the native nodes; six laws;
  pole       (P2c, reduced): hybrid2 with ε ∈ {0.24, …, 0.44};
  limb       (P2d): M2 vs M1, one r_c per model, μ nodes of ATLAS / MPS-ATLAS /
             Stagger, MPS-ATLAS Set 2 vs PHOENIX and Set 1 vs Set 2.

    .venv/bin/python studies/rev/p2_exact.py weighting|pole|limb
"""
import json
import sys

import numpy as np

import rc
from ldhybrid import bands, fit, laws, v1
from ldhybrid import lightcurve as lc
from ldhybrid.profile import Profile, band_profile, limb_cut

LAWS = ["quadratic", "power2", "hybrid2", "claret4", "hybrid4", "hybrid5"]
MEASURES = ("mu", "r", "area", "nodes")
EPS = (0.24, 0.28, 0.32, 0.36, 0.40, 0.44)
NODES = {
    "atlas17": np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.125,
                         0.1, 0.075, 0.05, 0.025, 0.01]),
    "mps24": np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.35, 0.3, 0.25, 0.22, 0.2, 0.17,
                       0.15, 0.12, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.03, 0.02, 0.01]),
    "stagger10": np.array([1.0, 0.9, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.01]),
}
_G = None
_C = {}


def G():
    global _G
    if _G is None:
        _G = lc.LCGrid()
    return _G


def cols(key, law):
    if key not in _C:
        _C[key] = G().linear_columns(law)
    return _C[key]


def con(name):
    return "simplex" if name in ("hybrid4", "hybrid5") else "physical"


def lin_err(name, law, w, T):
    return G().errors(*G().combine(*cols(name, law), w), *T)["max_ppm"]


def fit_err(name, pr, measure, T):
    rs = None if measure == "nodes" else v1.RESAMPLE
    if name == "power2":
        c, a = fit.fit_power2(pr, measure, rs).w
        return G().errors(*G().columns(lambda m: laws.power2(c, a, m)), *T)["max_ppm"], [c, a]
    L = laws.get(name)
    w = fit.fit_linear(L, pr, measure, rs, con(name)).w
    return lin_err(name, L, w, T), w.tolist()


def task_weighting(job):
    kind, key = job
    spec = rc.load(kind, key)
    out = []
    for band in rc.P2S_BANDS:
        pr = rc.profile(spec, band)
        if pr is None:
            continue
        T = G().columns(pr)
        row = dict(kind=kind, key=list(key), band=band, lc={})
        for m in MEASURES:
            for n in LAWS:
                row["lc"][f"{n}|{m}"] = fit_err(n, pr, m, T)[0]
        out.append(row)
    return out


def task_pole(job):
    kind, key = job
    spec = rc.load(kind, key)
    out = []
    for band in rc.P2S_BANDS:
        pr = rc.profile(spec, band)
        if pr is None:
            continue
        T = G().columns(pr)
        row = dict(kind=kind, key=list(key), band=band, lc={})
        for e in EPS:
            L = laws.hybrid2(float(e))
            w = fit.fit_linear(L, pr, "r", v1.RESAMPLE, "physical").w
            row["lc"][str(e)] = lin_err(f"h2_{e}", L, w, T)
        row["lc"]["quadratic"] = fit_err("quadratic", pr, "r", T)[0]
        out.append(row)
    return out


def task_limb(job):
    from ldhybrid.grids import mpsatlas, phoenix
    kind, key = job
    t, g, mh = key
    spec = phoenix.load(t, g, mh)
    lamall = np.array([spec.lam_edges[0], spec.lam_edges[-1]])
    _, rc_bol = limb_cut(spec.mu, spec.band(lamall, np.ones(2)), "M1")
    mkey = (float(mh), t, g)
    mps = {}
    if mkey in set(mpsatlas.catalogue("set2")):
        mps = {st: mpsatlas.load(st, *mkey) for st in ("set1", "set2")}
    out = []
    for band in rc.P2S_BANDS:
        lam, S = bands.response(band)
        base = band_profile(spec, lam, S, limb="M1")
        prs = {"M2": band_profile(spec, lam, S, limb="M2"),
               "rc_model": band_profile(spec, lam, S, limb="M1", r_c=rc_bol)}
        for k, nodes in NODES.items():
            prs[k] = Profile(nodes, base.interpolant()(nodes))
        for st, s in mps.items():
            prs[f"mps_{st}"] = band_profile(s, lam, S)
        row = dict(key=list(key), band=band, r_c=base.meta["r_c"], rc_bol=rc_bol, dlc={}, dw={})
        for n in ("hybrid4", "hybrid5", "claret4", "hybrid2"):
            L = laws.get(n)
            w0 = fit.fit_linear(L, base, "r", v1.RESAMPLE, con(n)).w
            F0 = G().combine(*cols(n, L), w0)
            W = {}
            for v, pr in prs.items():
                w = fit.fit_linear(L, pr, "r", v1.RESAMPLE, con(n)).w
                W[v] = w
                if v == "mps_set1":
                    continue
                row["dlc"][f"{n}|{v}"] = G().errors(*G().combine(*cols(n, L), w), *F0)["max_ppm"]
                row["dw"][f"{n}|{v}"] = (w - w0).tolist()
            if "mps_set1" in W:
                row["dlc"][f"{n}|set1_vs_set2"] = G().errors(*G().combine(*cols(n, L), W["mps_set1"]),
                                                            *G().combine(*cols(n, L), W["mps_set2"]))["max_ppm"]
            row["dw"][f"{n}|base"] = w0.tolist()
        out.append(row)
    return out


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "limb":
        jobs = [j for j in rc.p2s_models() if j[0] == "phoenix"]
    else:
        jobs = rc.p2s_models()
    task = dict(weighting=task_weighting, pole=task_pole, limb=task_limb)[what]
    res = rc.run_parallel(task, jobs, workers=5)
    (rc.HERE / f"p2x_{what}.json").write_text(json.dumps([r for rows in res for r in rows]))
    print("done", what)
