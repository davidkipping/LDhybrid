"""P2S study for the revision: D1 (null laws without poles), D4 (larger planets
and grazing transits), D5 (sensitivity to the interpolant of the 'exact'
profile), D7 (the limb cost of r-weighting) and D12 (where the light-curve
error comes from), on the P2 subset (34 PHOENIX + 16 MPS-ATLAS Set 2 models ×
8 bands), with production fits.

    .venv/bin/python studies/rev/p2s_study.py          # writes p2s_study.jsonl
"""
from __future__ import annotations

import json

import numpy as np
from scipy.interpolate import Akima1DInterpolator

import rc
from ldhybrid import lightcurve as lc
from ldhybrid import v1
from ldhybrid.profile import quadrature

LAWS = rc.LAWS6 + ["poly2", "poly4", "poly5"]
MU_SPLIT = 0.05
GEOM = {"p0.2": dict(p_grid=(0.2,), b_grid=(0.0, 0.3, 0.6, 0.9), n=60),
        "p0.3": dict(p_grid=(0.3,), b_grid=(0.0, 0.3, 0.6, 0.9), n=60),
        "p0.03": dict(p_grid=(0.03,), b_grid=(0.0, 0.3, 0.6, 0.9), n=60),
        "graze0.1": dict(p_grid=(0.1,), b_grid=(0.95,), n=60),
        "graze0.2": dict(p_grid=(0.2,), b_grid=(0.9,), n=60),
        "graze0.3": dict(p_grid=(0.3,), b_grid=(0.85,), n=60)}


def truths(pr):
    """Alternative 'exact' profiles through the same nodes."""
    mu0 = pr.mu[0]
    r_nodes = np.sqrt(np.clip(1 - pr.mu ** 2, 0, None))[::-1]
    y_r = pr.y[::-1]
    ak = Akima1DInterpolator(r_nodes, y_r)

    def below(mu, f):
        mu = np.asarray(mu, float)
        return np.where(mu < mu0, pr(mu), f(np.clip(mu, mu0, 1.0)))

    lin_mu = lambda mu: below(mu, lambda m: np.interp(m, pr.mu, pr.y))
    lin_r = lambda mu: below(mu, lambda m: np.interp(np.sqrt(1 - m * m), r_nodes, y_r))
    akima_r = lambda mu: below(mu, lambda m: ak(np.sqrt(1 - m * m)))
    return dict(lin_mu=lin_mu, lin_r=lin_r, akima_r=akima_r)


def wrms(f, pr, measure):
    mu, y, q = quadrature(pr, measure, v1.RESAMPLE)
    r = f(mu) - y
    return float(np.sqrt(np.sum(q * r * r) / np.sum(q)))


def task(job):
    kind, key = job
    spec = rc.load(kind, key)
    out = []
    for band in rc.P2S_BANDS:
        pr = rc.profile(spec, band)
        if pr is None:
            continue
        g = rc.lcgrid()
        T = g.columns(pr)
        alt = {k: g.columns(f) for k, f in truths(pr).items()}
        gG = {k: rc.lcgrid(**v) for k, v in GEOM.items()}
        TG = {k: gg.columns(pr) for k, gg in gG.items()}
        mu = np.linspace(pr.mu[0], 1.0, 4001)
        lim = mu < MU_SPLIT
        row = dict(kind=kind, key=list(key), band=band, mu_min=float(pr.mu[0]), laws={},
                   truth_diff={k: rc.score(f, T, g) for k, f in truths(pr).items()})
        for name in LAWS:
            f, w, F = rc.fit_law(name, pr, "r")
            fm, wm, Fm = rc.fit_law(name, pr, "mu")
            d = f(mu) - pr(mu)
            dm = fm(mu) - pr(mu)
            # D12: law below μ = 0.05 and truth above (limb part), and the reverse
            f_limb = lambda m, f=f: np.where(np.asarray(m) < MU_SPLIT, f(m), pr(m))
            f_in = lambda m, f=f: np.where(np.asarray(m) < MU_SPLIT, pr(m), f(m))
            e = dict(
                w=[float(x) for x in w], active=F.active,
                lc=rc.score(f, T, g), lc_mufit=rc.score(fm, T, g),
                **{f"lc_{k}": rc.score(f, A, g) for k, A in alt.items()},
                sig_mu=wrms(f, pr, "mu"), sig_r=wrms(f, pr, "r"),
                sig_mu_mufit=wrms(fm, pr, "mu"), sig_r_mufit=wrms(fm, pr, "r"),
                dI_limb=float(np.max(np.abs(d[lim]))) if lim.any() else np.nan,
                dI_in=float(np.max(np.abs(d[~lim]))),
                dI_limb_mufit=float(np.max(np.abs(dm[lim]))) if lim.any() else np.nan,
                dI_in_mufit=float(np.max(np.abs(dm[~lim]))),
                lc_limbpart=rc.score(f_limb, T, g), lc_inpart=rc.score(f_in, T, g),
                **{f"lc_{k}": rc.score(f, TG[k], gG[k]) for k in GEOM})
            row["laws"][name] = e
        out.append(row)
    return out


def main():
    jobs = rc.p2s_models()
    res = rc.run_parallel(task, jobs, workers=5)
    with open(rc.HERE / "p2s_study.jsonl", "w") as fh:
        for rows in res:
            for r in rows:
                fh.write(json.dumps(r) + "\n")
    print("done", sum(len(r) for r in res), "model-bands")


if __name__ == "__main__":
    main()
