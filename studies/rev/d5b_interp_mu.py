"""D5b: further interpolants of the 'exact' profile, in μ: Akima in μ and a
not-a-knot cubic spline in μ, against the baseline PCHIP in μ (P2S). Reports
the light-curve difference between the truths themselves and the re-scored
hybrid5 / four-term / hybrid4 production fits.

    .venv/bin/python studies/rev/d5b_interp_mu.py      # writes d5b.json
"""
import json

import numpy as np
from scipy.interpolate import Akima1DInterpolator, CubicSpline

import rc


def alts(pr):
    mu0 = pr.mu[0]
    ak = Akima1DInterpolator(pr.mu, pr.y)
    cs = CubicSpline(pr.mu, pr.y)
    wrap = lambda f: (lambda mu: np.where(np.asarray(mu) < mu0, pr(mu), f(np.clip(mu, mu0, 1.0))))
    return dict(akima_mu=wrap(ak), spline_mu=wrap(cs))


def task(job):
    kind, key = job
    spec = rc.load(kind, key)
    g = rc.lcgrid()
    out = []
    for band in rc.P2S_BANDS:
        pr = rc.profile(spec, band)
        if pr is None:
            continue
        T = g.columns(pr)
        A = {k: (f, g.columns(f)) for k, f in alts(pr).items()}
        row = dict(kind=kind, key=list(key), band=band,
                   truth_diff={k: rc.score(f, T, g) for k, (f, _) in A.items()}, laws={})
        for name in ("hybrid5", "hybrid4", "claret4", "power2", "quadratic"):
            f, _, _ = rc.fit_law(name, pr, "r")
            row["laws"][name] = dict(base=rc.score(f, T, g),
                                     **{k: rc.score(f, c, g) for k, (_, c) in A.items()})
        out.append(row)
    return out


if __name__ == "__main__":
    res = rc.run_parallel(task, rc.p2s_models(), workers=5)
    (rc.HERE / "d5b.json").write_text(json.dumps([r for rows in res for r in rows]))
    print("done")
