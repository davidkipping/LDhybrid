"""D12c: light-curve error of r-weighted fits that ignore the extreme limb
(fit only μ ≥ μ_min, μ_min = 0, 0.02, 0.05, 0.1), scored against the full
model light curve on the production grid. P2S.

    .venv/bin/python studies/rev/d12c_mu_min.py       # writes d12c.json
"""
import json

import rc
from ldhybrid import fit, laws, v1

MUMIN = (0.0, 0.02, 0.05, 0.1)
LAWS = ["hybrid5", "hybrid4", "claret4", "power2", "hybrid2", "quadratic"]


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
        row = dict(kind=kind, key=list(key), band=band, lc={})
        for n in LAWS:
            for m0 in MUMIN:
                if n == "power2":
                    c, a = fit.fit_power2(pr, "r", v1.RESAMPLE, mu_min=m0).w
                    f = lambda m, c=c, a=a: laws.power2(c, a, m)
                else:
                    L = laws.get(n)
                    w = fit.fit_linear(L, pr, "r", v1.RESAMPLE, rc.constraint_for(n), mu_min=m0).w
                    f = lambda m, L=L, w=w: L.intensity(w, m)
                row["lc"][f"{n}|{m0}"] = rc.score(f, T, g)
        out.append(row)
    return out


if __name__ == "__main__":
    res = rc.run_parallel(task, rc.p2s_models(), workers=5)
    (rc.HERE / "d12c.json").write_text(json.dumps([r for rows in res for r in rows]))
    print("done")
