"""T3 study (the 15 cases of the paper's Table 3): D1 null laws, D4 planet
size and grazing geometry (worst and rms |ΔF|), D12 where the error comes from.

    .venv/bin/python studies/rev/t3_study.py      # writes t3_study.json
"""
import json

import numpy as np

import rc
from ldhybrid import lightcurve as lc

LAWS = rc.LAWS6 + ["poly2", "poly4", "poly5"]
P = (0.03, 0.1, 0.2, 0.3)
MU_SPLIT = 0.05


def task(star):
    sname, (t, g) = star
    spec = rc.load_t3(t, g)
    out = {}
    for bkey, blab in rc.T3_BANDS:
        pr = rc.profile(spec, bkey)
        for name in LAWS:
            f, w, F = rc.fit_law(name, pr, "r")
            geo = {}
            for p in P:
                e = lc.lc_errors(f, pr, p_grid=(p,), b_grid=(0.0, 0.5, 1 - p / 2, 1.0))
                geo[str(p)] = dict(worst=e["worst"][0].tolist(), rms=e["rms"][0].tolist())
            f_limb = lambda m, f=f: np.where(np.asarray(m) < MU_SPLIT, f(m), pr(m))
            f_in = lambda m, f=f: np.where(np.asarray(m) < MU_SPLIT, pr(m), f(m))
            b0 = dict(p_grid=(0.1,), b_grid=(0.0,))
            out[f"{sname}|{blab}|{name}"] = dict(
                w=[float(x) for x in w], geo=geo,
                b0=lc.lc_errors(f, pr, **b0)["max_ppm"],
                limbpart=lc.lc_errors(f_limb, pr, **b0)["max_ppm"],
                inpart=lc.lc_errors(f_in, pr, **b0)["max_ppm"])
        print(sname, blab, flush=True)
    return out


if __name__ == "__main__":
    res = rc.run_parallel(task, rc.T3_STARS, workers=5)
    allr = {k: v for r in res for k, v in r.items()}
    (rc.HERE / "t3_study.json").write_text(json.dumps(allr))
    print("done")
