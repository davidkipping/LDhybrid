"""D12b: is the M-dwarf quadratic deficit caused by the steep limb? Fit each law
r-weighted to μ ≥ 0.05 only, then score the light curve with the model kept
below μ = 0.05 (so only the disc interior is compared). T3 cases, b = 0, p = 0.1.

    .venv/bin/python studies/rev/d12b_limb_excluded.py     # writes d12b.json
"""
import json

import numpy as np

import rc
from ldhybrid import fit, laws, v1
from ldhybrid import lightcurve as lc

MU = 0.05
LAWS = ["quadratic", "power2", "claret4", "hybrid2", "hybrid4", "hybrid5"]
out = {}
for sname, (t, g) in rc.T3_STARS:
    spec = rc.load_t3(t, g)
    for bkey, blab in rc.T3_BANDS:
        pr = rc.profile(spec, bkey)
        b0 = dict(p_grid=(0.1,), b_grid=(0.0,))
        for n in LAWS:
            if n == "power2":
                c, a = fit.fit_power2(pr, "r", v1.RESAMPLE, mu_min=MU).w
                f = lambda m, c=c, a=a: laws.power2(c, a, m)
            else:
                L = laws.get(n)
                w = fit.fit_linear(L, pr, "r", v1.RESAMPLE, rc.constraint_for(n), mu_min=MU).w
                f = lambda m, L=L, w=w: L.intensity(w, m)
            f_in = lambda m, f=f: np.where(np.asarray(m) < MU, pr(m), f(m))
            out[f"{sname}|{blab}|{n}"] = dict(interior=lc.lc_errors(f_in, pr, **b0)["max_ppm"],
                                              full=lc.lc_errors(f, pr, **b0)["max_ppm"])
(rc.HERE / "d12b.json").write_text(json.dumps(out))
for n in LAWS:
    v = [x for k, x in out.items() if k.endswith("|" + n)]
    vm = [x for k, x in out.items() if k.endswith("|" + n) and " M|" in k]
    print(n, "interior-only error, M dwarfs: %.1f–%.1f, all: median %.1f; full LC error of these fits: median %.0f"
          % (min(x["interior"] for x in vm), max(x["interior"] for x in vm),
             np.median([x["interior"] for x in v]), np.median([x["full"] for x in v])))
