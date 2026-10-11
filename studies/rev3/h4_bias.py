"""Equal-footing parameter bias (studies/rev2/e1b_bias.py set-up) for hybrid4
with candidate poles, against the nominal ladder and the four-term law.

    H4_EPS=0.0286,0.3858 .venv/bin/python studies/rev3/h4_bias.py  -> h4_bias_<tag>.json
"""
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / "rev2"))
sys.path.insert(0, str(HERE))
from h4common import hybrid4_law  # noqa: E402
from ldhybrid import laws  # noqa: E402

EPS = tuple(float(x) for x in os.environ["H4_EPS"].split(","))
laws.LAWS["hybrid4"] = lambda: hybrid4_law(EPS)          # patched in every worker on import

import numpy as np  # noqa: E402
import e1b_bias as E  # noqa: E402
import rc  # noqa: E402
from ldhybrid import lightcurve as lc  # noqa: E402
from ldhybrid.grids import mpsatlas  # noqa: E402

LAWS = ("hybrid4", "claret4")


def task(job):
    grid, sname, spec_key = job
    spec = rc.load_t3(*spec_key) if grid == "phoenix" else mpsatlas.load("set2", 0.0, *spec_key)
    out = {}
    for bkey, blab in rc.T3_BANDS:
        pr = rc.profile(spec, bkey)
        if pr is None:
            continue
        F_true = lc.light_curve(pr, E.z_of(E.B0, E.V0), E.P0)
        for name in LAWS:
            _, w_tab, _ = rc.fit_law(name, pr, "r")
            M = E.Model(name)
            res = {"fixed": E.summarise(*E.run_fixed(M, w_tab, F_true))}
            for mode in ("free", "prior"):
                res[mode] = E.summarise(*E.run_linear_free(M, name, np.array(w_tab, float), F_true, mode == "prior"))
            out[f"{grid}|{sname}|{blab}|{name}"] = res
    return out


def main():
    jobs = [("phoenix", n, k) for n, k in rc.T3_STARS]
    jobs += [("mps2", "early M", (3800, 5.0)), ("mps2", "K dwarf", (4600, 4.5)),
             ("mps2", "Sun analogue", (5800, 4.5))]
    res = {k: v for r in rc.run_parallel(task, jobs, workers=8) for k, v in r.items()}
    tag = "_".join(f"{e:.4g}" for e in EPS)
    (HERE / f"h4_bias_{tag}.json").write_text(json.dumps(res))
    for grid in ("phoenix", "mps2"):
        for name in LAWS:
            for mode in ("fixed", "free", "prior"):
                dp = np.abs([v[mode]["dp2_ppm"] for k, v in res.items()
                             if k.startswith(grid + "|") and k.endswith("|" + name)])
                print(f"{tag} {grid} {name} {mode}: |dp2| median {np.median(dp):.2f} max {dp.max():.1f}")


if __name__ == "__main__":
    main()
