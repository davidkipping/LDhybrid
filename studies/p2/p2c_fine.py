"""P2c follow-up: fine ε scan for hybrid2 around the coarse optimum, with
quadratic and power-2 scored on exactly the same model-bands (r-weighted).

    .venv/bin/python studies/p2/p2c_fine.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402
import p2c_hybrid2_pole as base  # noqa: E402

HERE = Path(__file__).parent
OUT = HERE / "p2c_fine.jsonl"
EPS = np.round(np.arange(0.30, 0.4601, 0.01), 3)
_COLS = {}


def task(job):
    from ldhybrid import bands, fit, laws
    from ldhybrid import lightcurve as lc
    from ldhybrid.profile import band_profile
    g = lc.LCGrid()
    spec = common.load(job["kind"], tuple(job["key"]))
    rows = []
    for b in base.bands_for(spec):
        pr = band_profile(spec, *bands.response(b), limb="M1")
        truth = g.columns(pr)
        worst = []
        for eps in EPS:
            if eps not in _COLS:
                _COLS[eps] = g.linear_columns(laws.hybrid2(float(eps)))
            f = fit.fit_linear(laws.hybrid2(float(eps)), pr, "r", 2000, "physical")
            worst.append(g.errors(*g.combine(*_COLS[eps], f.w), *truth)["max_ppm"])
        if "quad" not in _COLS:
            _COLS["quad"] = g.linear_columns(laws.quadratic())
        fq = fit.fit_linear(laws.quadratic(), pr, "r", 2000, "physical")
        quad = g.errors(*g.combine(*_COLS["quad"], fq.w), *truth)["max_ppm"]
        fp = fit.fit_power2(pr, "r", 2000)
        c, a = fp.w
        p2 = g.errors(*g.columns(lambda mu: laws.power2(c, a, mu)), *truth)["max_ppm"]
        rows.append(dict(job_id=job["job_id"], kind=job["kind"], band=b, worst=worst,
                         quadratic=quad, power2=p2))
    return rows


def report():
    rows = common.read(OUT)
    W = np.array([r["worst"] for r in rows])
    q = np.array([r["quadratic"] for r in rows])
    p2 = np.array([r["power2"] for r in rows])
    lines = ["# P2c (fine): hybrid2 pole vs quadratic and power-2", "",
             f"{len(rows)} model-bands; worst light-curve error per fit (ppm), r-weighted.", "",
             "| law | median | 90th | 99th | worst | beats quadratic |", "|---|---|---|---|---|---|"]
    for name, e in (("quadratic", q), ("power-2", p2)):
        lines.append(f"| {name} | {np.median(e):.1f} | {np.percentile(e, 90):.1f} | "
                     f"{np.percentile(e, 99):.1f} | {e.max():.1f} | |")
    for k, eps in enumerate(EPS):
        e = W[:, k]
        lines.append(f"| hybrid2 ε={eps:.2f} | {np.median(e):.1f} | {np.percentile(e, 90):.1f} | "
                     f"{np.percentile(e, 99):.1f} | {e.max():.1f} | {100 * np.mean(e < q):.0f}% |")
    (HERE / "p2c_fine.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    if "--report" not in sys.argv:
        common.run_parallel(task, base.jobs(), OUT)
    report()
