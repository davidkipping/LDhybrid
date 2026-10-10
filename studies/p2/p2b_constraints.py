"""P2b — physical constraints and flux conservation.

For each subset model × band and law, fit with
  none | exact physical region | simplex (hybrid4/5)   ×   fcm off/on,
r- and area-weighted (2000-pt resample), and record: light-curve error over
the (p, b) grid, whether the unconstrained fit is physical, whether it lies in
the simplex, and which constraints bind.

    .venv/bin/python studies/p2/p2b_constraints.py   (then --report)
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402

HERE = Path(__file__).parent
OUT = HERE / "p2b_constraints.jsonl"
LAWS = ["hybrid2", "hybrid4", "hybrid5", "claret4", "quadratic"]
MEASURES = ("r", "area")
_CACHE = {}


def task(job):
    from ldhybrid import bands, fit, laws
    from ldhybrid import lightcurve as lc
    from ldhybrid.profile import band_profile
    g = lc.LCGrid()
    spec = common.load(job["kind"], tuple(job["key"]))
    pr = band_profile(spec, *bands.response(job["band"]), limb="M1")
    truth = g.columns(pr)
    rows = []
    for name in LAWS:
        L = laws.get(name)
        if name not in _CACHE:
            _CACHE[name] = g.linear_columns(L)
        cons = ["none", "physical"] + (["simplex"] if name in ("hybrid4", "hybrid5") else [])
        for m in MEASURES:
            for con in cons:
                for fcm in (False, True):
                    f = fit.fit_linear(L, pr, m, 2000, con, fcm)
                    e = g.errors(*g.combine(*_CACHE[name], f.w), *truth)
                    rows.append(dict(job_id=job["job_id"], kind=job["kind"], key=job["key"],
                                     band=job["band"], law=name, measure=m, constraint=con,
                                     fcm=fcm, w=f.w.tolist(), lc_max_ppm=e["max_ppm"],
                                     lc_mean_rms_ppm=e["mean_rms_ppm"], flux_err=f.flux_err,
                                     physical=f.extra["physical"],
                                     in_simplex=bool(np.all(f.w >= -1e-12) and f.w.sum() <= 1 + 1e-12),
                                     active=f.active))
    return rows


def jobs():
    return [dict(job_id=f"{k}|{key}|{b}", kind=k, key=list(key), band=b)
            for k, key in common.models() for b in common.BANDS]


def report():
    rows = common.read(OUT)
    n = len({r["job_id"] for r in rows})
    lines = ["# P2b: constraints and flux conservation", "",
             f"{n} model-bands (PHOENIX M1 + MPS-ATLAS Set 2). Light-curve error (ppm), "
             "median / worst over model-bands.", ""]
    for m in MEASURES:
        lines += [f"## {m}-weighted", "", "| law | constraint | fcm | median | worst | "
                  "unconstrained fit physical | in simplex | binds |", "|---|---|---|---|---|---|---|---|"]
        for l in LAWS:
            for con in ("none", "physical", "simplex"):
                for fcm in (False, True):
                    S = [r for r in rows if r["law"] == l and r["measure"] == m
                         and r["constraint"] == con and r["fcm"] == fcm]
                    if not S:
                        continue
                    e = np.array([r["lc_max_ppm"] for r in S])
                    phys = np.mean([r["physical"] for r in S]) if con == "none" else np.nan
                    simp = np.mean([r["in_simplex"] for r in S]) if con == "none" else np.nan
                    binds = np.mean([bool([a for a in r["active"] if a != "fcm"]) for r in S])
                    lines.append(f"| {l} | {con} | {'yes' if fcm else 'no'} | {np.median(e):.2f} | "
                                 f"{e.max():.1f} | {'' if np.isnan(phys) else f'{100 * phys:.0f}%'} | "
                                 f"{'' if np.isnan(simp) else f'{100 * simp:.0f}%'} | {100 * binds:.0f}% |")
        lines.append("")
    (HERE / "p2b_constraints.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    if "--report" in sys.argv:
        report()
    else:
        common.run_parallel(task, jobs(), OUT)
        report()
