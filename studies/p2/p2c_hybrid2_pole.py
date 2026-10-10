"""P2c — re-tune hybrid2's single pole over the dwarf grid and all bands.

hybrid2: I/I(1) = 1 − w1 (1 − μ²) − w2 Π_ε(μ). For each ε on a log grid, every
sample profile is fitted (r- and area-weighted, exact physical triangle) and
scored in light-curve space; we report the median, 90th percentile and worst
error over all model-bands as functions of ε. The truth light curve of each
profile and the shape columns of each ε are computed once.

Sample: PHOENIX v3.0 (M1) Teff 2300–7500 step 400 × log g {4.0, 4.5, 5.0} ×
[M/H] {−1, 0, +0.5}, and MPS-ATLAS Set 2 Teff 3600–7200 step 400 × the same
log g and [M/H]; every registered band the models cover (MIRI excluded).

    .venv/bin/python studies/p2/p2c_hybrid2_pole.py   (then --report)
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402

HERE = Path(__file__).parent
OUT = HERE / "p2c_hybrid2_pole.jsonl"
EPS = np.round(np.geomspace(0.02, 2.0, 41), 5)
MEASURES = ("r", "area")


def sample():
    from ldhybrid.grids import mpsatlas, phoenix
    cat = {(e["teff"], e["logg"], e["mhs"]) for e in
           json.load(open(phoenix.REDUCED / "phoenix" / "catalogue.json"))}
    out = []
    for t in range(2300, 7600, 400):
        for g in (4.0, 4.5, 5.0):
            for mh in ("-1.0", "-0.0", "+0.5"):
                if (t, g, mh) in cat:
                    out.append(("phoenix", (t, g, mh)))
    mps = set(mpsatlas.catalogue("set2"))
    for t in range(3600, 7300, 400):
        for g in (4.0, 4.5, 5.0):
            for mh in (-1.0, 0.0, 0.5):
                if (mh, t, g) in mps:
                    out.append(("mps2", (mh, t, g)))
    return out


def bands_for(spec):
    from ldhybrid import bands
    return [b for b in bands.BANDS if b != "JWST_MIRI_LRS"
            and spec.coverage(*bands.response(b)) >= 0.999]


_COLS = {}


def _columns(g, eps):
    from ldhybrid import laws
    if eps not in _COLS:
        _COLS[eps] = g.linear_columns(laws.hybrid2(float(eps)))
    return _COLS[eps]


def task(job):
    from ldhybrid import bands, fit, laws
    from ldhybrid import lightcurve as lc
    from ldhybrid.profile import band_profile
    g = lc.LCGrid()
    spec = common.load(job["kind"], tuple(job["key"]))
    rows = []
    for b in bands_for(spec):
        pr = band_profile(spec, *bands.response(b), limb="M1")
        truth = g.columns(pr)
        for measure in MEASURES:
            worst = []
            for eps in EPS:
                law = laws.hybrid2(float(eps))
                f = fit.fit_linear(law, pr, measure, 2000, "physical")
                worst.append(g.errors(*g.combine(*_columns(g, eps), f.w), *truth)["max_ppm"])
            rows.append(dict(job_id=job["job_id"], kind=job["kind"], key=job["key"],
                             band=b, measure=measure, worst=worst))
    return rows


def jobs():
    return [dict(job_id=f"{k}|{key}", kind=k, key=list(key)) for k, key in sample()]


def report():
    rows = common.read(OUT)
    lines = ["# P2c: hybrid2 pole re-tune", "",
             f"{len({r['job_id'] for r in rows})} models, "
             f"{len(rows) // len(MEASURES)} model-bands. Worst light-curve error (ppm) per fit "
             "over the (p, b) grid; statistics over model-bands.", ""]
    best = {}
    for m in MEASURES:
        W = np.array([r["worst"] for r in rows if r["measure"] == m])
        med, p90, mx = np.median(W, 0), np.percentile(W, 90, 0), W.max(0)
        lines += [f"## {m}-weighted", "", "| ε | median | 90th | worst |", "|---|---|---|---|"]
        for k, e in enumerate(EPS):
            mark = " **" if k in (np.argmin(med), np.argmin(p90), np.argmin(mx)) else ""
            lines.append(f"| {e:.4g}{mark} | {med[k]:.1f} | {p90[k]:.1f} | {mx[k]:.1f} |")
        best[m] = dict(median=float(EPS[np.argmin(med)]), p90=float(EPS[np.argmin(p90)]),
                       worst=float(EPS[np.argmin(mx)]))
        k0 = np.argmin(np.abs(EPS - 0.208))
        lines += ["", f"Minimisers: median at ε = {best[m]['median']}, 90th at "
                  f"{best[m]['p90']}, worst at {best[m]['worst']}. Current ε = 0.208 "
                  f"(nearest grid {EPS[k0]}): median {med[k0]:.1f}, 90th {p90[k0]:.1f}, "
                  f"worst {mx[k0]:.1f}.", ""]
        for kind in ("phoenix", "mps2"):
            Wk = np.array([r["worst"] for r in rows if r["measure"] == m and r["kind"] == kind])
            lines.append(f"* {kind}: 90th-percentile minimiser ε = "
                         f"{EPS[np.argmin(np.percentile(Wk, 90, 0))]:.4g}")
        lines.append("")
    (HERE / "p2c_hybrid2_pole.md").write_text("\n".join(lines) + "\n")
    (HERE / "p2c_best.json").write_text(json.dumps(best, indent=1))
    print("\n".join(lines))


if __name__ == "__main__":
    if "--report" in sys.argv:
        report()
    else:
        common.run_parallel(task, jobs(), OUT)
        report()
