"""P2a — fit weighting. For every subset model × band × law, fit the profile
with each weighting and score in light-curve space over the (p, b) grid.

Weightings: measure ∈ {mu, r, area} × quadrature ∈ {native trapezoid, 2000-pt
PCHIP resample}, plus Claret/Kostogryz equal weights on native nodes.
Constraints fixed: simplex for hybrid4/5, exact physical region otherwise.

    .venv/bin/python studies/p2/p2a_weighting.py      (then  --report)
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402

HERE = Path(__file__).parent
OUT = HERE / "p2a_weighting.jsonl"
LAWS = ["quadratic", "power2", "hybrid2", "claret4", "hybrid4", "hybrid5"]
METHODS = [(m, rs) for m in ("mu", "r", "area") for rs in (None, 2000)] + [("nodes", None)]


def task(job):
    from ldhybrid import pipeline
    spec = common.load(job["kind"], tuple(job["key"]))
    rows = []
    for measure, rs in METHODS:
        for r in pipeline.fit_model(spec, job["band"], LAWS, limb="M1", measure=measure,
                                    resample=rs, constraint="simplex"):
            r.update(job_id=job["job_id"], kind=job["kind"], resample=rs)
            rows.append(r)
    return rows


def jobs():
    return [dict(job_id=f"{k}|{key}|{b}", kind=k, key=list(key), band=b)
            for k, key in common.models() for b in common.BANDS]


def report():
    rows = common.read(OUT)
    lines = ["# P2a: fit weighting", "",
             "Light-curve error (ppm) of profile fits, over the (p, b) grid "
             "(p ∈ {0.03, 0.1, 0.15}, b ∈ {0, …, 0.9}). Each cell: median over "
             "models × bands of the per-fit *worst* error / of the per-fit "
             "b-averaged rms. M1 limb cut for PHOENIX.", ""]
    for kind in ("phoenix", "mps2"):
        R = [r for r in rows if r["kind"] == kind]
        n = len({r["job_id"] for r in R})
        lines += [f"## {kind} ({n} model-bands)", "",
                  "| weighting | " + " | ".join(LAWS) + " |", "|---|" + "---|" * len(LAWS)]
        for m, rs in METHODS:
            S = [r for r in R if r["measure"] == m and r["resample"] == rs]
            cells = []
            for l in LAWS:
                e = np.array([r["lc_max_ppm"] for r in S if r["law"] == l])
                q = np.array([r["lc_mean_rms_ppm"] for r in S if r["law"] == l])
                cells.append(f"{np.median(e):.2f} / {np.median(q):.2f}")
            lab = f"{m}, {'native' if rs is None else 'resample'}" if m != "nodes" else "nodes (Claret)"
            lines.append(f"| {lab} | " + " | ".join(cells) + " |")
        lines += ["", "Worst case over all model-bands (ppm):", "",
                  "| weighting | " + " | ".join(LAWS) + " |", "|---|" + "---|" * len(LAWS)]
        for m, rs in METHODS:
            S = [r for r in R if r["measure"] == m and r["resample"] == rs]
            lab = f"{m}, {'native' if rs is None else 'resample'}" if m != "nodes" else "nodes (Claret)"
            lines.append(f"| {lab} | " + " | ".join(
                f"{max(r['lc_max_ppm'] for r in S if r['law'] == l):.1f}" for l in LAWS) + " |")
        lines.append("")
    (HERE / "p2a_weighting.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    if "--report" in sys.argv:
        report()
    else:
        common.run_parallel(task, jobs(), OUT)
        report()
