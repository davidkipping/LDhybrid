"""Re-run of the p2_exact weighting and limb studies for hybrid4 only, with the
curvature-simplex basis (v1.HYBRID4_BASIS, 2026-10-14), constraint = its simplex.

    .venv/bin/python studies/rev4/rv_p2.py weighting|limb
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / "rev"))
import p2_exact as P  # noqa: E402
import rc  # noqa: E402

P.LAWS = ["hybrid4", "claret4"]          # claret4 as an unchanged reference


def task_limb(job):
    # p2_exact.task_limb loops over a fixed tuple of laws; keep only hybrid4 rows
    rows = P.task_limb(job)
    for r in rows:
        r["dlc"] = {k: v for k, v in r["dlc"].items() if k.startswith("hybrid4|")}
        r["dw"] = {k: v for k, v in r["dw"].items() if k.startswith("hybrid4|")}
    return rows


if __name__ == "__main__":
    what = sys.argv[1]
    jobs = [j for j in rc.p2s_models() if j[0] == "phoenix"] if what == "limb" else rc.p2s_models()
    task = dict(weighting=P.task_weighting, limb=task_limb)[what]
    res = rc.run_parallel(task, jobs, workers=5)
    (HERE / f"rv_p2_{what}.json").write_text(json.dumps([r for rows in res for r in rows]))
    print("done", what)
