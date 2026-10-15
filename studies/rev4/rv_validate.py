"""Out-of-sample check of the hybrid4 basis change: the same 100 PHOENIX + 100
MPS-ATLAS Set 2 models (all passbands) as studies/rev3/h4_validate.py, with
hybrid4 fitted in the curvature-simplex basis (production) and, for comparison,
on the monomial-pole shapes constrained to the w-simplex (the previous tables).

    .venv/bin/python studies/rev4/rv_validate.py   -> rv_validate.md
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / "rev3"))
import h4_validate as V  # noqa: E402
from h4common import NOMINAL, err_h4, rc  # noqa: E402


def task(job):
    kind, key, bands = job
    spec = rc.load(kind, key)
    g = rc.lcgrid()
    out = []
    for band in bands:
        pr = rc.profile(spec, band)
        if pr is None:
            continue
        tr = g.columns(pr)
        r = {"hybrid4 curvature": rc.score(rc.fit_law("hybrid4", pr)[0], tr, g),
             "hybrid4 w-simplex": err_h4(NOMINAL, pr, tr, g)}
        r.update({n: rc.score(rc.fit_law(n, pr)[0], tr, g) for n in ("claret4", "hybrid5")})
        out.append(dict(kind=kind, band=band, **r))
    return out


def main():
    t0 = time.time()
    rows = [r for res in rc.run_parallel(task, V.jobs(), workers=8) for r in res]
    (HERE / "rv_validate.json").write_text(json.dumps(rows))
    names = ["hybrid4 curvature", "hybrid4 w-simplex", "claret4", "hybrid5"]
    L = ["# hybrid4 basis change, out of sample (100 PHOENIX + 100 MPS-ATLAS 2 models, all passbands)", "",
         "| law | all | PHOENIX | MPS-ATLAS 2 |", "|---|---|---|---|"]
    for n in names:
        cells = []
        for sel in (None, "phoenix", "mps2"):
            x = np.array([r[n] for r in rows if sel is None or r["kind"] == sel])
            cells.append(f"{np.median(x):.2f} / {np.percentile(x, 90):.1f} / {np.percentile(x, 99):.0f}")
        L.append(f"| {n} | " + " | ".join(cells) + " |")
    d = np.array([r["hybrid4 curvature"] - r["hybrid4 w-simplex"] for r in rows])
    L += ["", f"curvature − w-simplex per model-band: median {np.median(d):+.3f}, min {d.min():+.2f}, "
          f"max {d.max():+.2f} ppm; differs by > 0.01 ppm in {np.mean(np.abs(d) > 0.01) * 100:.2f}%.",
          "", f"(median / 90th / 99th, ppm; {len(rows)} model-bands; {(time.time() - t0) / 60:.1f} min)"]
    (HERE / "rv_validate.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
