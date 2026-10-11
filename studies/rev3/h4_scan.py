"""hybrid4 pole scan on the method study (P2S): 1-D ladder scale and 2-D (ε0, ε1).

    .venv/bin/python studies/rev3/h4_scan.py
"""
import json
import time

import numpy as np

from h4common import HERE, err_h4, pole_sets, rc

SETS = pole_sets()


def task(job):
    kind, key = job
    spec = rc.load(kind, key)
    g = rc.lcgrid()
    out = []
    for band in rc.P2S_BANDS:
        pr = rc.profile(spec, band)
        if pr is None:
            continue
        tr = g.columns(pr)
        errs = [err_h4(eps, pr, tr, g) for _, _, eps in SETS]
        ref = {n: rc.score(rc.fit_law(n, pr)[0], tr, g) for n in ("claret4", "hybrid5")}
        out.append(dict(kind=kind, key=[float(x) for x in key], band=band, errs=errs, ref=ref))
    return out


def main():
    t0 = time.time()
    rows = [r for res in rc.run_parallel(task, rc.p2s_models(), workers=10) for r in res]
    (HERE / "h4_scan.json").write_text(json.dumps(dict(
        sets=[[t, s, list(e)] for t, s, e in SETS], rows=rows)))
    print(f"{len(rows)} model-bands x {len(SETS)} pole sets in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
