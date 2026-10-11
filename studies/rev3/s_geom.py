"""D4 geometry (studies/rev/t3_study.py) for hybrid4 with the new poles.

    .venv/bin/python studies/rev3/s_geom.py
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / "rev"))
import rc  # noqa: E402
from ldhybrid import lightcurve as lc  # noqa: E402

P = (0.03, 0.1, 0.2, 0.3)


def task(star):
    sname, (t, g) = star
    spec = rc.load_t3(t, g)
    out = {}
    for bkey, blab in rc.T3_BANDS:
        pr = rc.profile(spec, bkey)
        f, w, _ = rc.fit_law("hybrid4", pr, "r")
        geo = {}
        for p in P:
            e = lc.lc_errors(f, pr, p_grid=(p,), b_grid=(0.0, 0.5, 1 - p / 2, 1.0))
            geo[str(p)] = dict(worst=e["worst"][0].tolist(), rms=e["rms"][0].tolist())
        out[f"{sname}|{blab}"] = dict(w=[float(x) for x in w], geo=geo)
    return out


def main():
    res = {k: v for r in rc.run_parallel(task, rc.T3_STARS, workers=5) for k, v in r.items()}
    (HERE / "s_geom.json").write_text(json.dumps(res))
    L = ["| p | b=0 | b=0.5 | b=1-p/2 | b=1 | rms b=0 |", "|---|---|---|---|---|---|"]
    for p in P:
        W = np.array([v["geo"][str(p)]["worst"] for v in res.values()])
        R = np.array([v["geo"][str(p)]["rms"] for v in res.values()])
        L.append(f"| {p} | " + " | ".join(f"{np.median(W[:, j]):.2g} ({W[:, j].max():.2g})" for j in range(4))
                 + f" | {np.median(R[:, 0]):.2g} |")
    ws = {k: [v["geo"][str(p)]["worst"][0] for p in P] for k, v in res.items()}
    slope = np.median([np.log(x[-1] / x[0]) / np.log(0.3 / 0.03) for x in ws.values()])
    L.append(f"\nlog-slope (p 0.03→0.3, b=0) median {slope:.2f}")
    (HERE / "s_geom.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
