"""C14: grazing transits relative to the transit depth (T3 cases).

For p = 0.1 and 0.2 at b = 0 and b = 1 − p/2: the depth of the model light
curve (largest blocked fraction) and each law's worst |ΔF| as a fraction of
that depth.

    .venv/bin/python studies/rev2/c14_grazing.py
"""
import json

import numpy as np

import rc
from ldhybrid import lightcurve as lc


def task(job):
    name, key = job
    spec = rc.load_t3(*key)
    rows = []
    for b, blab in rc.T3_BANDS:
        pr = rc.profile(spec, b)
        fits = {n: rc.fit_law(n, pr)[0] for n in rc.LAWS6}
        for p in (0.1, 0.2):
            for bb in (0.0, 1 - p / 2):
                z = lc.chord_z(p, bb, 240)
                Ft = lc.light_curve(pr, z, p)
                depth = float(1 - Ft.min()) * 1e6
                e = {n: float(np.max(np.abs(lc.light_curve(f, z, p) - Ft)) * 1e6) for n, f in fits.items()}
                rows.append(dict(case=name, band=blab, p=p, b=bb, depth_ppm=depth, err=e))
    return rows


def main():
    rows = [r for res in rc.run_parallel(task, rc.T3_STARS) for r in res]
    (rc.HERE / "c14_grazing.json").write_text(json.dumps(rows))
    L = ["## C14 — grazing transits relative to depth (T3, median over 15 cases)", "",
         "| p | b | depth (ppm) | " + " | ".join(f"{n} err/depth" for n in rc.LAWS6) + " |",
         "|---|---|---|" + "---|" * len(rc.LAWS6)]
    for p in (0.1, 0.2):
        for bb in (0.0, 1 - p / 2):
            sel = [r for r in rows if r["p"] == p and abs(r["b"] - bb) < 1e-9]
            d = np.median([r["depth_ppm"] for r in sel])
            fr = [np.median([r["err"][n] / r["depth_ppm"] for r in sel]) for n in rc.LAWS6]
            L.append(f"| {p} | {bb:.2f} | {d:.0f} | " + " | ".join(f"{1e2 * f:.3f}%" for f in fr) + " |")
    L.append("")
    L.append("Absolute worst errors (ppm, median): " + "; ".join(
        f"p={p} b={bb:.2f}: " + ", ".join(
            f"{n} {np.median([r['err'][n] for r in rows if r['p'] == p and abs(r['b'] - bb) < 1e-9]):.2g}"
            for n in rc.LAWS6) for p in (0.1, 0.2) for bb in (0.0, 1 - p / 2)))
    (rc.HERE / "c14_grazing.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
