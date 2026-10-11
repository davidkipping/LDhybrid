"""C8: sensitivity of hybrid4/hybrid5 to the pole ladder (ln ε_k scaled by
0.8 and 1.2). C11: effect of imposing flux conservation (FCM) on the
r-weighted fits of every linear law. Method study (P2S), production metric.

    .venv/bin/python studies/rev2/c8_c11_p2s.py
"""
import json

import numpy as np
from squishierplanet import laws as sl

import rc
from ldhybrid import fit, laws, v1
from ldhybrid.products import constraint_for

SCALES = (0.8, 1.0, 1.2)
LIN = ["hybrid5", "hybrid4", "hybrid2", "claret4", "quadratic"]


def ladder_law(name, s):
    eps = {"hybrid4": v1.HYBRID4_EPS, "hybrid5": v1.HYBRID5_EPS}[name]
    eps = [float(np.exp(s * np.log(e))) for e in eps]
    shapes = [sl.Law([1.0, 0.0, -1.0]), sl.Law([1.0, -2.0, 1.0])] + [sl.pole_shape(e) for e in eps]
    return laws._from_shapes(f"{name}x{s}", [f"w{j}" for j in range(len(shapes))], shapes), eps


def task(job):
    kind, key = job
    spec = rc.load(kind, key)
    rows = []
    g = rc.lcgrid()
    for band in rc.P2S_BANDS:
        pr = rc.profile(spec, band)
        if pr is None:
            continue
        tr = g.columns(pr)
        row = dict(kind=kind, key=str(key), band=band, ladder={}, fcm={})
        for name in ("hybrid4", "hybrid5"):
            for s in SCALES:
                L, _ = ladder_law(name, s)
                w = fit.fit_linear(L, pr, "r", v1.RESAMPLE, "simplex").w
                row["ladder"][f"{name}|{s}"] = g.errors(*g.columns(lambda m, L=L, w=w: L.intensity(w, m)), *tr)["max_ppm"]
        for name in LIN:
            L = laws.get(name)
            f0 = fit.fit_linear(L, pr, "r", v1.RESAMPLE, constraint_for(name), False)
            f1 = fit.fit_linear(L, pr, "r", v1.RESAMPLE, constraint_for(name), True)
            e0 = g.errors(*g.columns(lambda m: L.intensity(f0.w, m)), *tr)["max_ppm"]
            e1 = g.errors(*g.columns(lambda m: L.intensity(f1.w, m)), *tr)["max_ppm"]
            row["fcm"][name] = dict(dw=float(np.max(np.abs(f1.w - f0.w))), lc0=e0, lc1=e1,
                                    flux_err0=f0.flux_err)
        rows.append(row)
    return rows


def main():
    rows = [r for res in rc.run_parallel(task, rc.p2s_models()) for r in res]
    (rc.HERE / "c8_c11_p2s.json").write_text(json.dumps(rows))
    L = ["## C8 — pole-ladder sensitivity (ln ε_k × s; P2S, production metric, median / 90th ppm)", "",
         "| law | grid | s = 0.8 | s = 1 (nominal) | s = 1.2 |", "|---|---|---|---|---|"]
    for name in ("hybrid4", "hybrid5"):
        for kind in ("phoenix", "mps2"):
            sel = [r for r in rows if r["kind"] == kind]
            cells = []
            for s in SCALES:
                x = np.array([r["ladder"][f"{name}|{s}"] for r in sel])
                cells.append(f"{np.median(x):.2f} / {np.percentile(x, 90):.1f}")
            L.append(f"| {name} | {kind} | " + " | ".join(cells) + " |")
    L += ["", "Pole positions: " + "; ".join(
        f"{n} s={s}: " + ", ".join(f"{e:.4g}" for e in ladder_law(n, s)[1]) for n in ("hybrid4", "hybrid5") for s in SCALES), ""]
    L += ["## C11 — flux conservation (FCM) on P2S", "",
          "| law | grid | median / max |Δw| | median Δ lc (ppm) | max |Δ lc| | median |flux err| without FCM |",
          "|---|---|---|---|---|---|"]
    for name in LIN:
        for kind in ("phoenix", "mps2"):
            sel = [r["fcm"][name] for r in rows if r["kind"] == kind]
            dw = np.array([s["dw"] for s in sel])
            dl = np.array([s["lc1"] - s["lc0"] for s in sel])
            fe = np.array([abs(s["flux_err0"]) for s in sel])
            L.append(f"| {name} | {kind} | {np.median(dw):.2g} / {dw.max():.2g} | {np.median(dl):+.2f} | "
                     f"{np.abs(dl).max():.1f} | {np.median(fe):.1e} |")
    (rc.HERE / "c8_c11_p2s.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
