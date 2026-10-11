"""C8b: wider scan of the ladder scale s (ln ε_k × s) for hybrid4 and hybrid5
on P2S, to locate the optimum (production metric).

    .venv/bin/python studies/rev2/c8b_ladder_scan.py
"""
import json

import numpy as np

import rc
from c8_c11_p2s import ladder_law
from ldhybrid import fit, v1

SCALES = (0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1)


def task(job):
    kind, key = job
    spec = rc.load(kind, key)
    g = rc.lcgrid()
    rows = []
    for band in rc.P2S_BANDS:
        pr = rc.profile(spec, band)
        if pr is None:
            continue
        tr = g.columns(pr)
        row = dict(kind=kind, band=band)
        for name in ("hybrid4", "hybrid5"):
            for s in SCALES:
                L, _ = ladder_law(name, s)
                w = fit.fit_linear(L, pr, "r", v1.RESAMPLE, "simplex").w
                row[f"{name}|{s}"] = g.errors(*g.columns(lambda m, L=L, w=w: L.intensity(w, m)), *tr)["max_ppm"]
        rows.append(row)
    return rows


def main():
    rows = [r for res in rc.run_parallel(task, rc.p2s_models()) for r in res]
    (rc.HERE / "c8b_ladder_scan.json").write_text(json.dumps(rows))
    L = ["## C8b — ladder scale scan (P2S, production metric, median / 90th ppm)", "",
         "| law | grid | " + " | ".join(f"s={s}" for s in SCALES) + " |", "|---|---|" + "---|" * len(SCALES)]
    for name in ("hybrid4", "hybrid5"):
        for kind in ("phoenix", "mps2", "all"):
            sel = [r for r in rows if kind == "all" or r["kind"] == kind]
            L.append(f"| {name} | {kind} | " + " | ".join(
                f"{np.median([r[f'{name}|{s}'] for r in sel]):.2f} / {np.percentile([r[f'{name}|{s}'] for r in sel], 90):.1f}"
                for s in SCALES) + " |")
    L.append("")
    L.append("Pole positions: " + "; ".join(f"{n} s={s}: " + ", ".join(f"{e:.4g}" for e in ladder_law(n, s)[1])
                                            for n in ("hybrid4", "hybrid5") for s in SCALES))
    (rc.HERE / "c8b_ladder_scan.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
