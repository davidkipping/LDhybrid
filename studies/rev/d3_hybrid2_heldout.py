"""D3: held-out test of hybrid2's pole (ε = 0.36, tuned in P2c on 7,524
PHOENIX + MPS-ATLAS Set 2 model-bands).

Re-uses the P2c light-curve errors (studies/p2/p2c_fine.jsonl, ε = 0.30–0.46 in
steps of 0.01, r-weighted). Objective as in P2c: the 90th percentile of the
worst light-curve error over model-bands (ties to the 99th).

    .venv/bin/python studies/rev/d3_hybrid2_heldout.py
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
EPS = np.round(np.arange(0.30, 0.4601, 0.01), 3)
rows = [json.loads(l) for l in open(HERE.parent / "p2" / "p2c_fine.jsonl")]
W = np.array([r["worst"] for r in rows])
kind = np.array([r["kind"] for r in rows])


def stats(e):
    return np.median(e), np.percentile(e, 90), np.percentile(e, 99)


def tune(Wt):
    return int(np.argmin(np.percentile(Wt, 90, 0)))


out = ["## D3 — hybrid2 pole, held-out", "",
       f"{len(rows)} model-bands (PHOENIX {np.sum(kind == 'phoenix')}, MPS-ATLAS Set 2 "
       f"{np.sum(kind == 'mps2')}), ε ∈ 0.30–0.46, objective: 90th percentile of the worst "
       "light-curve error (ppm). 'held-out' = evaluated on the set not used for tuning.", "",
       "| tune on | evaluate on | ε chosen | median | 90th | 99th | at ε = 0.36 (median / 90th / 99th) |",
       "|---|---|---|---|---|---|---|"]
k36 = int(np.argmin(np.abs(EPS - 0.36)))
for a, b in (("phoenix", "mps2"), ("mps2", "phoenix"), ("all", "all")):
    ta = W if a == "all" else W[kind == a]
    tb = W if b == "all" else W[kind == b]
    k = tune(ta)
    m, p90, p99 = stats(tb[:, k])
    m0, q0, r0 = stats(tb[:, k36])
    out.append(f"| {a} | {b} | {EPS[k]:.2f} | {m:.1f} | {p90:.1f} | {p99:.1f} | "
               f"{m0:.1f} / {q0:.1f} / {r0:.1f} |")

rng = np.random.default_rng(1)
# split by model (not by model-band), so a model's bands stay together
models = np.array([r["job_id"] for r in rows])
um = np.unique(models)
ch, d90, dlt = [], [], []
for _ in range(20):
    half = set(rng.choice(um, len(um) // 2, replace=False))
    tr = np.array([m in half for m in models])
    k = tune(W[tr])
    ch.append(EPS[k])
    held = W[~tr]
    d90.append(np.percentile(held[:, k], 90))
    dlt.append(np.percentile(held[:, k], 90) - np.percentile(held[:, tune(held)], 90))
out += ["", f"Random half/half splits by model (20×): chosen ε = {np.min(ch):.2f}–{np.max(ch):.2f} "
        f"(median {np.median(ch):.2f}); held-out 90th percentile {np.median(d90):.1f} ppm "
        f"(range {np.min(d90):.1f}–{np.max(d90):.1f}); excess over the held-out optimum "
        f"median {np.median(dlt):.2f} ppm, max {np.max(dlt):.2f} ppm.", "",
        "Sensitivity, all 7,524 model-bands:", "", "| ε | median | 90th | 99th |", "|---|---|---|---|"]
for e in (0.30, 0.33, 0.36, 0.39, 0.42):
    k = int(np.argmin(np.abs(EPS - e)))
    m, p90, p99 = stats(W[:, k])
    out.append(f"| {e:.2f} | {m:.1f} | {p90:.1f} | {p99:.1f} |")
q = np.array([r["quadratic"] for r in rows])
out += ["", f"Quadratic on the same set: median {np.median(q):.1f}, 90th {np.percentile(q, 90):.1f}, "
        f"99th {np.percentile(q, 99):.1f}."]
for kd in ("phoenix", "mps2"):
    s = kind == kd
    out.append(f"* {kd}: hybrid2(0.36) better than quadratic in "
               f"{100 * np.mean(W[s, k36] < q[s]):.0f}% of model-bands; quadratic median "
               f"{np.median(q[s]):.1f}, hybrid2 {np.median(W[s, k36]):.1f}.")
txt = "\n".join(out) + "\n"
(HERE / "d3.md").write_text(txt)
print(txt)
