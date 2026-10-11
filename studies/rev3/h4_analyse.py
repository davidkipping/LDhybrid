"""Analyse the hybrid4 pole scan: in-sample optimum, held-out tests, candidates.

    .venv/bin/python studies/rev3/h4_analyse.py   -> h4_analysis.md, h4_candidates.json
"""
import json

import numpy as np

from h4common import HERE, NOMINAL

d = json.loads((HERE / "h4_scan.json").read_text())
sets, rows = d["sets"], d["rows"]
E = np.array([r["errs"] for r in rows])                       # (n_mb, n_sets)
kind = np.array([r["kind"] for r in rows])
key = [tuple(r["key"]) + (r["kind"],) for r in rows]
ref = {n: np.array([r["ref"][n] for r in rows]) for n in ("claret4", "hybrid5")}
i_nom = next(i for i, (t, s, e) in enumerate(sets) if t == "s" and abs(s - 1.0) < 1e-9) \
    if any(t == "s" and abs(s - 1.0) < 1e-9 for t, s, _ in sets) else None
nom_err = None


def stat(x):
    return np.median(x), np.percentile(x, 90)


def objective(sub):
    """Median plus a tenth of the 90th percentile (tail-aware), per pole set."""
    return np.median(E[sub], axis=0) + 0.1 * np.percentile(E[sub], 90, axis=0)


L = ["# hybrid4 pole study (method study, production metric, ppm: median / 90th)", ""]
allm = np.ones(len(rows), bool)
best = {}
for name, sub in (("all", allm), ("phoenix", kind == "phoenix"), ("mps2", kind == "mps2")):
    j = int(np.argmin(objective(sub)))
    best[name] = j
    L.append(f"* tuned on {name}: ε = ({sets[j][2][0]:.4g}, {sets[j][2][1]:.4g}) [{sets[j][0]} {sets[j][1]}], "
             f"in-sample {stat(E[sub, j])[0]:.2f} / {stat(E[sub, j])[1]:.1f}")
# 1-D scan table
L += ["", "## Ladder scale s (α = 1/s²)", "", "| s | ε | all | phoenix | mps2 |", "|---|---|---|---|---|"]
for j, (t, s, e) in enumerate(sets):
    if t != "s":
        continue
    L.append(f"| {s:.2f} | {e[0]:.4g}, {e[1]:.4g} | " + " | ".join(
        "{:.2f} / {:.1f}".format(*stat(E[m, j])) for m in (allm, kind == "phoenix", kind == "mps2")) + " |")
# held out by grid
L += ["", "## Held out by grid", ""]
for tr, te in (("phoenix", "mps2"), ("mps2", "phoenix")):
    j = best[tr]
    L.append(f"* tune {tr} → test {te}: ε = ({sets[j][2][0]:.4g}, {sets[j][2][1]:.4g}); "
             f"{te} error {stat(E[kind == te, j])[0]:.2f} / {stat(E[kind == te, j])[1]:.1f} "
             f"(best possible on {te}: {stat(E[kind == te, best[te]])[0]:.2f})")
# random half splits by model
rng = np.random.default_rng(0)
models = sorted(set(key[i][:3] + (key[i][3],) for i in range(len(rows))))
mid = np.array([models.index(key[i][:3] + (key[i][3],)) for i in range(len(rows))])
ex, chosen = [], []
for _ in range(40):
    half = rng.permutation(len(models))[: len(models) // 2]
    tr = np.isin(mid, half)
    j = int(np.argmin(objective(tr)))
    chosen.append(sets[j][2])
    jo = int(np.argmin(objective(~tr)))
    ex.append(np.median(E[~tr, j]) - np.median(E[~tr, jo]))
ch = np.array(chosen)
L += ["", f"* 40 random half-splits by model: chosen ε0 {np.percentile(ch[:, 0], [10, 50, 90]).round(4)}, "
      f"ε1 {np.percentile(ch[:, 1], [10, 50, 90]).round(3)}; held-out excess over the held-out optimum: "
      f"median {np.median(ex):.2f}, 90th {np.percentile(ex, 90):.2f} ppm"]
L += ["", "## References on the same model-bands", "",
      "* nominal hybrid4 ({:.4g}, {:.4g}): ".format(*NOMINAL) + "{:.2f} / {:.1f}".format(*stat(
          E[:, [j for j, (t, s, e) in enumerate(sets) if t == "s" and abs(s - 0.95) < 1e-9][0]]))
      + " (s = 0.95 shown; nominal from h4_validate)",
      "* four-term: all {:.2f} / {:.1f}, phoenix {:.2f} / {:.1f}, mps2 {:.2f} / {:.1f}".format(
          *stat(ref["claret4"]), *stat(ref["claret4"][kind == "phoenix"]), *stat(ref["claret4"][kind == "mps2"])),
      "* hybrid5: all {:.2f} / {:.1f}".format(*stat(ref["hybrid5"]))]
s08 = [j for j, (t, s, e) in enumerate(sets) if t == "s" and abs(s - 0.80) < 1e-9][0]
cands = dict(best_all=sets[best["all"]][2], s080=sets[s08][2],
             best_phoenix=sets[best["phoenix"]][2], best_mps2=sets[best["mps2"]][2])
(HERE / "h4_candidates.json").write_text(json.dumps(cands, indent=1))
(HERE / "h4_analysis.md").write_text("\n".join(L) + "\n")
print("\n".join(L))
