"""Out-of-sample validation of candidate hybrid4 poles: 100 PHOENIX and 100
MPS-ATLAS Set 2 models outside the method study, all their passbands.

    .venv/bin/python studies/rev3/h4_validate.py
"""
import json
import time

import h5py
import numpy as np

from h4common import HERE, NOMINAL, err_h4, rc
from ldhybrid.paths import DATA

CANDS = dict(nominal=NOMINAL, **json.loads((HERE / "h4_candidates.json").read_text()))


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
        r = {c: err_h4(tuple(e), pr, tr, g) for c, e in CANDS.items()}
        r.update({n: rc.score(rc.fit_law(n, pr)[0], tr, g) for n in ("claret4", "hybrid5")})
        out.append(dict(kind=kind, band=band, **r))
    return out


def jobs():
    rng = np.random.default_rng(42)
    p2s = {(k, tuple(float(x) for x in key)) for k, key in rc.p2s_models()}
    out = []
    for grid, kind in (("phoenix", "phoenix"), ("mps2", "mps2")):
        h = h5py.File(DATA / "products" / f"{grid}.h5")
        T, G, M = h["teff"][()], h["logg"][()], h["mh"][()]
        bands = [b.decode() for b in h["bands"][()]]
        idx = rng.permutation(len(T))
        n = 0
        for i in idx:
            key = (int(T[i]), float(G[i]), "-0.0" if M[i] == 0 else f"{M[i]:+.1f}") if kind == "phoenix" \
                else (float(M[i]), int(T[i]), float(G[i]))
            kk = (kind, tuple(float(x) for x in key))
            if kk in p2s:
                continue
            out.append((kind, key, bands))
            n += 1
            if n == 100:
                break
    return out


def main():
    t0 = time.time()
    rows = [r for res in rc.run_parallel(task, jobs(), workers=10) for r in res]
    (HERE / "h4_validate.json").write_text(json.dumps(rows))
    names = list(CANDS) + ["claret4", "hybrid5"]
    L = ["# hybrid4 candidate poles, out of sample (100 PHOENIX + 100 MPS-ATLAS 2 models, all passbands)", "",
         "| law / poles | ε | all | PHOENIX | MPS-ATLAS 2 |", "|---|---|---|---|---|"]
    for n in names:
        eps = ", ".join(f"{e:.4g}" for e in CANDS[n]) if n in CANDS else ""
        cells = []
        for sel in (None, "phoenix", "mps2"):
            x = np.array([r[n] for r in rows if sel is None or r["kind"] == sel])
            cells.append(f"{np.median(x):.2f} / {np.percentile(x, 90):.1f} / {np.percentile(x, 99):.0f}")
        L.append(f"| {n} | {eps} | " + " | ".join(cells) + " |")
    L += ["", f"(median / 90th / 99th, ppm; {len(rows)} model-bands; {(time.time() - t0) / 60:.1f} min)"]
    (HERE / "h4_validate.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
