"""D2 / D6 with the exact quadrature, recomputed from the model atmospheres (not
from the products): every PHOENIX model × its 38 bands, and a random sample of
300 MPS-ATLAS Set 2 models × 38 bands. Per row, r-weighted production fits of
the six laws, their light-curve error (production grid), whether the
constraint is active, power-2 refitted inside 0 ≤ c ≤ 1 when the LM fit is
unphysical, and the four-term law refitted without its constraint when the
constraint is active. Uniform-μ fits: constraint-active flags only.

    .venv/bin/python studies/rev/d2_d6_full.py        # writes d2_d6_full.jsonl
"""
import json

import numpy as np
from scipy.optimize import least_squares

import rc
from ldhybrid import bands, fit, laws, v1
from ldhybrid.profile import quadrature

LAWS = ["hybrid5", "hybrid4", "hybrid2", "claret4", "quadratic", "power2"]
BANDS = [b for b in bands.BANDS if b != "JWST_MIRI_LRS" and b not in getattr(bands, "DERIVED", {})]


def power2_bounded(pr, measure="r"):
    mu_q, y_q, q = quadrature(pr, measure, v1.RESAMPLE)
    sq = np.sqrt(q)
    res = least_squares(lambda p: sq * (laws.power2(p[0], p[1], mu_q) - y_q), (0.7, 0.6),
                        bounds=([0.0, 1e-6], [1.0, np.inf]), method="trf",
                        xtol=1e-14, ftol=1e-14, gtol=1e-14)
    return res.x


def task(job):
    kind, key = job
    spec = rc.load(kind, key)
    g = rc.lcgrid()
    out = []
    for band in BANDS:
        pr = rc.profile(spec, band)
        if pr is None:
            continue
        T = g.columns(pr)
        row = dict(kind=kind, key=list(key), band=band, lc={}, bound={}, bound_mu={})
        for n in LAWS:
            f, w, F = rc.fit_law(n, pr, "r")
            row["lc"][n] = rc.score(f, T, g)
            if n == "power2":
                c, a = w
                phys = bool(0 <= c <= 1 and a > 0)
                row["bound"][n] = not phys                 # 'unphysical' for power-2
                cm, am = fit.fit_power2(pr, "mu", v1.RESAMPLE).w
                row["bound_mu"][n] = not (0 <= cm <= 1 and am > 0)
                if not phys:
                    c1, a1 = power2_bounded(pr)
                    row["power2_c"] = [float(c), float(c1)]
                    row["lc"]["power2_con"] = rc.score(lambda m: laws.power2(c1, a1, m), T, g)
                continue
            row["bound"][n] = bool(F.active)
            _, _, Fm = rc.fit_law(n, pr, "mu")
            row["bound_mu"][n] = bool(Fm.active)
            if n == "claret4" and F.active:
                L = laws.get("claret4")
                w0 = fit.fit_linear(L, pr, "r", v1.RESAMPLE, "none").w
                row["lc"]["claret4_uncon"] = rc.score(lambda m: L.intensity(w0, m), T, g)
        out.append(row)
    return out


def jobs():
    import json as _j
    from ldhybrid.grids import mpsatlas
    from ldhybrid.paths import REDUCED
    cat = _j.load(open(REDUCED / "phoenix" / "catalogue.json"))
    ph = sorted(("phoenix", (e["teff"], e["logg"], e["mhs"])) for e in cat)
    mps = [m for m in mpsatlas.catalogue("set2") if 3500 <= m[1] <= 7500 and 3.5 <= m[2] <= 5.0 and -2.0 <= m[0] <= 1.0]
    rng = np.random.default_rng(5)
    pick = rng.choice(len(mps), 300, replace=False)
    return ph + [("mps2", mps[i]) for i in sorted(pick)]


if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor, as_completed
    path = rc.HERE / "d2_d6_full.jsonl"
    done = set()
    if path.exists():
        for l in open(path):
            r = json.loads(l)
            done.add((r["kind"], tuple(r["key"])))
    J = [j for j in jobs() if (j[0], tuple(j[1])) not in done]
    print(len(J), "models to do,", len(done), "done", flush=True)
    with ProcessPoolExecutor(5) as ex, open(path, "a") as fh:
        futs = [ex.submit(task, j) for j in J]
        for n, f in enumerate(as_completed(futs), 1):
            for r in f.result():
                fh.write(json.dumps(r) + "\n")
            fh.flush()
            if n % 50 == 0:
                print(f"  {n}/{len(J)}", flush=True)
    print("done")
