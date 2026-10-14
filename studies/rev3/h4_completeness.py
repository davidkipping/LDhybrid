"""Is the hybrid4 simplex complete? The simplex (w >= 0, sum w <= 1) is a
sufficient physical region; the exact region of Kipping's conditions
(I >= 0, dI/dmu >= 0) is larger. This script asks (A) how much larger, and
(B) whether any model-atmosphere profile prefers a physical hybrid4 fit
outside the simplex, and what that costs in light-curve accuracy.

(A) Geometry: the exact region P on a dense xi = mu^2 grid (fit._inequalities
    "physical"); its extent in each w_j (LPs), its volume relative to the
    simplex (Monte Carlo), and the most extreme profiles outside the simplex.
(B) Data: only fits where the simplex binds can move (otherwise the simplex
    optimum is the unconstrained optimum, which is already physical). For every
    bound hybrid4 row of the CDS tables, refit the band profile with the
    simplex (reproduces the table), the exact physical region and no
    constraint, and compare weights and light-curve errors.

    .venv/bin/python -W ignore studies/rev3/h4_completeness.py   -> h4_completeness.md
"""

import os
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from ldhybrid import fit, laws, v1  # noqa: E402
from ldhybrid.interpolate import LDTables  # noqa: E402

HERE = Path(__file__).parent
OUT = HERE / "h4_completeness.md"
CDS = ROOT / "data" / "release" / "cds"
PROF = ROOT / "data" / "release" / "zenodo"
N_LC_MU = 400                     # light curves for a random subset of the mu-method rows
RNG = np.random.default_rng(20261014)


# ---------------------------------------------------------------- (A) geometry
def region(name):
    L = laws.get(name)
    C, d = fit._inequalities(L, "physical")
    return L, C, d


def extent(C, d, m):
    lo, hi = [], []
    for j in range(m):
        c = np.zeros(m)
        c[j] = 1.0
        lo.append(linprog(c, A_ub=C, b_ub=d, bounds=[(None, None)] * m, method="highs").x)
        hi.append(linprog(-c, A_ub=C, b_ub=d, bounds=[(None, None)] * m, method="highs").x)
    return np.array(lo), np.array(hi)


def volume(C, d, lo, hi, n=2_000_000, chunk=50_000):
    m = len(lo)
    inside = in_simplex = 0
    for _ in range(n // chunk):
        W = lo + (hi - lo) * RNG.random((chunk, m))
        ok = np.all(W @ C.T <= d + 1e-12, axis=1)
        inside += ok.sum()
        in_simplex += np.sum(ok & np.all(W >= 0, axis=1))
    box = np.prod(hi - lo)
    return box * inside / n, box * in_simplex / n, inside


def geometry(name):
    L, C, d = region(name)
    m = L.m
    lo, hi = extent(C, d, m)
    wmin = np.array([lo[j][j] for j in range(m)])
    wmax = np.array([hi[j][j] for j in range(m)])
    vP, vS_mc, nin = volume(C, d, wmin, wmax)
    from math import factorial
    vS = 1.0 / factorial(m)
    return dict(L=L, wmin=wmin, wmax=wmax, argmin=lo, vP=vP, vS=vS, vS_mc=vS_mc, nin=nin)


# ---------------------------------------------------------------- (B) data
_G = {}


def _lc():
    if "g" not in _G:
        from ldhybrid import lightcurve as lc
        g = lc.LCGrid()
        _G["g"] = g
        _G["cols"] = g.linear_columns(laws.get("hybrid4"))
    return _G["g"], _G["cols"]


def task(job):
    from ldhybrid.profile import Profile
    mu, y, measure, w_tab, do_lc = job
    L = laws.get("hybrid4")
    pr = Profile(mu, y)
    out = {}
    for con in ("simplex", "physical", "none"):
        f = fit.fit_linear(L, pr, measure, v1.RESAMPLE, con, v1.FCM)
        out[con] = dict(w=np.asarray(f.w), rms=f.rms, phys=bool(f.extra["physical"]))
    out["dtab"] = float(np.max(np.abs(out["simplex"]["w"] - w_tab)))
    if do_lc:
        g, (one, cols) = _lc()
        truth = g.columns(pr)
        for con in ("simplex", "physical"):
            out[con]["lc"] = g.errors(*g.combine(one, cols, out[con]["w"]), *truth)["max_ppm"]
    return out


def bound_jobs():
    import h5py
    t = LDTables(CDS, cache=False)
    d = t.table("phoenix", "hybrid4")
    with h5py.File(PROF / "ldhybrid_v1.0_band_profiles_phoenix.h5", "r") as h:
        P = {k: h[k][()] for k in ("teff", "logg", "mh", "mu", "y")}
        bands = [b.decode() for b in h["bands"][()]]
    key = {(round(a, 3), round(b, 3), round(c, 3)): i
           for i, (a, b, c) in enumerate(zip(P["teff"], P["logg"], P["mh"]))}
    jobs, meta = [], []
    for method in ("r", "mu"):
        rows = np.flatnonzero((d["method"] == method) & (d["Bound"] == 1))
        lc_rows = set(rows) if method == "r" else set(RNG.choice(rows, min(N_LC_MU, len(rows)), replace=False))
        for r in rows:
            i = key[(round(d["teff"][r], 3), round(d["logg"][r], 3), round(d["mh"][r], 3))]
            j = bands.index(d["band"][r])
            ok = np.isfinite(P["mu"][i, j])
            jobs.append((P["mu"][i, j][ok].astype(float), P["y"][i, j][ok].astype(float),
                         method, d["coeffs"][r], r in lc_rows))
            meta.append((method, d["teff"][r], d["logg"][r], d["mh"][r], d["band"][r]))
    return jobs, meta


def check_monte_carlo_count(n):
    return "" if n > 20_000 else f" (only {n} Monte Carlo points inside; volume uncertain)"


def main():
    lines = ["# Completeness of the hybrid4 simplex", "",
             "Generated by `studies/rev3/h4_completeness.py`.", ""]

    # (A)
    lines += ["## A. The exact physical region versus the simplex", ""]
    for name in ("hybrid4", "hybrid5"):
        G = geometry(name)
        lines += [f"**{name}**", "",
                  "| | " + " | ".join(f"w{j + 1}" for j in range(G['L'].m)) + " |",
                  "|---|" + "---|" * G["L"].m,
                  "| min over exact region | " + " | ".join(f"{v:+.4f}" for v in G["wmin"]) + " |",
                  "| max over exact region | " + " | ".join(f"{v:+.4f}" for v in G["wmax"]) + " |",
                  "| simplex range | " + " | ".join("0 to 1" for _ in G["wmin"]) + " |", "",
                  f"Volume: exact region {G['vP']:.4g}, simplex {G['vS']:.4g} "
                  f"(Monte Carlo {G['vS_mc']:.4g}); the simplex is "
                  f"{100 * G['vS'] / G['vP']:.1f}% of the exact region"
                  f"{check_monte_carlo_count(G['nin'])}.", ""]
        for j in range(G["L"].m):
            w = G["argmin"][j]
            lines.append(f"- most negative w{j + 1}: w = ({', '.join(f'{v:+.4f}' for v in w)}), "
                         f"I(0) = {1 - w.sum():.4f}")
        lines.append("")

    # (B)
    jobs, meta = bound_jobs()
    nw = min(4, os.cpu_count() or 1)
    with ProcessPoolExecutor(nw) as ex:
        res = list(ex.map(task, jobs, chunksize=20))
    meta = np.array(meta, dtype=object)
    lines += ["## B. Tabulated fits where the simplex binds (PHOENIX only; MPS-ATLAS never binds)", "",
              "| method | bound rows | refit = table (max abs dw) | exact-region fit outside simplex | "
              "unconstrained fit physical | max outside (most negative w) | rms I change (median / max) |",
              "|---|---|---|---|---|---|---|"]
    for method in ("r", "mu"):
        k = [i for i, mm in enumerate(meta[:, 0]) if mm == method]
        if not k:
            continue
        R = [res[i] for i in k]
        dtab = max(r["dtab"] for r in R)
        Wp = np.array([r["physical"]["w"] for r in R])
        outside = np.min(Wp, axis=1) < -1e-8
        free_phys = np.mean([r["none"]["phys"] for r in R])
        drms = np.array([r["simplex"]["rms"] - r["physical"]["rms"] for r in R])
        lines.append(f"| {method} | {len(R)} | {dtab:.1e} | {outside.sum()} ({100 * outside.mean():.1f}%) | "
                     f"{100 * free_phys:.1f}% | {Wp.min():+.4f} | "
                     f"{np.median(drms):.2e} / {drms.max():.2e} |")
    lines.append("")
    lines += ["Light-curve error (max ppm over the production (p, b) grid), simplex vs exact region:", "",
              "| method | rows with LC | simplex median / max | exact region median / max | "
              "improvement median / max |", "|---|---|---|---|---|"]
    for method in ("r", "mu"):
        R = [res[i] for i, mm in enumerate(meta[:, 0]) if mm == method and "lc" in res[i]["simplex"]]
        if not R:
            continue
        s = np.array([r["simplex"]["lc"] for r in R])
        p = np.array([r["physical"]["lc"] for r in R])
        lines.append(f"| {method} | {len(R)} | {np.median(s):.2f} / {s.max():.2f} | "
                     f"{np.median(p):.2f} / {p.max():.2f} | {np.median(s - p):.2f} / {(s - p).max():.2f} |")
    lines.append("")
    # where the outside fits live
    lines += ["Where the exact-region fits leave the simplex (rows outside, by Teff and band):", ""]
    for method in ("r", "mu"):
        k = [i for i, mm in enumerate(meta[:, 0]) if mm == method and np.min(res[i]["physical"]["w"]) < -1e-8]
        if not k:
            lines.append(f"- {method}: none")
            continue
        T = np.array([meta[i, 1] for i in k], float)
        from collections import Counter
        bands = Counter(meta[i, 4] for i in k).most_common(6)
        neg = Counter(int(np.argmin(res[i]["physical"]["w"])) + 1 for i in k)
        lines.append(f"- {method}: {len(k)} rows, Teff {T.min():.0f}-{T.max():.0f} K "
                     f"(median {np.median(T):.0f}); negative weight: "
                     + ", ".join(f"w{j} x{c}" for j, c in sorted(neg.items()))
                     + "; bands: " + ", ".join(f"{b} ({c})" for b, c in bands))
    OUT.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
