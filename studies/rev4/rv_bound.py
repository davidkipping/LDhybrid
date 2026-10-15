"""Tabulated hybrid4 fits (curvature-simplex basis) where the constraint binds:
refit each with the exact region of Kipping's conditions I-II alone (no
curvature conditions) and compare. Fits where the constraint does not bind are
the unconstrained optimum, so they are unchanged by any larger region.

    .venv/bin/python studies/rev4/rv_bound.py   -> rv_bound.md
"""
import os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import h5py
import numpy as np

from ldhybrid import fit, laws, v1
from ldhybrid.interpolate import LDTables
from ldhybrid.paths import DATA

HERE = Path(__file__).parent
CDS = DATA / "release" / "cds"
PROF = DATA / "release" / "zenodo" / "ldhybrid_v1.0_band_profiles_{}.h5"
N_LC_MU = 300
RNG = np.random.default_rng(20261014)
_G = {}


def task(job):
    from ldhybrid.profile import Profile
    mu, y, measure, w_tab, do_lc = job
    L = laws.get("hybrid4")
    pr = Profile(mu, y)
    out = {c: np.asarray(fit.fit_linear(L, pr, measure, v1.RESAMPLE, c, v1.FCM).w)
           for c in ("simplex", "physical")}
    out["dtab"] = float(np.max(np.abs(out["simplex"] - w_tab)))
    if do_lc:
        if "g" not in _G:
            from ldhybrid import lightcurve as lc
            _G["g"] = lc.LCGrid()
            _G["cols"] = _G["g"].linear_columns(L)
        g, (one, cols) = _G["g"], _G["cols"]
        truth = g.columns(pr)
        for c in ("simplex", "physical"):
            out[c + "_lc"] = g.errors(*g.combine(one, cols, out[c]), *truth)["max_ppm"]
    return out


def jobs(grid):
    d = LDTables(CDS, cache=False).table(grid, "hybrid4")
    with h5py.File(str(PROF).format(grid), "r") as h:
        P = {k: h[k][()] for k in ("teff", "logg", "mh", "mu", "y")}
        bands = [b.decode() for b in h["bands"][()]]
    key = {(round(a, 3), round(b, 3), round(c, 3)): i
           for i, (a, b, c) in enumerate(zip(P["teff"], P["logg"], P["mh"]))}
    J, M = [], []
    for method in ("r", "mu"):
        rows = np.flatnonzero((d["method"] == method) & (d["Bound"] == 1))
        lc_rows = set(rows) if method == "r" else set(RNG.choice(rows, min(N_LC_MU, len(rows)), replace=False)) if len(rows) else set()
        for r in rows:
            i = key[(round(d["teff"][r], 3), round(d["logg"][r], 3), round(d["mh"][r], 3))]
            j = bands.index(d["band"][r])
            ok = np.isfinite(P["mu"][i, j])
            J.append((P["mu"][i, j][ok].astype(float), P["y"][i, j][ok].astype(float), method,
                      d["coeffs"][r], r in lc_rows))
            M.append((grid, method, d["teff"][r], d["band"][r]))
    return J, M


def main():
    L = ["# hybrid4 (curvature-simplex basis): fits where the constraint binds", "",
         "| grid | method | rows | bound | refit = table (max dw) | Kipping-region fit leaves the region | "
         "LC simplex median / max | LC Kipping region median / max |", "|---|---|---|---|---|---|---|---|"]
    t = LDTables(CDS, cache=False)
    notes = []
    for grid in ("phoenix", "mps1", "mps2"):
        d = t.table(grid, "hybrid4")
        J, M = jobs(grid)
        with ProcessPoolExecutor(min(6, os.cpu_count() or 1)) as ex:
            R = list(ex.map(task, J, chunksize=20)) if J else []
        for method in ("r", "mu"):
            n = int(np.sum(d["method"] == method))
            k = [i for i, m in enumerate(M) if m[1] == method]
            if not k:
                L.append(f"| {grid} | {method} | {n} | 0 | | | | |")
                continue
            Rk = [R[i] for i in k]
            leaves = [i for i in k if np.min(R[i]["physical"]) < -1e-8 or np.sum(R[i]["physical"]) > 1 + 1e-8]
            lc = [r for r in Rk if "simplex_lc" in r]
            s = np.array([r["simplex_lc"] for r in lc]) if lc else np.array([np.nan])
            p = np.array([r["physical_lc"] for r in lc]) if lc else np.array([np.nan])
            L.append(f"| {grid} | {method} | {n} | {len(k)} ({100 * len(k) / n:.2f}%) | "
                     f"{max(r['dtab'] for r in Rk):.1e} | {len(leaves)} | "
                     f"{np.nanmedian(s):.1f} / {np.nanmax(s):.1f} | {np.nanmedian(p):.1f} / {np.nanmax(p):.1f} |")
            if leaves:
                T = np.array([M[i][2] for i in leaves], float)
                bands = Counter(M[i][3] for i in leaves).most_common(6)
                notes.append(f"- {grid} {method}: {len(leaves)} rows, Teff {T.min():.0f}-{T.max():.0f} K; bands: "
                             + ", ".join(f"{b} ({c})" for b, c in bands))
    L += ["", "LC: max ppm over the production (p, b) grid (all bound r rows; a random "
          f"{N_LC_MU} of the bound mu rows).", ""] + notes
    (HERE / "rv_bound.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
