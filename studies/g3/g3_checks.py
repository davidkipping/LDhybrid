"""Gate G3 checks on the P3 products, and the systematic-spread table.

    .venv/bin/python studies/g3/g3_checks.py

1. Physicality: every hybrid row physical, or flagged (counts per grid/law).
2. Interpolation: leave-one-out along Teff (±100 K), log g (±0.5) and [M/H]
   (±0.5) in the PHOENIX headline tables: light-curve error (ppm, production
   grid) of the linearly interpolated weights against the tabulated ones.
3. Coefficient-spectra closure: band weights from the spectra vs the direct
   same-settings band fit, for a sample of models × bands.
4. Systematic spread: at grid points shared by PHOENIX, MPS-ATLAS Set 1 and
   Set 2, the s.d. of each weight across the three and the largest pairwise
   light-curve difference; written to data/products/systematics.csv.gz.
   Gate: the interpolation error must be well below this spread.
"""

import gzip
import json
from pathlib import Path

import h5py
import numpy as np

from ldhybrid import bands, laws, spectra
from ldhybrid import lightcurve as lc
from ldhybrid.grids import phoenix
from ldhybrid.paths import DATA
from ldhybrid.products import PROD_LC

HERE = Path(__file__).parent
P = DATA / "products"
HYB = ("hybrid2", "hybrid4", "hybrid5")
G = lc.LCGrid(**PROD_LC)
COLS = {n: G.linear_columns(laws.get(n)) for n in HYB}


def lcdiff(name, w1, w2):
    one, cols = COLS[name]
    return G.errors(*G.combine(one, cols, w1), *G.combine(one, cols, w2))["max_ppm"]


def load(grid):
    h = h5py.File(P / f"{grid}.h5")
    return h, list(zip(h["teff"][()], h["logg"][()], np.round(h["mh"][()], 2))), \
        [b.decode() for b in h["bands"][()]]


def physicality(report):
    report["physical"] = {}
    for grid in ("phoenix", "mps1", "mps2"):
        if not (P / f"{grid}.h5").exists():
            continue
        h, _, _ = load(grid)
        cov = h["covered"][()]
        for name in HYB + ("power2", "quadratic", "claret4"):
            for m in ("r", "mu"):
                ph = h[f"{name}|{m}|physical"][()][cov]
                bd = h[f"{name}|{m}|bound"][()][cov]
                report["physical"][f"{grid}|{name}|{m}"] = dict(
                    rows=int(cov.sum()), unphysical=int((~ph).sum()), bound=int(bd.sum()))


def interpolation(report, n_max=4000, seed=0):
    h, keys, bnames = load("phoenix")
    idx = {k: i for i, k in enumerate(keys)}
    rng = np.random.default_rng(seed)
    out = {}
    for axis, step in (("teff", 100), ("logg", 0.5), ("mh", 0.5)):
        res = {n: [] for n in HYB}
        trial = []
        for k in keys:
            lo, hi = list(k), list(k)
            a = ("teff", "logg", "mh").index(axis)
            lo[a] = round(k[a] - step, 2)
            hi[a] = round(k[a] + step, 2)
            if tuple(lo) in idx and tuple(hi) in idx:
                trial.append((idx[k], idx[tuple(lo)], idx[tuple(hi)]))
        rng.shuffle(trial)
        for i, a, b in trial[: n_max // len(bnames) + 1]:
            for j in range(len(bnames)):
                if not h["covered"][i, j]:
                    continue
                for n in HYB:
                    W = h[f"{n}|r|w"]
                    w_int = 0.5 * (W[a, j] + W[b, j])
                    res[n].append(lcdiff(n, W[i, j], w_int))
        out[axis] = {n: dict(median=float(np.median(v)), p90=float(np.percentile(v, 90)),
                             n=len(v)) for n, v in res.items()}
    report["interpolation_ppm"] = out


def closure(report, n_models=25, seed=1):
    rng = np.random.default_rng(seed)
    cat = sorted((P / "spectra" / "phoenix").glob("*.npz"))
    worst, nphys, ntot = 0.0, 0, 0
    test_bands = ["TESS", "HST_WFC3_G141", "HST_STIS_G430L", "JWST_NIRSpec_PRISM",
                  "JWST_NIRSpec_G395H", "JWST_NIRISS_SOSS_o1"]
    for f in rng.choice(cat, n_models, replace=False):
        z = np.load(f)
        cs = spectra.CoeffSpectrum(z["lam_edges"], z["I1"].astype(float),
                                   {n: z[f"a_{n}"].astype(float) for n in HYB}, float(z["r_c"]), {})
        spec = phoenix.load(int(z["teff"]), float(z["logg"]),
                            "-0.0" if float(z["mh"]) == 0 else f"{float(z['mh']):+.1f}")
        for b in test_bands:
            lam, S = bands.response(b)
            for n in HYB:
                w, ph = cs.band_coefficients(lam, S, n)
                ref = spectra.band_fit_same_settings(spec, lam, S, n, cs.r_c)
                worst = max(worst, float(np.max(np.abs(w - ref))))
                nphys += ph
                ntot += 1
    report["spectra_closure"] = dict(max_abs_dw=worst, physical=f"{nphys}/{ntot}",
                                     note="float32 storage of a_kj and I1")


def systematics(report):
    grids = {g: load(g) for g in ("phoenix", "mps1", "mps2")}
    hP, kP, bP = grids["phoenix"]
    maps = {g: {k: i for i, k in enumerate(grids[g][1])} for g in grids}
    common = [k for k in kP if all(k in maps[g] for g in ("mps1", "mps2"))]
    bidx = {g: {b: j for j, b in enumerate(grids[g][2])} for g in grids}
    rows, spreads = [], {n: [] for n in HYB}
    for k in common:
        for b in bP:
            if not all(grids[g][0]["covered"][maps[g][k], bidx[g][b]] for g in grids):
                continue
            for n in HYB:
                ws = [grids[g][0][f"{n}|r|w"][maps[g][k], bidx[g][b]] for g in grids]
                sd = np.std(ws, axis=0, ddof=1)
                dl = max(lcdiff(n, ws[i], ws[j]) for i in range(3) for j in range(i + 1, 3))
                spreads[n].append(dl)
                rows.append([f"{k[0]:.0f}", f"{k[1]:.2f}", f"{k[2]:+.2f}", b, n] +
                            [f"{x:.5f}" for x in sd] + [f"{dl:.2f}"])
    with gzip.open(P / "systematics.csv.gz", "wt") as f:
        f.write("teff,logg,mh,band,law,sd_w1,sd_w2,sd_w3,sd_w4,sd_w5,max_pairwise_lc_ppm\n")
        for r in rows:
            sd = r[5:-1] + [""] * (5 - len(r[5:-1]))
            f.write(",".join(r[:5] + sd + [r[-1]]) + "\n")
    report["systematics"] = dict(points=len(common), rows=len(rows),
                                 max_pairwise_lc_ppm={n: dict(median=float(np.median(v)),
                                                              p90=float(np.percentile(v, 90)))
                                                      for n, v in spreads.items()})


if __name__ == "__main__":
    rep = {}
    for step in (physicality, interpolation, closure, systematics):
        try:
            step(rep)
            print(step.__name__, "ok", flush=True)
        except Exception as e:                     # report and continue
            rep[step.__name__ + "_error"] = repr(e)
            print(step.__name__, "FAILED", e, flush=True)
    (HERE / "g3_report.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1))
