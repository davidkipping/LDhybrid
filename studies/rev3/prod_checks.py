"""Checks after the hybrid4 pole change (FORK P):
(a) narrow-bin physicality from the coefficient spectra (D14 method, same
    150 models, seed 3), hybrid4 (hybrid2/5 for reference);
(b) coefficient spectra vs band tables in light-curve terms (25 random PHOENIX
    models, seed 1; TESS, WFC3 G141, NIRSpec PRISM, NIRISS SOSS o1).

    .venv/bin/python studies/rev3/prod_checks.py   -> prod_checks.md/json
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / "rev"))
import d14_narrow_bins as D  # noqa: E402

from ldhybrid import api, bands, laws  # noqa: E402
from ldhybrid import lightcurve as lc  # noqa: E402
from ldhybrid.grids import phoenix  # noqa: E402
from ldhybrid.paths import DATA  # noqa: E402
from ldhybrid.products import PROD_LC  # noqa: E402
from ldhybrid.profile import band_profile  # noqa: E402
from ldhybrid.specint import SpecInt  # noqa: E402

LAWS = ("hybrid2", "hybrid4", "hybrid5")


def narrow():
    files = sorted(D.SPEC.glob("*.npz"))
    files = [f for f in files if ".tmp" not in f.name]
    rng = np.random.default_rng(3)
    pick = rng.choice(len(files), 150, replace=False)
    z0 = np.load(files[0])
    edges0 = z0["lam_edges"]
    dummy = SpecInt(edges0, np.array([0.0, 1.0]), np.zeros((len(edges0) - 1, 2)))
    Wm = {}
    for R in D.RES:
        e = D.bins(R)
        M = np.array([dummy.bin_throughput(np.array([a - 1e-3, a, b, b + 1e-3]), np.array([0.0, 1, 1, 0]))
                      for a, b in zip(e[:-1], e[1:])])
        Wm[R] = (M, e)
    G = {n: laws.get(n).dT(D.X) for n in LAWS}
    cls = lambda t: "M" if t < 3900 else ("K" if t < 5200 else "FG")
    agg = {}
    for i in pick:
        z = np.load(files[i])
        c = cls(float(z["teff"]))
        for R, (M, e) in Wm.items():
            I1 = M @ z["I1"].astype(float)
            for n in LAWS:
                w = (M @ z[f"a_{n}"].astype(float)) / I1[:, None]
                mono = np.all(-(w @ G[n].T) >= -1e-10, axis=1)
                pos = 1.0 - w.sum(1) >= -1e-10
                bad = ~(mono & pos)
                a = agg.setdefault((R, n, c), dict(bins=0, bad=0, models=0, models_bad=0, lam=[]))
                a["bins"] += len(w)
                a["bad"] += int(bad.sum())
                a["models"] += 1
                a["models_bad"] += int(bad.any())
                a["lam"] += ((e[:-1][bad] + e[1:][bad]) / 2).tolist()
    return agg


def spec_vs_table():
    G = lc.LCGrid(**PROD_LC)
    h = h5py.File(DATA / "products" / "phoenix.h5")
    B = [b.decode() for b in h["bands"][()]]
    rng = np.random.default_rng(1)
    idx = rng.choice(len(h["teff"]), 25, replace=False)
    res = {l: [] for l in LAWS}
    for k in idx:
        t, g, mh = int(h["teff"][k]), float(h["logg"][k]), float(h["mh"][k])
        mhs = "-0.0" if mh == 0 else f"{mh:+.1f}"
        cs = api._spectrum(t, g, mhs)
        spec = phoenix.load(t, g, mhs)
        for b in ("TESS", "HST_WFC3_G141", "JWST_NIRSpec_PRISM", "JWST_NIRISS_SOSS_o1"):
            j = B.index(b)
            pr = band_profile(spec, *bands.response(b), limb="M1")
            tr = G.columns(pr)
            for l in LAWS:
                L = laws.get(l)
                w_tab = h[f"{l}|r|w"][k, j]
                w_sp, _ = cs.band_coefficients(*bands.response(b), l)
                res[l].append([G.errors(*G.columns(lambda m: L.intensity(w_sp, m)),
                                        *G.columns(lambda m: L.intensity(w_tab, m)))["max_ppm"],
                               G.errors(*G.columns(lambda m: L.intensity(w_sp, m)), *tr)["max_ppm"],
                               float(h[f"{l}|r|lc_ppm"][k, j])])
    return res


def main():
    api._spectrum.cache_clear()
    agg = narrow()
    svt = spec_vs_table()
    L = ["# FORK P checks after the hybrid4 pole change", "",
         "## (a) narrow bins summed from the coefficient spectra: % unphysical bins (% models with ≥ 1)", "",
         "| R | law | M | K | FG | 90% of bad bins within (µm) |", "|---|---|---|---|---|---|"]
    for R in D.RES:
        for n in LAWS:
            cells, lam = [], []
            for c in ("M", "K", "FG"):
                a = agg[(R, n, c)]
                cells.append(f"{100 * a['bad'] / a['bins']:.2f} ({100 * a['models_bad'] / a['models']:.0f})")
                lam += a["lam"]
            rng = (f"{np.percentile(lam, 5) / 1e4:.2f}–{np.percentile(lam, 95) / 1e4:.2f}" if lam else "—")
            L.append(f"| {R} | {n} | " + " | ".join(cells) + f" | {rng} |")
    L += ["", "## (b) coefficient spectra vs band tables (25 PHOENIX models × 4 bands), ppm median / 90th / max", "",
          "| law | spectra vs table | spectra vs model | table vs model |", "|---|---|---|---|"]
    for l, v in svt.items():
        v = np.array(v)
        L.append(f"| {l} | " + " | ".join(
            f"{np.median(v[:, i]):.2f} / {np.percentile(v[:, i], 90):.1f} / {v[:, i].max():.1f}" for i in range(3)) + " |")
    (HERE / "prod_checks.md").write_text("\n".join(L) + "\n")
    (HERE / "prod_checks.json").write_text(json.dumps(dict(
        narrow={f"{R}|{n}|{c}": {k: v for k, v in a.items() if k != "lam"} for (R, n, c), a in agg.items()},
        spec_vs_table=svt)))
    print("\n".join(L))


if __name__ == "__main__":
    main()
