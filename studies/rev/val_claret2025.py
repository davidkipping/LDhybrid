"""Revision D8: comparison with Claret, Hauschildt & Torres (2025, J/A+A/699/A97),
four-term and power-2 coefficients for 11 JWST passbands from new
high-resolution spherical PHOENIX models (solar composition, ξ = 1 km/s).

Their intensities are not public, so this is a consistency check (as G1(c) for
Claret & Southworth 2023), not a reproduction. For every PHOENIX-ACES v3.0
model at [M/H] = 0 that is also in their table (241 models) and the nine
overlapping passbands we compute:

  recipe      their recipe as described in the ReadMe: limb cut at our M1 edge,
              60 000 points uniform in μ (PCHIP), unweighted, unconstrained
              least squares (four-term) / Levenberg–Marquardt (power-2);
  recipe_mc   the same with the limb placed at their tabulated μ_crit;
  tab_mu/r    our tabulated uniform-μ and r-weighted coefficients.

Light-curve differences (ppm): worst |ΔF| on the production grid (p = 0.1,
b = 0, 0.3, 0.6, 0.9) between pairs of laws, and against our model profile.

    .venv/bin/python studies/rev/val_claret2025.py
"""

import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import h5py
import numpy as np

from ldhybrid import bands, fit, laws
from ldhybrid import lightcurve as lc
from ldhybrid.grids import phoenix
from ldhybrid.paths import DATA, RAW
from ldhybrid.products import PROD_LC
from ldhybrid.profile import band_profile

HERE = Path(__file__).parent
V = RAW / "vizier" / "699_A97"
BANDS = {"F322W2": "JWST_NIRCam_F322W2", "F444W": "JWST_NIRCam_F444W",
         "SOSS1": "JWST_NIRISS_SOSS_o1", "SOSS2": "JWST_NIRISS_SOSS_o2",
         "G235H": "JWST_NIRSpec_G235H", "G235M": "JWST_NIRSpec_G235M",
         "G395H": "JWST_NIRSpec_G395H", "G395M": "JWST_NIRSpec_G395M",
         "PRISM": "JWST_NIRSpec_PRISM"}
C4 = laws.get("claret4")
N_RES = 60000


def read_claret():
    out = {}
    for line in open(V / "table2.dat"):
        g, t, f = float(line[0:4]), int(line[5:9]), line[10:16].strip()
        out.setdefault((t, g, f), {})["a"] = [float(line[17 + 11 * k:27 + 11 * k]) for k in range(4)]
        out[(t, g, f)]["mucrit4"] = float(line[61:71])
    for line in open(V / "table3.dat"):
        g, t, f = float(line[0:4]), int(line[5:9]), line[10:16].strip()
        out[(t, g, f)]["p2"] = [float(line[17:28]), float(line[29:40])]
        out[(t, g, f)]["mucrit2"] = float(line[41:51])
    return out


def tabulated():
    h = h5py.File(DATA / "products" / "phoenix.h5", "r")
    keys = np.c_[h["teff"][:], h["logg"][:], h["mh"][:]]
    bl = [b.decode() for b in h["bands"][:]]
    idx = {(int(t), float(g)): i for i, (t, g, m) in enumerate(keys) if m == 0}
    T = {k: h[k][:] for k in ("claret4|r|w", "claret4|mu|w", "power2|r|w", "power2|mu|w", "r_c")}
    return idx, bl, T


def f4(a):
    return lambda mu: C4.intensity(a, mu)


def fp2(cg):
    return lambda mu: laws.power2(cg[0], cg[1], mu)


def lcd(f, g):
    return lc.lc_errors(f, g, **PROD_LC)["max_ppm"]


def job(args):
    (t, g), cl, tab = args
    s = phoenix.load(t, g, "-0.0")
    rows = []
    for f, b in BANDS.items():
        c = cl.get((t, g, f))
        if c is None:
            continue
        lam, S = bands.response(b)
        pr = band_profile(s, lam, S, limb="M1")
        mucri = float(np.sqrt(1 - pr.meta["r_c"] ** 2))
        # their μ_crit imposed as the edge (same rescaling)
        pr_mc = band_profile(s, lam, S, limb="M1", r_c=float(np.sqrt(1 - c["mucrit4"] ** 2)))
        A, P = f4(c["a"]), fp2(c["p2"])
        tb = tab[b]
        r = dict(teff=t, logg=g, band=f, mucri_ours=mucri, mucri_claret=c["mucrit4"],
                 a_cl=c["a"], p_cl=c["p2"], a_mu=tb["c4mu"], a_r=tb["c4r"], p_mu=tb["p2mu"], p_r=tb["p2r"],
                 lc4={}, lcp2={}, var={})
        for lim, prof in (("M1", pr), ("mc", pr_mc)):
            for meas in ("mu", "r", "area"):
                a = fit.fit_linear(C4, prof, meas, N_RES, "none").w
                p = fit.fit_power2(prof, meas, N_RES).w
                k = f"{lim}_{meas}"
                r["var"][k] = dict(a=a.tolist(), p=p.tolist())
                r["lc4"]["cl_" + k] = lcd(A, f4(a))
                r["lcp2"]["cl_" + k] = lcd(P, fp2(p))
                if lim == "M1":
                    r["lc4"][k + "_truth"] = lcd(f4(a), pr)
                    r["lcp2"][k + "_truth"] = lcd(fp2(p), pr)
        r["a_rec"], r["p_rec"] = r["var"]["M1_mu"]["a"], r["var"]["M1_mu"]["p"]
        for k, f_ in (("cl_mu", f4(tb["c4mu"])), ("cl_r", f4(tb["c4r"])), ("cl_truth", pr)):
            r["lc4"][k] = lcd(A, f_)
        for k, f_ in (("cl_mu", fp2(tb["p2mu"])), ("cl_r", fp2(tb["p2r"])), ("cl_truth", pr)):
            r["lcp2"][k] = lcd(P, f_)
        r["lc4"]["tab_r_truth"] = lcd(f4(tb["c4r"]), pr)
        r["lcp2"]["tab_r_truth"] = lcd(fp2(tb["p2r"]), pr)
        rows.append(r)
    return rows


def main():
    cl = read_claret()
    idx, bl, T = tabulated()
    jobs = []
    for (t, g), i in sorted(idx.items()):
        if not any((t, g, f) in cl for f in BANDS):
            continue
        tab = {}
        for f, b in BANDS.items():
            j = bl.index(b)
            tab[b] = dict(c4r=T["claret4|r|w"][i, j].tolist(), c4mu=T["claret4|mu|w"][i, j].tolist(),
                          p2r=T["power2|r|w"][i, j].tolist(), p2mu=T["power2|mu|w"][i, j].tolist())
        jobs.append(((t, g), {k: v for k, v in cl.items() if k[:2] == (t, g)}, tab))
    print(len(jobs), "models", flush=True)
    rows = []
    with ProcessPoolExecutor(4) as ex:
        for n, r in enumerate(ex.map(job, jobs, chunksize=4), 1):
            rows += r
            if n % 20 == 0:
                print(n, flush=True)
    (HERE / "val_claret2025.json").write_text(json.dumps(rows))
    summarise(rows)


def summarise(rows):
    q = lambda x: f"{np.median(x):.3g} / {np.percentile(x, 90):.3g} / {np.max(x):.3g}"
    L = ["# D8: Claret, Hauschildt & Torres (2025) JWST tables vs our PHOENIX-ACES v3.0", "",
         f"{len(rows)} model-passbands ({len({(r['teff'], r['logg']) for r in rows})} models at [M/H] = 0, "
         f"9 passbands). Entries: median / 90th pct / max.", ""]
    dm = np.array([r["mucri_ours"] - r["mucri_claret"] for r in rows])
    L += [f"* μ_cri (ours M1 − theirs): median {np.median(dm):+.4f}, median |Δ| {np.median(abs(dm)):.4f}, "
          f"90th |Δ| {np.percentile(abs(dm), 90):.4f}", ""]
    da = np.array([np.subtract(r["a_rec"], r["a_cl"]) for r in rows])
    dp = np.array([np.subtract(r["p_rec"], r["p_cl"]) for r in rows])
    L += ["## Coefficients, our recipe fit − theirs (|Δ|)", "",
          "| coef | median | 90th | max |", "|---|---|---|---|"]
    for k, n in enumerate(["a1", "a2", "a3", "a4"]):
        x = abs(da[:, k])
        L.append(f"| {n} | {np.median(x):.3f} | {np.percentile(x, 90):.3f} | {x.max():.3f} |")
    for k, n in enumerate(["c (g)", "α (h)"]):
        x = abs(dp[:, k])
        L.append(f"| {n} | {np.median(x):.3f} | {np.percentile(x, 90):.3f} | {x.max():.3f} |")
    L += ["", "## Light-curve differences (ppm, p = 0.1)", "",
          "| comparison | four-term | power-2 |", "|---|---|---|"]
    lab = {"cl_M1_mu": "theirs vs our fit, uniform μ (stated recipe), our M1 edge",
           "cl_mc_mu": "theirs vs our fit, uniform μ, their μ_crit",
           "cl_M1_r": "theirs vs our fit, uniform r, M1", "cl_mc_r": "theirs vs our fit, uniform r, their μ_crit",
           "cl_M1_area": "theirs vs our fit, uniform area, M1",
           "cl_mu": "theirs vs our tabulated uniform-μ (constrained)",
           "cl_r": "theirs vs our tabulated r-weighted", "cl_truth": "theirs vs our model profile",
           "M1_mu_truth": "our uniform-μ fit vs our model profile",
           "M1_r_truth": "our uniform-r fit vs our model profile",
           "tab_r_truth": "our tabulated r-weighted vs our model profile"}
    for k, t in lab.items():
        L.append(f"| {t} | {q([r['lc4'][k] for r in rows])} | {q([r['lcp2'][k] for r in rows])} |")
    L += ["", "By temperature (their coefficients vs our model profile, median ppm):", ""]
    for lo, hi in ((2400, 3500), (3500, 5000), (5000, 7400)):
        x = [r["lc4"]["cl_truth"] for r in rows if lo <= r["teff"] < hi]
        y = [r["lcp2"]["cl_truth"] for r in rows if lo <= r["teff"] < hi]
        L.append(f"* {lo}–{hi} K: four-term {np.median(x):.0f}, power-2 {np.median(y):.0f} (n = {len(x)})")
    (HERE / "val_claret2025.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
