"""Revision statistics from the rebuilt production files (data/products/*.h5).

    .venv/bin/python studies/rev/production_stats.py [--only corrections]

a. light-curve error per law × grid (Table 2), plus MIRI LRS and NRS1/NRS2 subsets;
b. constraint-active ('bound') fractions; power-2 unphysical fractions;
c. chromatic stellar edge r_c across JWST bands (PHOENIX);
d. Δw corrections added to MPS-ATLAS Set 2: fraction leaving the physical region;
e. passbands per grid and CDS record counts.

Writes studies/rev/production_stats.{json,md}.
"""

import csv
import gzip
import json
import sys
from pathlib import Path

import h5py
import numpy as np

from ldhybrid import laws, predict
from ldhybrid.paths import DATA
from ldhybrid.products import LAWS, METHODS

HERE = Path(__file__).parent
P = DATA / "products"
GRIDS = ("phoenix", "mps1", "mps2")
HYB = ("hybrid2", "hybrid4", "hybrid5")
NRS = ("JWST_NIRSpec_G395H_NRS1", "JWST_NIRSpec_G395H_NRS2")
OUT = HERE / "production_stats.json"


def h5(grid):
    return h5py.File(P / f"{grid}.h5")


def bands_of(h):
    return [b.decode() for b in h["bands"][()]]


def stat(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return dict(n=int(x.size), median=float(np.median(x)), p99=float(np.percentile(x, 99)),
                max=float(np.max(x)))


def lc_stats():
    out = {}
    for g in GRIDS:
        with h5(g) as h:
            bl, cov = bands_of(h), h["covered"][()]
            sub = {"all": np.ones(len(bl), bool),
                   "NRS1+NRS2": np.isin(bl, NRS),
                   "MIRI LRS": np.array([b == "JWST_MIRI_LRS" for b in bl])}
            for name, sel in sub.items():
                if not sel.any() or not cov[:, sel].any():
                    continue
                for law in LAWS:
                    for m in METHODS:
                        x = h[f"{law}|{m}|lc_ppm"][()][:, sel][cov[:, sel]]
                        out[f"{g}|{name}|{law}|{m}"] = stat(x)
    return out


def bound_stats():
    out = {}
    for g in GRIDS:
        with h5(g) as h:
            cov = h["covered"][()]
            for law in LAWS:
                for m in METHODS:
                    b = h[f"{law}|{m}|bound"][()][cov]
                    ph = h[f"{law}|{m}|physical"][()][cov]
                    out[f"{g}|{law}|{m}"] = dict(n=int(b.size), bound=float(b.mean()),
                                                 unphysical=float(1 - ph.mean()))
    return out


JWST = ("JWST_NIRISS", "JWST_NIRSpec", "JWST_NIRCam")


def edge_stats():
    with h5("phoenix") as h:
        bl, cov, rc = bands_of(h), h["covered"][()], h["r_c"][()]
    rc = np.where(cov, rc, np.nan)
    out = {}
    for name, sel in {"JWST": np.array([b.startswith(JWST) for b in bl]),
                      "all": np.ones(len(bl), bool),
                      "NRS1-NRS2": np.isin(bl, NRS)}.items():
        r = rc[:, sel]
        d = (np.nanmax(r, 1) - np.nanmin(r, 1)) / np.nanmean(r, 1)
        out[name] = dict(n_bands=int(sel.sum()), median=float(np.nanmedian(d)),
                         p90=float(np.nanpercentile(d, 90)), max=float(np.nanmax(d)),
                         depth_ppm_1pct=dict(median=float(2 * np.nanmedian(d) * 1e4),
                                             p90=float(2 * np.nanpercentile(d, 90) * 1e4),
                                             max=float(2 * np.nanmax(d) * 1e4)))
    return out


def _physical(law, W, tol=1e-10):
    """Exact positive–monotonic check for many weight vectors at once."""
    L = laws.get(law)
    x = np.linspace(0.0, 1.0, 2001)
    G = -L.dT(x)                                         # dI/dx = G w
    ok = 1.0 - W.sum(1) >= -tol
    for i in range(0, len(W), 20000):
        ok[i:i + 20000] &= np.min(G @ W[i:i + 20000].T, axis=0) >= -tol
    return ok


def _simplex(W, tol=1e-10):
    return np.all(W >= -tol, 1) & (W.sum(1) <= 1 + tol)


def correction_stats():
    rows = list(csv.DictReader(open(P / "corrections.csv")))
    out = {}
    with h5("mps2") as h:
        bl, cov = bands_of(h), h["covered"][()]
        T, G_, M = h["teff"][()], h["logg"][()], h["mh"][()]
        W = {law: h[f"{law}|r|w"][()] for law in HYB}
    muram_T = sorted({float(r["teff"]) for r in rows if r["source"] == "MURaM"
                      and r["kind"] == "magnetic" and float(r["mh"]) == 0.0})
    inT = (T >= muram_T[0]) & (T <= muram_T[-1])
    cbands = sorted({r["band"] for r in rows})
    # magnetic: predict.magnetic (Δw linear in Teff between the [M/H] = 0 stars)
    for state in ("ssd", "B100G", "B200G", "B300G"):
        for law in HYB:
            n = bad = bad_s = 0
            for b in cbands:
                if b not in bl:
                    continue
                j = bl.index(b)
                sel = cov[:, j] & inT
                Ts = np.unique(T[sel])
                dw = {t: predict.magnetic(law, b, t, state) for t in Ts}
                if any(v is None for v in dw.values()):
                    continue
                Wc = W[law][sel, j] + np.array([dw[t] for t in T[sel]])
                ok = _physical(law, Wc)
                n += len(Wc)
                bad += int((~ok).sum())
                if law != "hybrid2":
                    bad_s += int((~_simplex(Wc)).sum())
            out[f"MURaM magnetic {state}|{law}"] = dict(
                n=n, unphysical=bad / n, outside_simplex=(bad_s / n if law != "hybrid2" else None),
                applied="MPS2 grid points with Teff in MURaM range, all log g and [M/H]")
    # 3D: each row's own 1D base (MPS2 nearest for MURaM, Kurucz for Stagger), and
    # for Stagger also the MPS2 point at the same (Teff, log g, [M/H])
    key = {(round(t), round(g, 2), round(m, 2)): i for i, (t, g, m) in enumerate(zip(T, G_, M))}
    for src in ("MURaM", "Stagger"):
        for law in HYB:
            m_ = laws.get(law).m
            sub = [r for r in rows if r["source"] == src and r["kind"] == "3d" and r["law"] == law]
            w0 = np.array([[float(r[f"w{k + 1}"]) for k in range(m_)] for r in sub])
            dw = np.array([[float(r[f"dw{k + 1}"]) for k in range(m_)] for r in sub])
            ok = _physical(law, w0 + dw)
            res = dict(n=len(sub), unphysical_own_base=float(1 - ok.mean()))
            if law != "hybrid2":
                res["outside_simplex_own_base"] = float(1 - _simplex(w0 + dw).mean())
            if src == "Stagger":
                Wm = []
                for r, d in zip(sub, dw):
                    i = key.get((round(float(r["teff"])), round(float(r["logg"]), 2),
                                 round(float(r["mh"]), 2)))
                    if i is not None and r["band"] in bl and cov[i, bl.index(r["band"])]:
                        Wm.append(W[law][i, bl.index(r["band"])] + d)
                Wm = np.array(Wm)
                res["n_on_mps2"] = len(Wm)
                res["unphysical_on_mps2"] = float(1 - _physical(law, Wm).mean())
                if law != "hybrid2":
                    res["outside_simplex_on_mps2"] = float(1 - _simplex(Wm).mean())
            out[f"{src} 3D|{law}"] = res
    out["_conventions"] = (
        "Δw = fit(I_base · C) − fit(I_base), added to the 1D weights (I_corrected ≈ "
        "I_1D − Σ Δw_j T_j). Magnetic: C = I_B/I_hydro (MURaM, same code with and "
        "without field), I_base = MPS-ATLAS Set 2 band profile at the nearest grid point. "
        "MURaM 3D: C = I_hydro/I_MPS2 (MPS-ATLAS Set 2 = the code that synthesised the "
        "MURaM spectra; nearest grid point within 150 K). Stagger 3D: "
        "Δw = fit(Stagger) − fit(ATLAS9/Kurucz) at identical (Teff, log g, [M/H]) "
        "from the ExoTiC-LD bundle. MURaM ratios held constant below μ = 0.1.")
    return out


def cds_stats():
    out = {}
    for g in GRIDS:
        with h5(g) as h:
            out[f"bands|{g}"] = bands_of(h)
    cds = DATA / "release" / "cds"
    for f in sorted(cds.glob("*.dat.gz")):
        with gzip.open(f, "rt") as fh:
            out[f"records|{f.name[:-3]}"] = sum(1 for _ in fh)
    return out


def md(R):
    L = ["# Revision statistics from the rebuilt production files", "",
         "Generated by `production_stats.py`; r-weighted unless stated.", ""]
    if "lc" in R:
        L += ["## a. Light-curve error (ppm), median / 99th / max", "",
              "| grid | subset | " + " | ".join(LAWS) + " |", "|---|---|" + "---|" * len(LAWS)]
        for g in GRIDS:
            for sub in ("all", "NRS1+NRS2", "MIRI LRS"):
                k0 = f"{g}|{sub}|{LAWS[0]}|r"
                if k0 not in R["lc"]:
                    continue
                cells = []
                for law in LAWS:
                    s = R["lc"][f"{g}|{sub}|{law}|r"]
                    cells.append(f"{s['median']:.3g} / {s['p99']:.3g} / {s['max']:.3g}")
                L.append(f"| {g} | {sub} (n={R['lc'][k0]['n']}) | " + " | ".join(cells) + " |")
        L.append("")
    if "bound" in R:
        L += ["## b. Constraint active / unphysical (fraction of rows)", "",
              "| grid | law | bound r | bound μ | unphysical r | unphysical μ |",
              "|---|---|---|---|---|---|"]
        for g in GRIDS:
            for law in LAWS:
                a, b = R["bound"][f"{g}|{law}|r"], R["bound"][f"{g}|{law}|mu"]
                L.append(f"| {g} | {law} | {a['bound']:.4f} | {b['bound']:.4f} | "
                         f"{a['unphysical']:.4f} | {b['unphysical']:.4f} |")
        L.append("")
    if "edge" in R:
        L += ["## c. Chromatic stellar edge (PHOENIX, per model, (max−min)/mean of r_c)", "",
              "| bands | n | median | 90th | max | depth error for a 1% transit (ppm): med / 90th / max |",
              "|---|---|---|---|---|---|"]
        for k, s in R["edge"].items():
            d = s["depth_ppm_1pct"]
            L.append(f"| {k} | {s['n_bands']} | {s['median']:.2e} | {s['p90']:.2e} | {s['max']:.2e} | "
                     f"{d['median']:.2f} / {d['p90']:.2f} / {d['max']:.2f} |")
        L.append("")
    if "corr" in R:
        L += ["## d. Δw corrections on MPS-ATLAS Set 2: fraction leaving the physical region", ""]
        for k, s in R["corr"].items():
            if k != "_conventions":
                L.append(f"* {k}: " + ", ".join(f"{a}={b:.4g}" if isinstance(b, float) else f"{a}={b}"
                                               for a, b in s.items() if b is not None))
        L += ["", R["corr"]["_conventions"], ""]
    if "cds" in R:
        L += ["## e. Passbands and CDS records", ""]
        for k, v in R["cds"].items():
            if k.startswith("bands|"):
                L.append(f"* {k[6:]} ({len(v)}): {', '.join(v)}")
        L.append("")
        for k, v in R["cds"].items():
            if k.startswith("records|"):
                L.append(f"* {k[8:]}: {v}")
        L.append("")
    return "\n".join(L)


def main():
    R = json.loads(OUT.read_text()) if OUT.exists() else {}
    if "--only" in sys.argv:
        parts = sys.argv[sys.argv.index("--only") + 1].split(",")
    else:
        parts = ["lc", "bound", "edge", "corr", "cds"]
    fns = dict(lc=lc_stats, bound=bound_stats, edge=edge_stats, corr=correction_stats, cds=cds_stats)
    for p in parts:
        print("computing", p, flush=True)
        R[p] = fns[p]()
    OUT.write_text(json.dumps(R, indent=1))
    (HERE / "production_stats.md").write_text(md(R))
    print(md(R))


if __name__ == "__main__":
    main()
