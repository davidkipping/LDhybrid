"""C10: size and assumptions of the MURaM Δw corrections.

(a) Light-curve size of the magnetic corrections as applied by the code
    (predict.magnetic: Δw linear in Teff between the [M/H] = 0 MURaM stars,
    clamped) to MPS-ATLAS Set 2 hybrid5 weights at grid points within the
    MURaM Teff range: |ΔF| between w and w + Δw, production metric, by
    passband class.
(b) The ratio C(μ) is held constant below μ = 0.1 (no MURaM data there).
    Alternative: continue it linearly in μ from its two innermost nodes.
    Light-curve difference between the two Δw (hybrid5) for the six MURaM
    stars in several passbands.
(c) The MURaM stars available.

    .venv/bin/python studies/rev2/c10_corrections.py
"""
import json
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import PchipInterpolator

import rc
from ldhybrid import bands as B
from ldhybrid import fit, laws, predict, v1
from ldhybrid.grids import mpsatlas, muram
from ldhybrid.products import constraint_for
from ldhybrid.profile import Profile, band_profile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import p3_corrections as P3  # noqa: E402

STATES = ("ssd", "B100G", "B200G")
BANDS_B = ["Kepler", "TESS", "HST_WFC3_G141", "JWST_NIRSpec_PRISM", "JWST_NIRSpec_G395H", "2MASS_H"]


def band_class(b):
    lam, S = B.response(b)
    lc_ = np.sum(lam * S) / np.sum(S)
    return "optical (<1 µm)" if lc_ < 1e4 else ("near-IR (1–2.5 µm)" if lc_ < 2.5e4 else "IR (>2.5 µm)")


def part_a(n=3000, seed=0):
    h = h5py.File(ROOT / "data" / "products" / "mps2.h5")
    hb = [b.decode() for b in h["bands"][()]]
    T = h["teff"][()]
    rows = predict._corrections()
    corr_bands = sorted({r["band"] for r in rows if r["source"] == "MURaM"})
    Tm = sorted({float(r["teff"]) for r in rows if r["source"] == "MURaM" and float(r["mh"]) == 0})
    ok = np.flatnonzero((T >= Tm[0]) & (T <= Tm[-1]))
    rng = np.random.default_rng(seed)
    L = laws.get("hybrid5")
    g = rc.lcgrid()
    one, cols = g.linear_columns(L)
    out = []
    for _ in range(n):
        i = ok[rng.integers(len(ok))]
        b = corr_bands[rng.integers(len(corr_bands))]
        if b not in hb or not h["covered"][i, hb.index(b)]:
            continue
        w = h["hybrid5|r|w"][i, hb.index(b)]
        r = dict(teff=float(T[i]), band=b, cls=band_class(b))
        for s in STATES:
            dw = predict.magnetic("hybrid5", b, float(T[i]), s)
            r[s] = g.errors(*g.combine(one, cols, w + dw), *g.combine(one, cols, w))["max_ppm"]
        out.append(r)
    return out, Tm


def ratio_linear_below(num, den, mu_cut=0.1):
    mu = num.mu
    r = num.y / den.interpolant()(mu)
    f = PchipInterpolator(mu, r)
    m0, m1 = mu[0], mu[1]
    s = (r[1] - r[0]) / (m1 - m0)
    return lambda m: np.where(m < m0, r[0] + s * (m - m0), f(np.clip(m, m0, 1.0)))


def part_b_task(star):
    par = muram.params(star)
    I = {s: muram.load(star, s) for s in ("hydro",) + STATES}
    mkey = P3.nearest_mps(par["teff"]["hydro"], par["logg"], par["mh"])
    base_spec = mpsatlas.load("set2", *mkey)
    L = laws.get("hybrid5")
    g = rc.lcgrid()
    out = []
    for b in BANDS_B:
        lam, S = B.response(b)
        if min(I["hydro"].coverage(lam, S), base_spec.coverage(lam, S)) < 0.999:
            continue
        base = band_profile(base_spec, lam, S)
        prof = {s: band_profile(I[s], lam, S) for s in I}
        w0 = P3.fitw("hybrid5", base)
        for s in STATES:
            dwc = P3.fitw("hybrid5", P3.corrected(base, P3.ratio(prof[s], prof["hydro"]))) - w0
            dwl = P3.fitw("hybrid5", P3.corrected(base, ratio_linear_below(prof[s], prof["hydro"]))) - w0
            f = lambda w: (lambda m: L.intensity(w, m))
            d_lc = g.errors(*g.columns(f(w0 + dwl)), *g.columns(f(w0 + dwc)))["max_ppm"]
            size = g.errors(*g.columns(f(w0 + dwc)), *g.columns(f(w0)))["max_ppm"]
            out.append(dict(star=star, teff=par["teff"]["hydro"], band=b, state=s, size_ppm=size,
                            d_extrap_ppm=d_lc, dw_max=float(np.max(np.abs(dwl - dwc))),
                            mu_min=float(prof["hydro"].mu[0])))
    return out


def main():
    a, Tm = part_a()
    stars = muram.stars()
    b = [r for res in rc.run_parallel(part_b_task, stars) for r in res]
    c = [dict(star=s, **{k: v for k, v in muram.params(s).items() if k != "star"}) for s in stars]
    (rc.HERE / "c10_corrections.json").write_text(json.dumps(dict(a=a, b=b, c=c), default=str))
    L = ["## C10 — MURaM corrections", "",
         f"(a) Light-curve size of the correction as applied to MPS-ATLAS Set 2 hybrid5 weights "
         f"(grid points {Tm[0]:.0f}–{Tm[-1]:.0f} K, all log g and [M/H]; {len(a)} point-bands), ppm, median / 90th / max:", "",
         "| class | " + " | ".join(STATES) + " |", "|---|" + "---|" * len(STATES)]
    for cls in sorted({r["cls"] for r in a}) + ["all"]:
        sel = [r for r in a if cls == "all" or r["cls"] == cls]
        L.append(f"| {cls} ({len(sel)}) | " + " | ".join(
            f"{np.median([r[s] for r in sel]):.0f} / {np.percentile([r[s] for r in sel], 90):.0f} / "
            f"{np.max([r[s] for r in sel]):.0f}" for s in STATES) + " |")
    L += ["", "(b) Constant vs linearly extrapolated ratio below μ = 0.1 (MURaM stars, hybrid5): "
          "light-curve difference between the two Δw, ppm (median / max), and the size of the correction:", "",
          "| state | Δ (extrapolation) median / max | correction size median / max | max |Δ Δw| |", "|---|---|---|---|"]
    for s in STATES:
        sel = [r for r in b if r["state"] == s]
        L.append(f"| {s} | {np.median([r['d_extrap_ppm'] for r in sel]):.1f} / {np.max([r['d_extrap_ppm'] for r in sel]):.1f} | "
                 f"{np.median([r['size_ppm'] for r in sel]):.0f} / {np.max([r['size_ppm'] for r in sel]):.0f} | "
                 f"{np.max([r['dw_max'] for r in sel]):.3f} |")
    L += ["", "(c) MURaM stars: " + "; ".join(
        f"{r['star']} Teff(hydro) {r['teff']['hydro']:.0f} K, log g {r['logg']}, [M/H] {r['mh']}" for r in c),
          f"Magnetic Δw are interpolated in Teff between the [M/H] = 0 stars ({', '.join(f'{t:.0f}' for t in Tm)} K), "
          "clamped at the ends, and applied at every log g and [M/H] (predict.magnetic)."]
    (rc.HERE / "c10_corrections.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
