"""P3: 3D and magnetic correction layer, as weight shifts Δw to add to 1D rows.

    .venv/bin/python scripts/p3_corrections.py

Corrections are computed in profile space and refitted with the v1 method:

  * magnetic (MURaM): C_B(μ) = I_B / I_hydro for B ∈ {ssd, B100G, B200G,
    B300G}, the same simulation code with and without surface fields;
  * 3D (MURaM): C(μ) = I_hydro / I_1D with I_1D = MPS-ATLAS Set 2 (the code
    used to synthesise the MURaM spectra) at the nearest grid point;
  * 3D (Stagger): Δw = fit(Stagger) − fit(Kurucz) at identical (Teff, log g,
    [M/H]) from the ExoTiC-LD bundle (both plane-parallel, μ ≥ 0.01).

For MURaM the ratio is a PCHIP through its 10 nodes (μ = 0.1–1) and is held
constant below μ = 0.1, where MURaM has no data (flagged
`extrapolated_below_0.1`). Δw = fit(I_base · C) − fit(I_base), with I_base the
MPS-ATLAS Set 2 band profile. Each row also gives the light-curve change the
correction implies (ppm; production LC grid).

Output: data/products/corrections.csv
"""

import csv

import numpy as np
from scipy.interpolate import PchipInterpolator

from ldhybrid import bands as B
from ldhybrid import fit, laws, v1
from ldhybrid import lightcurve as lc
from ldhybrid.grids import exoticld, mpsatlas, muram
from ldhybrid.paths import DATA
from ldhybrid.products import PROD_LC, constraint_for
from ldhybrid.profile import Profile, band_profile

LAWS = ("hybrid2", "hybrid4", "hybrid5")
BANDS = [b for b in B.BANDS if b != "JWST_MIRI_LRS"]
G = lc.LCGrid(**PROD_LC)
COLS = {n: G.linear_columns(laws.get(n)) for n in LAWS}


def fitw(name, pr):
    return fit.fit_linear(laws.get(name), pr, v1.MEASURE, v1.RESAMPLE, constraint_for(name), v1.FCM).w


def corrected(base: Profile, C):
    return Profile(base.mu, base.y * C(base.mu), dict(base.meta))


def ratio(num: Profile, den: Profile):
    mu = num.mu
    r = num.y / den.interpolant()(mu)
    f = PchipInterpolator(mu, r)
    return lambda m: f(np.clip(m, mu[0], 1.0))


def nearest_mps(teff, logg, mh):
    cat = mpsatlas.catalogue("set2")
    return min(cat, key=lambda c: (abs(c[0] - mh) * 10, abs(c[1] - teff) / 100 + abs(c[2] - logg) * 5))


def row(source, kind, state, meta, band, name, w0, dw, extra=""):
    one, cols = COLS[name]
    dlc = G.errors(*G.combine(one, cols, w0 + dw), *G.combine(one, cols, w0))["max_ppm"]
    return dict(source=source, kind=kind, state=state, star=meta.get("star", ""),
                teff=round(meta["teff"]), logg=meta["logg"], mh=meta["mh"], band=band, law=name,
                **{f"w{j + 1}": f"{x:.6f}" for j, x in enumerate(w0)},
                **{f"dw{j + 1}": f"{x:+.6f}" for j, x in enumerate(dw)},
                dlc_ppm=f"{dlc:.2f}", flags=extra)


def main():
    rows = []
    for star in muram.stars():
        par = muram.params(star)
        I = {s: muram.load(star, s) for s in muram.STATES}
        mkey = nearest_mps(par["teff"]["hydro"], par["logg"], par["mh"])
        have_1d = abs(mkey[1] - par["teff"]["hydro"]) <= 150
        base_spec = mpsatlas.load("set2", *mkey)
        print(star, "1D partner", mkey if have_1d else "none", flush=True)
        for b in BANDS:
            lam, S = B.response(b)
            if min(I["hydro"].coverage(lam, S), base_spec.coverage(lam, S)) < 0.999:
                continue
            base = band_profile(base_spec, lam, S)
            prof = {s: band_profile(I[s], lam, S) for s in muram.STATES}
            for name in LAWS:
                w0 = fitw(name, base)
                meta = dict(par, teff=par["teff"]["hydro"])
                for s in muram.STATES[1:]:
                    dw = fitw(name, corrected(base, ratio(prof[s], prof["hydro"]))) - w0
                    rows.append(row("MURaM", "magnetic", s, meta, b, name, w0, dw,
                                    "extrapolated_below_0.1"))
                if have_1d:
                    dw = fitw(name, corrected(base, ratio(prof["hydro"], base))) - w0
                    rows.append(row("MURaM", "3d", "hydro", meta, b, name, w0, dw,
                                    f"extrapolated_below_0.1;1D=MPS2{mkey}"))
    # Stagger vs Kurucz
    for p in sorted((exoticld.ROOT / "stagger").glob("MH*/teff*/logg*/stagger_spectra.dat")):
        mh, t, g = float(p.parts[-4][2:]), int(p.parts[-3][4:]), float(p.parts[-2][4:])
        if not (g >= 3.5 and mh >= -2.0 and exoticld.path("kurucz", mh, t, g).exists()):
            continue
        st, ku = exoticld.load("stagger", mh, t, g), exoticld.load("kurucz", mh, t, g)
        print("Stagger/Kurucz", mh, t, g, flush=True)
        for b in BANDS:
            lam, S = B.response(b)
            if min(st.coverage(lam, S), ku.coverage(lam, S)) < 0.999:
                continue
            pk, ps = band_profile(ku, lam, S), band_profile(st, lam, S)
            for name in LAWS:
                w0 = fitw(name, pk)
                rows.append(row("Stagger", "3d", "hydro", dict(teff=t, logg=g, mh=mh), b, name,
                                w0, fitw(name, ps) - w0, "1D=Kurucz"))
    out = DATA / "products" / "corrections.csv"
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(out, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=keys)
        wr.writeheader()
        wr.writerows(rows)
    print(f"{len(rows)} rows -> {out}")


if __name__ == "__main__":
    main()
