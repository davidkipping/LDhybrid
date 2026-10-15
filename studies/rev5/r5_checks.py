"""Round-3 referee checks (computations behind the round-3 revision).

Items (numbering as in r5_checks.md):
  1  Sun: laws fitted to the Neckel & Labs (1994) polynomials, hybrid4 with the
     adopted poles and curvature basis (p4a_sun part (i)).
  2  Neckel & Labs polynomials vs the curvature conditions (I'' <= 0,
     I''' >= 0 in x = mu^2), analytic derivatives.
  3  Prior implied by uniform q for hybrid4 (q-map) and hybrid5 (stick-breaking);
     tabulated hybrid4 fits in NIRSpec G395H.
  4  Tabulated hybrid5 r-weighted fits vs the curvature conditions.
  5  Grey LTE (Eddington) monochromatic profiles vs the curvature conditions.
  6  Narrow bins from the PHOENIX hybrid4 coefficient spectra: outside the
     curvature simplex vs violating conditions I-II (0.1-5.3 um).
  7  3D-corrected hybrid4 sets (MURaM hydro, Stagger) vs the curvature simplex.
  8  The PHOENIX r-weighted hybrid4 rows bound on w3 = 0 vs the hybrid5
     limb-brightening rows.
  9  Bound uniform-mu PHOENIX hybrid4 rows: sigma_mu with the curvature simplex
     vs the region of conditions I-II alone (random sample).
  10 Fisher best-measured hybrid4 (and hybrid5) combinations.
  11 Chord-velocity bias in the bias test (from the stored json).
  12 MURaM stars and passband coverage of corrections.dat.

    .venv/bin/python studies/rev5/r5_checks.py [items...]   -> r5_checks.md
Each item caches its result in studies/rev5/r5_cache/; the report is rebuilt
from all cached items.
"""
import json
import os
import sys
import warnings
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent
CACHE = HERE / "r5_cache"
CACHE.mkdir(exist_ok=True)
sys.path[:0] = [str(ROOT / "studies" / "p4"), str(ROOT / "studies" / "rev")]

from ldhybrid import fit, laws, qspace, v1  # noqa: E402
from ldhybrid.interpolate import LDTables  # noqa: E402
from ldhybrid.paths import DATA, RAW  # noqa: E402

CDS = DATA / "release" / "cds"
PROF = DATA / "release" / "zenodo" / "ldhybrid_v1.0_band_profiles_{}.h5"
_T = {}


def tables():
    if "t" not in _T:
        _T["t"] = LDTables(CDS)
    return _T["t"]


def pct(x, q):
    return float(np.percentile(x, q))


# --- curvature derivatives of the hybrid laws in x = mu^2 ---------------------
MU = np.unique(np.r_[0.0, np.geomspace(1e-9, 1e-2, 1000), np.linspace(1e-2, 1, 3001)])
XI = MU * MU


def curv_matrices(law):
    """d2T/dx2 and d3T/dx3 of the law's shapes on XI (rows), as in
    studies/rev3/h4_kipping_check.py."""
    eps = v1.HYBRID4_EPS if law == "hybrid4" else v1.HYBRID5_EPS
    a = [1 / (e ** -2 - (1 + e) ** -2) for e in eps]
    z = 0 * XI
    D2 = np.column_stack([-2 + z, 2 + z] + [6 * ak * (XI + e) ** -4 for ak, e in zip(a, eps)])
    D3 = np.column_stack([z, z] + [-24 * ak * (XI + e) ** -5 for ak, e in zip(a, eps)])
    if law == "hybrid4":
        D2, D3 = D2 @ v1.HYBRID4_BASIS.T, D3 @ v1.HYBRID4_BASIS.T
    return D2, D3


def curvature_violations(law, W, tol=1e-9, chunk=2000):
    """Per row: violates I'' <= 0 / I''' >= 0 somewhere, and the largest mu of a
    violation (relative tolerance, as the Kipping check)."""
    D2, D3 = curv_matrices(law)
    v2, v3, mu2, mu3 = [], [], [], []
    for i in range(0, len(W), chunk):
        w = W[i:i + chunk]
        d2, d3 = -(w @ D2.T), -(w @ D3.T)
        b2 = d2 > tol * np.maximum(1, np.abs(d2))
        b3 = d3 < -tol * np.maximum(1, np.abs(d3))
        v2.append(b2.any(1))
        v3.append(b3.any(1))
        mu2.append(np.where(b2, MU, -1).max(1))
        mu3.append(np.where(b3, MU, -1).max(1))
    return tuple(np.concatenate(x) for x in (v2, v3, mu2, mu3))


# --- items ---------------------------------------------------------------------
def item1():
    import pandas as pd
    import p4a_sun as S
    from ldhybrid.profile import Profile
    nl = pd.read_csv(RAW / "validation" / "neckel_labs" / "neckel_labs_1994_tableI.csv", comment="#")
    rows = []
    for _, r in nl.iterrows():
        a = np.array([r[f"a{k}"] for k in range(6)])
        obs = Profile(S.MU, np.polynomial.polynomial.polyval(S.MU, a))
        rows.append(dict(lam=float(r.wavelength_nm),
                         **{n: S.law_fit_error(n, obs) for n in ("hybrid4", "hybrid5", "claret4")}))
    out = {}
    for lo, hi in ((300, 450), (450, 700), (700, 1100)):
        sel = [x for x in rows if lo <= x["lam"] < hi]
        out[f"{lo}-{hi}"] = {n: [float(np.median([x[n] for x in sel])), float(max(x[n] for x in sel)), len(sel)]
                             for n in ("hybrid4", "hybrid5", "claret4")}
    allv = {n: [float(np.median([x[n] for x in rows])), float(max(x[n] for x in rows))]
            for n in ("hybrid4", "hybrid5", "claret4")}
    return dict(ranges=out, all=allv, rows=rows)


def item2():
    import pandas as pd
    nl = pd.read_csv(RAW / "validation" / "neckel_labs" / "neckel_labs_1994_tableI.csv", comment="#")
    mu = np.r_[np.geomspace(1e-4, 1e-2, 200), np.linspace(1e-2, 1, 5001)[1:]]
    P = np.polynomial.polynomial
    res = []
    for _, r in nl.iterrows():
        a = np.array([r[f"a{k}"] for k in range(6)])
        f1, f2, f3 = (P.polyval(mu, P.polyder(a, k)) for k in (1, 2, 3))
        I2 = (mu * f2 - f1) / (4 * mu ** 3)                       # d2I/dx2
        I3 = (mu ** 2 * f3 - 3 * mu * f2 + 3 * f1) / (8 * mu ** 5)  # d3I/dx3
        # relative size: compare with the magnitude of the terms
        s2 = (np.abs(mu * f2) + np.abs(f1)) / (4 * mu ** 3)
        s3 = (np.abs(mu ** 2 * f3) + 3 * np.abs(mu * f2) + 3 * np.abs(f1)) / (8 * mu ** 5)
        b2, b3 = I2 > 1e-12 * s2, I3 < -1e-12 * s3
        res.append(dict(lam=float(r.wavelength_nm), mono=bool(np.all(f1 >= 0)),
                        v2=bool(b2.any()), v3=bool(b3.any()),
                        v2_mu=[float(mu[b2].min()), float(mu[b2].max())] if b2.any() else None,
                        v3_mu=[float(mu[b3].min()), float(mu[b3].max())] if b3.any() else None,
                        v2_rel=float((I2[b2] / s2[b2]).max()) if b2.any() else 0.0,
                        v3_rel=float((-I3[b3] / s3[b3]).max()) if b3.any() else 0.0,
                        I3_min=float(I3[b3].min()) if b3.any() else 0.0))
    return dict(n=len(res), both_ok=sum(not (x["v2"] or x["v3"]) for x in res),
                bad=[x for x in res if x["v2"] or x["v3"]], nonmono=[x["lam"] for x in res if not x["mono"]])


def _prior_stats(L, W):
    mu = np.array([0.0, 1 / 3, 2 / 3])
    I = 1.0 - W @ L.T(mu).T
    I0, h1 = I[:, 0], I[:, 2]
    h2 = I[:, 2] - I[:, 1]
    return dict(I0=[pct(I0, 5), pct(I0, 50), pct(I0, 95)], P05=float(np.mean(I0 > 0.5)),
                P07=float(np.mean(I0 > 0.7)), h1=[pct(h1, 5), pct(h1, 95)], h2=[pct(h2, 5), pct(h2, 95)],
                I0_90=[pct(I0, 5), pct(I0, 95)])


def item3():
    rng = np.random.default_rng(20261015)
    n = 200000
    W4 = qspace.hybrid4_from_q(rng.random((n, 4)))
    q = rng.random((n, 5))
    W5 = np.empty((n, 5))
    stick = np.ones(n)
    for j in range(5):                       # w_j = (1 - q_j^{1/(n+1-j)}) prod_{i<j} q_i^{1/(n+1-i)}
        r = q[:, j] ** (1.0 / (5 - j))
        W5[:, j] = stick * (1 - r)
        stick = stick * r
    out = dict(hybrid4=_prior_stats(laws.get("hybrid4"), W4), hybrid5=_prior_stats(laws.get("hybrid5"), W5))
    # tabulated fits
    t = tables()
    L4 = laws.get("hybrid4")
    tab = {}
    for grid, (te, lg, mh) in (("phoenix", (3000, 5.0, 0.0)), ("mps2", (5600, 4.5, 0.0))):
        d = t.table(grid, "hybrid4")
        for band in ("JWST_NIRSpec_G395H", "JWST_NIRSpec_PRISM", "JWST_MIRI_LRS", "TESS"):
            k = np.flatnonzero((d["teff"] == te) & (d["logg"] == lg) & (d["mh"] == mh) & (d["band"] == band)
                               & (d["method"] == "r"))
            if not len(k):
                continue
            w = d["coeffs"][k[0]]
            I = L4.intensity(w, np.array([0.0, 1 / 3, 2 / 3]))
            tab[f"{grid} {te} {lg} {mh} {band}"] = dict(I0=float(I[0]), h1=float(I[2]), h2=float(I[2] - I[1]),
                                                       q=d["q"][k[0]].tolist(), w=w.tolist())
    out["tables"] = tab
    # where the tabulated r fits sit in q1, by band (PHOENIX + MPS2)
    q1 = {}
    for grid in ("phoenix", "mps2"):
        d = t.table(grid, "hybrid4")
        m = d["method"] == "r"
        for band in ("TESS", "HST_WFC3_G141", "JWST_NIRSpec_PRISM", "JWST_NIRSpec_G395H", "JWST_MIRI_LRS"):
            k = m & (d["band"] == band)
            if k.any():
                I0 = 1 - d["coeffs"][k].sum(1)
                q1[f"{grid} {band}"] = dict(q1_med=float(np.median(d["q"][k, 0])), q1_min=float(d["q"][k, 0].min()),
                                            I0_med=float(np.median(I0)), I0_max=float(I0.max()),
                                            frac_I0_gt_prior95=float(np.mean(I0 > out["hybrid4"]["I0"][2])))
    out["q1_by_band"] = q1
    return out


def item4():
    t = tables()
    res = {}
    rng = np.random.default_rng(4)
    for grid, nsamp in (("phoenix", None), ("mps2", 60000)):
        d = t.table(grid, "hybrid5")
        idx = np.flatnonzero(d["method"] == "r")
        if nsamp:
            idx = rng.choice(idx, nsamp, replace=False)
        W = d["coeffs"][idx]
        v2, v3, m2, m3 = curvature_violations("hybrid5", W)
        any_ = v2 | v3
        # violations only at mu < 0.1 / 0.05?
        mx = np.maximum(m2, m3)
        res[grid] = dict(rows=int(len(idx)), sampled=bool(nsamp), v2=float(v2.mean()), v3=float(v3.mean()),
                         any=float(any_.mean()),
                         any_mu_ge_005=float(np.mean(any_ & (mx >= 0.05))),
                         any_mu_ge_01=float(np.mean(any_ & (mx >= 0.1))),
                         maxmu_pcts=[pct(mx[any_], q) for q in (50, 90, 100)] if any_.any() else None,
                         bands_top=Counter(d["band"][idx][any_]).most_common(6),
                         teff_range=[float(d["teff"][idx][any_].min()), float(d["teff"][idx][any_].max())]
                         if any_.any() else None)
    return res


def _jet_mul(a, b):
    """Truncated (order 3) product of Taylor coefficient arrays."""
    return np.array([sum(a[i] * b[k - i] for i in range(k + 1)) for k in range(4)])


def _jet_apply(a, f):
    """f(a) for a Taylor series a, given f = (f(a0), f'(a0), f''(a0), f'''(a0))."""
    d = np.array([0.0 * a[0], a[1], a[2], a[3]])
    d2 = _jet_mul(d, d)
    d3 = _jet_mul(d2, d)
    out = f[1] * d + f[2] / 2 * d2 + f[3] / 6 * d3
    out[0] = out[0] + f[0]
    return out


def _grey_B_derivs(tau, teff, c2):
    """B(tau) = 1/(exp(c2/T) - 1), T^4 = 3/4 Teff^4 (tau + 2/3), and its first
    three tau-derivatives, by Taylor (jet) arithmetic."""
    A = (0.75 * teff ** 4) ** 0.25
    x = tau + 2.0 / 3.0
    p = 0.25
    T = _jet_apply(np.array([x, 1.0, 0.0, 0.0]),
                   (A * x ** p, A * p * x ** (p - 1), A * p * (p - 1) * x ** (p - 2),
                    A * p * (p - 1) * (p - 2) * x ** (p - 3)))
    t0 = T[0]
    u = _jet_apply(T, (c2 / t0, -c2 / t0 ** 2, 2 * c2 / t0 ** 3, -6 * c2 / t0 ** 4))
    y0 = np.exp(u[0])
    y = _jet_apply(u, (y0, y0, y0, y0))
    m = y0 - 1.0
    B = _jet_apply(y, (1 / m, -1 / m ** 2, 2 / m ** 3, -6 / m ** 4))
    return B * np.array([1.0, 1.0, 2.0, 6.0])             # k! * coefficient = k-th derivative


def _grey_case(args):
    from scipy.integrate import quad
    teff, lam_nm, mus = args
    c2 = 1.438776877e-2 / (lam_nm * 1e-9)                       # hc/(lambda k) in K
    scale = float(_grey_B_derivs(2.0 / 3.0, teff, c2)[0])

    def f(mu, k):                       # d^k/dmu^k of I(mu) = int B(mu s) e^-s ds
        g = lambda s: s ** k * _grey_B_derivs(mu * s, teff, c2)[k] / scale * np.exp(-s)
        return quad(g, 0, np.inf, epsabs=0, epsrel=1e-11, limit=400)[0]

    F = np.array([[f(m, k) for k in range(4)] for m in mus])
    f1, f2, f3 = F[:, 1], F[:, 2], F[:, 3]
    I2 = (mus * f2 - f1) / (4 * mus ** 3)
    I3 = (mus ** 2 * f3 - 3 * mus * f2 + 3 * f1) / (8 * mus ** 5)
    s2 = (np.abs(mus * f2) + np.abs(f1)) / (4 * mus ** 3)
    s3 = (np.abs(mus ** 2 * f3) + 3 * np.abs(mus * f2) + 3 * np.abs(f1)) / (8 * mus ** 5)
    b2, b3 = I2 > 1e-8 * s2, I3 < -1e-8 * s3
    rng_ = lambda b: [float(mus[b].min()), float(mus[b].max())] if b.any() else None
    return dict(teff=teff, lam_nm=lam_nm, mono=bool(np.all(f1 >= 0)), v2=rng_(b2), v3=rng_(b3),
                I0_over_I1=float(F[0, 0] / F[-1, 0]))


def item5():
    mus = np.r_[np.geomspace(1e-3, 0.05, 40), np.linspace(0.05, 1, 191)[1:]]
    jobs = [(te, lam, mus) for te in (3000, 4000, 5800) for lam in (200, 300, 500, 1000)]
    with ProcessPoolExecutor(min(6, os.cpu_count() or 1)) as ex:
        return list(ex.map(_grey_case, jobs))


def item6():
    import scipy.sparse as sps
    import d14_narrow_bins as D
    from ldhybrid.specint import SpecInt
    files = sorted(f for f in D.SPEC.glob("*.npz") if ".tmp" not in f.name)
    rng = np.random.default_rng(3)
    pick = rng.choice(len(files), 150, replace=False)
    edges0 = np.load(files[0])["lam_edges"]
    dummy = SpecInt(edges0, np.array([0.0, 1.0]), np.zeros((len(edges0) - 1, 2)))
    Wm = {}
    for R in D.RES:
        e = D.bins(R, lo=1000.0, hi=53000.0)
        M = sps.vstack([sps.csr_matrix(dummy.bin_throughput(np.array([a - 1e-3, a, b, b + 1e-3]),
                                                            np.array([0.0, 1, 1, 0])))
                        for a, b in zip(e[:-1], e[1:])]).tocsr()
        Wm[R] = (M, e)
    L = laws.get("hybrid4")
    G = L.dT(D.X)
    cls = lambda t: "M" if t < 3900 else ("K" if t < 5200 else "FG")
    agg = {}
    for i in pick:
        z = np.load(files[i])
        c = cls(float(z["teff"]))
        for R, (M, e) in Wm.items():
            I1 = M @ z["I1"].astype(float)
            w = (M @ z["a_hybrid4"].astype(float)) / I1[:, None]
            mono = np.all(-(w @ G.T) >= -1e-10, axis=1)
            pos = 1.0 - w.sum(1) >= -1e-10
            bad12 = ~(mono & pos)
            badS = (w < -1e-9).any(1) | (w.sum(1) > 1 + 1e-9)
            lamc = (e[:-1] + e[1:]) / 2 / 1e4
            for reg, sel in (("all", lamc > 0), ("<0.3", lamc < 0.3), ("0.3-0.6", (lamc >= 0.3) & (lamc < 0.6)),
                             (">=0.6", lamc >= 0.6)):
                a = agg.setdefault(f"{R}|{c}|{reg}", dict(bins=0, bad12=0, badS=0, badS_only=0, lamS=[]))
                a["bins"] += int(sel.sum())
                a["bad12"] += int((bad12 & sel).sum())
                a["badS"] += int((badS & sel).sum())
                a["badS_only"] += int((badS & ~bad12 & sel).sum())
                if reg == "all":
                    a["lamS"] += lamc[badS].tolist()
    for a in agg.values():
        lam = a.pop("lamS", None)
        if lam is not None and len(lam):
            a["lamS_5_50_95"] = [pct(lam, 5), pct(lam, 50), pct(lam, 95)]
            a["lamS_frac_lt03"] = float(np.mean(np.array(lam) < 0.3))
    # basis sanity: summed spectra vs the band table (TESS), for 3 picked models
    from ldhybrid import api, bands
    t = tables()
    d = t.table("phoenix", "hybrid4")
    san = []
    for i in pick[:3]:
        z = np.load(files[i])
        te, lg, mh = int(z["teff"]), float(z["logg"]), float(z["mh"])
        mhs = "-0.0" if mh == 0 else f"{mh:+.1f}"
        ws, _ = api._spectrum(te, lg, mhs).band_coefficients(*bands.response("TESS"), "hybrid4")
        k = np.flatnonzero((d["teff"] == te) & (d["logg"] == lg) & (d["mh"] == mh) & (d["band"] == "TESS")
                           & (d["method"] == "r"))
        if len(k):
            san.append(dict(model=[te, lg, mh], spectra=np.round(ws, 4).tolist(),
                            table=np.round(d["coeffs"][k[0]], 4).tolist()))
    return dict(agg=agg, sanity=san, models=150)


def item7():
    t = tables()
    d = t.table("mps2", "hybrid4")
    m = d["method"] == "r"
    key = {(d["teff"][i], d["logg"][i], d["mh"][i], d["band"][i]): i for i in np.flatnonzero(m)}
    out = {}
    Lw = laws.get("hybrid4")
    for src in ("MURaM", "Stagger"):
        rows = [r for r in t.corrections() if r["Source"] == src and r["Kind"] == "3d" and r["Law"] == "hybrid4"]
        W, miss, approx = [], 0, 0
        for r in rows:
            if src == "MURaM":
                b = r["Base"].split()              # MPS2 mh teff logg
                k = (float(b[2]), float(b[3]), float(b[1]), r["Band"])
            else:
                k = (float(r["Teff"]), float(r["logg"]), float(r["mh"]), r["Band"])
            i = key.get(k)
            if i is None:                          # nearest grid point (same band)
                c = np.flatnonzero(m & (d["band"] == r["Band"]))
                if not len(c):
                    miss += 1
                    continue
                dist = (np.abs(d["teff"][c] - k[0]) / 100 + np.abs(d["logg"][c] - k[1]) / 0.5
                        + np.abs(d["mh"][c] - k[2]) / 0.1)
                i = c[np.argmin(dist)]
                approx += 1
            dw = np.array([float(r[f"dw{j + 1}"]) for j in range(4)])
            W.append(d["coeffs"][i] + dw)
        W = np.array(W)
        outS = (W < -1e-9).any(1) | (W.sum(1) > 1 + 1e-9)
        phys = np.array([fit.is_physical(Lw, w, tol=1e-9) for w in W])
        out[src] = dict(sets=len(W), matched_exact=len(W) - approx, nearest=approx, missing=miss,
                        outside_simplex=int(outS.sum()), violate_I_II=int((~phys).sum()),
                        min_w=float(W.min()), max_sum=float(W.sum(1).max()),
                        worst=W[outS][np.argsort(W[outS].min(1))[:3]].round(4).tolist() if outS.any() else [])
        if outS.any():
            out[src]["which_w_negative"] = Counter(int(np.argmin(w)) + 1 for w in W[outS]).most_common()
    return out


def item8():
    t = tables()
    d4, d5 = t.table("phoenix", "hybrid4"), t.table("phoenix", "hybrid5")
    m4 = (d4["method"] == "r") & (d4["Bound"] == 1)
    on_w3 = m4 & (d4["coeffs"][:, 2] == 0)
    on_sum = m4 & (d4["coeffs"].sum(1) > 1 - 1e-6)
    key5 = {(d5["teff"][i], d5["logg"][i], d5["mh"][i], d5["band"][i]): i
            for i in np.flatnonzero(d5["method"] == "r")}
    idx = np.flatnonzero(on_w3)
    j5 = np.array([key5[(d4["teff"][i], d4["logg"][i], d4["mh"][i], d4["band"][i])] for i in idx])
    b5 = d5["Bound"][j5] == 1
    e5 = d5["LCppm"][j5] > 15
    # the 57 hybrid5 rows with LC > 15 ppm
    all57 = np.flatnonzero((d5["method"] == "r") & (d5["LCppm"] > 15))
    k57 = {(d5["teff"][i], d5["logg"][i], d5["mh"][i], d5["band"][i]) for i in all57}
    k123 = {(d4["teff"][i], d4["logg"][i], d4["mh"][i], d4["band"][i]) for i in idx}
    # model profile non-monotonic?
    with h5py.File(str(PROF).format("phoenix"), "r") as h:
        P = {k: h[k][()] for k in ("teff", "logg", "mh", "mu", "y")}
        bands = [b.decode() for b in h["bands"][()]]
    pk = {(round(a, 3), round(b, 3), round(c, 3)): i for i, (a, b, c) in enumerate(zip(P["teff"], P["logg"], P["mh"]))}
    nonmono, upto, rise = [], [], []
    for i in idx:
        a = pk[(round(d4["teff"][i], 3), round(d4["logg"][i], 3), round(d4["mh"][i], 3))]
        j = bands.index(d4["band"][i])
        mu, y = P["mu"][a, j], P["y"][a, j]
        ok = np.isfinite(mu)
        mu, y = mu[ok], y[ok]
        o = np.argsort(mu)
        mu, y = mu[o], y[o]
        dec = np.diff(y) < -1e-12
        nonmono.append(bool(dec.any()))
        if dec.any():
            upto.append(float(mu[1:][dec].max()))
            rise.append(float(y[:np.searchsorted(mu, mu[1:][dec].max()) + 1].max() - y[np.searchsorted(mu, mu[1:][dec].max())]))
    mh = d4["mh"][idx]
    return dict(bound_r=int(m4.sum()), on_w3=int(on_w3.sum()), on_sum1=int(on_sum.sum()),
                both=int((on_w3 & on_sum).sum()),
                h5_bound=int(b5.sum()), h5_lc_gt15=int(e5.sum()), h5_bound_or_lc=int((b5 | e5).sum()),
                n57=len(all57), overlap_57=len(k57 & k123),
                nonmono=int(np.sum(nonmono)), nonmono_upto_mu_median=float(np.median(upto)) if upto else None,
                mh_counts=Counter(mh.tolist()).most_common(), teff=[float(d4["teff"][idx].min()), float(d4["teff"][idx].max())],
                h4_lc=[float(np.median(d4["LCppm"][idx])), float(d4["LCppm"][idx].max())],
                h5_lc=[float(np.median(d5["LCppm"][j5])), float(d5["LCppm"][j5].max())])


def _refit_mu(job):
    from ldhybrid.profile import Profile
    mu, y, w_tab = job
    L = laws.get("hybrid4")
    pr = Profile(mu, y)
    o = {}
    for c in ("simplex", "physical"):
        F = fit.fit_linear(L, pr, "mu", v1.RESAMPLE, c, v1.FCM)
        o[c] = (float(F.rms), float(F.max_abs), np.asarray(F.w).tolist())
    o["dtab"] = float(np.max(np.abs(np.array(o["simplex"][2]) - w_tab)))
    return o


def item9():
    t = tables()
    d = t.table("phoenix", "hybrid4")
    rows = np.flatnonzero((d["method"] == "mu") & (d["Bound"] == 1))
    rng = np.random.default_rng(9)
    pick = rng.choice(rows, min(1000, len(rows)), replace=False)
    with h5py.File(str(PROF).format("phoenix"), "r") as h:
        P = {k: h[k][()] for k in ("teff", "logg", "mh", "mu", "y")}
        bands = [b.decode() for b in h["bands"][()]]
    pk = {(round(a, 3), round(b, 3), round(c, 3)): i for i, (a, b, c) in enumerate(zip(P["teff"], P["logg"], P["mh"]))}
    jobs = []
    for r in pick:
        a = pk[(round(d["teff"][r], 3), round(d["logg"][r], 3), round(d["mh"][r], 3))]
        j = bands.index(d["band"][r])
        ok = np.isfinite(P["mu"][a, j])
        jobs.append((P["mu"][a, j][ok].astype(float), P["y"][a, j][ok].astype(float), d["coeffs"][r]))
    with ProcessPoolExecutor(min(6, os.cpu_count() or 1)) as ex:
        R = list(ex.map(_refit_mu, jobs, chunksize=20))
    s = np.array([x["simplex"][0] for x in R])
    p = np.array([x["physical"][0] for x in R])
    ms = np.array([x["simplex"][1] for x in R])
    mp = np.array([x["physical"][1] for x in R])
    rel = (s - p) / s
    return dict(total_bound=int(len(rows)), sample=int(len(pick)), max_dtab=float(max(x["dtab"] for x in R)),
                sigma_simplex=[float(np.median(s)), float(s.max())], sigma_I_II=[float(np.median(p)), float(p.max())],
                dsigma_abs=[float(np.median(s - p)), pct(s - p, 90), float((s - p).max())],
                dsigma_rel=[float(np.median(rel)), pct(rel, 90), float(rel.max())],
                maxabs_simplex=[float(np.median(ms)), float(ms.max())], maxabs_I_II=[float(np.median(mp)), float(mp.max())],
                tabulated_sigma_median=float(np.median(d["sigma"][pick])))


def _dirs(rows, law):
    out = {}
    for case in ("JWST-like", "TESS-like"):
        sel = [r for r in rows if r["law"] == law and r["case"] == case]
        if not sel:
            continue
        D1 = []
        for r in sel:
            v = np.array(json.loads(r["dirs"]) if isinstance(r["dirs"], str) else r["dirs"])[0]
            if v[np.argmax(np.abs(v))] < 0:
                v = -v
            D1.append(v)
        D1 = np.array(D1)
        ref = np.median(D1, 0)
        D1 = np.where((D1 @ ref)[:, None] < 0, -D1, D1)
        med = np.median(D1, 0)
        out[case] = dict(median=(med / np.linalg.norm(med)).round(3).tolist(),
                         min=D1.min(0).round(2).tolist(), max=D1.max(0).round(2).tolist(), n=len(sel))
    return out


def item10():
    r4 = json.load(open(ROOT / "studies" / "rev4" / "rv_fisher.json"))
    r2 = json.load(open(ROOT / "studies" / "rev2" / "c15_fisher.json"))
    return dict(hybrid4=_dirs(r4, "hybrid4"), hybrid5_rev2=_dirs(r2, "hybrid5"),
                why_nan="studies/rev4/rv_fisher.py ran hybrid4 only (rv_fisher.json has no hybrid5 rows); "
                        "its hybrid5 summary loop takes the median of an empty selection")


def item11():
    e = json.load(open(ROOT / "studies" / "rev2" / "e1b_bias.json"))
    r = json.load(open(ROOT / "studies" / "rev4" / "rv_bias.json"))
    out = {}
    for grid in ("phoenix", "mps2"):
        for law in ("hybrid5", "hybrid4", "claret4", "power2", "hybrid2", "quadratic"):
            src = r if law in ("hybrid4", "claret4") else e
            rows = {k: v for k, v in src.items() if k.startswith(grid + "|") and k.endswith("|" + law)}
            if not rows:
                continue
            o = {}
            for mode in ("fixed", "free", "prior"):
                dv = np.abs([v[mode]["dv"] for v in rows.values() if mode in v])
                dp = np.abs([v[mode]["dp2_ppm"] for v in rows.values() if mode in v])
                if len(dv):
                    o[mode] = dict(dv=[float(np.median(dv)), float(dv.max())], dp2=[float(np.median(dp)), float(dp.max())])
            out[f"{grid}|{law}"] = dict(n=len(rows), src="rev4/rv_bias" if src is r else "rev2/e1b_bias", **o)
    return out


def item12():
    t = tables()
    C = t.corrections()
    out = {}
    for (src, kind), grp in sorted(Counter((r["Source"], r["Kind"]) for r in C).items()):
        rows = [r for r in C if r["Source"] == src and r["Kind"] == kind]
        stars = sorted({(int(r["Teff"]), float(r["logg"]), float(r["mh"])) for r in rows})
        out[f"{src}|{kind}"] = dict(rows=grp, stars=stars if src == "MURaM" else len(stars),
                                    bands=len({r["Band"] for r in rows}))
    allb = set(t.bands("mps2")) | set(t.bands("phoenix"))
    corr_b = {r["Band"] for r in C if r["Source"] == "MURaM" and r["Kind"] == "magnetic"}
    corr_3d = {r["Band"] for r in C if r["Kind"] == "3d" and r["Source"] == "MURaM"}
    st_b = {r["Band"] for r in C if r["Source"] == "Stagger"}
    out["bands_total"] = len(allb)
    out["no_magnetic"] = sorted(allb - corr_b)
    out["no_muram3d"] = sorted(allb - corr_3d)
    out["no_stagger"] = sorted(allb - st_b)
    return out


ITEMS = {i: globals()[f"item{i}"] for i in range(1, 13)}


# --- report ----------------------------------------------------------------------
def f3(x):
    return f"{x:.3g}"


def report(R):
    L = ["# Round-3 revision checks", "", "Generated by `studies/rev5/r5_checks.py`. hybrid4: curvature basis, "
         "alpha = 25/16 poles.", ""]
    if 1 in R:
        L += ["## 1. Sun: laws fitted to the Neckel & Labs polynomials, worst LC error (ppm), median / max", "",
              "| range | N | hybrid4 | hybrid5 | four-term |", "|---|---|---|---|---|"]
        for k, v in R[1]["ranges"].items():
            L.append(f"| {k} nm | {v['hybrid4'][2]} | " + " | ".join(
                f"{v[n][0]:.2f} / {v[n][1]:.2f}" for n in ("hybrid4", "hybrid5", "claret4")) + " |")
        a = R[1]["all"]
        L += [f"| all | 30 | " + " | ".join(f"{a[n][0]:.2f} / {a[n][1]:.2f}" for n in ("hybrid4", "hybrid5", "claret4")) + " |",
              "", "(p4a_sun.md, old alpha = 1 poles: hybrid4 8.3/10.3/9.1 median by range.)", ""]
    if 2 in R:
        r = R[2]
        L += ["## 2. Neckel & Labs polynomials vs the curvature conditions", "",
              f"{r['both_ok']} of {r['n']} wavelengths satisfy I'' <= 0 and I''' >= 0 on 1e-4 <= mu <= 1; "
              f"non-monotonic: {r['nonmono'] or 'none'}.", ""]
        for x in r["bad"]:
            L.append(f"- {x['lam']:.1f} nm: I''>0 at mu {x['v2_mu']} (rel {x['v2_rel']:.2g}); "
                     f"I'''<0 at mu {x['v3_mu']} (rel {x['v3_rel']:.2g}, min I''' {x['I3_min']:.3g})")
        L.append("")
    if 3 in R:
        r = R[3]
        L += ["## 3. Prior implied by uniform q (2e5 draws)", "",
              "| law | I(0) 5/50/95% | P(I(0)>0.5) | P(I(0)>0.7) | h'1 90% | h'2 90% |", "|---|---|---|---|---|---|"]
        for n in ("hybrid4", "hybrid5"):
            s = r[n]
            L.append(f"| {n} | {s['I0'][0]:.3f} / {s['I0'][1]:.3f} / {s['I0'][2]:.3f} | {100 * s['P05']:.1f}% | "
                     f"{100 * s['P07']:.2f}% | {s['h1'][0]:.3f}–{s['h1'][1]:.3f} | {s['h2'][0]:.3f}–{s['h2'][1]:.3f} |")
        L += ["", "Tabulated hybrid4 r fits:", ""]
        for k, v in r["tables"].items():
            L.append(f"- {k}: I(0) = {v['I0']:.3f}, h'1 = {v['h1']:.3f}, h'2 = {v['h2']:.3f}, q1 = {v['q'][0]:.2e}")
        L += ["", "Tabulated hybrid4 r fits by band (q1 median/min; I(0) median/max; fraction with I(0) above the prior 95th pct):", ""]
        for k, v in r["q1_by_band"].items():
            L.append(f"- {k}: q1 {v['q1_med']:.3g} / {v['q1_min']:.2g}; I(0) {v['I0_med']:.3f} / {v['I0_max']:.3f}; "
                     f"{100 * v['frac_I0_gt_prior95']:.0f}%")
        L.append("")
    if 4 in R:
        L += ["## 4. Tabulated hybrid5 r fits vs the curvature conditions (dense mu grid, rel. tol 1e-9)", "",
              "| grid | rows | I''>0 | I'''<0 | either | either at mu>=0.05 | either at mu>=0.1 | largest mu of a violation 50/90/100% |",
              "|---|---|---|---|---|---|---|---|"]
        for g, v in R[4].items():
            mp = "/".join(f"{x:.3f}" for x in v["maxmu_pcts"]) if v["maxmu_pcts"] else "—"
            L.append(f"| {g}{' (sample)' if v['sampled'] else ''} | {v['rows']} | {100 * v['v2']:.2f}% | {100 * v['v3']:.2f}% | "
                     f"{100 * v['any']:.2f}% | {100 * v['any_mu_ge_005']:.2f}% | {100 * v['any_mu_ge_01']:.3f}% | {mp} |")
        for g, v in R[4].items():
            L.append(f"- {g}: Teff of violating rows {v['teff_range']}; bands {v['bands_top']}")
        L.append("")
    if 5 in R:
        L += ["## 5. Grey LTE (Eddington) monochromatic profiles", "",
              "| Teff | lambda | I(0)/I(1) | I''>0 at mu in | I'''<0 at mu in |", "|---|---|---|---|---|"]
        for v in R[5]:
            fm = lambda x: f"{x[0]:.3f}–{x[1]:.3f}" if x else "never"
            L.append(f"| {v['teff']} | {v['lam_nm']} nm | {v['I0_over_I1']:.3f} | {fm(v['v2'])} | {fm(v['v3'])} |")
        L.append("")
    if 6 in R:
        r = R[6]
        L += [f"## 6. Narrow bins from the PHOENIX hybrid4 coefficient spectra ({r['models']} models, seed 3; 0.1–5.3 um)", "",
              "% of bins violating conditions I–II / outside the curvature simplex (of which still I–II-physical)", "",
              "| R | class | all | < 0.3 um | 0.3–0.6 um | >= 0.6 um | outside-simplex bins: 5/50/95% lambda (um); share < 0.3 um |",
              "|---|---|---|---|---|---|---|"]
        for R_ in (100, 300, 1000):
            for c in ("M", "K", "FG"):
                cells = []
                for reg in ("all", "<0.3", "0.3-0.6", ">=0.6"):
                    a = r["agg"][f"{R_}|{c}|{reg}"]
                    n = max(a["bins"], 1)
                    cells.append(f"{100 * a['bad12'] / n:.2f} / {100 * a['badS'] / n:.2f} ({100 * a['badS_only'] / n:.2f})")
                a = r["agg"][f"{R_}|{c}|all"]
                lam = (" / ".join(f"{x:.2f}" for x in a["lamS_5_50_95"]) + f"; {100 * a['lamS_frac_lt03']:.0f}%"
                       if "lamS_5_50_95" in a else "—")
                L.append(f"| {R_} | {c} | " + " | ".join(cells) + f" | {lam} |")
        L += ["", "Basis sanity (TESS, summed spectra vs band table): " +
              "; ".join(f"{s['model']}: {s['spectra']} vs {s['table']}" for s in r["sanity"]), ""]
    if 7 in R:
        L += ["## 7. 3D-corrected hybrid4 sets (MPS-ATLAS Set 2 r weights + Δw)", ""]
        for k, v in R[7].items():
            L.append(f"- {k}: {v['sets']} sets ({v['matched_exact']} at the exact grid point, {v['nearest']} nearest, "
                     f"{v['missing']} unmatched); outside the curvature simplex {v['outside_simplex']}; violating I–II "
                     f"{v['violate_I_II']}; min w {v['min_w']:.4f}, max Σw {v['max_sum']:.4f}"
                     + (f"; negative weight index {v.get('which_w_negative')}; worst {v['worst']}" if v["outside_simplex"] else ""))
        L.append("")
    if 8 in R:
        r = R[8]
        L += ["## 8. PHOENIX r-weighted hybrid4 rows bound on w3 = 0", "",
              f"- bound r rows {r['bound_r']}: on w3 = 0 {r['on_w3']}, on Σw = 1 {r['on_sum1']} (both {r['both']})",
              f"- of the {r['on_w3']}: hybrid5 constraint binds in {r['h5_bound']}, hybrid5 LC > 15 ppm in {r['h5_lc_gt15']}, either {r['h5_bound_or_lc']}",
              f"- overlap with the {r['n57']} hybrid5 rows with LC > 15 ppm: {r['overlap_57']}",
              f"- model profile non-monotonic (brightens towards the limb): {r['nonmono']} of {r['on_w3']}; "
              f"brightening up to mu = {r['nonmono_upto_mu_median']} (median)",
              f"- [M/H] counts {r['mh_counts']}; Teff {r['teff']}; hybrid4 LC median/max {r['h4_lc']}; hybrid5 LC {r['h5_lc']}", ""]
    if 9 in R:
        r = R[9]
        L += [f"## 9. Bound uniform-mu PHOENIX hybrid4 rows: curvature simplex vs conditions I–II only ({r['sample']} of {r['total_bound']})", "",
              f"- refit reproduces the table to {r['max_dtab']:.1e}",
              f"- sigma_mu: curvature simplex median {r['sigma_simplex'][0]:.3g} (max {r['sigma_simplex'][1]:.3g}); "
              f"I–II only {r['sigma_I_II'][0]:.3g} (max {r['sigma_I_II'][1]:.3g})",
              f"- reduction in sigma_mu: absolute median {r['dsigma_abs'][0]:.2g}, 90th {r['dsigma_abs'][1]:.2g}, max {r['dsigma_abs'][2]:.2g}; "
              f"relative median {100 * r['dsigma_rel'][0]:.2f}%, 90th {100 * r['dsigma_rel'][1]:.1f}%, max {100 * r['dsigma_rel'][2]:.1f}%",
              f"- max |ΔI/I(1)|: simplex median {r['maxabs_simplex'][0]:.3g} (max {r['maxabs_simplex'][1]:.3g}); "
              f"I–II {r['maxabs_I_II'][0]:.3g} (max {r['maxabs_I_II'][1]:.3g})", ""]
    if 10 in R:
        r = R[10]
        L += ["## 10. Best-measured combination (unit vector, free geometry; median over cases, sign-aligned)", ""]
        for law in ("hybrid4", "hybrid5_rev2"):
            for case, v in r[law].items():
                L.append(f"- {law} {case} (n={v['n']}): {v['median']} (component ranges {v['min']} – {v['max']})")
        L += [f"- rv_fisher.md nan: {r['why_nan']}", ""]
    if 11 in R:
        L += ["## 11. Bias test: |Δv/v| and |Δ(p²)| (ppm), median / max", "",
              "| grid, law | n | source | fixed |Δv/v| | free | prior | fixed |Δp²| | free | prior |", "|---|---|---|---|---|---|---|---|---|"]
        for k, v in R[11].items():
            L.append(f"| {k} | {v['n']} | {v['src']} | " + " | ".join(
                f"{v[m]['dv'][0]:.1e} / {v[m]['dv'][1]:.1e}" if m in v else "—" for m in ("fixed", "free", "prior")) + " | "
                + " | ".join(f"{v[m]['dp2'][0]:.2f} / {v[m]['dp2'][1]:.1f}" if m in v else "—" for m in ("fixed", "free", "prior")) + " |")
        L.append("")
    if 12 in R:
        r = R[12]
        L += ["## 12. corrections.dat", ""]
        for k, v in r.items():
            if "|" in k:
                L.append(f"- {k}: {v['rows']} rows, {v['bands']} bands, stars {v['stars']}")
        L += [f"- passbands (union of grids): {r['bands_total']}",
              f"- no magnetic Δw: {r['no_magnetic']}", f"- no MURaM 3D Δw: {r['no_muram3d']}",
              f"- no Stagger Δw: {r['no_stagger']}", ""]
    return "\n".join(L) + "\n"


def main(argv):
    want = [int(a) for a in argv] or list(ITEMS)
    for i in want:
        print(f"item {i} ...", flush=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = ITEMS[i]()
        (CACHE / f"item{i}.json").write_text(json.dumps(res, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))
    R = {i: json.loads((CACHE / f"item{i}.json").read_text()) for i in ITEMS if (CACHE / f"item{i}.json").exists()}
    txt = report(R)
    (HERE / "r5_checks.md").write_text(txt)
    print(txt)


if __name__ == "__main__":
    main(sys.argv[1:])
