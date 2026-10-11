"""Revision E2: checks of the user-code additions (ldhybrid/interpolate.py,
ldhybrid/api.py): prior(), condition(), observables(), refit.

  1. prior() with σ → 0 returns the interpolated coefficients, cov → 0;
  2. condition() shrinks the variance along the observable's gradient and
     moves its mean by gCg/(gCg + σ²)·δ; observables() equal direct evaluation;
  3. refit at grid models reproduces the tabulated coefficients;
  4. leave-one-out across a 200 K Teff gap (PHOENIX, 200 trials): light-curve
     error (p = 0.1 production grid) of interpolated coefficients and of the
     refit of interpolated profiles, against the left-out model's profile.

    .venv/bin/python studies/rev/val_code_checks.py
"""

import json
import sys
from pathlib import Path

import h5py
import numpy as np

from ldhybrid import interpolate as I
from ldhybrid import laws
from ldhybrid import lightcurve as lc
from ldhybrid.paths import DATA
from ldhybrid.products import PROD_LC
from ldhybrid.profile import Profile

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parents[0] / "p4"))
import p4d_sing as P4D  # noqa: E402

LAWS = ("hybrid5", "hybrid4", "hybrid2", "claret4", "quadratic", "power2")


def law_fn(law, c):
    return lambda mu: I.intensity(law, c, mu)


def main():
    rng = np.random.default_rng(1)
    L = ["# E2: user-code checks", ""]
    out = {}

    # 1. prior at σ → 0 (CDS reader, mps2)
    t = I.LDTables(DATA / "release" / "cds", profiles=DATA / "products" / "profiles")
    d1 = []
    for teff, logg, mh, band in ((5600, 4.4, 0.1, "TESS"), (4300, 4.6, -0.3, "JWST_NIRSpec_PRISM"),
                                 (6200, 4.2, 0.2, "HST_WFC3_G141")):
        for law in ("hybrid5", "claret4"):
            p = t.prior(law, teff, 1e-6, logg, 1e-8, mh, 1e-8, band, "mps2", n=200)
            c = np.array(list(t.coefficients(law, teff, logg, mh, band, "mps2")["coeffs"].values()))
            d1.append((np.max(np.abs(p["mean"] - c)), np.sqrt(np.max(np.diag(p["cov"])))))
    d1 = np.array(d1)
    L.append(f"1. prior(σ → 0) − coefficients(): max |Δmean| {d1[:, 0].max():.1e}, "
             f"max sd {d1[:, 1].max():.1e} (6 cases, mps2 CDS tables)")
    p = t.prior("hybrid5", 5600, 100, 4.4, 0.1, 0.1, 0.1, "TESS", "mps2", systematic=True)
    L.append(f"   example hybrid5 TESS (5600±100, 4.4±0.1, 0.1±0.1), + systematic over {p['grids']}: "
             f"mean {np.round(p['mean'], 4).tolist()}, sd {np.round(np.sqrt(np.diag(p['cov'])), 4).tolist()}")

    # 2. condition and observables
    m2, C2 = I.condition(p["mean"], p["cov"], "hybrid5", "h1", 0.008, 0.006)
    g = I.observables("hybrid5")["h1"][1]
    s0, s1 = g @ p["cov"] @ g, g @ C2 @ g
    shift = g @ (m2 - p["mean"])
    L.append(f"2. condition(h1: δ = +0.008, σ = 0.006): sd(h1) {np.sqrt(s0):.4f} → {np.sqrt(s1):.4f} "
             f"(expected {np.sqrt(s0 * 0.006 ** 2 / (s0 + 0.006 ** 2)):.4f}); h1 shift {shift:+.4f} "
             f"(expected {s0 / (s0 + 0.006 ** 2) * 0.008:+.4f})")
    dev = 0.0
    mu = P4D.MUU
    for law in ("hybrid5", "hybrid4", "hybrid2", "claret4", "quadratic"):
        for _ in range(20):
            c = rng.uniform(0, 0.3, len(I.COEFFS[law]))
            o = I.observables(law, c)
            y = lambda m: I.intensity(law, c, m)
            ell, delta = P4D.quad_ell_delta(mu, y(mu))
            dev = max(dev, abs(o["h1"] - y(2 / 3)), abs(o["h2"] - y(2 / 3) + y(1 / 3)),
                      abs(o["ell"] - ell), abs(o["delta"] - delta))
    L.append(f"   observables() vs direct evaluation (h1, h2; ℓ, δ as p4d_sing): max |Δ| {dev:.1e}")

    # 3. refit at the models
    H = h5py.File(DATA / "products" / "phoenix.h5", "r")
    hb = [b.decode() for b in H["bands"][()]]
    hk = {(float(a), float(b), float(c)): i for i, (a, b, c) in
          enumerate(zip(H["teff"][()], H["logg"][()], H["mh"][()]))}
    P = t._profiles("phoenix")
    pk = {(float(a), float(b), float(c)): i for i, (a, b, c) in enumerate(zip(P["teff"], P["logg"], P["mh"]))}
    keys = list(hk)
    res3 = {law: [] for law in LAWS}
    for _ in range(40):
        k = keys[rng.integers(len(keys))]
        band = hb[rng.integers(len(hb))]
        if not H["covered"][hk[k], hb.index(band)]:
            continue
        for law in LAWS:
            r = t.refit(law, *k, band, "phoenix")
            w = H[f"{law}|r|w"][hk[k], hb.index(band)]
            res3[law].append((np.max(np.abs(np.array(list(r["coeffs"].values())) - w)),
                              bool(H[f"{law}|r|bound"][hk[k], hb.index(band)])))
    L.append("3. refit at PHOENIX models vs tabulated (40 random model-bands): " + ", ".join(
        f"{law} max |Δ| {max(x for x, b in v if not b):.1e}" for law, v in res3.items()))
    nb = sum(b for v in res3.values() for _, b in v)
    if nb:
        L.append(f"   rows with an active constraint ({nb}): max |Δ| "
                 f"{max(x for v in res3.values() for x, b in v if b):.1e}")

    # 4. leave-one-out across a 200 K gap
    trials = []
    while len(trials) < 200:
        k = keys[rng.integers(len(keys))]
        lo, hi = (k[0] - 100, k[1], k[2]), (k[0] + 100, k[1], k[2])
        if k[0] > 7000 or lo not in hk or hi not in hk or k not in pk:
            continue
        band = hb[rng.integers(len(hb))]
        j = hb.index(band)
        if not all(H["covered"][hk[x], j] for x in (k, lo, hi)):
            continue
        trials.append((k, lo, hi, band))
    jP = {b: P["bands"].index(b) for b in hb}
    rows = []
    for k, lo, hi, band in trials:
        jh, jp = hb.index(band), jP[band]
        nodes = []
        for x in (k, lo, hi):
            m = np.isfinite(P["mu"][pk[x], jp])
            nodes.append((P["mu"][pk[x], jp][m].astype(float), P["y"][pk[x], jp][m].astype(float)))
        truth = Profile(*nodes[0])
        from scipy.interpolate import PchipInterpolator
        mu_min = max(n[0][0] for n in nodes[1:])
        mq, qq = I.quadrature_r(np.concatenate([n[0][n[0] >= mu_min] for n in nodes[1:]] + [[mu_min]]))
        y1, y2 = (PchipInterpolator(*n)(mq) for n in nodes[1:])
        row = dict(teff=k[0], logg=k[1], mh=k[2], band=band)
        for law in ("hybrid5", "claret4", "power2"):
            wi = 0.5 * (H[f"{law}|r|w"][hk[lo], jh] + H[f"{law}|r|w"][hk[hi], jh])
            wr = I.fit_profile(law, mq, 0.5 * (y1 + y2), qq)
            row[law] = dict(coef=lc.lc_errors(law_fn(law, wi), truth, **PROD_LC)["max_ppm"],
                            refit=lc.lc_errors(law_fn(law, wr), truth, **PROD_LC)["max_ppm"],
                            dcoef=float(np.max(np.abs(wi - wr))))
        rows.append(row)
    L.append("4. leave-one-out across 200 K in Teff (PHOENIX, 200 trials), light-curve error vs the "
             "left-out profile, median / 90th pct (ppm):")
    for law in ("hybrid5", "claret4", "power2"):
        a = np.array([r[law]["coef"] for r in rows])
        b = np.array([r[law]["refit"] for r in rows])
        dc = np.array([r[law]["dcoef"] for r in rows])
        L.append(f"   {law}: interpolated coefficients {np.median(a):.1f} / {np.percentile(a, 90):.1f}; "
                 f"refit of interpolated profiles {np.median(b):.1f} / {np.percentile(b, 90):.1f}; "
                 f"max |Δw| between the two {dc.max():.1e}")
    L += ["", "For the linear laws an unconstrained r-weighted fit is a fixed linear map of the "
          "profile, so refitting an interpolated profile equals interpolating the coefficients "
          "wherever no constraint binds; the interpolation error is that of the profiles themselves. "
          "Refitting differs only for power-2 (nonlinear) and where a constraint is active."]
    out["loo"] = rows
    (HERE / "val_code_checks.json").write_text(json.dumps(out, default=float))
    (HERE / "val_code_checks.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
