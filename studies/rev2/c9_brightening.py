"""C9: the metal-poor cool PHOENIX rows whose hybrid5 light-curve error exceeds
15 ppm (production). Does the limb brightening survive a different edge?

For each row: the post-cut profile with the production M1 edge, with M2
(midway to the next node inward), and with the edge moved one and two nodes
inward. Reported per variant: whether the profile is non-monotonic (brightens
towards the limb), where, whether the hybrid5 constraint binds, and the
hybrid5 light-curve error against that variant's own profile.

    .venv/bin/python studies/rev2/c9_brightening.py
"""
import json
from pathlib import Path

import h5py
import numpy as np

import rc
from ldhybrid import bands
from ldhybrid.grids import phoenix
from ldhybrid.profile import band_profile, limb_cut

ROOT = Path(__file__).resolve().parents[2]


def rows_to_check():
    h = h5py.File(ROOT / "data" / "products" / "phoenix.h5")
    e = h["hybrid5|r|lc_ppm"][()]
    cov = h["covered"][()]
    B = [b.decode() for b in h["bands"][()]]
    i, j = np.where(cov & (e > 15))
    return [(float(h["teff"][a]), float(h["logg"][a]), float(h["mh"][a]), B[b], float(e[a, b]))
            for a, b in zip(i, j)]


def analyse(pr):
    y = pr.y
    dec = np.diff(y) < -1e-12                     # μ ascending: y must not decrease
    where = float(pr.mu[1:][dec].max()) if dec.any() else None
    I, w, f = rc.fit_law("hybrid5", pr)
    err = rc.score(I, rc.lcgrid().columns(pr))
    return dict(nonmono=bool(dec.any()), nonmono_upto_mu=where, bound=bool(f.active), lc=err)


def task(row):
    t, g, mh, band, e0 = row
    spec = phoenix.load(int(t), g, "-0.0" if mh == 0 else f"{mh:+.1f}")
    lam, S = bands.response(band)
    I = spec.band(lam, S)
    out = dict(teff=t, logg=g, mh=mh, band=band, lc_prod=e0)
    o = np.argsort(spec.mu)[::-1]
    r = np.sqrt(np.clip(1 - spec.mu[o] ** 2, 0, None))
    _, rc1 = limb_cut(spec.mu, I, "M1")
    k = int(np.argmin(np.abs(r - rc1)))
    variants = {"M1": None, "M2": "M2", "in1": float(r[k - 1]), "in2": float(r[k - 2])}
    for name, v in variants.items():
        if v is None:
            pr = band_profile(spec, lam, S, limb="M1")
        elif v == "M2":
            pr = band_profile(spec, lam, S, limb="M2")
        else:
            pr = band_profile(spec, lam, S, limb="M1", r_c=v)
        out[name] = analyse(pr)
    return out


def main():
    rows = rows_to_check()
    res = rc.run_parallel(task, rows)
    (rc.HERE / "c9_brightening.json").write_text(json.dumps(res))
    L = [f"## C9 — metal-poor limb brightening ({len(res)} PHOENIX hybrid5 rows with lc > 15 ppm)", ""]
    for v in ("M1", "M2", "in1", "in2"):
        nm = sum(r[v]["nonmono"] for r in res)
        bd = sum(r[v]["bound"] for r in res)
        lcs = np.array([r[v]["lc"] for r in res])
        mus = [r[v]["nonmono_upto_mu"] for r in res if r[v]["nonmono"]]
        L.append(f"* {v}: non-monotonic {nm}/{len(res)} (brightening at μ ≤ {max(mus) if mus else 0:.3f}), "
                 f"constraint binds {bd}/{len(res)}, hybrid5 lc median {np.median(lcs):.1f}, max {lcs.max():.1f} ppm")
    (rc.HERE / "c9_brightening.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
