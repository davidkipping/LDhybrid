"""C1: light-curve effect of the limb-rescaling convention for spherical models.

Same PHOENIX band intensities and the same edge, two rescalings:
  geom  production M1: r_c at the node of max |dI/dr|, μ' = √(1 − (r/r_c)²);
  lin   Sing et al. (2026): μ_cri at max dI/dμ (the same node),
        μ' = (μ − μ_cri)/(1 − μ_cri);
  lin4  as lin, after removing the four nodes inward of μ_cri.
Both profiles use the same PCHIP interpolant (Profile). Score: worst |ΔF|
between their light curves on the production grid (p = 0.1,
b ∈ {0, 0.3, 0.6, 0.9}); T3 and the PHOENIX part of P2S.

    .venv/bin/python studies/rev2/c1_rescale.py
"""
import json

import numpy as np

import rc
from ldhybrid import bands
from ldhybrid.profile import Profile, band_profile


def lin_profiles(spec, band):
    lam, S = bands.response(band)
    I = spec.band(lam, S)
    o = np.argsort(spec.mu)
    mu_ph, I = spec.mu[o], I[o]
    d = np.gradient(I, mu_ph)
    k = int(np.argmax(d))
    mc = mu_ph[k]
    out = []
    for drop in (0, 4):
        keep = np.r_[k, np.arange(k + 1 + drop, len(mu_ph))]
        out.append(Profile((mu_ph[keep] - mc) / (1 - mc), I[keep]))
    return out, float(mc)


def one(spec, band):
    lam, S = bands.response(band)
    if spec.coverage(lam, S) < 0.999:
        return None
    geo = band_profile(spec, lam, S, limb="M1")
    (lin, lin4), mc = lin_profiles(spec, band)
    g = rc.lcgrid()
    tg = g.columns(geo)
    mu_geo = float(np.sqrt(1 - geo.meta["r_c"] ** 2))
    return dict(band=band, mu_cri_geo=mu_geo, mu_cri_lin=mc,
                lin=g.errors(*g.columns(lin), *tg)["max_ppm"],
                lin4=g.errors(*g.columns(lin4), *tg)["max_ppm"])


def task(job):
    kind, key, bandlist, label = job
    spec = rc.load_t3(*key) if kind == "t3" else rc.load(kind, key)
    rows = []
    for b in bandlist:
        r = one(spec, b)
        if r:
            r.update(case=label, kind=kind)
            rows.append(r)
    return rows


def main():
    jobs = [("t3", key, [b for b, _ in rc.T3_BANDS], name) for name, key in rc.T3_STARS]
    jobs += [("phoenix", key, rc.P2S_BANDS, str(key)) for kind, key in rc.p2s_models() if kind == "phoenix"]
    rows = [r for res in rc.run_parallel(task, jobs) for r in res]
    (rc.HERE / "c1_rescale.json").write_text(json.dumps(rows))
    L = ["## C1 — limb rescaling (geometric vs Sing et al. linear), light-curve difference (ppm)", ""]
    for kind, lab in (("t3", "T3 (15 cases)"), ("phoenix", "P2S PHOENIX (272)")):
        sel = [r for r in rows if r["kind"] == kind]
        for k in ("lin", "lin4"):
            x = np.array([r[k] for r in sel])
            L.append(f"* {lab}, {k}: median {np.median(x):.0f}, 90th {np.percentile(x, 90):.0f}, "
                     f"min {x.min():.0f}, max {x.max():.0f}")
        dm = np.array([abs(r["mu_cri_geo"] - r["mu_cri_lin"]) for r in sel])
        L.append(f"  same edge node: max |Δμ_cri| = {dm.max():.2g}")
    L.append("")
    L.append("T3 per case (lin): " + "; ".join(f"{r['case']} {r['band']} {r['lin']:.0f}"
                                                for r in rows if r["kind"] == "t3"))
    (rc.HERE / "c1_rescale.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
