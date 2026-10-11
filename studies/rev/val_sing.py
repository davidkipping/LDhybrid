"""Revision D11: the JWST ℓ comparison with Sing et al. (2026), made explicit.

1. Sing et al.'s own per-star Δℓ (PHOENIX, MPS1, Stagger; μ_min = 0, 0.1,
   0.2) recovered by exact vector extraction from their offset-vs-μ_min
   figures. These include WASP-12, whose per-wavelength ℓ spectrum is not
   published (no figure), so we can show how much its omission matters.
2. PHOENIX ℓ under three limb rescalings, per wavelength bin, at μ_min = 0.2:
     geom     our M1 cut (node of max |dI/dr|), μ' = √(1 − (r/r_c)²);
     lin      Sing et al.: μ_cri at max dI/dμ_ph, μ = (μ_ph − μ_cri)/(1 − μ_cri),
              intensities interpolated linearly on a uniform μ grid;
     lin4     as lin, after removing the four nodes inward of μ_cri.
   ℓ = 1 − (c1 + c2) of the quadratic least-squares fit, uniform in μ over
   μ ≥ 0.2 (161 points), as in p4d_sing.py.
3. Statistics of the per-star offsets from p4d_sing_rows.csv (and from 2):
   mean, star-to-star s.d., standard error of the mean, mean/SEM, median.
4. Kepler/TESS (Maxted 2023, b < 0.65): sample sizes and weighted means.

    .venv/bin/python studies/rev/val_sing.py
"""

import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pymupdf

from ldhybrid.grids import phoenix
from ldhybrid.paths import RAW, REDUCED
from ldhybrid.profile import limb_cut

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "studies" / "p4"))
import p4d_sing as P4D  # noqa: E402
import p4b_maxted as P4B  # noqa: E402

D = RAW / "validation" / "sing2026"
FIG = str(D / "raw" / "figures" / "Limb_Darkening_{}_offset_vs_mumin_raw_and_offset_chi2.pdf")
MUU = P4D.MUU


# ------------------------------------------------------------ 1. extraction
def extract(model):
    page = pymupdf.open(FIG.format(model))[0]
    words = page.get_text("words")
    # Δℓ panel tick labels (x ≈ 28, y < 168): minus signs are paths, so the
    # labels below "0.00" are negative
    lab = sorted([(0.5 * (w[1] + w[3]), float(w[4])) for w in words
                  if abs(w[0] - 28.1) < 1 and 0.5 * (w[1] + w[3]) < 168 and w[4][0].isdigit()])
    y0 = [y for y, v in lab if v == 0][0]
    vals = np.array([v if y <= y0 else -v for y, v in lab])
    ys = np.array([y for y, _ in lab])
    a, b = np.polyfit(ys, vals, 1)
    names = [(0.5 * (w[1] + w[3]), w[0], w[4]) for w in words if w[4].startswith(("WASP", "Hat"))]
    dr = page.get_drawings()
    # legend handles: short horizontal coloured segments (20 pt) just left of a name
    colour = {}
    for x in dr:
        it = x["items"]
        if len(it) == 1 and it[0][0] == "l" and x.get("color") and x["color"] != (0, 0, 0):
            p, q = it[0][1], it[0][2]
            if abs(p.y - q.y) < 0.01 and abs(abs(q.x - p.x) - 20) < 0.5:
                y, xr = p.y, max(p.x, q.x)
                c = [n for yy, xx, n in names if abs(yy - y) < 3 and 0 < xx - xr < 20]
                if c:
                    colour[tuple(round(v, 4) for v in x["color"])] = c[0].replace("Hat", "HAT")
    out = {}
    for x in dr:
        it = x["items"]
        col = tuple(round(v, 4) for v in x["color"]) if x.get("color") else None
        if col in colour and len(it) == 2 and all(i[0] == "l" for i in it):
            pts = [it[0][1], it[0][2], it[1][2]]
            if max(p.y for p in pts) < 168:                     # top panel (Δℓ)
                out[colour[col]] = [a * p.y + b for p in pts]   # μ_min = 0, 0.1, 0.2
    return out


# ------------------------------------------------------------ 2. rescalings
def quad_ell(mu, y):
    return P4D.quad_ell_delta(mu, y)[0]


def ell_three(mu_ph, I):
    o = np.argsort(mu_ph)
    mu_ph, I = mu_ph[o], I[o]
    pr, _ = limb_cut(mu_ph, I, "M1")
    geom = quad_ell(MUU, pr.interpolant()(MUU))
    d = np.gradient(I, mu_ph)
    k = int(np.argmax(d))
    mc = mu_ph[k]
    res = []
    for drop in (0, 4):
        keep = np.r_[k, np.arange(k + 1 + drop, len(mu_ph))]
        mu = (mu_ph[keep] - mc) / (1 - mc)
        y = I[keep] / I[-1]
        res.append(quad_ell(MUU, np.interp(MUU, mu, y)))
    return geom, res[0], res[1], float(mc)


def star_job(args):
    star, teff, logg, feh, lams = args
    cat = [(e["teff"], e["logg"], float(e["mhs"])) for e in
           json.load(open(REDUCED / "phoenix" / "catalogue.json"))]
    corners = P4D.interp_models(cat, teff, logg, feh)
    specs = [(phoenix.load(t, g, "-0.0" if m == 0 else f"{m:+.1f}"), w) for (t, g, m), w in corners]
    rows = []
    for lam, th in zip(lams, P4D.bins_for(np.array(lams))):
        if min(s.coverage(*th) for s, _ in specs) < 0.999:
            continue
        acc = np.zeros(3)
        mcs = []
        for s, w in specs:
            I = s.band(*th)
            g_, l_, l4_, mc = ell_three(s.mu, I)
            acc += w * np.array([g_, l_, l4_])
            _, rc = limb_cut(*[a[np.argsort(s.mu)] for a in (s.mu, I)], "M1")
            mcs.append((mc, float(np.sqrt(1 - rc ** 2))))
        rows.append(dict(star=star, lam=lam, geom=acc[0], lin=acc[1], lin4=acc[2],
                         mucri_lin=float(np.mean([m[0] for m in mcs])),
                         mucri_M1=float(np.mean([m[1] for m in mcs]))))
    return rows


# ------------------------------------------------------------ 3. statistics
def per_star(R, col_dl="dl"):
    out = []
    for star, T in R.groupby("star"):
        w = 1 / T.ell_err ** 2
        out.append((star, float(T.teff.iloc[0]), float(np.sum(w * T[col_dl]) / np.sum(w)),
                    float(np.sqrt(1 / np.sum(w))), len(T)))
    return out


def summary(per):
    x = np.array([p[2] for p in per])
    sd = x.std(ddof=1)
    sem = sd / np.sqrt(len(x))
    return dict(n=len(x), mean=float(x.mean()), sd=float(sd), sem=float(sem),
                signif=float(x.mean() / sem), median=float(np.median(x)))


def main():
    L = ["# D11: JWST ℓ comparison with Sing et al. (2026), μ_min = 0.2", ""]
    # 1. their per-star values
    ext = {m: extract(m) for m in ("phoenix", "mps1", "stagger")}
    ref = pd.read_csv(D / "sing2026_model_offsets.csv", comment="#")
    ref = ref[ref.mu_min == 0.2].set_index("model")
    L += ["## 1. Sing et al.'s per-star Δℓ (vector-extracted from their offset-vs-μ_min figures)", "",
          "| model | " + " | ".join(sorted(ext["mps1"])) + " | mean 7 | mean w/o WASP-12 | Table 3 |",
          "|---|" + "---|" * (len(ext["mps1"]) + 3)]
    for m, e in ext.items():
        v = {k: x[2] for k, x in e.items()}
        a7 = np.mean(list(v.values()))
        a6 = np.mean([x for k, x in v.items() if k != "WASP-12"])
        L.append(f"| {m} | " + " | ".join(f"{v[k]:+.3f}" for k in sorted(v)) +
                 f" | {a7:+.3f} | {a6:+.3f} | {ref.loc[m.upper(), 'dl']:+.3f} |")
    L += ["", "(Table 3 means differ from the plain mean of the plotted per-star values; "
          "their averaging is not stated.)", ""]

    # 2. rescalings for PHOENIX
    obs = pd.read_csv(D / "sing2026_lambda_bins.csv", comment="#")
    par = pd.read_csv(D / "sing2026_targets_wlc_fits.csv", comment="#")
    par["star"] = par.star_planet.str.replace(" b", "").str.replace(" A", "")
    jobs = []
    for star, O in obs.groupby("star"):
        p = par[par.star == star].iloc[0]
        jobs.append((star, p.Teff_K, p.logg_cgs, p.feh_dex, sorted(O.wavelength_um.values)))
    rows = []
    with ProcessPoolExecutor(4) as ex:
        for r in ex.map(star_job, jobs):
            rows += r
    Rr = pd.DataFrame(rows).merge(obs[["star", "wavelength_um", "ell", "ell_err"]],
                                  left_on=["star", "lam"], right_on=["star", "wavelength_um"])
    Rr["teff"] = Rr.star.map(dict(zip(par.star, par.Teff_K)))
    L += ["## 2. PHOENIX ℓ under three limb rescalings (μ_min = 0.2)", "",
          "| star | μ_cri M1 (ours) | μ_cri max dI/dμ (Sing) | Δℓ geom | Δℓ lin | Δℓ lin4 | lin − geom |",
          "|---|---|---|---|---|---|---|"]
    stat2 = {}
    for case in ("geom", "lin", "lin4"):
        Rr["dl_" + case] = Rr.ell - Rr[case]
        stat2[case] = per_star(Rr, "dl_" + case)
    for i, (star, *_r) in enumerate(stat2["geom"]):
        T = Rr[Rr.star == star]
        g_, l_, l4_ = (stat2[c][i][2] for c in ("geom", "lin", "lin4"))
        L.append(f"| {star} | {T.mucri_M1.median():.4f} | {T.mucri_lin.median():.4f} | {g_:+.3f} | "
                 f"{l_:+.3f} | {l4_:+.3f} | {l_ - g_:+.3f} |")
    for c in ("geom", "lin", "lin4"):
        s = summary(stat2[c])
        L.append(f"* {c}: mean Δℓ {s['mean']:+.3f} ± {s['sd']:.3f} (s.d.), SEM {s['sem']:.3f}")
    shift = np.mean([stat2["lin"][i][2] - stat2["geom"][i][2] for i in range(len(stat2["geom"]))])
    L += [f"* mean shift lin − geom: {shift:+.3f}; lin4 − lin: "
          f"{np.mean([stat2['lin4'][i][2] - stat2['lin'][i][2] for i in range(len(stat2['lin']))]):+.4f}", ""]

    # 3. statistics of the per-star offsets
    R = pd.read_csv(ROOT / "studies" / "p4" / "p4d_sing_rows.csv")
    pts = []
    L += ["## 3. Per-star offset statistics (observed − model)", "",
          "| case | n | mean | s.d. | SEM | mean/SEM | median |", "|---|---|---|---|---|---|---|"]
    allstats = {}
    for case, S in R.groupby("case"):
        per = per_star(S)
        s = summary(per)
        allstats[case] = dict(s, per_star=per)
        L.append(f"| {case} | {s['n']} | {s['mean']:+.3f} | {s['sd']:.3f} | {s['sem']:.3f} | "
                 f"{s['signif']:.1f} | {s['median']:+.3f} |")
        pts += [dict(case=case, star=st, teff=t, dl=d, dl_err=e, n_bins=n) for st, t, d, e, n in per]
    for c in ("geom", "lin", "lin4"):
        s = summary(stat2[c])
        allstats[f"phoenix|recipe_{c}"] = dict(s, per_star=stat2[c])
        L.append(f"| phoenix recipe, {c} rescaling | {s['n']} | {s['mean']:+.3f} | {s['sd']:.3f} | "
                 f"{s['sem']:.3f} | {s['signif']:.1f} | {s['median']:+.3f} |")
        pts += [dict(case=f"phoenix|recipe_{c}", star=st, teff=t, dl=d, dl_err=e, n_bins=n)
                for st, t, d, e, n in stat2[c]]
    pd.DataFrame(pts).to_csv(HERE / "val_sing_points.csv", index=False)
    L += ["", "dl_err in val_sing_points.csv is the statistical error of the per-star weighted "
          "mean (1/√Σσ⁻²), from the extracted per-bin ℓ errors; it ignores correlations between "
          "bins and is far smaller than the star-to-star scatter.", ""]

    # our recipe vs their per-star values
    L += ["Our recipe vs their plotted per-star Δℓ (μ_min = 0.2):", ""]
    for ours, theirs in (("mps1|recipe", "mps1"), ("phoenix|recipe", "phoenix")):
        per = {p[0]: p[2] for p in allstats[ours]["per_star"]}
        d = [per[k] - ext[theirs][k][2] for k in per if k in ext[theirs]]
        L.append(f"* {ours} − Sing {theirs}: " + ", ".join(
            f"{k} {per[k] - ext[theirs][k][2]:+.3f}" for k in sorted(per) if k in ext[theirs]) +
            f"; mean {np.mean(d):+.3f}")
    per = {p[0]: p[2] for p in allstats["phoenix|recipe_lin4"]["per_star"]}
    d = [per[k] - ext["phoenix"][k][2] for k in per]
    L.append(f"* phoenix recipe with Sing's lin4 rescaling − Sing phoenix: mean {np.mean(d):+.3f}, "
             f"max |Δ| {np.max(np.abs(d)):.3f}")

    # 4. Maxted
    M = pd.read_csv(ROOT / "studies" / "p4" / "p4b_maxted_rows.csv")
    M = M[M.b.abs() < 0.65]
    mp = []
    L += ["", "## 4. Kepler/TESS (Maxted 2023), b < 0.65", ""]
    for (band, grid, law, state), T in M.groupby(["band", "grid", "law", "state"]):
        m1, e1, s1 = P4B.wmean(T.d1, T.e_h1)
        m2, e2, s2 = P4B.wmean(T.d2, T.e_h2)
        mp.append(dict(band=band, grid=grid, law=law, state=state, n=len(T), dh1=m1, dh1_err=e1,
                       dh1_ext=s1, dh2=m2, dh2_err=e2, dh2_ext=s2))
    MP = pd.DataFrame(mp)
    MP.to_csv(HERE / "val_maxted_points.csv", index=False)
    n = M.drop_duplicates(["star", "band"]).groupby("band").size()
    L.append(f"* hosts with b < 0.65: Kepler {n.get('Kepler', 0)}, TESS {n.get('TESS', 0)} "
             f"(of {M.star.nunique()} stars; 43 in Maxted's Table 2)")
    for _, r in MP[(MP.grid == "mps2") & (MP.law == "hybrid5")].iterrows():
        L.append(f"* {r.band} mps2 hybrid5 {r.state}: Δh'1 {r.dh1:+.4f} ± {r.dh1_err:.4f} "
                 f"({abs(r.dh1) / r.dh1_err:.1f}σ), Δh'2 {r.dh2:+.4f} ± {r.dh2_err:.4f} "
                 f"({abs(r.dh2) / r.dh2_err:.1f}σ); σ_ext {r.dh1_ext:.4f}, {r.dh2_ext:.4f}")
    (HERE / "val_sing.md").write_text("\n".join(L) + "\n")
    (HERE / "val_sing.json").write_text(json.dumps(dict(extracted=ext, stats=allstats), indent=1, default=str))
    print("\n".join(L))


if __name__ == "__main__":
    main()
