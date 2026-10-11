"""P4d — JWST: Sing et al. (2026) per-bin empirical (ℓ, δ) for six FGK hosts
(NIRISS/SOSS, NIRSpec/PRISM), extracted from their μ_min = 0.2 figures.

Definitions (Sing et al.): quadratic I/I(1) = 1 − c1(1−μ) − c2(1−μ)²,
ℓ = I(0) = 1 − c1 − c2, δ = c2/4; model values from a quadratic fitted to the
model intensities with μ ≥ μ_min = 0.2.

Cases, per wavelength bin (top-hat between bin midpoints):
  recipe       (headline) quadratic fitted directly to the grid's band profile
               (μ ≥ 0.2; native nodes for MPS-ATLAS, uniform μ for PHOENIX
               after the M1 cut and geometric rescaling): their method,
               checked against their per-star values;
  recipe+S     MPS-ATLAS Set 2 profile multiplied by the MURaM S/hydro
               centre-to-limb ratio (S = B100G, ssd; linear in Teff between
               the bracketing [M/H] = 0 MURaM stars; ratio held constant below
               μ = 0.1, as in scripts/p3_corrections.py), then the recipe;
  recipe_lin   PHOENIX with Sing et al.'s linear rescaling: μ_cri at max
               dI/dμ, μ' = (μ − μ_cri)/(1 − μ_cri), then the recipe;
  hybrid5      (reference) quadratic (μ ≥ 0.2, uniform μ) fitted to our
               hybrid5 law, itself fitted to the band profile (v1 method);
  hybrid5+B100G  the same on the B100G-corrected MPS-ATLAS Set 2 profile.
Model values: trilinear interpolation in (Teff, log g, [M/H]) between grid
models. Offsets observed − model; per star the inverse-variance weighted
mean over wavelength bins (weights 1/σ_ℓ² of the measured bins; correlations
between bins ignored), with error √(1/Σw); then the plain mean over stars,
the star-to-star s.d., the standard error of the mean (s.d./√n), the mean in
units of the SEM, and the median (Sing et al. Table 3 convention for the mean).
Per-star values for the figure: studies/rev2/code_sing_points.csv.

    .venv/bin/python studies/p4/p4d_sing.py
"""

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

from ldhybrid import fit, laws, v1
from ldhybrid.grids import mpsatlas, muram, phoenix
from ldhybrid.paths import RAW, REDUCED
from ldhybrid.products import constraint_for
from ldhybrid.profile import Profile, band_profile

HERE = Path(__file__).parent
D = RAW / "validation" / "sing2026"
MU_MIN = 0.2
MUU = np.linspace(MU_MIN, 1.0, 161)
L5 = laws.get("hybrid5")


def quad_ell_delta(mu, y):
    A = np.column_stack([1 - mu, (1 - mu) ** 2])
    c1, c2 = np.linalg.lstsq(A, 1 - y, rcond=None)[0]
    return 1 - c1 - c2, c2 / 4


def bins_for(lam_um):
    lam = np.sort(lam_um) * 1e4
    mid = 0.5 * (lam[1:] + lam[:-1])
    edges = np.concatenate([[lam[0] - (mid[0] - lam[0])], mid, [lam[-1] + (lam[-1] - mid[-1])]])
    return [([a - 1e-3, a, b, b + 1e-3], [0, 1, 1, 0]) for a, b in zip(edges[:-1], edges[1:])]


def corners(values, x):
    v = np.array(sorted(set(values)))
    if x <= v[0]:
        return [(v[0], 1.0)]
    if x >= v[-1]:
        return [(v[-1], 1.0)]
    i = np.searchsorted(v, x) - 1
    t = (x - v[i]) / (v[i + 1] - v[i])
    return [(v[i], 1 - t), (v[i + 1], t)]


def interp_models(catalogue, teff, logg, mh):
    """[(key, weight)] over the trilinear corners present in the catalogue."""
    cat = set(catalogue)
    out = []
    for (t, a), (g, b), (m, c) in itertools.product(
            corners([k[0] for k in catalogue], teff), corners([k[1] for k in catalogue], logg),
            corners([k[2] for k in catalogue], mh)):
        if (t, g, m) in cat and a * b * c > 0:
            out.append(((t, g, m), a * b * c))
    s = sum(w for _, w in out)
    return [(k, w / s) for k, w in out]


def stats(spec, th, spherical, muram_ratio=None, hybrid=True):
    """(ℓ, δ) by the recipe and via hybrid5 for one model and bin."""
    pr = band_profile(spec, *th, limb="M1") if spherical else band_profile(spec, *th)
    if muram_ratio is not None:
        pr = Profile(pr.mu, pr.y * muram_ratio(pr.mu), pr.meta)
    if spherical:
        rec = quad_ell_delta(MUU, pr.interpolant()(MUU))
    else:
        k = pr.mu >= MU_MIN - 1e-9
        rec = quad_ell_delta(pr.mu[k], pr.y[k])
    if not hybrid:
        return np.array(rec), None
    w = fit.fit_linear(L5, pr, v1.MEASURE, v1.RESAMPLE, constraint_for("hybrid5")).w
    hyb = quad_ell_delta(MUU, L5.intensity(w, MUU))
    return np.array(rec), np.array(hyb)


def ell_linear(mu_ph, I):
    """ℓ with Sing et al.'s linear rescaling: μ_cri at max dI/dμ of the raw
    profile, μ' = (μ − μ_cri)/(1 − μ_cri), quadratic uniform in μ' ≥ 0.2."""
    o = np.argsort(mu_ph)
    mu_ph, I = mu_ph[o], I[o]
    k = int(np.argmax(np.gradient(I, mu_ph)))
    mc = mu_ph[k]
    mu = (mu_ph[k:] - mc) / (1 - mc)
    return quad_ell_delta(MUU, np.interp(MUU, mu, I[k:] / I[-1]))


def muram_ratio_fn(teff, th, state="B100G", cache={}):
    """MURaM state/hydro ratio for a bin, linear in Teff between bracketing
    [M/H] = 0 stars, constant below μ = 0.1."""
    st = [(muram.params(s)["teff"]["hydro"], s) for s in muram.stars() if s.endswith("MH_00")]
    st.sort()
    T = np.array([t for t, _ in st])
    i = int(np.clip(np.searchsorted(T, teff) - 1, 0, len(T) - 2))
    t = float(np.clip((teff - T[i]) / (T[i + 1] - T[i]), 0, 1))
    fs = []
    for s in (st[i][1], st[i + 1][1]):
        if (s, state) not in cache:
            cache[(s, state)] = {k: muram.load(s, k) for k in ("hydro", state)}
        ph = band_profile(cache[(s, state)]["hydro"], *th)
        pb = band_profile(cache[(s, state)][state], *th)
        fs.append(PchipInterpolator(ph.mu, pb.y / ph.y))
    return lambda mu: (1 - t) * fs[0](np.clip(mu, 0.1, 1)) + t * fs[1](np.clip(mu, 0.1, 1))


def main():
    obs = pd.read_csv(D / "sing2026_lambda_bins.csv", comment="#")
    par = pd.read_csv(D / "sing2026_targets_wlc_fits.csv", comment="#")
    par["star"] = par.star_planet.str.replace(" b", "").str.replace(" A", "")
    ph_cat = [(e["teff"], e["logg"], float(e["mhs"])) for e in
              json.load(open(REDUCED / "phoenix" / "catalogue.json"))]
    mps_cat = [(t, g, mh) for mh, t, g in mpsatlas.catalogue("set2")]
    rows = []
    for star, O in obs.groupby("star"):
        p = par[par.star == star].iloc[0]
        O = O.sort_values("wavelength_um")
        ths = bins_for(O.wavelength_um.values)
        print(star, p.Teff_K, p.logg_cgs, p.feh_dex, len(O), flush=True)
        specs = {
            "phoenix": [(phoenix.load(t, g, "-0.0" if m == 0 else f"{m:+.1f}"), w, True)
                        for (t, g, m), w in interp_models(ph_cat, p.Teff_K, p.logg_cgs, p.feh_dex)],
            "mps1": [(mpsatlas.load("set1", m, t, g), w, False)
                     for (t, g, m), w in interp_models(mps_cat, p.Teff_K, p.logg_cgs, p.feh_dex)],
            "mps2": [(mpsatlas.load("set2", m, t, g), w, False)
                     for (t, g, m), w in interp_models(mps_cat, p.Teff_K, p.logg_cgs, p.feh_dex)],
        }
        for j, (th, (_, o)) in enumerate(zip(ths, O.iterrows())):
            rec = {}
            for grid, lst in specs.items():
                if min(s.coverage(*th) for s, _, _ in lst) < 0.999:
                    continue
                r = [stats(s, th, sph) for s, w, sph in lst]
                rec[f"{grid}|recipe"] = sum(w * x[0] for (_, w, _), x in zip(lst, r))
                rec[f"{grid}|hybrid5"] = sum(w * x[1] for (_, w, _), x in zip(lst, r))
                if grid == "phoenix":
                    rec["phoenix|recipe_lin"] = sum(
                        w * np.array(ell_linear(s.mu, s.band(*th))) for s, w, _ in lst)
                if grid == "mps2" and th[0][1] >= 2000.0 and th[0][2] <= 99800.0:
                    for state in ("B100G", "ssd"):
                        C = muram_ratio_fn(p.Teff_K, th, state)
                        r = [stats(s, th, sph, C, hybrid=state == "B100G") for s, w, sph in lst]
                        rec[f"mps2|recipe+{state}"] = sum(w * x[0] for (_, w, _), x in zip(lst, r))
                        if state == "B100G":
                            rec["mps2|hybrid5+B100G"] = sum(w * x[1] for (_, w, _), x in zip(lst, r))
            for k, (l, d) in rec.items():
                rows.append(dict(star=star, teff=p.Teff_K, lam=o.wavelength_um, case=k,
                                 ell_obs=o.ell, ell_err=o.ell_err, delta_obs=o.delta,
                                 delta_err=o.delta_err, ell_mod=l, delta_mod=d))
    R = pd.DataFrame(rows)
    R["dl"], R["dd"] = R.ell_obs - R.ell_mod, R.delta_obs - R.delta_mod
    R.to_csv(HERE / "p4d_sing_rows.csv", index=False)

    ref = pd.read_csv(D / "sing2026_model_offsets.csv", comment="#")
    ref = ref[ref.mu_min == 0.2].set_index("model")
    lines = ["# P4d: JWST per-bin (ℓ, δ), Sing et al. (2026), μ_min = 0.2", "",
             "Per star: inverse-variance (1/σ_ℓ² per bin) weighted mean of observed − model "
             "over wavelength bins, error √(1/Σw) (bin correlations ignored). Over stars: "
             "plain mean, s.d., SEM = s.d./√n, mean/SEM, median. Headline route: recipe "
             "(quadratic fitted to the model profile at μ ≥ 0.2).", "",
             "| case | n | mean Δℓ | s.d. | SEM | mean/SEM | median | Δδ (mean ± s.d.) | per star Δℓ |",
             "|---|---|---|---|---|---|---|---|---|"]
    summary = {}
    points = []
    order = ["mps2|recipe", "mps2|recipe+B100G", "mps2|recipe+ssd", "mps1|recipe",
             "phoenix|recipe", "phoenix|recipe_lin", "mps2|hybrid5", "mps2|hybrid5+B100G",
             "mps1|hybrid5", "phoenix|hybrid5"]
    groups = dict(list(R.groupby("case")))
    for case in [c for c in order if c in groups] + [c for c in groups if c not in order]:
        S = groups[case]
        per = []
        for star, T in S.groupby("star"):
            wl, wd = 1 / T.ell_err ** 2, 1 / T.delta_err ** 2
            per.append((star, float(np.sum(wl * T.dl) / np.sum(wl)),
                        float(np.sum(wd * T.dd) / np.sum(wd)), float(np.sqrt(1 / np.sum(wl))),
                        float(T.teff.iloc[0]), len(T)))
            points.append(dict(case=case, star=star, teff=float(T.teff.iloc[0]), dl=per[-1][1],
                               dl_err=per[-1][3], n_bins=len(T)))
        dl = np.array([x[1] for x in per])
        dd = np.array([x[2] for x in per])
        sd = dl.std(ddof=1)
        sem = sd / np.sqrt(len(dl))
        summary[case] = dict(dl=float(dl.mean()), dl_sd=float(sd), dl_sem=float(sem),
                             signif=float(dl.mean() / sem), dl_median=float(np.median(dl)),
                             dd=float(dd.mean()), dd_sd=float(dd.std(ddof=1)), per_star=per)
        lines.append(f"| {case} | {len(per)} | {dl.mean():+.3f} | {sd:.3f} | {sem:.3f} | "
                     f"{dl.mean() / sem:.1f} | {np.median(dl):+.3f} | "
                     f"{dd.mean():+.3f} ± {dd.std(ddof=1):.3f} | " +
                     ", ".join(f"{x[0]} {x[1]:+.3f}" for x in per) + " |")
    out = HERE.parent / "rev2"
    out.mkdir(exist_ok=True)
    pd.DataFrame(points).to_csv(out / "code_sing_points.csv", index=False)
    lines += ["", "Sing et al. (2026) Table 3, μ_min = 0.2 (7 stars incl. WASP-12):", ""]
    for m in ("PHOENIX", "MPS1", "MPS2", "STAGGER"):
        r = ref.loc[m]
        lines.append(f"* {m}: Δℓ = {r.dl:+.3f} ± {r.dl_sig:.3f}, Δδ = {r.ddelta:+.3f} ± {r.ddelta_sig:.3f}")
    (HERE / "p4d_sing.md").write_text("\n".join(lines) + "\n")
    (HERE / "p4d_sing.json").write_text(json.dumps(summary, indent=1, default=str))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
