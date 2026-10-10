"""P4a — the Sun: Neckel & Labs (1994) continuum centre-to-limb profiles at
30 wavelengths (303–1099 nm), I(μ)/I(1) = Σ a_k μ^k.

(i) Can the laws represent the real Sun? Fit each law to each observed
    profile (v1 method) and score in light-curve space (production grid).
(ii) Do the models predict it? Model profiles at each wavelength (1-nm
    photon-weighted top-hat) for PHOENIX v3.0 (5700/5800 K, log g 4.5,
    interpolated to 5772 K), MPS-ATLAS Set 1/2 (log g 4.4), and MPS-ATLAS
    Set 2 corrected to MURaM G2 (3D hydro; + small-scale dynamo; + 100 G):
    I_1D × (I_MURaM / I_MPS2), the ratio held flat below μ = 0.1.
    Score: max |Δ I/I(1)| over μ ∈ [0.1, 1], and the light-curve difference
    (ppm) between model and observed profiles.

Caveat: NL measured line-free continuum windows; the models' 1-nm bands
include lines, which matters most in the blue/UV.

    .venv/bin/python studies/p4/p4a_sun.py
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

from ldhybrid import fit, laws
from ldhybrid import lightcurve as lc
from ldhybrid.grids import mpsatlas, muram, phoenix
from ldhybrid.paths import RAW
from ldhybrid.products import PROD_LC, constraint_for
from ldhybrid.profile import Profile, band_profile

HERE = Path(__file__).parent
G = lc.LCGrid(**PROD_LC)
MU = np.linspace(0.0, 1.0, 401)
LAWS = ("quadratic", "power2", "hybrid2", "claret4", "hybrid4", "hybrid5")


def tophat(lam_nm, width_nm=1.0):
    c = lam_nm * 10
    h = width_nm * 5
    return [c - h - 1e-3, c - h, c + h, c + h + 1e-3], [0, 1, 1, 0]


def law_fit_error(name, pr):
    if name == "power2":
        c, a = fit.fit_power2(pr, "r", 2000).w
        model = G.columns(lambda m: laws.power2(c, a, m))
    else:
        L = laws.get(name)
        w = fit.fit_linear(L, pr, "r", 2000, constraint_for(name)).w
        model = G.combine(*G.linear_columns(L), w) if name not in _C else G.combine(*_C[name], w)
    return G.errors(*model, *G.columns(pr))["max_ppm"]


_C = {n: G.linear_columns(laws.get(n)) for n in LAWS if n != "power2"}


def main():
    nl = pd.read_csv(RAW / "validation" / "neckel_labs" / "neckel_labs_1994_tableI.csv", comment="#")
    ph = {t: phoenix.load(t, 4.5, "-0.0") for t in (5700, 5800)}
    m1 = {t: mpsatlas.load("set1", 0.0, t, 4.4) for t in (5700, 5800)}
    m2 = {t: mpsatlas.load("set2", 0.0, t, 4.4) for t in (5700, 5800)}
    mu_states = {s: muram.load("G2_MH_00", s) for s in ("hydro", "ssd", "B100G")}
    f = (5772 - 5700) / 100.0
    rows = []
    for _, r in nl.iterrows():
        lam = r.wavelength_nm
        a = np.array([r[f"a{k}"] for k in range(6)])
        obs = Profile(MU, np.polynomial.polynomial.polyval(MU, a))
        # (i) representation
        rep = {n: law_fit_error(n, obs) for n in LAWS}
        # (ii) models, interpolated to 5772 K on a common μ grid
        th = tophat(lam)

        def grid_profile(models, spherical):
            ys = []
            for t in (5700, 5800):
                p = band_profile(models[t], *th, limb="M1") if spherical else band_profile(models[t], *th)
                ys.append(p.interpolant()(np.clip(MU, p.mu[0], 1)))
            return Profile(MU, (1 - f) * ys[0] + f * ys[1])

        model = {"phoenix": grid_profile(ph, True), "mps1": grid_profile(m1, False),
                 "mps2": grid_profile(m2, False)}
        mps2_558 = band_profile(m2[5800], *th)
        for s, spec in mu_states.items():
            pm = band_profile(spec, *th)
            ratio = PchipInterpolator(pm.mu, pm.y / mps2_558.interpolant()(pm.mu))
            C = ratio(np.clip(MU, pm.mu[0], 1.0))
            model[f"muram_{s}"] = Profile(MU, model["mps2"].y * C)
        sel = MU >= 0.1
        cmp = {k: dict(max_dI=float(np.max(np.abs(p.y[sel] - obs.y[sel]))),
                       lc_ppm=G.errors(*G.columns(p), *G.columns(obs))["max_ppm"])
               for k, p in model.items()}
        rows.append(dict(lam_nm=lam, representation_ppm=rep, models=cmp))
        print(f"{lam:7.1f} nm  hybrid5 fit {rep['hybrid5']:.2f} ppm | model LC ppm: " +
              " ".join(f"{k}={v['lc_ppm']:.0f}" for k, v in cmp.items()), flush=True)
    (HERE / "p4a_sun.json").write_text(json.dumps(rows, indent=1))
    lines = ["# P4a: the Sun (Neckel & Labs 1994)", "",
             "## (i) Representation: law fitted to the observed profile, worst LC error (ppm)", "",
             "| range | " + " | ".join(LAWS) + " |", "|---|" + "---|" * len(LAWS)]
    for lo, hi in ((300, 450), (450, 700), (700, 1100)):
        S = [x for x in rows if lo <= x["lam_nm"] < hi]
        lines.append(f"| {lo}–{hi} nm (N={len(S)}) | " + " | ".join(
            f"{np.median([x['representation_ppm'][n] for x in S]):.2f} / "
            f"{max(x['representation_ppm'][n] for x in S):.2f}" for n in LAWS) + " |")
    keys = list(rows[0]["models"])
    lines += ["", "(median / max over wavelengths)", "",
              "## (ii) Models vs observed: LC difference (ppm), median / max over wavelengths; "
              "max |ΔI/I(1)| (μ ≥ 0.1) median", "", "| range | " + " | ".join(keys) + " |",
              "|---|" + "---|" * len(keys)]
    for lo, hi in ((300, 450), (450, 700), (700, 1100)):
        S = [x for x in rows if lo <= x["lam_nm"] < hi]
        lines.append(f"| {lo}–{hi} nm | " + " | ".join(
            f"{np.median([x['models'][k]['lc_ppm'] for x in S]):.0f} / "
            f"{max(x['models'][k]['lc_ppm'] for x in S):.0f} "
            f"({np.median([x['models'][k]['max_dI'] for x in S]):.3f})" for k in keys) + " |")
    (HERE / "p4a_sun.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
