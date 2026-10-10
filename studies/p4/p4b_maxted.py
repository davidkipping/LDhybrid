"""P4b — Kepler and TESS: Maxted (2023) empirical h'1 = I(2/3),
h'2 = I(2/3) − I(1/3) for 43 transiting-planet hosts, against predictions from
the v1.0 tables (direct profile evaluation, as Maxted compared tables).

Predictions: grids PHOENIX, MPS-ATLAS Set 1, Set 2; laws hybrid5 (headline)
and Claret 4-term (sanity check: our MPS-ATLAS offsets should reproduce
Maxted's, which used Kostogryz et al. 2022 coefficients); hybrid5 also with
the MURaM magnetic corrections (ssd, B100G, B200G) added.

Offsets: observed − predicted, error-weighted means with an external scatter
σ_ext fitted so that χ²/ν = 1 (Maxted's approach), for all stars and for
b < 0.65 (his Table 3 cut).

    .venv/bin/python studies/p4/p4b_maxted.py
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from ldhybrid import laws, predict
from ldhybrid.paths import RAW

HERE = Path(__file__).parent
DATA = RAW / "validation" / "maxted2023"
GRIDS = ("phoenix", "mps1", "mps2")
CASES = [("claret4", None), ("hybrid5", None), ("hybrid5", "ssd"), ("hybrid5", "B100G"),
         ("hybrid5", "B200G")]


def wmean(d, e):
    """Weighted mean with external scatter s: Σ (d − m)² / (e² + s²) = N − 1."""
    d, e = np.asarray(d), np.asarray(e)
    def chi(s):
        w = 1 / (e ** 2 + s ** 2)
        m = np.sum(w * d) / np.sum(w)
        return np.sum(w * (d - m) ** 2) - (len(d) - 1)
    s = 0.0 if chi(0.0) <= 0 else brentq(chi, 0.0, 1.0)
    w = 1 / (e ** 2 + s ** 2)
    return float(np.sum(w * d) / np.sum(w)), float(np.sqrt(1 / np.sum(w))), float(s)


def main():
    df = pd.read_csv(DATA / "maxted2023_results.csv", comment="#").dropna(subset=["Teff_K", "h1p"])
    rows = []
    for _, s in df.iterrows():
        band = s["instrument"]
        for grid in GRIDS:
            for law, state in CASES:
                w = predict.weights(grid, law, band, s["Teff_K"], s["logg_cgs"], s["feh_dex"])
                if np.any(np.isnan(w)):
                    continue
                if state:
                    dw = predict.magnetic(law, band, s["Teff_K"], state)
                    if dw is None:
                        continue
                    w = w + dw
                L = laws.get(law)
                h1, h2 = predict.h_prime(lambda mu: L.intensity(w, mu))
                rows.append(dict(star=s["star"], band=band, grid=grid, law=law,
                                 state=state or "none", teff=s["Teff_K"], b=s["b"],
                                 obs_h1=s["h1p"], obs_h2=s["h2p"], e_h1=s["h1p_err"],
                                 e_h2=s["h2p_err"], pred_h1=h1, pred_h2=h2))
    R = pd.DataFrame(rows)
    R["d1"], R["d2"] = R.obs_h1 - R.pred_h1, R.obs_h2 - R.pred_h2
    R.to_csv(HERE / "p4b_maxted_rows.csv", index=False)

    ref = pd.read_csv(DATA / "maxted2023_model_offsets.csv", comment="#")
    out, lines = [], ["# P4b: Kepler/TESS h'1, h'2 offsets (observed − predicted), Maxted (2023) sample",
                      "", "Weighted mean ± error (σ_ext), b < 0.65 as in Maxted's Table 3.", "",
                      "| band | grid | law | correction | N | Δh'1 | Δh'2 |", "|---|---|---|---|---|---|---|"]
    for band in ("Kepler", "TESS"):
        for grid in GRIDS:
            for law, state in CASES:
                S = R[(R.band == band) & (R.grid == grid) & (R.law == law) &
                      (R.state == (state or "none")) & (R.b < 0.65)]
                if len(S) < 3:
                    continue
                m1 = wmean(S.d1, S.e_h1)
                m2 = wmean(S.d2, S.e_h2)
                out.append(dict(band=band, grid=grid, law=law, state=state or "none", N=len(S),
                                dh1=m1, dh2=m2))
                lines.append(f"| {band} | {grid} | {law} | {state or '—'} | {len(S)} | "
                             f"{m1[0]:+.4f} ± {m1[1]:.4f} ({m1[2]:.4f}) | "
                             f"{m2[0]:+.4f} ± {m2[1]:.4f} ({m2[2]:.4f}) |")
    lines += ["", "Maxted (2023) Table 3 for reference (MPS-ATLAS, Kostogryz et al. 2022 coefficients):", ""]
    for _, r in ref[ref.model.str.contains("MPS|PHOENIX|Stagger")].iterrows():
        lines.append(f"* {r.band} {r.model} {r.notes if isinstance(r.notes, str) else ''}: "
                     f"Δh'1 = {r.mean_dh1p:+.3f} ± {r.mean_dh1p_err:.3f}, "
                     f"Δh'2 = {r.mean_dh2p:+.3f} ± {r.mean_dh2p_err:.3f} (N = {r.N})")
    (HERE / "p4b_maxted.md").write_text("\n".join(lines) + "\n")
    (HERE / "p4b_maxted.json").write_text(json.dumps(out, indent=1))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
