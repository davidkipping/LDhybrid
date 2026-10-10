"""P4c — TESS: Patel & Espinoza (2022) empirical quadratic u1, u2 for 176
systems, against SPAM predictions from the v1.0 tables: the quadratic law a
light-curve fit would return for the hybrid5 profile at each system's k and b
(Howarth 2011), like for like with their light-curve fits.

Cases: PHOENIX, MPS-ATLAS Set 1/2; hybrid5 with and without the MURaM B100G
correction. Offsets are reported as predicted − empirical (their sign
convention, Table 4), error-weighted with a fitted external scatter.
Also: our PHOENIX SPAM vs their C17 PHOENIX r-method SPAM values.

    .venv/bin/python studies/p4/p4c_patel.py
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from p4b_maxted import wmean  # noqa: E402

from ldhybrid import laws, predict
from ldhybrid.paths import RAW

HERE = Path(__file__).parent
DATA = RAW / "validation" / "patel_espinoza2022"
GRIDS = ("phoenix", "mps1", "mps2")


def main():
    df = pd.read_csv(DATA / "patel_espinoza2022_tess_ldcs.csv", comment="#")
    df = df.dropna(subset=["Teff_K", "logg_cgs", "u1", "u2", "k", "b"])
    df["feh"] = df["feh_dex"].fillna(0.0)
    L = laws.get("hybrid5")
    rows = []
    for _, s in df.iterrows():
        for grid in GRIDS:
            w0 = predict.weights(grid, "hybrid5", "TESS", s.Teff_K, s.logg_cgs, s.feh)
            if np.any(np.isnan(w0)):
                continue
            for state in ("none", "B100G"):
                w = w0 if state == "none" else w0 + predict.magnetic("hybrid5", "TESS", s.Teff_K, state)
                u = predict.spam(lambda mu: L.intensity(w, mu), float(s.k), float(min(s.b, 0.95)))
                rows.append(dict(planet=s.planet, grid=grid, state=state, teff=s.Teff_K,
                                 feh_missing=bool(np.isnan(s.feh_dex)), b=s.b,
                                 u1=s.u1, u2=s.u2, e1=0.5 * (s.u1_errp + s.u1_errm),
                                 e2=0.5 * (s.u2_errp + s.u2_errm), p1=u[0], p2=u[1],
                                 c17r_1=s.u1_C17_PHOENIX_r_SPAM, c17r_2=s.u2_C17_PHOENIX_r_SPAM))
    R = pd.DataFrame(rows)
    R["d1"], R["d2"] = R.p1 - R.u1, R.p2 - R.u2
    R.to_csv(HERE / "p4c_patel_rows.csv", index=False)
    lines = ["# P4c: TESS quadratic coefficients, Patel & Espinoza (2022) sample",
             "", "SPAM predictions from hybrid5 profiles at each system's (k, b); "
             "offsets predicted − empirical (their convention), weighted mean ± error (σ_ext).", "",
             "| grid | correction | N | Δu1 | Δu2 |", "|---|---|---|---|---|"]
    out = []
    for grid in GRIDS:
        for state in ("none", "B100G"):
            S = R[(R.grid == grid) & (R.state == state)]
            m1, m2 = wmean(S.d1, S.e1), wmean(S.d2, S.e2)
            out.append(dict(grid=grid, state=state, N=len(S), du1=m1, du2=m2))
            lines.append(f"| {grid} | {state} | {len(S)} | {m1[0]:+.3f} ± {m1[1]:.3f} ({m1[2]:.3f}) | "
                         f"{m2[0]:+.3f} ± {m2[1]:.3f} ({m2[2]:.3f}) |")
    S = R[(R.grid == "phoenix") & (R.state == "none")].dropna(subset=["c17r_1"])
    lines += ["", f"Our PHOENIX-ACES v3.0 SPAM minus their C17 PHOENIX r-method SPAM "
              f"({len(S)} systems): median Δu1 = {np.median(S.p1 - S.c17r_1):+.3f}, "
              f"Δu2 = {np.median(S.p2 - S.c17r_2):+.3f}.",
              "", "Patel & Espinoza (2022) Table 4, SPAM, for reference: C17 PHOENIX r-method "
              "Δu1 = −0.020 ± 0.018, Δu2 = +0.042 ± 0.022; C17 ATLAS −0.084, +0.093."]
    for lo, hi in ((2500, 4000), (4000, 5000), (5000, 6000), (6000, 7500)):
        S = R[(R.grid == "phoenix") & (R.state == "none") & (R.teff >= lo) & (R.teff < hi)]
        if len(S) >= 3:
            lines.append(f"* PHOENIX, {lo}–{hi} K (N = {len(S)}): Δu1 = {wmean(S.d1, S.e1)[0]:+.3f}, "
                         f"Δu2 = {wmean(S.d2, S.e2)[0]:+.3f}")
    (HERE / "p4c_patel.md").write_text("\n".join(lines) + "\n")
    (HERE / "p4c_patel.json").write_text(json.dumps(out, indent=1))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
