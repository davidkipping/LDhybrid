"""Rebuild data/products/empirical_offsets.csv (CDS offsets.dat) with the
correction state, the route of the model value, the passbands, and the
correlation of the per-star h'1 and h'2 offsets.

  Kepler/TESS (Maxted 2023, b < 0.65): hybrid5 predictions (studies/p4/p4b);
      offset = weighted mean ± error, scatter = fitted external scatter,
      rho = correlation of the per-star (h'1, h'2) offsets.
  JWST (Sing et al. 2026, μ_min = 0.2): per-star inverse-variance means over
      bins (studies/p4/p4d); offset = mean over stars, offset_err = SEM,
      scatter = star-to-star s.d. Two routes: hybrid5 (statistic of the
      tabulated hybrid5 fit; use to condition weights) and profile (statistic
      of the model profile; the like-for-like comparison with Sing et al.).

    .venv/bin/python studies/rev2/code_offsets.py
"""

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "products" / "empirical_offsets.csv"
P4 = ROOT / "studies" / "p4"


def main():
    old = [r for r in csv.DictReader(l for l in open(OUT) if not l.startswith("#"))]
    rows_b = pd.read_csv(P4 / "p4b_maxted_rows.csv")
    rows_b = rows_b[(rows_b.law == "hybrid5") & (rows_b.b < 0.65)]
    state = {"none": "none", "MURaM_B100G": "B100G"}
    out = []
    for r in old:
        if r["regime"] not in ("Kepler", "TESS"):
            continue
        T = rows_b[(rows_b.band == r["regime"]) & (rows_b.grid == r["grid"]) &
                   (rows_b.state == state[r["correction"]])]
        rho = float(np.corrcoef(T.d1, T.d2)[0, 1]) if len(T) > 2 else ""
        out.append(dict(regime=r["regime"], bands=r["regime"], grid=r["grid"],
                        correction=r["correction"], route="hybrid5", statistic=r["statistic"],
                        offset=r["offset"], offset_err=r["offset_err"], scatter=r["scatter"],
                        rho_h1h2=f"{rho:+.3f}" if rho != "" else "", N=r["N"],
                        source="Maxted 2023 (b<0.65)"))
    sj = json.loads((P4 / "p4d_sing.json").read_text())
    cases = [("mps2|hybrid5", "mps2", "none", "hybrid5"),
             ("mps2|hybrid5+B100G", "mps2", "MURaM_B100G", "hybrid5"),
             ("mps1|hybrid5", "mps1", "none", "hybrid5"),
             ("phoenix|hybrid5", "phoenix", "none", "hybrid5"),
             ("mps2|recipe", "mps2", "none", "profile"),
             ("mps2|recipe+B100G", "mps2", "MURaM_B100G", "profile"),
             ("mps2|recipe+ssd", "mps2", "MURaM_ssd", "profile"),
             ("mps1|recipe", "mps1", "none", "profile"),
             ("phoenix|recipe", "phoenix", "none", "profile")]
    for case, grid, corr, route in cases:
        s = sj[case]
        n = len(s["per_star"])
        for stat, mean, sd in (("ell", s["dl"], s["dl_sd"]), ("delta", s["dd"], s["dd_sd"])):
            out.append(dict(regime="JWST_0.86-5.5um", bands="NIRISS SOSS, NIRSpec PRISM",
                            grid=grid, correction=corr, route=route, statistic=stat,
                            offset=f"{mean:+.4f}", offset_err=f"{sd / np.sqrt(n):.4f}",
                            scatter=f"{sd:.4f}", rho_h1h2="", N=n,
                            source="Sing et al. 2026 (per-star)"))
    head = ["# LDhybrid v1.0 empirical offsets (observed - model); rebuilt by "
            "studies/rev2/code_offsets.py",
            "# h1p = I(2/3)/I(1), h2p = h1p - I(1/3)/I(1) (Maxted 2023); ell = 1-(c1+c2), "
            "delta = c2/4 of a quadratic fitted uniformly in mu >= 0.2 (Sing et al. 2026)",
            "# route: hybrid5 = statistic of the tabulated hybrid5 fit (condition weights "
            "with these); profile = statistic of the model profile (comparison with Sing et al.)",
            "# correction: the state of the predictions; apply only to a prior in the same "
            "state. scatter = star-to-star (1-sigma prior width for one star); rho = "
            "correlation of the per-star h1p, h2p offsets"]
    cols = list(out[0])
    with open(OUT, "w") as f:
        f.write("\n".join(head) + "\n")
        w = csv.DictWriter(f, cols)
        w.writeheader()
        w.writerows(out)
    print(open(OUT).read())


if __name__ == "__main__":
    main()
