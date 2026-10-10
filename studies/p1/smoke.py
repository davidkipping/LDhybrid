"""P1 end-to-end smoke run: PHOENIX v3.0 dwarf sequence (P2 subset, solar
metallicity) × key bands × all laws, headline-candidate method (r-weighted,
physical constraints), M1 limb. Writes smoke.json and smoke.md.

    .venv/bin/python studies/p1/smoke.py
"""

import json
from pathlib import Path

import numpy as np

from ldhybrid import pipeline
from ldhybrid.grids import phoenix

HERE = Path(__file__).parent
BANDS = ["TESS", "Kepler", "HST_WFC3_G141", "HST_STIS_G430L", "JWST_NIRISS_SOSS_o1",
         "JWST_NIRSpec_G395H", "JWST_NIRSpec_PRISM"]
STARS = [(t, 5.0 if t <= 3800 else 4.5) for t in
         (2600, 3000, 3400, 3800, 4200, 4600, 5000, 5400, 5800, 6200, 6600, 7000, 7400)]
LAWS = ["quadratic", "power2", "hybrid2", "claret4", "hybrid4", "hybrid5"]


def main():
    rows = []
    for teff, logg in STARS:
        s = phoenix.load(teff, logg, "-0.0")
        for b in BANDS:
            rows += pipeline.fit_model(s, b, LAWS, limb="M1", measure="r",
                                       constraint="simplex")
        print(teff, "done", flush=True)
    (HERE / "smoke.json").write_text(json.dumps(rows, indent=0))

    lines = ["# P1 smoke run: worst light-curve error (ppm) of profile fits",
             "",
             "PHOENIX-ACES v3.0, [M/H] = 0, log g 5.0 (≤ 3800 K) / 4.5; M1 limb cut; "
             "r-weighted fit; hybrid4/5 on the simplex, other linear laws on their exact "
             "physical region; error = max |ΔF| over p ∈ {0.03, 0.1, 0.15} × "
             "b ∈ {0, …, 0.9} against the exact light curve of the tabulated profile.",
             "", "| Teff | band | " + " | ".join(LAWS) + " |",
             "|---|---|" + "---|" * len(LAWS)]
    for teff, logg in STARS:
        for b in BANDS:
            r = {x["law"]: x for x in rows if x["teff"] == teff and x["band"] == b}
            lines.append(f"| {teff} | {b} | " +
                         " | ".join(f"{r[l]['lc_max_ppm']:.1f}" for l in LAWS) + " |")
    lines += ["", "| law | median | 90th pct | worst | unphysical fits |", "|---|---|---|---|---|"]
    for l in LAWS:
        e = np.array([x["lc_max_ppm"] for x in rows if x["law"] == l])
        bad = sum(1 for x in rows if x["law"] == l and x["physical"] is False)
        lines.append(f"| {l} | {np.median(e):.2f} | {np.percentile(e, 90):.2f} | "
                     f"{e.max():.2f} | {bad}/{len(e)} |")
    (HERE / "smoke.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[-len(LAWS) - 2:]))


if __name__ == "__main__":
    main()
