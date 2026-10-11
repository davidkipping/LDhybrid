"""C15: the P2f Fisher analysis (studies/p2/p2f_fisher.py) rerun with the v1
(exact-quadrature) fits, with the set-up stated and the best-measured
eigen-directions of hybrid5 reported.

Set-up (from P2f): circular orbit, P = 3 d, a/R★ = 10, one transit, sampled
over ±0.75 T14; JWST-like σ = 100 ppm per 1-min point, TESS-like σ = 1000 ppm
per 2-min point; p ∈ {0.1, 0.15}, b ∈ {0.3, 0.7}; free geometry p, b, a/R★, t0;
a weak prior σ = 1 per weight. "Constrained" = posterior s.d. < 0.05 along an
eigen-direction of the weight covariance (unit-norm direction in w).

    .venv/bin/python studies/rev2/c15_fisher.py
"""
import importlib.util
import json
from pathlib import Path

import numpy as np

import rc
from ldhybrid import bands, fit, laws, v1
from ldhybrid.grids import phoenix
from ldhybrid.profile import band_profile

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("p2f", ROOT / "studies" / "p2" / "p2f_fisher.py")
P2F = importlib.util.module_from_spec(spec)
spec.loader.exec_module(P2F)


def task(sb):
    key, band = sb
    pr = band_profile(phoenix.load(*key), *bands.response(band), limb="M1")
    out = []
    for name in ("hybrid4", "hybrid5"):
        L = laws.get(name)
        w = fit.fit_linear(L, pr, "r", v1.RESAMPLE, "simplex").w
        for case in P2F.CASES:
            for p, b in P2F.GEOM:
                Fm = P2F.fisher(L, w, p, b, case)
                a = P2F.analyse(Fm, L.m)
                Cfree = np.linalg.inv(Fm + np.pad(np.eye(L.m), ((0, Fm.shape[0] - L.m),) * 2))[:L.m, :L.m]
                ev, V = np.linalg.eigh(Cfree)
                out.append(dict(star=key, band=band, law=name, case=case["name"], p=p, b=b,
                                sd_free=np.sqrt(ev).tolist(), n_free=int(np.sum(np.sqrt(ev) < 0.05)),
                                n_fixed=a["fixed"]["n_lt_005"], sd_fixed=a["fixed"]["sd"],
                                dirs=[(V[:, k] * np.sign(V[np.argmax(np.abs(V[:, k])), k])).tolist() for k in range(2)]))
    return out


def main():
    rows = [r for res in rc.run_parallel(task, P2F.STARS, workers=4) for r in res]
    (rc.HERE / "c15_fisher.json").write_text(json.dumps(rows, default=str))
    L = ["## C15 — Fisher analysis with the exact-quadrature fits", "",
         "Set-up: P = 3 d, a/R★ = 10, one transit sampled over ±0.75 T14; JWST-like 100 ppm per 1-min point; "
         "TESS-like 1000 ppm per 2-min point; p ∈ {0.1, 0.15}, b ∈ {0.3, 0.7}; geometry (p, b, a/R★, t0) free; "
         "4 star-bands (5800 K G141 and TESS, 3400 K PRISM, 4600 K SOSS o1).", "",
         "| law | case | combinations with s.d. < 0.05, free geometry (min–max over cases) | median s.d. of 1st / 2nd / 3rd direction |",
         "|---|---|---|---|"]
    for name in ("hybrid4", "hybrid5"):
        for case in ("JWST-like", "TESS-like"):
            sel = [r for r in rows if r["law"] == name and r["case"] == case]
            n = [r["n_free"] for r in sel]
            sd = np.median([r["sd_free"][:3] for r in sel], axis=0)
            L.append(f"| {name} | {case} | {min(n)}–{max(n)} | " + " / ".join(f"{x:.3f}" for x in sd) + " |")
    L += ["", "hybrid5 best-measured direction (free geometry), components (w1…w5), median over cases:"]
    for case in ("JWST-like", "TESS-like"):
        sel = [r for r in rows if r["law"] == "hybrid5" and r["case"] == case]
        d1 = np.median([r["dirs"][0] for r in sel], axis=0)
        d2 = np.median([r["dirs"][1] for r in sel], axis=0)
        L.append(f"* {case}: 1st {np.round(d1, 2).tolist()}; 2nd {np.round(d2, 2).tolist()}")
    (rc.HERE / "c15_fisher.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
