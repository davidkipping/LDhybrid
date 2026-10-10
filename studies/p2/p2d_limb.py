"""P2d — limb sensitivity. How much do the coefficients (and the light curves
they imply) move with choices made at the extreme limb, compared with the
spread between model grids?

Baseline: PHOENIX v3.0, M1 cut per band, fit on a 2000-pt resample, r-weighted
(MEASURE), simplex (hybrid4/5) or exact region. Variants, each refitted:
  limb:    M2, peak (parabolic) cut
  r_c:     one radius per model, from the photon-bolometric profile, for all bands
  nodes:   the M1 profile sampled only at the μ nodes of ATLAS (17), MPS-ATLAS
           (24) or Stagger (10), fitted on those nodes (what a coarser grid gives)
  grid:    MPS-ATLAS Set 2 at the same (Teff, log g, [M/H]) — plane-parallel,
           different code: the inter-grid spread
  set:     MPS-ATLAS Set 1 vs Set 2 (same code, different convection/abundances)

For each variant: |Δw| per weight, and the light-curve difference between the
two fitted laws, max over the (p, b) grid (ppm) — the consequence for a user.

    .venv/bin/python studies/p2/p2d_limb.py   (then --report)
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import common  # noqa: E402

HERE = Path(__file__).parent
OUT = HERE / "p2d_limb.jsonl"
MEASURE = "r"
LAWS = ["hybrid4", "hybrid5", "claret4", "hybrid2"]
NODES = {
    "atlas17": np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.125,
                         0.1, 0.075, 0.05, 0.025, 0.01]),
    "mps24": np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.35, 0.3, 0.25, 0.22, 0.2, 0.17,
                       0.15, 0.12, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.03, 0.02, 0.01]),
    "stagger10": np.array([1.0, 0.9, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.01]),
}


def pairs():
    """PHOENIX P2-subset models at [M/H] = 0 with an MPS-ATLAS match."""
    from ldhybrid.grids import mpsatlas, phoenix
    mps = set(mpsatlas.catalogue("set2"))
    out = []
    for t, g, mh in phoenix.P2_SUBSET:
        m = (float(mh), t, g)
        out.append(((t, g, mh), m if m in mps else None))
    return out


def task(job):
    from ldhybrid import bands, fit, laws
    from ldhybrid import lightcurve as lc
    from ldhybrid.grids import mpsatlas, phoenix
    from ldhybrid.profile import Profile, band_profile, limb_cut

    t, g, mh = job["key"]
    spec = phoenix.load(t, g, mh)
    lamall = np.array([spec.lam_edges[0], spec.lam_edges[-1]])
    _, rc_bol = limb_cut(spec.mu, spec.band(lamall, np.ones(2)), "M1")
    grid = lc.LCGrid()
    lam, S = bands.response(job["band"])

    profiles = {"base": band_profile(spec, lam, S, limb="M1"),
                "M2": band_profile(spec, lam, S, limb="M2"),
                "peak": band_profile(spec, lam, S, limb="peak"),
                "rc_model": band_profile(spec, lam, S, limb="M1", r_c=rc_bol)}
    base = profiles["base"]
    for k, nodes in NODES.items():
        profiles[k] = Profile(nodes, base.interpolant()(nodes))
    if job["mps"] is not None:
        for st in ("set1", "set2"):
            s = mpsatlas.load(st, *job["mps"])
            profiles[f"mps_{st}"] = band_profile(s, lam, S)

    def fitw(name, pr, native):
        L = laws.get(name)
        con = "simplex" if name in ("hybrid4", "hybrid5") else "physical"
        return L, fit.fit_linear(L, pr, MEASURE, None if native else 2000, con).w

    rows = []
    for name in LAWS:
        L, w0 = fitw(name, base, False)
        one, cols = grid.linear_columns(L) if name not in _CACHE else _CACHE[name]
        _CACHE[name] = (one, cols)
        dF0, F0 = grid.combine(one, cols, w0)
        out = dict(job_id=job["job_id"], key=job["key"], band=job["band"], law=name,
                   w_base=w0.tolist(), r_c=base.meta["r_c"], rc_bol=rc_bol, var={})
        for v, pr in profiles.items():
            if v == "base":
                continue
            _, w = fitw(name, pr, v in NODES)
            dF, F = grid.combine(one, cols, w)
            out["var"][v] = dict(w=w.tolist(), dw=(w - w0).tolist(),
                                 dlc_ppm=grid.errors(dF, F, dF0, F0)["max_ppm"],
                                 r_c=pr.meta.get("r_c"))
        rows.append(out)
    return rows


_CACHE = {}


def jobs():
    return [dict(job_id=f"{key}|{b}", key=list(key), mps=list(m) if m else None, band=b)
            for key, m in pairs() for b in common.BANDS]


def report():
    rows = common.read(OUT)
    variants = ["M2", "peak", "rc_model", "atlas17", "mps24", "stagger10", "mps_set2", "mps_set1"]
    lines = ["# P2d: limb sensitivity", "",
             f"PHOENIX v3.0 P2 subset × {len(common.BANDS)} bands; baseline M1 per band, "
             f"{MEASURE}-weighted, 2000-pt resample. Entries: median / 90th percentile of "
             "the light-curve difference (ppm, max over the (p, b) grid) between the law "
             "fitted to the variant and to the baseline. 'mps_set1' compares MPS Set 1 "
             "with Set 2 (not with PHOENIX).", "",
             "| variant | " + " | ".join(LAWS) + " |", "|---|" + "---|" * len(LAWS)]
    for v in variants:
        cells = []
        for l in LAWS:
            if v == "mps_set1":
                x = [_cross(r) for r in rows if r["law"] == l and "mps_set1" in r["var"]]
            else:
                x = [r["var"][v]["dlc_ppm"] for r in rows if r["law"] == l and v in r["var"]]
            cells.append(f"{np.median(x):.2f} / {np.percentile(x, 90):.2f}" if x else "—")
        lines.append(f"| {v} | " + " | ".join(cells) + " |")
    lines += ["", "## Weight shifts, hybrid5 (median |Δw_j| over model-bands)", "",
              "| variant | " + " | ".join(f"w{j + 1}" for j in range(5)) + " |",
              "|---|" + "---|" * 5]
    for v in variants[:-1]:
        D = np.array([np.abs(r["var"][v]["dw"]) for r in rows
                      if r["law"] == "hybrid5" and v in r["var"]])
        lines.append(f"| {v} | " + " | ".join(f"{x:.4f}" for x in np.median(D, 0)) + " |")
    W = np.array([r["w_base"] for r in rows if r["law"] == "hybrid5"])
    lines.append("| (baseline median w) | " + " | ".join(f"{x:.4f}" for x in np.median(W, 0)) + " |")
    (HERE / "p2d_limb.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def _cross(r):
    """MPS Set 1 vs Set 2 light-curve difference: via the stored weights."""
    from ldhybrid import laws
    from ldhybrid import lightcurve as lc
    global _G
    try:
        _G
    except NameError:
        _G = lc.LCGrid()
    if r["law"] not in _CACHE:
        _CACHE[r["law"]] = _G.linear_columns(laws.get(r["law"]))
    one, cols = _CACHE[r["law"]]
    a = _G.combine(one, cols, r["var"]["mps_set1"]["w"])
    b = _G.combine(one, cols, r["var"]["mps_set2"]["w"])
    return _G.errors(*a, *b)["max_ppm"]


if __name__ == "__main__":
    if "--report" in sys.argv:
        report()
    else:
        common.run_parallel(task, jobs(), OUT)
        report()
