"""Assemble results.md (analysis section) from the revision-study outputs.

    .venv/bin/python studies/rev/summarise.py
"""
import json
from pathlib import Path

import numpy as np

H = Path(__file__).parent
L6 = ["hybrid5", "hybrid4", "claret4", "power2", "hybrid2", "quadratic"]
NULL = ["poly2", "poly4", "poly5"]
NAMES = dict(hybrid5="hybrid5", hybrid4="hybrid4", claret4="four-term", power2="power-2",
             hybrid2="hybrid2", quadratic="quadratic", poly2="poly2 (2 terms)",
             poly4="poly4 (4 terms)", poly5="poly5 (5 terms)")


def g(v):
    return f"{v:.2f}" if v < 1 else f"{v:.1f}" if v < 10 else f"{v:.0f}"


def s3(x):
    x = np.asarray(x, float)
    return f"{g(np.median(x))} / {g(np.percentile(x, 99))} / {g(x.max())}"


out = ["# Revision analyses (studies/rev)", "",
       "Scripts in this directory; `rc.py` holds the shared set-up. P2S = the P2 subset "
       "(34 PHOENIX v3.0 + 16 MPS-ATLAS Set 2 models × 8 bands: TESS, Kepler, STIS G430L, WFC3 G141, "
       "NIRISS SOSS o1, NIRSpec G235H/G395H/PRISM). T3 = the 15 cases of the paper's Table 3. "
       "All fits use the exact quadrature (v1.RESAMPLE = 'exact': Gauss–Legendre on every interval "
       "between profile nodes), which replaced the 2000-point uniform-r resample on 2026-10-10. "
       "Unless stated: production fits (v1.0), light-curve (LC) error = worst |ΔF| in ppm at "
       "p = 0.1, b ∈ {0, 0.3, 0.6, 0.9} against the exact light curve of the tabulated profile "
       "(PCHIP in μ); statistics are median / 99th percentile / maximum.", ""]

# ---------------------------------------------------------------- D1
R = [json.loads(l) for l in open(H / "p2s_study.jsonl")]
T3 = json.load(open(H / "t3_study.json"))
cases = sorted({k.rsplit("|", 1)[0] for k in T3}, key=lambda s: list(T3).index(s + "|hybrid5"))
out += ["## D1 — null laws: even polynomials in μ² without poles (`p2s_study.py`, `t3_study.py`)", "",
        "poly*n* = 1 − Σ_{k≤n} w_k (1 − μ^{2k}), fitted r-weighted inside its exact physical region.", "",
        "| law (parameters) | PHOENIX P2S (272) | MPS-ATLAS 2 P2S (128) | T3 range (b = 0) |", "|---|---|---|---|"]
npar = dict(hybrid5=5, hybrid4=4, claret4=4, power2=2, hybrid2=2, quadratic=2, poly2=2, poly4=4, poly5=5)
for l in L6 + NULL:
    ph = [r["laws"][l]["lc"] for r in R if r["kind"] == "phoenix"]
    mp = [r["laws"][l]["lc"] for r in R if r["kind"] == "mps2"]
    t3 = [T3[c + "|" + l]["b0"] for c in cases]
    out.append(f"| {NAMES[l]} ({npar[l]}) | {s3(ph)} | {s3(mp)} | {g(min(t3))}–{g(max(t3))} |")
med = lambda l, k: np.median([r["laws"][l]["lc"] for r in R if r["kind"] == k])
out += ["", "Adding the poles at fixed parameter count (P2S medians, PHOENIX; MPS-ATLAS 2): " + "; ".join(
    f"{a} {g(med(a, 'phoenix'))} → {b} {g(med(b, 'phoenix'))} ({g(med(a, 'mps2'))} → {g(med(b, 'mps2'))})"
    for a, b in (("poly2", "hybrid2"), ("poly4", "hybrid4"), ("poly5", "hybrid5"))) + ".", ""]

# ---------------------------------------------------------------- D2 / D6
if (H / "d2_d6_full.jsonl").exists():
    FULL = [json.loads(l) for l in open(H / "d2_d6_full.jsonl")]
    ph = [r for r in FULL if r["kind"] == "phoenix"]
    mp = [r for r in FULL if r["kind"] == "mps2"]
    nm = len({tuple(r["key"]) for r in mp})
    out += ["## Full-grid fit quality with the exact quadrature (`d2_d6_full.py`)", "",
            f"Recomputed from the atmospheres (not from the products): all PHOENIX models × 38 bands "
            f"({len(ph)} rows) and {nm} random MPS-ATLAS Set 2 models × 38 bands ({len(mp)} rows); "
            "r-weighted production fits, production light-curve grid. Median / 99th / max (ppm).", "",
            "| law | PHOENIX | MPS-ATLAS 2 (sample) |", "|---|---|---|"]
    for l in L6:
        out.append(f"| {NAMES[l]} | {s3([r['lc'][l] for r in ph])} | {s3([r['lc'][l] for r in mp])} |")
    for kind, rr in (("PHOENIX", ph), ("MPS 2", mp)):
        con = [r["lc"].get("power2_con", r["lc"]["power2"]) for r in rr]
        out.append(f"| power-2, unphysical rows refitted in 0 ≤ c ≤ 1 ({kind}) | " +
                   (f"{s3(con)} |  |" if kind == "PHOENIX" else f" | {s3(con)} |"))
    out += ["", "## D6 — how often the constraints are active (`d2_d6_full.py`)", "",
            "% of rows with an active constraint (hybrids: a simplex face, or for hybrid2 the triangle; "
            "quadratic/four-term: their exact positive–monotonic region). Power-2 (fitted without "
            "constraint): % of rows unphysical.", "",
            "| law | PHOENIX r | PHOENIX μ | MPS 2 r | MPS 2 μ |", "|---|---|---|---|---|"]
    for l in ["hybrid5", "hybrid4", "hybrid2", "claret4", "quadratic", "power2"]:
        cells = [f"{100 * np.mean([r[k][l] for r in rr]):.2f}" for rr in (ph, mp) for k in ("bound", "bound_mu")]
        out.append(f"| {NAMES[l]}{' (unphysical)' if l == 'power2' else ''} | " + " | ".join(cells) + " |")
    out += ["", "## D2 — power-2 inside its physical region; D6b — what the four-term constraint costs", "",
            "| rows | n | LC unconstrained | LC constrained | increase |", "|---|---|---|---|---|"]
    for kind, rr in (("PHOENIX", ph), ("MPS 2", mp)):
        u = [r for r in rr if "power2_con" in r["lc"]]
        if u:
            a0 = np.array([r["lc"]["power2"] for r in u]); a1 = np.array([r["lc"]["power2_con"] for r in u])
            out.append(f"| power-2, {kind}, LM fit unphysical | {len(u)} | {s3(a0)} | {s3(a1)} | {s3(a1 - a0)} |")
        u = [r for r in rr if "claret4_uncon" in r["lc"]]
        if u:
            a0 = np.array([r["lc"]["claret4_uncon"] for r in u]); a1 = np.array([r["lc"]["claret4"] for r in u])
            out.append(f"| four-term, {kind}, constraint active | {len(u)} | {s3(a0)} | {s3(a1)} | "
                       f"{g(np.median(a1 - a0))} / {g(np.percentile(a1 - a0, 99))} / {g((a1 - a0).max())} (min {(a1 - a0).min():.1f}) |")
    cs = [r["power2_c"][0] for r in ph if "power2_c" in r]
    out += ["", f"Refitted power-2 rows all end on c = 1; their LM fits have c = {min(cs):.3f}–{max(cs):.3f} "
            f"(median {np.median(cs):.3f}, 90th percentile {np.percentile(cs, 90):.2f}; the largest have α → 0, a log-like law). Four-term: a constrained fit can score better in light-curve terms, "
            "because least squares in intensity is not the light-curve metric.", ""]

# ---------------------------------------------------------------- D3
out += [open(H / "d3.md").read().replace(
    "## D3 — hybrid2 pole, held-out",
    "## D3 — hybrid2 pole, held-out (`d3_hybrid2_heldout.py`; P2c data, 2000-point resample)"), ""]
PX = json.load(open(H / "p2x_pole.json"))
out += ["### hybrid2 pole re-scanned with the exact quadrature (`p2_exact.py pole`)", "",
        "P2S, P2 light-curve metric (p ∈ {0.03, 0.1, 0.15}, b ∈ {0, …, 0.9}); median / 90th / 99th (ppm).", "",
        "| ε | PHOENIX | MPS-ATLAS 2 | both |", "|---|---|---|---|"]
def s390(x):
    x = np.asarray(x, float)
    return f"{g(np.median(x))} / {g(np.percentile(x, 90))} / {g(np.percentile(x, 99))}"
for e in ["0.24", "0.28", "0.32", "0.36", "0.4", "0.44", "quadratic"]:
    out.append(f"| {e} | {s390([r['lc'][e] for r in PX if r['kind'] == 'phoenix'])} | "
               f"{s390([r['lc'][e] for r in PX if r['kind'] == 'mps2'])} | {s390([r['lc'][e] for r in PX])} |")
out.append("")

# ---------------------------------------------------------------- D4
out += ["## D4 — planet size and grazing geometry (`t3_study.py`, `p2s_study.py`)", "",
        "T3 cases (PHOENIX, b as listed; 'graze' = b = 1 − p/2). Median over the 15 cases of the worst "
        "|ΔF| (ppm); in brackets the maximum.", "",
        "| law | p | b = 0 | b = 0.5 | b = 1 − p/2 | b = 1 | rms (b = 0) |", "|---|---|---|---|---|---|---|"]
for l in L6 + NULL:
    for pp in ("0.03", "0.1", "0.2", "0.3"):
        W = np.array([T3[c + "|" + l]["geo"][pp]["worst"] for c in cases])
        Rm = np.array([T3[c + "|" + l]["geo"][pp]["rms"] for c in cases])
        cells = [f"{g(np.median(W[:, j]))} ({g(W[:, j].max())})" for j in range(4)]
        out.append(f"| {NAMES[l] if pp == '0.03' else ''} | {pp} | " + " | ".join(cells) +
                   f" | {g(np.median(Rm[:, 0]))} |")
# empirical exponent
out += ["", "Scaling with p (T3, b = 0, median over cases of the per-case log-slope between p = 0.03 "
        "and 0.3): " + ", ".join(
            f"{NAMES[l]} {np.median([np.log(T3[c+'|'+l]['geo']['0.3']['worst'][0] / T3[c+'|'+l]['geo']['0.03']['worst'][0]) / np.log(10) for c in cases]):.2f}"
            for l in L6) + ". The error grows roughly as p^1.1–1.5, not p².", "",
        "P2S (all 400 model-bands; b ∈ {0, 0.3, 0.6, 0.9} unless 'graze'): median / 99th percentile.", "",
        "| law | grid | p = 0.03 | p = 0.1 | p = 0.2 | p = 0.3 | graze p = 0.1 | graze 0.2 | graze 0.3 |",
        "|---|---|---|---|---|---|---|---|---|"]
for l in L6:
    for kind in ("phoenix", "mps2"):
        rr = [r["laws"][l] for r in R if r["kind"] == kind]
        f = lambda k: f"{g(np.median([x[k] for x in rr]))} / {g(np.percentile([x[k] for x in rr], 99))}"
        out.append(f"| {NAMES[l] if kind == 'phoenix' else ''} | {kind} | {f('lc_p0.03')} | {f('lc')} | "
                   f"{f('lc_p0.2')} | {f('lc_p0.3')} | {f('lc_graze0.1')} | {f('lc_graze0.2')} | {f('lc_graze0.3')} |")
out += ["", "A grazing transit never gives a larger error than a full transit of the same planet "
        "in these samples: the worst case for every law is the limb crossing of a full transit.", ""]

# ---------------------------------------------------------------- D5
out += ["## D5 — dependence on the interpolant of the 'exact' profile (`p2s_study.py`, `d5b_interp_mu.py`)", "",
        "First block: light-curve difference between alternative truths and the baseline (PCHIP in μ) "
        "through the same nodes. Second: the production fits re-scored against each truth.", "",
        "| interpolant | PHOENIX: truth vs PCHIP | MPS 2: truth vs PCHIP |", "|---|---|---|"]
for k, lab in (("lin_mu", "linear in μ"), ("lin_r", "linear in r"), ("akima_r", "Akima in r")):
    out.append(f"| {lab} | {s3([r['truth_diff'][k] for r in R if r['kind'] == 'phoenix'])} | "
               f"{s3([r['truth_diff'][k] for r in R if r['kind'] == 'mps2'])} |")
B5 = H / "d5b.json"
if B5.exists():
    D5 = json.load(open(B5))
    for k, lab in (("akima_mu", "Akima in μ"), ("spline_mu", "cubic spline in μ")):
        out.append(f"| {lab} | {s3([r['truth_diff'][k] for r in D5 if r['kind'] == 'phoenix'])} | "
                   f"{s3([r['truth_diff'][k] for r in D5 if r['kind'] == 'mps2'])} |")
out += ["", "| law | grid | PCHIP μ (baseline) | linear μ | Akima μ | spline μ | Akima r | linear r |",
        "|---|---|---|---|---|---|---|---|"]
for l in ["hybrid5", "hybrid4", "claret4", "power2", "quadratic"]:
    for kind in ("phoenix", "mps2"):
        rr = [r["laws"][l] for r in R if r["kind"] == kind]
        f = lambda k: s3([x[k] for x in rr])
        extra = ["–", "–"]
        if B5.exists():
            r5 = [r["laws"][l] for r in D5 if r["kind"] == kind]
            extra = [s3([x["akima_mu"] for x in r5]), s3([x["spline_mu"] for x in r5])]
        out.append(f"| {NAMES[l] if kind == 'phoenix' else ''} | {kind} | {f('lc')} | {f('lc_lin_mu')} | "
                   f"{extra[0]} | {extra[1]} | {f('lc_akima_r')} | {f('lc_lin_r')} |")
out += ["", "Linear interpolation in r is a poor interpolant near the disc centre, where the nodes are "
        "sparse in r and I is quadratic in r; it is shown only as an extreme case.", ""]

# ---------------------------------------------------------------- D7
out += ["## D7 — what r-weighting costs at the limb (`p2s_study.py`)", "",
        "σ_μ = rms residual uniform in μ (Claret's merit function), σ_r = rms uniform in r; "
        "max |ΔI| in I/I(1) for μ < 0.05 and μ ≥ 0.05; medians over P2S. 'r' = r-weighted fit "
        "(tabulated), 'μ' = uniform-μ fit (secondary table).", "",
        "| law | grid | σ_μ r / μ | σ_r r / μ | max|ΔI| μ<0.05 r / μ | max|ΔI| μ≥0.05 r / μ | LC error r / μ |",
        "|---|---|---|---|---|---|---|"]
for l in L6:
    for kind in ("phoenix", "mps2"):
        rr = [r["laws"][l] for r in R if r["kind"] == kind]
        m = lambda k: np.nanmedian([x[k] for x in rr])
        e = lambda v: f"{v:.1e}"
        out.append(f"| {NAMES[l] if kind == 'phoenix' else ''} | {kind} | {e(m('sig_mu'))} / {e(m('sig_mu_mufit'))} | "
                   f"{e(m('sig_r'))} / {e(m('sig_r_mufit'))} | {m('dI_limb'):.3f} / {m('dI_limb_mufit'):.3f} | "
                   f"{m('dI_in'):.4f} / {m('dI_in_mufit'):.4f} | {g(m('lc'))} / {g(m('lc_mufit'))} |")
out.append("")

# ---------------------------------------------------------------- D12
out += ["## D12 — where the light-curve error comes from (`t3_study.py`, `p2s_study.py`)", "",
        "The profile error ΔI = I_law − I_model is split at μ = 0.05; each part is scored on its own "
        "(law below 0.05 and model above, and the reverse). Ratios to the full error (median; T3 at "
        "b = 0, P2S on the production grid).", "",
        "| law | T3: μ<0.05 part | T3: μ≥0.05 part | P2S PHOENIX: μ<0.05 | μ≥0.05 | P2S MPS 2: μ<0.05 | μ≥0.05 |",
        "|---|---|---|---|---|---|---|"]
for l in L6 + NULL:
    a = [T3[c + "|" + l]["limbpart"] / T3[c + "|" + l]["b0"] for c in cases]
    b = [T3[c + "|" + l]["inpart"] / T3[c + "|" + l]["b0"] for c in cases]
    cells = []
    for kind in ("phoenix", "mps2"):
        rr = [r["laws"][l] for r in R if r["kind"] == kind]
        cells += [f"{np.median([x['lc_limbpart'] / x['lc'] for x in rr]):.2f}",
                  f"{np.median([x['lc_inpart'] / x['lc'] for x in rr]):.2f}"]
    out.append(f"| {NAMES[l]} | {np.median(a):.2f} | {np.median(b):.2f} | " + " | ".join(cells) + " |")
qd = [(T3[c + "|quadratic"]["b0"], T3[c + "|quadratic"]["limbpart"], T3[c + "|quadratic"]["inpart"])
      for c in cases[:9]]
out += ["", "M dwarfs (T3, quadratic): full 72–114 ppm; the part from μ < 0.05 alone gives "
        f"{min(x[1] for x in qd):.0f}–{max(x[1] for x in qd):.0f} ppm, the part from μ ≥ 0.05 "
        f"{min(x[2] for x in qd):.0f}–{max(x[2] for x in qd):.0f} ppm. The quadratic deficit lies in the "
        "disc interior: a law that cannot bend at the limb is pulled away from the profile over the "
        "whole disc. hybrid5 has the largest share of its (small) error at μ < 0.05.", ""]

# ---------------------------------------------------------------- D13
D13 = json.load(open(H / "d13.json"))
out += ["## D13 — photon vs energy integration (`d13_photon_energy.py`)", "",
        "Five T3 stars (PHOENIX). 'profile' = light-curve difference between the energy- and photon-"
        "weighted model profiles; the fitted laws carry the same difference (to within a few ppm). max |Δw|: "
        "largest coefficient shift.", "",
        "| band | profile LC difference (ppm), min–max over stars | max |Δw| hybrid5 | four-term | quadratic | power-2 |",
        "|---|---|---|---|---|---|"]
for b in ["TESS", "Kepler", "HST_WFC3_G141", "JWST_NIRSpec_PRISM", "2MASS_J"]:
    es = [v for k, v in D13.items() if k.endswith("|" + b)]
    pl = [e["profile_lc"] for e in es]
    dw = {n: max(np.max(np.abs(e[n]["dw"])) for e in es) for n in ("hybrid5", "claret4", "quadratic", "power2")}
    out.append(f"| {b} | {g(min(pl))}–{g(max(pl))} | {dw['hybrid5']:.3f} | {dw['claret4']:.3f} | "
               f"{dw['quadratic']:.3f} | {dw['power2']:.3f} |")
out.append("")

# ---------------------------------------------------------------- D14
D14 = json.load(open(H / "d14.json"))
cls = lambda t: "M (<3900 K)" if t < 3900 else ("K (3900–5200 K)" if t < 5200 else "FG (>5200 K)")
out += ["## D14 — narrow bins summed from the coefficient spectra (`d14_narrow_bins.py`)", "",
        f"{len(D14)} random PHOENIX models; top-hat bins at resolving power R over 0.6–5.3 µm. "
        "Percentage of bins whose summed coefficients are unphysical (positivity or monotonicity "
        "fails); in brackets the percentage of models with at least one such bin.", "",
        "| R | law | M (<3900 K) | K (3900–5200 K) | FG (>5200 K) | wavelengths of unphysical bins |",
        "|---|---|---|---|---|---|"]
for R_ in ("100", "300", "1000"):
    for n in ("hybrid2", "hybrid4", "hybrid5"):
        cells = []
        for c in ("M (<3900 K)", "K (3900–5200 K)", "FG (>5200 K)"):
            rr = [r["res"][R_][n] for r in D14 if cls(r["teff"]) == c]
            cells.append(f"{100 * sum(x['unphysical'] for x in rr) / sum(x['n'] for x in rr):.2f} "
                         f"({100 * np.mean([x['unphysical'] > 0 for x in rr]):.0f})")
        lam = np.concatenate([r["res"][R_][n]["lam_bad"] for r in D14]) / 1e4
        where = f"{np.percentile(lam, 5):.2f}–{np.percentile(lam, 95):.2f} µm (90%)" if lam.size else "–"
        out.append(f"| {R_} | {n} | " + " | ".join(cells) + f" | {where} |")
nM = sum(cls(r["teff"]) == "M (<3900 K)" for r in D14)
out += ["", f"Models: {nM} M, {sum(cls(r['teff']) == 'K (3900–5200 K)' for r in D14)} K, "
        f"{sum(cls(r['teff']) == 'FG (>5200 K)' for r in D14)} FG. Unphysical narrow bins occur only below "
        "~2.5 µm, in cool stars (molecular bands); none for FG stars.", ""]

# ---------------------------------------------------------------- P2 exact
W = json.load(open(H / "p2x_weighting.json"))
out += ["## Table 1 analogue — fit weighting with the exact quadrature (`p2_exact.py weighting`)", "",
        "P2S; P2 light-curve metric (p ∈ {0.03, 0.1, 0.15}, b ∈ {0, …, 0.9}); median / worst (ppm). "
        "'nodes' = equal weights on the native μ nodes (Claret 2000; Kostogryz et al. 2022).", ""]
for kind, lab in (("phoenix", "PHOENIX (34 × 8)"), ("mps2", "MPS-ATLAS Set 2 (16 × 8)")):
    rr = [r for r in W if r["kind"] == kind]
    out += [f"**{lab}**", "", "| weighting | " + " | ".join(NAMES[l] for l in
            ["hybrid5", "hybrid4", "claret4", "power2", "hybrid2", "quadratic"]) + " |", "|---|" + "---|" * 6]
    for m, ml in (("nodes", "equal per native node"), ("mu", "uniform in μ"), ("r", "uniform in r"),
                  ("area", "uniform in area")):
        cells = []
        for l in ["hybrid5", "hybrid4", "claret4", "power2", "hybrid2", "quadratic"]:
            x = [r["lc"][f"{l}|{m}"] for r in rr]
            cells.append(f"{g(np.median(x))} / {g(np.max(x))}")
        out.append(f"| {ml} | " + " | ".join(cells) + " |")
    out.append("")
    fac = []
    for l in ["hybrid5", "claret4", "power2", "quadratic"]:
        mu_ = np.array([r["lc"][f"{l}|mu"] for r in rr]); r_ = np.array([r["lc"][f"{l}|r"] for r in rr])
        fac.append(f"{NAMES[l]} ×{np.median(mu_) / np.median(r_):.1f} (median), ×{mu_.max() / r_.max():.1f} (worst)")
    out += ["Uniform-μ → uniform-r improvement: " + "; ".join(fac) + ".", ""]

LB = json.load(open(H / "p2x_limb.json"))
out += ["## Limb sensitivity with the exact quadrature (`p2_exact.py limb`)", "",
        "PHOENIX P2S (34 × 8); baseline M1 per band, r-weighted. Light-curve difference (ppm, P2 "
        "metric) between the law fitted to the variant and to the baseline: median / 90th "
        "percentile. 'MPS-ATLAS 2 vs PHOENIX' and 'Set 1 vs Set 2' where the MPS-ATLAS grid has the model.", "",
        "| variant | hybrid5 | hybrid4 | four-term | hybrid2 |", "|---|---|---|---|---|"]
for v, vl in (("M2", "M2 instead of M1"), ("rc_model", "one r_c per model"), ("atlas17", "μ nodes of ATLAS (17)"),
              ("mps24", "μ nodes of MPS-ATLAS (24)"), ("stagger10", "μ nodes of Stagger (10)"),
              ("mps_set2", "MPS-ATLAS 2 vs PHOENIX"), ("set1_vs_set2", "MPS-ATLAS Set 1 vs Set 2")):
    cells = []
    for l in ("hybrid5", "hybrid4", "claret4", "hybrid2"):
        x = [r["dlc"][f"{l}|{v}"] for r in LB if f"{l}|{v}" in r["dlc"]]
        cells.append(f"{g(np.median(x))} / {g(np.percentile(x, 90))}" if x else "–")
    out.append(f"| {vl} | " + " | ".join(cells) + " |")
D = np.array([np.abs(r["dw"]["hybrid5|M2"]) for r in LB])
Dn = np.array([np.abs(r["dw"]["hybrid5|mps24"]) for r in LB])
B0 = np.array([r["dw"]["hybrid5|base"] for r in LB])
out += ["", "hybrid5 median |Δw_j| (w1…w5): M2 " + ", ".join(f"{x:.3f}" for x in np.median(D, 0)) +
        "; MPS-ATLAS nodes " + ", ".join(f"{x:.3f}" for x in np.median(Dn, 0)) +
        "; baseline median w " + ", ".join(f"{x:.3f}" for x in np.median(B0, 0)) + ".", ""]

# ---------------------------------------------------------------- D12b / D12c
B = json.load(open(H / "d12b.json"))
out += ["## D12b — is the M-dwarf quadratic deficit caused by the steep limb? (`d12b_limb_excluded.py`)", "",
        "Each law fitted (r-weighted) to μ ≥ 0.05 only; 'interior' = light-curve error with the model "
        "profile kept below μ = 0.05 (only the disc interior compared); 'full' = against the full model. "
        "T3, b = 0, p = 0.1. Range over cases (ppm).", "",
        "| law | interior, M dwarfs | interior, K/G | full, all |", "|---|---|---|---|"]
for l in ["quadratic", "power2", "claret4", "hybrid2", "hybrid4", "hybrid5"]:
    vm = [x for k, x in B.items() if k.endswith("|" + l) and " M|" in k]
    vk = [x for k, x in B.items() if k.endswith("|" + l) and " M|" not in k]
    va = [x for k, x in B.items() if k.endswith("|" + l)]
    out.append(f"| {NAMES[l]} | {g(min(x['interior'] for x in vm))}–{g(max(x['interior'] for x in vm))} | "
               f"{g(min(x['interior'] for x in vk))}–{g(max(x['interior'] for x in vk))} | "
               f"{g(min(x['full'] for x in va))}–{g(max(x['full'] for x in va))} |")
C = json.load(open(H / "d12c.json"))
out += ["", "## D12c — fitting without the extreme limb (`d12c_mu_min.py`)", "",
        "r-weighted fits to μ ≥ μ_min only, scored against the full model light curve (production grid). "
        "P2S; median / 99th (ppm).", "",
        "| law | grid | μ_min = 0 (tabulated) | 0.02 | 0.05 | 0.1 |", "|---|---|---|---|---|---|"]
for l in L6:
    for kind in ("phoenix", "mps2"):
        rr = [r for r in C if r["kind"] == kind]
        cells = [f"{g(np.median([r['lc'][f'{l}|{m}'] for r in rr]))} / "
                 f"{g(np.percentile([r['lc'][f'{l}|{m}'] for r in rr], 99))}" for m in ("0.0", "0.02", "0.05", "0.1")]
        out.append(f"| {NAMES[l] if kind == 'phoenix' else ''} | {kind} | " + " | ".join(cells) + " |")
out.append("")

# ---------------------------------------------------------------- E1
E = H / "e1_bias.json"
if E.exists():
    E1 = json.load(open(E))
    out += ["## E1 — parameter bias, noiseless (`e1_bias.py`)", "",
            "T3 stars and bands; truth p = 0.1, b = 0.3, chord velocity v = 1 (600 points). "
            "Fit of (p, b, v) with coefficients fixed at the tabulated values, or free "
            "(quadratic/four-term unconstrained; hybrid2 triangle and hybrid4/5 simplex via q; "
            "power-2 with 0 ≤ c ≤ 1). Δ(p²) in ppm of a 10,000 ppm depth.", "",
            "| law | mode | Δ(p²) ppm: median |·| / max |·| | Δb median |·| / max | Δv/v max |·| | rms resid. ppm (median / max) | chromatic spread of Δ(p²) per star (max over stars) |",
            "|---|---|---|---|---|---|---|"]
    for l in L6:
        for mode in ("fixed", "free"):
            v = [E1[c + "|" + l][mode] for c in cases]
            dp = np.abs([x["dp2_ppm"] for x in v])
            db = np.abs([x["db"] for x in v])
            dv = np.abs([x["dv"] for x in v])
            rm = [x["rms_ppm"] for x in v]
            spread = []
            for st in ("late M", "mid M", "early M", "K dwarf", "Sun analogue"):
                s = [E1[c + "|" + l][mode]["dp2_ppm"] for c in cases if c.startswith(st + "|")]
                spread.append(max(s) - min(s))
            out.append(f"| {NAMES[l] if mode == 'fixed' else ''} | {mode} | {g(np.median(dp))} / {g(dp.max())} | "
                       f"{np.median(db):.4f} / {db.max():.4f} | {dv.max():.4f} | {g(np.median(rm))} / {g(max(rm))} | "
                       f"{g(max(spread))} |")
    out += ["", "Per case (Δ(p²) ppm, fixed / free):", "",
            "| case | " + " | ".join(NAMES[l] for l in L6) + " |", "|---|" + "---|" * len(L6)]
    for c in cases:
        out.append(f"| {c} | " + " | ".join(
            f"{E1[c+'|'+l]['fixed']['dp2_ppm']:+.1f} / {E1[c+'|'+l]['free']['dp2_ppm']:+.1f}" for l in L6) + " |")
    out.append("")

txt = "\n".join(out) + "\n"
(H / "results.md").write_text(txt)
print(txt)
