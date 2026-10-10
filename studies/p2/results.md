# P2: method study — results and recommendations for Gate G2

Run 2026-10-09/10. Subset: 34 PHOENIX v3.0 models (spherical, M1 cut) and 16
MPS-ATLAS Set 2 models (plane-parallel) × 8 bands (TESS, Kepler, STIS G430L,
WFC3 G141, NIRISS SOSS o1, NIRSpec G235H/G395H/PRISM); P2c used 198 models ×
all 37 covered bands (7,524 model-bands). "LC error" = worst |ΔF| over
p ∈ {0.03, 0.1, 0.15} × b ∈ {0, …, 0.9} between the fitted law and the exact
light curve of the tabulated profile. Per-study tables: `p2a_weighting.md`,
`p2b_constraints.md`, `p2c_hybrid2_pole.md`, `p2c_fine.md`, `p2d_limb.md`,
`p2f_fisher.md`.

## P2a — weighting: r (or area) wins by 3–10×

Median / worst LC error (ppm), PHOENIX, 272 model-bands:

| weighting | hybrid5 | hybrid4 | Claret 4-term | power-2 | quadratic |
|---|---|---|---|---|---|
| equal weights on native nodes (Claret, Kostogryz) | 43 / 333 | 127 / 349 | 75 / 204 | 247 / 582 | 289 / 914 |
| uniform in μ | 13–23 / 33–58 | 14–25 / 70 | 44–55 / 123 | 57–63 / 152 | 182–187 / 446 |
| **uniform in r** | **3.3–4.9 / 10–11** | **6.5–7.7 / 23–25** | 6–13 / 23–33 | 15–22 / 72–83 | 62–64 / 177–186 |
| uniform in area (r²) | 3.8–3.9 / 11–12 | 6.2–8.7 / 26–29 | 7–12 / 28–31 | 15–19 / 75–84 | 58–65 / 168–188 |

(Ranges: native-node trapezoid vs 2000-pt PCHIP resample.) On MPS-ATLAS,
hybrid5 ties Claret 4-term (1.4 / 9 ppm, r-weighted); hybrid4 is weaker there
(11 / 22 ppm). r and area are statistically indistinguishable; r matches
ExoTETHyS. **Recommend: r-weighted, 2000-pt PCHIP resample** (behaves the same
on coarse and dense grids).

## P2b — constraints never bind; flux conservation buys nothing

Over 400 model-bands, every unconstrained hybrid fit was already physical
and inside the simplex (hybrid4/5) or triangle (hybrid2): constraints bound in
0% of fits, so imposing them is free. Flux conservation changes the median LC
error by −0.3 to +0.5 ppm (hybrids) and worsens quadratic. **Recommend: impose
the simplex (hybrid4/5) and the exact triangle (hybrid2) as guarantees; no
FCM in the headline tables** (keep the flux error as a column).

## P2c — hybrid2's pole: ε = 0.36

| law (7,524 model-bands) | median | 90th | 99th | worst |
|---|---|---|---|---|
| quadratic | 67 | 123 | 179 | 479 |
| hybrid2, ε = 0.208 (J-band tuned) | 87 | 192 | — | 581 |
| **hybrid2, ε = 0.36** | **57** | **93** | **122** | **259** |
| power-2 | 17 | 41 | 72 | 145 |

The 90th and 99th percentiles are minimised at ε = 0.35–0.36 for both PHOENIX
and MPS-ATLAS separately, and r and area weighting agree. At 0.36 hybrid2
beats quadratic at every percentile and in 73% of model-bands; it remains
well behind power-2 (which is not elementary for an ellipse). Above ε ≈ 0.45
its advantage collapses (41% at 0.46). **Recommend ε = 0.36.** Note: this
changes a constant shared with squishierplanet (`HYBRID2_EPS = 0.208`) and
MetalPlanet v0.10.x.

hybrid4/5's analytic ladder (α_eff = 1) needs no change: hybrid5 matches or
beats Claret 4-term everywhere tested.

## P2d — limb choices are small; model choice is not

Median / 90th percentile LC difference (ppm) between laws fitted to a variant
and to the baseline (PHOENIX, M1 per band):

| variant | hybrid4 | hybrid5 | Claret 4-term |
|---|---|---|---|
| M2 cut | 1.6 / 2.4 | 1.4 / 2.3 | 3.1 / 4.8 |
| parabolic peak | 0.7 / 1.5 | 0.6 / 1.4 | 1.3 / 3.1 |
| one r_c per model (all bands) | 2.6 / 5.8 | 2.4 / 5.6 | 4.9 / 11.3 |
| μ nodes of ATLAS (17) | 4.7 / 7.1 | 4.4 / 7.0 | 6.8 / 14.3 |
| μ nodes of MPS-ATLAS (24) | 5.6 / 8.1 | 4.6 / 7.8 | 8.1 / 16.4 |
| μ nodes of Stagger (10) | 3.6 / 5.5 | 2.2 / 4.6 | 4.8 / 10.8 |
| **MPS-ATLAS Set 2 instead of PHOENIX** | **172 / 447** | **172 / 447** | **172 / 446** |
| MPS-ATLAS Set 1 vs Set 2 | 32 / 128 | 32 / 128 | 32 / 128 |

* Every limb and sampling choice moves the light curves by ≲ 8 ppm (90th);
  the choice of model atmosphere moves them by ~170 ppm (median). The
  systematic floor is set by the grids, not by our method.
* **hybrid5 is robust in light-curve terms** (its sensitivities are the
  smallest of all laws), **but its innermost-pole weight w3 is not**: it
  shifts by a median 0.03 (M2), 0.04–0.10 (coarse μ nodes) against a
  baseline 0.25 and an inter-grid shift of 0.16. The light curve barely
  depends on w3, so it is weakly determined (see P2f). Tabulate it, flag it,
  and give its systematic spread.
* A single r_c per model costs 2–6 ppm: acceptable for coefficient spectra,
  where it preserves linear passband integration.

## P2e — custom radiative transfer: no-go on sampling grounds

μ-sampling effects (≲ 8 ppm) are ≪ the inter-grid spread (~170 ppm), so dense-μ
radiative transfer is not needed for accuracy. The remaining reason is
coverage: MIRI LRS for M dwarfs (no public grid beyond 5.5 µm below 3500 K).
**Recommend: defer**, pending the NewEra request.

## P2f — a transit constrains 1–2 combinations

Posterior s.d. along Fisher eigen-directions (median over stars/geometries),
geometry free:

| | hybrid4 | hybrid5 |
|---|---|---|
| JWST-like (100 ppm/min, one transit) | 0.003, 0.03, 0.25, 0.85 | 0.003, 0.03, 0.19, 0.62, 1.0 |
| TESS-like (1000 ppm/2 min) | 0.05, 0.40, 0.91, 1.0 | 0.05, 0.39, 0.91, 1.0, 1.0 |

One (TESS) to two (JWST) combinations are measured below 0.05; the rest are
set by the prior. The best-measured direction for hybrid5 is mostly the
polynomial weights (w1, w2) plus the outer pole (w5) — not the inner pole.
So theoretical tables matter most as *priors* on the unconstrained
combinations, which argues for publishing the systematic spread σ_w (P3).

## G2 freeze (v1.0) — adopted 2026-10-10

Adopted as recommended below, except item 4: MIRI M dwarfs are dropped from
v1.0 (MIRI for FGK via Kurucz only). Frozen in `ldhybrid/v1.py`; upstream
note for squishierplanet/MetalPlanet in `docs/upstream/hybrid2_pole_036.md`.

### Recommendations

1. Laws: hybrid2 with **ε = 0.36**; hybrid4/5 with the analytic ladder
   (unchanged).
2. Headline method: **r-weighted** least squares on a 2000-pt PCHIP resample;
   simplex (hybrid4/5) / exact triangle (hybrid2); no FCM. Secondary table:
   uniform-μ, for comparison with Claret.
3. Limb: **M1 per band** for band tables (r_c tabulated); **one r_c per
   model** for coefficient spectra.
4. Custom RT: deferred (MIRI-only motivation).
