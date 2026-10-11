# Revision analyses (studies/rev)

Scripts in this directory; `rc.py` holds the shared set-up. P2S = the P2 subset (34 PHOENIX v3.0 + 16 MPS-ATLAS Set 2 models × 8 bands: TESS, Kepler, STIS G430L, WFC3 G141, NIRISS SOSS o1, NIRSpec G235H/G395H/PRISM). T3 = the 15 cases of the paper's Table 3. All fits use the exact quadrature (v1.RESAMPLE = 'exact': Gauss–Legendre on every interval between profile nodes), which replaced the 2000-point uniform-r resample on 2026-10-10. Unless stated: production fits (v1.0), light-curve (LC) error = worst |ΔF| in ppm at p = 0.1, b ∈ {0, 0.3, 0.6, 0.9} against the exact light curve of the tabulated profile (PCHIP in μ); statistics are median / 99th percentile / maximum.

## D1 — null laws: even polynomials in μ² without poles (`p2s_study.py`, `t3_study.py`)

poly*n* = 1 − Σ_{k≤n} w_k (1 − μ^{2k}), fitted r-weighted inside its exact physical region.

| law (parameters) | PHOENIX P2S (272) | MPS-ATLAS 2 P2S (128) | T3 range (b = 0) |
|---|---|---|---|
| hybrid5 (5) | 1.2 / 5.3 / 5.4 | 0.79 / 4.2 / 4.9 | 0.42–2.2 |
| hybrid4 (4) | 4.3 / 13 / 14 | 7.1 / 12 / 13 | 3.8–15 |
| four-term (4) | 3.7 / 10 / 11 | 0.84 / 3.8 / 5.2 | 0.16–5.4 |
| power-2 (2) | 8.4 / 50 / 53 | 6.5 / 22 / 40 | 0.89–15 |
| hybrid2 (2) | 41 / 95 / 107 | 43 / 69 / 70 | 28–97 |
| quadratic (2) | 42 / 115 / 125 | 51 / 80 / 83 | 25–113 |
| poly2 (2 terms) (2) | 193 / 377 / 388 | 194 / 299 / 331 | 164–375 |
| poly4 (4 terms) (4) | 44 / 83 / 87 | 44 / 63 / 65 | 37–87 |
| poly5 (5 terms) (5) | 24 / 47 / 49 | 24 / 34 / 35 | 21–49 |

Adding the poles at fixed parameter count (P2S medians, PHOENIX; MPS-ATLAS 2): poly2 193 → hybrid2 41 (194 → 43); poly4 44 → hybrid4 4.3 (44 → 7.1); poly5 24 → hybrid5 1.2 (24 → 0.79).

## Full-grid fit quality with the exact quadrature (`d2_d6_full.py`)

Recomputed from the atmospheres (not from the products): all PHOENIX models × 38 bands (63916 rows) and 300 random MPS-ATLAS Set 2 models × 38 bands (11400 rows); r-weighted production fits, production light-curve grid. Median / 99th / max (ppm).

| law | PHOENIX | MPS-ATLAS 2 (sample) |
|---|---|---|
| hybrid5 | 1.3 / 5.6 / 51 | 0.94 / 4.3 / 6.8 |
| hybrid4 | 4.7 / 17 / 62 | 7.0 / 14 / 17 |
| four-term | 3.4 / 12 / 75 | 1.2 / 5.7 / 7.4 |
| power-2 | 8.4 / 61 / 189 | 7.0 / 36 / 108 |
| hybrid2 | 43 / 113 / 344 | 43 / 84 / 105 |
| quadratic | 42 / 133 / 507 | 51 / 100 / 127 |
| power-2, unphysical rows refitted in 0 ≤ c ≤ 1 (PHOENIX) | 11 / 100 / 357 |  |
| power-2, unphysical rows refitted in 0 ≤ c ≤ 1 (MPS 2) |  | 7.2 / 38 / 108 |

## D6 — how often the constraints are active (`d2_d6_full.py`)

% of rows with an active constraint (hybrids: a simplex face, or for hybrid2 the triangle; quadratic/four-term: their exact positive–monotonic region). Power-2 (fitted without constraint): % of rows unphysical.

| law | PHOENIX r | PHOENIX μ | MPS 2 r | MPS 2 μ |
|---|---|---|---|---|
| hybrid5 | 0.54 | 8.74 | 0.00 | 0.00 |
| hybrid4 | 0.31 | 2.00 | 0.00 | 0.00 |
| hybrid2 | 0.06 | 0.20 | 0.00 | 0.00 |
| four-term | 6.38 | 5.78 | 1.32 | 5.77 |
| quadratic | 0.63 | 11.82 | 0.00 | 1.33 |
| power-2 (unphysical) | 21.79 | 27.93 | 3.93 | 3.11 |

## D2 — power-2 inside its physical region; D6b — what the four-term constraint costs

| rows | n | LC unconstrained | LC constrained | increase |
|---|---|---|---|---|
| power-2, PHOENIX, LM fit unphysical | 13926 | 4.0 / 83 / 189 | 17 / 180 / 357 | 13 / 113 / 208 |
| four-term, PHOENIX, constraint active | 4075 | 1.8 / 9.8 / 17 | 4.6 / 44 / 75 | 2.2 / 40 / 68 (min -5.8) |
| power-2, MPS 2, LM fit unphysical | 448 | 2.0 / 27 / 33 | 10 / 56 / 67 | 7.6 / 40 / 53 |
| four-term, MPS 2, constraint active | 151 | 2.6 / 4.9 / 5.0 | 2.8 / 6.0 / 6.4 | 0.03 / 5.1 / 5.2 (min -2.0) |

Refitted power-2 rows all end on c = 1; their LM fits have c = 1.000–76.276 (median 1.178, 90th percentile 1.56; the largest have α → 0, a log-like law). Four-term: a constrained fit can score better in light-curve terms, because least squares in intensity is not the light-curve metric.

## D3 — hybrid2 pole, held-out (`d3_hybrid2_heldout.py`; P2c data, 2000-point resample)

7524 model-bands (PHOENIX 4104, MPS-ATLAS Set 2 3420), ε ∈ 0.30–0.46, objective: 90th percentile of the worst light-curve error (ppm). 'held-out' = evaluated on the set not used for tuning.

| tune on | evaluate on | ε chosen | median | 90th | 99th | at ε = 0.36 (median / 90th / 99th) |
|---|---|---|---|---|---|---|
| phoenix | mps2 | 0.34 | 50.1 | 80.9 | 115.8 | 51.5 / 77.8 / 104.1 |
| mps2 | phoenix | 0.39 | 65.3 | 106.3 | 142.4 | 64.8 / 100.0 / 136.1 |
| all | all | 0.35 | 56.5 | 92.6 | 126.9 | 57.2 / 92.8 / 122.4 |

Random half/half splits by model (20×): chosen ε = 0.35–0.38 (median 0.35); held-out 90th percentile 93.3 ppm (range 88.1–96.1); excess over the held-out optimum median 0.47 ppm, max 3.87 ppm.

Sensitivity, all 7,524 model-bands:

| ε | median | 90th | 99th |
|---|---|---|---|
| 0.30 | 58.1 | 108.7 | 161.8 |
| 0.33 | 56.6 | 95.6 | 138.5 |
| 0.36 | 57.2 | 92.8 | 122.4 |
| 0.39 | 60.0 | 96.0 | 134.8 |
| 0.42 | 64.5 | 105.6 | 153.7 |

Quadratic on the same set: median 67.4, 90th 123.0, 99th 178.9.
* phoenix: hybrid2(0.36) better than quadratic in 68% of model-bands; quadratic median 62.5, hybrid2 64.8.
* mps2: hybrid2(0.36) better than quadratic in 79% of model-bands; quadratic median 70.6, hybrid2 51.5.


### hybrid2 pole re-scanned with the exact quadrature (`p2_exact.py pole`)

P2S, P2 light-curve metric (p ∈ {0.03, 0.1, 0.15}, b ∈ {0, …, 0.9}); median / 90th / 99th (ppm).

| ε | PHOENIX | MPS-ATLAS 2 | both |
|---|---|---|---|
| 0.24 | 62 / 160 / 207 | 58 / 137 / 189 | 61 / 147 / 204 |
| 0.28 | 55 / 129 / 168 | 49 / 107 / 155 | 53 / 119 / 168 |
| 0.32 | 54 / 102 / 139 | 46 / 85 / 128 | 50 / 97 / 138 |
| 0.36 | 56 / 93 / 115 | 49 / 72 / 106 | 52 / 90 / 115 |
| 0.4 | 57 / 99 / 123 | 59 / 71 / 92 | 58 / 94 / 121 |
| 0.44 | 62 / 115 / 150 | 69 / 80 / 105 | 64 / 108 / 147 |
| quadratic | 60 / 123 / 162 | 74 / 89 / 114 | 66 / 117 / 159 |

## D4 — planet size and grazing geometry (`t3_study.py`, `p2s_study.py`)

T3 cases (PHOENIX, b as listed; 'graze' = b = 1 − p/2). Median over the 15 cases of the worst |ΔF| (ppm); in brackets the maximum.

| law | p | b = 0 | b = 0.5 | b = 1 − p/2 | b = 1 | rms (b = 0) |
|---|---|---|---|---|---|---|
| hybrid5 | 0.03 | 0.30 (0.64) | 0.30 (0.64) | 0.18 (0.39) | 0.18 (0.39) | 0.08 |
|  | 0.1 | 0.84 (2.2) | 0.84 (2.2) | 0.40 (0.87) | 0.34 (0.74) | 0.50 |
|  | 0.2 | 2.6 (6.6) | 2.4 (5.1) | 0.58 (1.3) | 0.55 (1.1) | 1.4 |
|  | 0.3 | 4.7 (12) | 2.3 (4.8) | 1.7 (3.5) | 0.65 (1.5) | 2.6 |
| hybrid4 | 0.03 | 1.6 (3.7) | 1.6 (3.7) | 0.86 (1.9) | 0.80 (1.9) | 0.43 |
|  | 0.1 | 6.5 (15) | 6.5 (15) | 3.4 (7.7) | 1.8 (5.0) | 2.8 |
|  | 0.2 | 13 (28) | 13 (28) | 4.6 (11) | 4.6 (11) | 6.7 |
|  | 0.3 | 20 (42) | 9.9 (22) | 9.5 (20) | 5.4 (13) | 11 |
| four-term | 0.03 | 0.55 (0.96) | 0.56 (0.94) | 0.19 (0.33) | 0.19 (0.33) | 0.19 |
|  | 0.1 | 3.2 (5.4) | 3.2 (5.4) | 1.2 (2.0) | 1.1 (2.0) | 1.5 |
|  | 0.2 | 7.5 (12) | 6.7 (11) | 3.2 (5.7) | 1.7 (2.9) | 4.3 |
|  | 0.3 | 14 (22) | 7.8 (14) | 6.6 (11) | 3.6 (6.4) | 7.3 |
| power-2 | 0.03 | 0.80 (2.6) | 0.81 (2.7) | 0.66 (1.3) | 0.38 (0.86) | 0.23 |
|  | 0.1 | 4.5 (15) | 4.5 (15) | 3.4 (12) | 2.3 (8.9) | 1.7 |
|  | 0.2 | 8.4 (35) | 8.4 (35) | 6.2 (26) | 6.2 (23) | 4.7 |
|  | 0.3 | 15 (64) | 15 (64) | 8.0 (34) | 8.0 (34) | 7.8 |
| hybrid2 | 0.03 | 13 (21) | 13 (21) | 9.5 (14) | 9.5 (14) | 3.6 |
|  | 0.1 | 65 (97) | 66 (97) | 41 (60) | 20 (31) | 22 |
|  | 0.2 | 98 (145) | 98 (145) | 74 (111) | 70 (104) | 52 |
|  | 0.3 | 172 (256) | 172 (256) | 92 (137) | 92 (137) | 76 |
| quadratic | 0.03 | 10 (15) | 11 (16) | 8.1 (12) | 8.1 (12) | 3.8 |
|  | 0.1 | 76 (113) | 76 (113) | 41 (60) | 18 (27) | 30 |
|  | 0.2 | 146 (217) | 146 (217) | 113 (167) | 85 (125) | 85 |
|  | 0.3 | 273 (404) | 271 (400) | 145 (214) | 144 (212) | 144 |
| poly2 (2 terms) | 0.03 | 36 (64) | 36 (65) | 36 (65) | 33 (59) | 11 |
|  | 0.1 | 203 (375) | 203 (375) | 103 (184) | 103 (184) | 86 |
|  | 0.2 | 445 (805) | 446 (806) | 240 (444) | 157 (277) | 233 |
|  | 0.3 | 679 (1222) | 504 (929) | 368 (666) | 246 (454) | 382 |
| poly4 (4 terms) | 0.03 | 14 (24) | 14 (24) | 8.8 (15) | 8.8 (15) | 4.0 |
|  | 0.1 | 51 (87) | 51 (87) | 27 (47) | 19 (34) | 25 |
|  | 0.2 | 97 (188) | 92 (167) | 34 (61) | 34 (61) | 55 |
|  | 0.3 | 213 (419) | 76 (145) | 41 (72) | 40 (72) | 92 |
| poly5 (5 terms) | 0.03 | 10 (17) | 11 (17) | 5.0 (8.1) | 5.0 (8.1) | 2.9 |
|  | 0.1 | 29 (49) | 29 (49) | 13 (22) | 13 (22) | 17 |
|  | 0.2 | 71 (120) | 45 (76) | 21 (35) | 17 (29) | 36 |
|  | 0.3 | 185 (314) | 55 (92) | 25 (41) | 25 (41) | 64 |

Scaling with p (T3, b = 0, median over cases of the per-case log-slope between p = 0.03 and 0.3): hybrid5 1.22, hybrid4 1.07, four-term 1.38, power-2 1.31, hybrid2 1.08, quadratic 1.42. The error grows roughly as p^1.1–1.5, not p².

P2S (all 400 model-bands; b ∈ {0, 0.3, 0.6, 0.9} unless 'graze'): median / 99th percentile.

| law | grid | p = 0.03 | p = 0.1 | p = 0.2 | p = 0.3 | graze p = 0.1 | graze 0.2 | graze 0.3 |
|---|---|---|---|---|---|---|---|---|
| hybrid5 | phoenix | 0.25 / 0.64 | 1.2 / 5.3 | 2.4 / 15 | 4.3 / 27 | 0.42 / 1.4 | 0.80 / 3.3 | 1.6 / 8.1 |
|  | mps2 | 0.17 / 0.44 | 0.79 / 4.2 | 1.7 / 12 | 2.9 / 22 | 0.22 / 1.00 | 0.43 / 2.6 | 1.0 / 6.6 |
| hybrid4 | phoenix | 1.1 / 3.4 | 4.3 / 13 | 7.5 / 26 | 11 / 39 | 2.0 / 6.9 | 2.9 / 9.4 | 4.6 / 19 |
|  | mps2 | 1.9 / 3.0 | 7.1 / 12 | 13 / 23 | 21 / 35 | 3.6 / 6.2 | 4.9 / 8.4 | 9.7 / 17 |
| four-term | phoenix | 0.62 / 1.9 | 3.7 / 10 | 8.6 / 23 | 16 / 42 | 1.3 / 3.8 | 3.7 / 11 | 7.5 / 20 |
|  | mps2 | 0.13 / 0.65 | 0.84 / 3.8 | 1.8 / 9.1 | 3.3 / 17 | 0.33 / 1.5 | 0.62 / 2.8 | 1.4 / 7.7 |
| power-2 | phoenix | 1.4 / 7.8 | 8.4 / 50 | 19 / 113 | 34 / 207 | 5.6 / 37 | 14 / 84 | 18 / 109 |
|  | mps2 | 1.5 / 4.0 | 6.5 / 22 | 16 / 46 | 27 / 84 | 5.5 / 16 | 11 / 34 | 15 / 44 |
| hybrid2 | phoenix | 9.5 / 21 | 41 / 95 | 76 / 168 | 117 / 307 | 24 / 61 | 42 / 114 | 61 / 152 |
|  | mps2 | 9.9 / 14 | 43 / 69 | 72 / 157 | 111 / 285 | 27 / 42 | 46 / 91 | 58 / 141 |
| quadratic | phoenix | 5.7 / 15 | 42 / 115 | 83 / 221 | 156 / 412 | 22 / 59 | 64 / 171 | 83 / 218 |
|  | mps2 | 6.8 / 11 | 51 / 80 | 100 / 154 | 187 / 287 | 25 / 41 | 77 / 119 | 99 / 152 |

A grazing transit never gives a larger error than a full transit of the same planet in these samples: the worst case for every law is the limb crossing of a full transit.

## D5 — dependence on the interpolant of the 'exact' profile (`p2s_study.py`, `d5b_interp_mu.py`)

First block: light-curve difference between alternative truths and the baseline (PCHIP in μ) through the same nodes. Second: the production fits re-scored against each truth.

| interpolant | PHOENIX: truth vs PCHIP | MPS 2: truth vs PCHIP |
|---|---|---|
| linear in μ | 2.7 / 7.4 / 8.0 | 2.5 / 4.1 / 4.2 |
| linear in r | 30 / 147 / 190 | 54 / 230 / 247 |
| Akima in r | 5.6 / 27 / 34 | 6.4 / 37 / 42 |
| Akima in μ | 0.41 / 1.1 / 1.2 | 0.12 / 0.18 / 0.21 |
| cubic spline in μ | 0.47 / 1.3 / 1.4 | 0.19 / 0.30 / 0.32 |

| law | grid | PCHIP μ (baseline) | linear μ | Akima μ | spline μ | Akima r | linear r |
|---|---|---|---|---|---|---|---|
| hybrid5 | phoenix | 1.2 / 5.3 / 5.4 | 3.1 / 8.1 / 8.7 | 1.3 / 5.2 / 5.3 | 1.4 / 5.2 / 5.4 | 5.9 / 28 / 36 | 30 / 148 / 192 |
|  | mps2 | 0.79 / 4.2 / 4.9 | 2.6 / 4.6 / 4.9 | 0.82 / 4.3 / 4.9 | 0.91 / 4.2 / 4.9 | 7.0 / 39 / 43 | 55 / 232 / 249 |
| hybrid4 | phoenix | 4.3 / 13 / 14 | 3.9 / 12 / 13 | 4.4 / 13 / 15 | 4.4 / 13 / 15 | 7.0 / 28 / 36 | 30 / 148 / 193 |
|  | mps2 | 7.1 / 12 / 13 | 8.0 / 13 / 14 | 7.1 / 12 / 13 | 7.1 / 12 / 13 | 7.6 / 39 / 43 | 53 / 232 / 248 |
| four-term | phoenix | 3.7 / 10 / 11 | 5.2 / 11 / 12 | 3.5 / 10.0 / 11 | 3.5 / 9.9 / 11 | 6.6 / 27 / 34 | 31 / 147 / 191 |
|  | mps2 | 0.84 / 3.8 / 5.2 | 2.6 / 4.9 / 6.3 | 0.85 / 3.8 / 5.2 | 0.83 / 3.9 / 5.3 | 6.3 / 38 / 42 | 54 / 230 / 247 |
| power-2 | phoenix | 8.4 / 50 / 53 | 9.9 / 50 / 54 | 8.6 / 50 / 54 | 8.7 / 50 / 54 | 10 / 51 / 54 | 31 / 145 / 191 |
|  | mps2 | 6.5 / 22 / 40 | 6.9 / 23 / 39 | 6.5 / 22 / 40 | 6.5 / 22 / 40 | 8.5 / 39 / 43 | 55 / 228 / 246 |
| quadratic | phoenix | 42 / 115 / 125 | 42 / 113 / 124 | 42 / 114 / 125 | 42 / 115 / 125 | 41 / 113 / 124 | 60 / 146 / 191 |
|  | mps2 | 51 / 80 / 83 | 50 / 79 / 81 | 51 / 81 / 83 | 51 / 81 / 83 | 52 / 81 / 83 | 67 / 227 / 244 |

Linear interpolation in r is a poor interpolant near the disc centre, where the nodes are sparse in r and I is quadratic in r; it is shown only as an extreme case.

## D7 — what r-weighting costs at the limb (`p2s_study.py`)

σ_μ = rms residual uniform in μ (Claret's merit function), σ_r = rms uniform in r; max |ΔI| in I/I(1) for μ < 0.05 and μ ≥ 0.05; medians over P2S. 'r' = r-weighted fit (tabulated), 'μ' = uniform-μ fit (secondary table).

| law | grid | σ_μ r / μ | σ_r r / μ | max|ΔI| μ<0.05 r / μ | max|ΔI| μ≥0.05 r / μ | LC error r / μ |
|---|---|---|---|---|---|---|
| hybrid5 | phoenix | 8.9e-03 / 7.5e-03 | 7.5e-04 / 1.0e-03 | 0.109 / 0.069 | 0.0038 / 0.0099 | 1.2 / 7.6 |
|  | mps2 | 7.9e-04 / 6.7e-04 | 1.8e-04 / 2.3e-04 | 0.005 / 0.003 | 0.0033 / 0.0026 | 0.79 / 2.2 |
| hybrid4 | phoenix | 1.7e-02 / 1.5e-02 | 1.4e-03 / 2.1e-03 | 0.195 / 0.157 | 0.0117 / 0.0317 | 4.3 / 15 |
|  | mps2 | 2.5e-03 / 2.4e-03 | 8.2e-04 / 9.5e-04 | 0.021 / 0.023 | 0.0077 / 0.0066 | 7.1 / 13 |
| four-term | phoenix | 9.5e-03 / 6.0e-03 | 7.5e-04 / 2.3e-03 | 0.098 / 0.058 | 0.0039 / 0.0063 | 3.7 / 33 |
|  | mps2 | 4.5e-04 / 3.1e-04 | 9.2e-05 / 1.6e-04 | 0.005 / 0.003 | 0.0007 / 0.0008 | 0.84 / 2.1 |
| power-2 | phoenix | 1.1e-02 / 9.1e-03 | 1.3e-03 / 2.5e-03 | 0.140 / 0.095 | 0.0068 / 0.0144 | 8.4 / 31 |
|  | mps2 | 3.2e-03 / 2.1e-03 | 6.9e-04 / 1.2e-03 | 0.024 / 0.011 | 0.0063 / 0.0038 | 6.5 / 11 |
| hybrid2 | phoenix | 3.8e-02 / 3.3e-02 | 5.1e-03 / 9.8e-03 | 0.331 / 0.287 | 0.0596 / 0.0338 | 41 / 139 |
|  | mps2 | 2.3e-02 / 1.9e-02 | 5.0e-03 / 7.8e-03 | 0.132 / 0.104 | 0.0664 / 0.0378 | 43 / 107 |
| quadratic | phoenix | 3.1e-02 / 2.5e-02 | 3.9e-03 / 8.9e-03 | 0.282 / 0.230 | 0.0361 / 0.0214 | 42 / 111 |
|  | mps2 | 1.5e-02 / 1.0e-02 | 3.8e-03 / 6.3e-03 | 0.087 / 0.058 | 0.0409 / 0.0160 | 51 / 83 |

## D12 — where the light-curve error comes from (`t3_study.py`, `p2s_study.py`)

The profile error ΔI = I_law − I_model is split at μ = 0.05; each part is scored on its own (law below 0.05 and model above, and the reverse). Ratios to the full error (median; T3 at b = 0, P2S on the production grid).

| law | T3: μ<0.05 part | T3: μ≥0.05 part | P2S PHOENIX: μ<0.05 | μ≥0.05 | P2S MPS 2: μ<0.05 | μ≥0.05 |
|---|---|---|---|---|---|---|
| hybrid5 | 0.24 | 1.03 | 0.36 | 1.02 | 0.17 | 1.01 |
| hybrid4 | 0.04 | 1.00 | 0.05 | 1.00 | 0.05 | 0.99 |
| four-term | 0.04 | 1.01 | 0.04 | 1.01 | 0.09 | 0.97 |
| power-2 | 0.10 | 1.04 | 0.09 | 1.01 | 0.14 | 1.08 |
| hybrid2 | 0.18 | 0.95 | 0.20 | 0.95 | 0.19 | 0.95 |
| quadratic | 0.12 | 0.96 | 0.12 | 0.96 | 0.10 | 0.97 |
| poly2 (2 terms) | 0.08 | 0.98 | 0.08 | 0.98 | 0.08 | 0.98 |
| poly4 (4 terms) | 0.21 | 1.07 | 0.22 | 1.07 | 0.20 | 1.06 |
| poly5 (5 terms) | 0.31 | 1.10 | 0.32 | 1.10 | 0.29 | 1.09 |

M dwarfs (T3, quadratic): full 72–114 ppm; the part from μ < 0.05 alone gives 8–13 ppm, the part from μ ≥ 0.05 67–109 ppm. The quadratic deficit lies in the disc interior: a law that cannot bend at the limb is pulled away from the profile over the whole disc. hybrid5 has the largest share of its (small) error at μ < 0.05.

## D13 — photon vs energy integration (`d13_photon_energy.py`)

Five T3 stars (PHOENIX). 'profile' = light-curve difference between the energy- and photon-weighted model profiles; the fitted laws carry the same difference (to within a few ppm). max |Δw|: largest coefficient shift.

| band | profile LC difference (ppm), min–max over stars | max |Δw| hybrid5 | four-term | quadratic | power-2 |
|---|---|---|---|---|---|
| TESS | 44–69 | 0.007 | 0.095 | 0.016 | 0.018 |
| Kepler | 69–118 | 0.009 | 0.177 | 0.027 | 0.028 |
| HST_WFC3_G141 | 1.3–22 | 0.009 | 0.087 | 0.008 | 0.018 |
| JWST_NIRSpec_PRISM | 122–258 | 0.023 | 0.329 | 0.065 | 0.082 |
| 2MASS_J | 1.4–2.7 | 0.001 | 0.007 | 0.001 | 0.004 |

## D14 — narrow bins summed from the coefficient spectra (`d14_narrow_bins.py`)

150 random PHOENIX models; top-hat bins at resolving power R over 0.6–5.3 µm. Percentage of bins whose summed coefficients are unphysical (positivity or monotonicity fails); in brackets the percentage of models with at least one such bin.

| R | law | M (<3900 K) | K (3900–5200 K) | FG (>5200 K) | wavelengths of unphysical bins |
|---|---|---|---|---|---|
| 100 | hybrid2 | 0.25 (13) | 0.00 (0) | 0.00 (0) | 0.61–0.71 µm (90%) |
| 100 | hybrid4 | 0.51 (15) | 0.00 (0) | 0.00 (0) | 0.61–0.82 µm (90%) |
| 100 | hybrid5 | 3.64 (40) | 0.74 (13) | 0.00 (0) | 0.62–1.79 µm (90%) |
| 300 | hybrid2 | 0.26 (13) | 0.00 (0) | 0.00 (0) | 0.61–0.73 µm (90%) |
| 300 | hybrid4 | 0.56 (19) | 0.00 (0) | 0.00 (0) | 0.61–0.85 µm (90%) |
| 300 | hybrid5 | 3.59 (43) | 0.75 (13) | 0.00 (0) | 0.62–1.79 µm (90%) |
| 1000 | hybrid2 | 0.27 (17) | 0.00 (0) | 0.00 (0) | 0.61–0.73 µm (90%) |
| 1000 | hybrid4 | 0.62 (23) | 0.00 (0) | 0.00 (0) | 0.61–0.85 µm (90%) |
| 1000 | hybrid5 | 3.68 (43) | 0.78 (19) | 0.00 (0) | 0.62–1.80 µm (90%) |

Models: 47 M, 31 K, 72 FG. Unphysical narrow bins occur only below ~2.5 µm, in cool stars (molecular bands); none for FG stars.

## Table 1 analogue — fit weighting with the exact quadrature (`p2_exact.py weighting`)

P2S; P2 light-curve metric (p ∈ {0.03, 0.1, 0.15}, b ∈ {0, …, 0.9}); median / worst (ppm). 'nodes' = equal weights on the native μ nodes (Claret 2000; Kostogryz et al. 2022).

**PHOENIX (34 × 8)**

| weighting | hybrid5 | hybrid4 | four-term | power-2 | hybrid2 | quadratic |
|---|---|---|---|---|---|---|
| equal per native node | 43 / 333 | 127 / 348 | 75 / 204 | 247 / 582 | 465 / 1270 | 289 / 914 |
| uniform in μ | 13 / 33 | 25 / 70 | 55 / 123 | 56 / 166 | 213 / 461 | 182 / 422 |
| uniform in r | 1.9 / 9.7 | 6.5 / 24 | 6.3 / 20 | 12 / 73 | 56 / 146 | 60 / 176 |
| uniform in area | 2.3 / 12 | 7.1 / 27 | 7.3 / 23 | 12 / 75 | 58 / 133 | 58 / 168 |

Uniform-μ → uniform-r improvement: hybrid5 ×6.5 (median), ×3.4 (worst); four-term ×8.8 (median), ×6.3 (worst); power-2 ×4.6 (median), ×2.3 (worst); quadratic ×3.0 (median), ×2.4 (worst).

**MPS-ATLAS Set 2 (16 × 8)**

| weighting | hybrid5 | hybrid4 | four-term | power-2 | hybrid2 | quadratic |
|---|---|---|---|---|---|---|
| equal per native node | 8.1 / 23 | 15 / 84 | 8.7 / 33 | 60 / 326 | 326 / 499 | 256 / 456 |
| uniform in μ | 4.0 / 16 | 20 / 50 | 3.7 / 16 | 24 / 102 | 170 / 264 | 140 / 247 |
| uniform in r | 1.4 / 8.8 | 12 / 22 | 1.3 / 9.2 | 9.8 / 52 | 49 / 111 | 74 / 116 |
| uniform in area | 1.7 / 10 | 13 / 25 | 1.5 / 10.0 | 10 / 54 | 54 / 101 | 70 / 111 |

Uniform-μ → uniform-r improvement: hybrid5 ×2.9 (median), ×1.8 (worst); four-term ×2.7 (median), ×1.8 (worst); power-2 ×2.4 (median), ×1.9 (worst); quadratic ×1.9 (median), ×2.1 (worst).

## Limb sensitivity with the exact quadrature (`p2_exact.py limb`)

PHOENIX P2S (34 × 8); baseline M1 per band, r-weighted. Light-curve difference (ppm, P2 metric) between the law fitted to the variant and to the baseline: median / 90th percentile. 'MPS-ATLAS 2 vs PHOENIX' and 'Set 1 vs Set 2' where the MPS-ATLAS grid has the model.

| variant | hybrid5 | hybrid4 | four-term | hybrid2 |
|---|---|---|---|---|
| M2 instead of M1 | 0.86 / 1.4 | 0.85 / 1.5 | 1.6 / 2.7 | 1.00 / 1.8 |
| one r_c per model | 1.2 / 4.8 | 1.2 / 4.8 | 2.5 / 7.9 | 1.4 / 5.6 |
| μ nodes of ATLAS (17) | 0.55 / 1.1 | 0.66 / 1.2 | 1.8 / 3.3 | 1.6 / 2.3 |
| μ nodes of MPS-ATLAS (24) | 0.86 / 1.6 | 0.79 / 1.4 | 2.2 / 4.1 | 1.8 / 2.3 |
| μ nodes of Stagger (10) | 2.2 / 3.4 | 2.7 / 4.1 | 3.5 / 5.6 | 1.4 / 3.1 |
| MPS-ATLAS 2 vs PHOENIX | 172 / 446 | 172 / 447 | 172 / 446 | 173 / 454 |
| MPS-ATLAS Set 1 vs Set 2 | 32 / 128 | 32 / 128 | 32 / 128 | 32 / 128 |

hybrid5 median |Δw_j| (w1…w5): M2 0.000, 0.001, 0.012, 0.001, 0.001; MPS-ATLAS nodes 0.000, 0.001, 0.017, 0.003, 0.002; baseline median w 0.051, 0.090, 0.139, 0.099, 0.174.

## D12b — is the M-dwarf quadratic deficit caused by the steep limb? (`d12b_limb_excluded.py`)

Each law fitted (r-weighted) to μ ≥ 0.05 only; 'interior' = light-curve error with the model profile kept below μ = 0.05 (only the disc interior compared); 'full' = against the full model. T3, b = 0, p = 0.1. Range over cases (ppm).

| law | interior, M dwarfs | interior, K/G | full, all |
|---|---|---|---|
| quadratic | 49–86 | 15–36 | 16–89 |
| power-2 | 0.24–7.1 | 3.4–13 | 0.11–14 |
| four-term | 0.26–1.2 | 0.22–0.47 | 0.31–1.8 |
| hybrid2 | 33–58 | 15–33 | 21–64 |
| hybrid4 | 8.3–15 | 7.2–9.4 | 9.6–22 |
| hybrid5 | 0.54–2.5 | 1.1–1.9 | 0.92–6.9 |

## D12c — fitting without the extreme limb (`d12c_mu_min.py`)

r-weighted fits to μ ≥ μ_min only, scored against the full model light curve (production grid). P2S; median / 99th (ppm).

| law | grid | μ_min = 0 (tabulated) | 0.02 | 0.05 | 0.1 |
|---|---|---|---|---|---|
| hybrid5 | phoenix | 1.2 / 5.3 | 1.4 / 5.4 | 4.0 / 14 | 4.8 / 14 |
|  | mps2 | 0.79 / 4.2 | 0.76 / 4.1 | 1.5 / 4.3 | 3.1 / 9.8 |
| hybrid4 | phoenix | 4.3 / 13 | 5.0 / 14 | 10 / 24 | 11 / 51 |
|  | mps2 | 7.1 / 12 | 7.2 / 12 | 6.8 / 11 | 9.8 / 17 |
| four-term | phoenix | 3.7 / 10 | 1.3 / 4.8 | 1.3 / 4.5 | 1.3 / 4.6 |
|  | mps2 | 0.84 / 3.8 | 0.70 / 3.7 | 0.55 / 3.0 | 0.64 / 3.6 |
| power-2 | phoenix | 8.4 / 50 | 8.1 / 50 | 7.4 / 44 | 7.3 / 44 |
|  | mps2 | 6.5 / 22 | 6.0 / 22 | 5.1 / 20 | 4.6 / 17 |
| hybrid2 | phoenix | 41 / 95 | 39 / 91 | 32 / 70 | 32 / 77 |
|  | mps2 | 43 / 69 | 42 / 68 | 37 / 64 | 29 / 63 |
| quadratic | phoenix | 42 / 115 | 40 / 112 | 32 / 90 | 31 / 86 |
|  | mps2 | 51 / 80 | 50 / 79 | 48 / 73 | 42 / 64 |

## E1 — parameter bias, noiseless (`e1_bias.py`)

T3 stars and bands; truth p = 0.1, b = 0.3, chord velocity v = 1 (600 points). Fit of (p, b, v) with coefficients fixed at the tabulated values, or free (quadratic/four-term unconstrained; hybrid2 triangle and hybrid4/5 simplex via q; power-2 with 0 ≤ c ≤ 1). Δ(p²) in ppm of a 10,000 ppm depth.

| law | mode | Δ(p²) ppm: median |·| / max |·| | Δb median |·| / max | Δv/v max |·| | rms resid. ppm (median / max) | chromatic spread of Δ(p²) per star (max over stars) |
|---|---|---|---|---|---|---|
| hybrid5 | fixed | 0.08 / 0.29 | 0.0001 / 0.0004 | 0.0001 | 0.49 / 1.2 | 0.54 |
|  | free | 1.9 / 6.0 | 0.0000 / 0.0003 | 0.0002 | 0.12 / 0.30 | 4.3 |
| hybrid4 | fixed | 0.49 / 1.4 | 0.0012 / 0.0018 | 0.0006 | 2.6 / 5.7 | 0.88 |
|  | free | 10 / 34 | 0.0002 / 0.0007 | 0.0016 | 0.36 / 0.80 | 24 |
| four-term | fixed | 0.54 / 1.1 | 0.0001 / 0.0003 | 0.0001 | 1.5 / 2.4 | 1.3 |
|  | free | 1.3 / 1.8 | 0.0000 / 0.0000 | 0.0001 | 0.03 / 0.05 | 3.2 |
| power-2 | fixed | 2.1 / 8.9 | 0.0010 / 0.0051 | 0.0015 | 1.2 / 4.4 | 10 |
|  | free | 8.3 / 58 | 0.0019 / 0.0140 | 0.0027 | 1.3 / 6.4 | 49 |
| hybrid2 | fixed | 22 / 60 | 0.0226 / 0.0368 | 0.0105 | 14 / 20 | 29 |
|  | free | 29 / 50 | 0.0025 / 0.0104 | 0.0041 | 3.6 / 8.2 | 32 |
| quadratic | fixed | 31 / 97 | 0.0268 / 0.0534 | 0.0147 | 22 / 30 | 53 |
|  | free | 41 / 113 | 0.0079 / 0.0161 | 0.0007 | 4.8 / 7.0 | 57 |

Per case (Δ(p²) ppm, fixed / free):

| case | hybrid5 | hybrid4 | four-term | power-2 | hybrid2 | quadratic |
|---|---|---|---|---|---|---|
| late M|TESS | -0.2 / -2.4 | -1.1 / -10.9 | +0.7 / -1.3 | +8.6 / +58.0 | -60.2 / +49.1 | -96.7 / +112.8 |
| late M|WFC3 G141 | -0.0 / +1.9 | -0.4 / -33.8 | -0.3 / +1.8 | +2.1 / +29.6 | -41.7 / +50.0 | -57.6 / +71.1 |
| late M|NIRSpec PRISM | -0.1 / +1.4 | -0.3 / -17.7 | -0.6 / +1.6 | +1.1 / +8.7 | -31.4 / +38.6 | -44.1 / +55.5 |
| mid M|TESS | -0.0 / -0.3 | -1.0 / -34.5 | +0.0 / +1.1 | +1.9 / +6.9 | -44.3 / +44.1 | -64.1 / +74.3 |
| mid M|WFC3 G141 | -0.0 / +1.0 | -0.1 / -14.1 | -0.4 / +1.6 | +0.6 / +10.4 | -26.6 / +35.7 | -37.0 / +48.6 |
| mid M|NIRSpec PRISM | -0.1 / +0.8 | -0.1 / -10.5 | -0.6 / +1.3 | +0.4 / +0.2 | -21.7 / +29.5 | -30.8 / +40.9 |
| early M|TESS | -0.1 / +0.1 | -1.0 / -21.7 | -0.0 / +0.1 | +2.1 / -2.9 | -34.9 / +33.7 | -51.6 / +60.7 |
| early M|WFC3 G141 | -0.0 / +0.5 | -0.2 / -10.6 | -0.5 / +1.4 | -1.0 / +1.3 | -22.4 / +29.5 | -31.7 / +41.4 |
| early M|NIRSpec PRISM | -0.1 / +0.8 | -0.1 / -9.2 | -0.5 / +1.0 | -0.2 / -0.1 | -19.5 / +26.2 | -27.9 / +37.2 |
| K dwarf|TESS | +0.3 / -5.9 | -1.4 / -2.3 | -0.1 / +1.0 | -2.6 / -1.7 | -11.6 / -21.0 | -21.5 / +19.8 |
| K dwarf|WFC3 G141 | -0.3 / -4.5 | -0.7 / -6.0 | -1.1 / +1.5 | +7.6 / -12.9 | -17.2 / +11.1 | -30.6 / +35.7 |
| K dwarf|NIRSpec PRISM | -0.2 / -4.1 | -0.7 / -5.4 | -0.7 / +1.1 | +1.7 / -4.0 | -15.5 / +11.0 | -25.9 / +30.6 |
| Sun analogue|TESS | -0.1 / -6.0 | -1.1 / -2.1 | -0.3 / +1.1 | +8.9 / -14.9 | -7.0 / -11.2 | -18.9 / +19.2 |
| Sun analogue|WFC3 G141 | -0.2 / -4.7 | -0.5 / -4.1 | -0.9 / +1.5 | +8.6 / -14.9 | -12.3 / +10.7 | -22.5 / +27.9 |
| Sun analogue|NIRSpec PRISM | -0.2 / -5.1 | -0.5 / -3.4 | -0.6 / +1.1 | +4.7 / -8.3 | -10.2 / +7.4 | -18.9 / +23.0 |

