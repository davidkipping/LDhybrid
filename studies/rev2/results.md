# Round-2 revision computations (studies/rev2)

All fits with v1 settings (exact quadrature). Production light-curve metric (p = 0.1, b ∈ {0, 0.3, 0.6, 0.9}) unless stated. Scripts alongside; JSON outputs per item.

## C1 — limb rescaling (geometric vs Sing et al. linear), light-curve difference (ppm)

* T3 (15 cases), lin: median 99, 90th 131, min 63, max 151
* T3 (15 cases), lin4: median 99, 90th 131, min 63, max 151
  same edge node: max |Δμ_cri| = 0.00051
* P2S PHOENIX (272), lin: median 100, 90th 212, min 30, max 434
* P2S PHOENIX (272), lin4: median 100, 90th 212, min 30, max 434
  same edge node: max |Δμ_cri| = 0.0011

T3 per case (lin): late M TESS 129; late M HST_WFC3_G141 88; late M JWST_NIRSpec_PRISM 72; mid M TESS 108; mid M HST_WFC3_G141 72; mid M JWST_NIRSpec_PRISM 63; early M TESS 102; early M HST_WFC3_G141 70; early M JWST_NIRSpec_PRISM 64; K dwarf TESS 151; K dwarf HST_WFC3_G141 113; K dwarf JWST_NIRSpec_PRISM 100; Sun analogue TESS 132; Sun analogue HST_WFC3_G141 99; Sun analogue JWST_NIRSpec_PRISM 91

## C5 — power-2 physicality and uniform-μ constraint fractions

**phoenix**: power-2 unphysical r 21.5%, μ 28.8%
* by Teff (r): 2300-3000 K 74%, 3000-3500 K 69%, 3500-4000 K 30%, 4000-99999 K 1%
* c of unphysical rows (r): median 1.181, 90th 1.58, max 76.3, c > 2 in 3.4% (n = 14447)
* < 3500 K by band (r): 2MASS_H 96%, JWST_NIRSpec_G140M_F100LP 95%, HST_WFC3_G141 95%, JWST_NIRSpec_G140H_F100LP 95%, 2MASS_J 94%, JWST_NIRISS_SOSS_o1 94%, JWST_NIRSpec_G140H_F070LP 93%, JWST_NIRSpec_G140M_F070LP 93%, Sloan_z 93%, HST_WFC3_G102 92%, Cousins_I 82%, JWST_NIRSpec_PRISM 81%, Sloan_i 81%, TESS 80%, JWST_NIRISS_SOSS_o2 77%, NGTS 74%, HST_STIS_G750L 74%, 2MASS_Ks 73%, Kepler 71%, Sloan_g 70%, HST_STIS_G430L 70%, CHEOPS 70%, JWST_NIRSpec_G395H_NRS2 69%, JWST_NIRCam_F444W 69%, Johnson_B 67%, HST_WFC3_G280n1 66%, HST_WFC3_G280p1 66%, Spitzer_IRAC1_full 63%, Spitzer_IRAC1 62%, Johnson_V 61%, JWST_NIRSpec_G235H 58%, Spitzer_IRAC2_full 57%, Spitzer_IRAC2 56%, JWST_NIRSpec_G235M 56%, Cousins_R 50%, Sloan_r 50%, JWST_NIRSpec_G395H 48%, JWST_NIRSpec_G395M 48%, JWST_NIRCam_F322W2 40%, JWST_NIRSpec_G395H_NRS1 37%
* bound fractions (r): hybrid5 0.51%, hybrid4 0.30%, hybrid2 0.06%, claret4 6.06%, quadratic 0.60%
* by Teff (mu): 2300-3000 K 69%, 3000-3500 K 71%, 3500-4000 K 47%, 4000-99999 K 10%
* c of unphysical rows (mu): median 1.199, 90th 15.61, max 591.4, c > 2 in 22.6% (n = 19408)
* bound fractions (mu): hybrid5 10.09%, hybrid4 2.28%, hybrid2 0.24%, claret4 6.02%, quadratic 12.62%

**mps1**: power-2 unphysical r 3.8%, μ 2.9%
* by Teff (r): 2300-3000 K 0%, 3000-3500 K 0%, 3500-4000 K 29%, 4000-99999 K 0%
* c of unphysical rows (r): median 1.139, 90th 2.27, max 117.4, c > 2 in 12.3% (n = 26323)
* bound fractions (r): hybrid5 0.00%, hybrid4 0.00%, hybrid2 0.00%, claret4 1.09%, quadratic 0.00%
* by Teff (mu): 2300-3000 K 0%, 3000-3500 K 0%, 3500-4000 K 23%, 4000-99999 K 0%
* c of unphysical rows (mu): median 1.117, 90th 1.90, max 111.7, c > 2 in 9.0% (n = 20387)
* bound fractions (mu): hybrid5 0.00%, hybrid4 0.00%, hybrid2 0.00%, claret4 4.26%, quadratic 1.30%

**mps2**: power-2 unphysical r 3.4%, μ 2.7%
* by Teff (r): 2300-3000 K 0%, 3000-3500 K 0%, 3500-4000 K 27%, 4000-99999 K 0%
* c of unphysical rows (r): median 1.137, 90th 2.34, max 106.2, c > 2 in 13.0% (n = 23989)
* bound fractions (r): hybrid5 0.00%, hybrid4 0.00%, hybrid2 0.00%, claret4 1.37%, quadratic 0.00%
* by Teff (mu): 2300-3000 K 0%, 3000-3500 K 0%, 3500-4000 K 22%, 4000-99999 K 0%
* c of unphysical rows (mu): median 1.111, 90th 1.92, max 113.6, c > 2 in 9.2% (n = 19125)
* bound fractions (mu): hybrid5 0.00%, hybrid4 0.00%, hybrid2 0.00%, claret4 5.01%, quadratic 1.40%

T3 power-2 (tabulated, r): late M|TESS c=1.131 UNPHYSICAL; late M|WFC3 G141 c=1.483 UNPHYSICAL; late M|NIRSpec PRISM c=1.150 UNPHYSICAL; mid M|TESS c=1.030 UNPHYSICAL; mid M|WFC3 G141 c=1.335 UNPHYSICAL; mid M|NIRSpec PRISM c=1.013 UNPHYSICAL; early M|TESS c=0.866; early M|WFC3 G141 c=0.914; early M|NIRSpec PRISM c=0.823; K dwarf|TESS c=0.669; K dwarf|WFC3 G141 c=0.679; K dwarf|NIRSpec PRISM c=0.595; Sun analogue|TESS c=0.630; Sun analogue|WFC3 G141 c=0.614; Sun analogue|NIRSpec PRISM c=0.535

## C6 — interpolation error vs grid spacing (production metric, ppm)

PHOENIX Teff (native step h = 100 K), 1500 node-bands:

| law | gap 2h (200 K) med / 90th | gap 3h (300 K) med / 90th | k (median ratio) | native h (extrapolated) med / 90th |
|---|---|---|---|---|
| hybrid5 | 4.8 / 26.6 | 8.4 / 42.9 | 1.34 (90th: 1.18) | 1.9 / 11.8 |
| claret4 | 6.7 / 26.7 | 9.3 / 42.7 | 0.82 (90th: 1.16) | 3.8 / 12.0 |
| power2 | 14.6 / 48.0 | 19.3 / 72.5 | 0.69 (90th: 1.02) | 9.0 / 23.8 |

MPS-ATLAS Set 2 [M/H] (native step 0.05 dex in −1…+0.5), 600 node-bands:

| law | gap 0.1 dex | gap 0.2 dex | gap 0.3 dex | k (0.1→0.2) |
|---|---|---|---|---|
| hybrid5 | 1.10 / 3.4 | 1.32 / 4.0 | 1.54 / 5.3 | 0.26 |
| claret4 | 1.42 / 3.5 | 1.57 / 4.1 | 1.77 / 5.5 | 0.14 |
| power2 | 7.10 / 17.4 | 7.18 / 17.6 | 7.44 / 19.1 | 0.02 |

Pure interpolation error (interpolated vs tabulated coefficients at the node), median / 90th:

| law | grid | gap 2h/0.1 | gap 3h/0.2 | (MPS 0.3) | k | native (extrap.) |
|---|---|---|---|---|---|---|
| hybrid5 | PHOENIX Teff | 4.2 / 26.5 | 8.2 / 42.9 | – | 1.67 | 1.3 / 11.6 |
| claret4 | PHOENIX Teff | 4.1 / 26.5 | 8.3 / 43.0 | – | 1.70 | 1.3 / 11.5 |
| power2 | PHOENIX Teff | 5.2 / 38.2 | 9.4 / 67.5 | – | 1.46 | 1.9 / 14.4 |
| hybrid5 | MPS [M/H] | 0.30 / 1.8 | 0.57 / 3.2 | 0.86 / 5.3 | 0.91 | (native = 0.1 gap) |
| claret4 | MPS [M/H] | 0.30 / 1.8 | 0.57 / 3.2 | 0.86 / 5.3 | 0.91 | (native = 0.1 gap) |
| power2 | MPS [M/H] | 0.29 / 1.9 | 0.53 / 3.7 | 0.83 / 7.3 | 0.88 | (native = 0.1 gap) |

Tabulated law error at the same PHOENIX nodes (median): hybrid5 1.3, claret4 3.4, power2 7.1

## C7 — convergence of the light-curve error (ppm)

45 model-bands (T3 + 30 PHOENIX P2S); change relative to the production quadrature (N_GL = 48, n_sub = 4, 60 points per half chord).

| variant | max |Δ lc error| over laws & cases | median |Δ| | max |Δ model curve| (b = 0) |
|---|---|---|---|
| GL96 | 0.004 | 0.0000 | 1.15e-02 |
| nsub16 | 0.004 | 0.0000 | 9.90e-03 |
| GL96+nsub16 | 0.004 | 0.0000 | 1.01e-02 |
| chord240 | 0.239 | 0.0017 | 0.00e+00 |
| chord960 | 0.246 | 0.0015 | 0.00e+00 |
| all | 0.246 | 0.0022 | 1.01e-02 |

Largest chord-sampling change: 0.25 ppm ((4600, 4.5, '+0.5') Kepler power2, production 12.28 ppm).
Relative change (all refinements): median 0.02%, max 2.0%.

## C8 — pole-ladder sensitivity (ln ε_k × s; P2S, production metric, median / 90th ppm)

| law | grid | s = 0.8 | s = 1 (nominal) | s = 1.2 |
|---|---|---|---|---|
| hybrid4 | phoenix | 2.11 / 3.5 | 4.31 / 10.2 | 7.52 / 16.4 |
| hybrid4 | mps2 | 1.51 / 3.6 | 7.09 / 10.1 | 10.15 / 13.0 |
| hybrid5 | phoenix | 0.97 / 2.6 | 1.17 / 3.2 | 2.33 / 4.2 |
| hybrid5 | mps2 | 1.22 / 2.5 | 0.79 / 2.6 | 2.79 / 4.7 |

Pole positions: hybrid4 s=0.8: 0.0286, 0.3858; hybrid4 s=1.0: 0.01176, 0.3041; hybrid4 s=1.2: 0.004837, 0.2397; hybrid5 s=0.8: 0.005788, 0.07808, 0.4683; hybrid5 s=1.0: 0.001597, 0.04128, 0.3874; hybrid5 s=1.2: 0.0004404, 0.02182, 0.3205

## C11 — flux conservation (FCM) on P2S

| law | grid | median / max |Δw| | median Δ lc (ppm) | max |Δ lc| | median |flux err| without FCM |
|---|---|---|---|---|---|
| hybrid5 | phoenix | 1.4e-05 / 7.2e-05 | +0.01 | 0.1 | 5.7e-06 |
| hybrid5 | mps2 | 8.2e-06 / 6.1e-05 | +0.00 | 0.1 | 3.4e-06 |
| hybrid4 | phoenix | 2.2e-05 / 9.2e-05 | +0.02 | 0.1 | 1.3e-05 |
| hybrid4 | mps2 | 4.7e-05 / 8e-05 | +0.06 | 0.1 | 2.8e-05 |
| hybrid2 | phoenix | 0.00047 / 0.0017 | +1.14 | 6.7 | 2.0e-04 |
| hybrid2 | mps2 | 0.00047 / 0.0013 | +1.12 | 4.2 | 2.1e-04 |
| claret4 | phoenix | 0.0034 / 0.011 | +0.08 | 0.2 | 2.9e-05 |
| claret4 | mps2 | 0.00065 / 0.0044 | +0.01 | 0.1 | 5.6e-06 |
| quadratic | phoenix | 0.002 / 0.0056 | +2.47 | 7.6 | 4.3e-04 |
| quadratic | mps2 | 0.0024 / 0.0038 | +2.98 | 4.9 | 5.2e-04 |

## C8b — ladder scale scan (P2S, production metric, median / 90th ppm)

| law | grid | s=0.5 | s=0.6 | s=0.7 | s=0.8 | s=0.9 | s=1.0 | s=1.1 |
|---|---|---|---|---|---|---|---|---|
| hybrid4 | phoenix | 12.21 / 18.6 | 8.76 / 13.8 | 5.14 / 8.0 | 2.11 / 3.5 | 2.15 / 5.7 | 4.31 / 10.2 | 6.17 / 13.9 |
| hybrid4 | mps2 | 10.74 / 12.8 | 6.67 / 8.3 | 2.76 / 4.6 | 1.51 / 3.6 | 4.30 / 7.1 | 7.09 / 10.1 | 9.10 / 12.0 |
| hybrid4 | all | 11.62 / 17.6 | 8.09 / 12.7 | 4.54 / 7.3 | 1.86 / 3.5 | 2.65 / 6.1 | 5.00 / 10.2 | 6.98 / 12.9 |
| hybrid5 | phoenix | 2.44 / 4.6 | 2.12 / 3.8 | 1.18 / 2.1 | 0.97 / 2.6 | 1.23 / 3.0 | 1.17 / 3.2 | 1.43 / 3.6 |
| hybrid5 | mps2 | 2.32 / 2.8 | 0.94 / 1.2 | 0.79 / 1.9 | 1.22 / 2.5 | 1.00 / 2.6 | 0.79 / 2.6 | 1.54 / 3.3 |
| hybrid5 | all | 2.36 / 4.3 | 1.69 / 3.6 | 1.10 / 2.0 | 1.11 / 2.6 | 1.11 / 2.9 | 1.06 / 3.1 | 1.44 / 3.4 |

Pole positions: hybrid4 s=0.5: 0.1085, 0.5514; hybrid4 s=0.6: 0.06955, 0.4895; hybrid4 s=0.7: 0.0446, 0.4346; hybrid4 s=0.8: 0.0286, 0.3858; hybrid4 s=0.9: 0.01834, 0.3425; hybrid4 s=1.0: 0.01176, 0.3041; hybrid4 s=1.1: 0.007543, 0.27; hybrid5 s=0.5: 0.03996, 0.2032, 0.6224; hybrid5 s=0.6: 0.02098, 0.1477, 0.5661; hybrid5 s=0.7: 0.01102, 0.1074, 0.5149; hybrid5 s=0.8: 0.005788, 0.07808, 0.4683; hybrid5 s=0.9: 0.00304, 0.05677, 0.426; hybrid5 s=1.0: 0.001597, 0.04128, 0.3874; hybrid5 s=1.1: 0.0008385, 0.03001, 0.3524

## C9 — metal-poor limb brightening (57 PHOENIX hybrid5 rows with lc > 15 ppm)

* M1: non-monotonic 57/57 (brightening at μ ≤ 0.135), constraint binds 57/57, hybrid5 lc median 23.2, max 50.7 ppm
* M2: non-monotonic 57/57 (brightening at μ ≤ 0.135), constraint binds 57/57, hybrid5 lc median 23.1, max 50.7 ppm
* in1: non-monotonic 57/57 (brightening at μ ≤ 0.134), constraint binds 57/57, hybrid5 lc median 23.0, max 50.6 ppm
* in2: non-monotonic 57/57 (brightening at μ ≤ 0.134), constraint binds 57/57, hybrid5 lc median 22.9, max 50.5 ppm

Profile detail (M1): the non-monotonic part runs from μ ≈ 0.016 (median lowest
μ of a decrease) to μ ≈ 0.11 (median; 10–90%: 0.044–0.13), with a rise towards
the limb of 0.083 (median; max 0.11) in I/I(1). It is the same for every edge
choice, so it is a feature of the model intensities inside the disc, not of
where the edge falls.

## C10 — MURaM corrections

(a) Light-curve size of the correction as applied to MPS-ATLAS Set 2 hybrid5 weights (grid points 3196–6807 K, all log g and [M/H]; 3000 point-bands), ppm, median / 90th / max:

| class | ssd | B100G | B200G |
|---|---|---|---|
| IR (>2.5 µm) (840) | 3 / 4 / 11 | 13 / 26 / 78 | 30 / 54 / 166 |
| near-IR (1–2.5 µm) (1035) | 5 / 11 / 15 | 38 / 76 / 100 | 81 / 157 / 210 |
| optical (<1 µm) (1125) | 13 / 25 / 42 | 87 / 199 / 287 | 186 / 401 / 603 |
| all (3000) | 6 / 18 / 42 | 34 / 133 / 287 | 78 / 273 / 603 |

(b) Constant vs linearly extrapolated ratio below μ = 0.1 (MURaM stars, hybrid5): light-curve difference between the two Δw, ppm (median / max), and the size of the correction:

| state | Δ (extrapolation) median / max | correction size median / max | max |Δ Δw| |
|---|---|---|---|
| ssd | 0.0 / 0.2 | 5 / 25 | 0.002 |
| B100G | 0.1 / 0.6 | 27 / 185 | 0.006 |
| B200G | 0.1 / 0.7 | 66 / 393 | 0.008 |

(c) MURaM stars: F3_MH_00 Teff(hydro) 6807 K, log g 4.0, [M/H] 0.0; G2_MH_00 Teff(hydro) 5778 K, log g 4.438, [M/H] 0.0; G2_MH_m10 Teff(hydro) 5750 K, log g 4.438, [M/H] -1.0; G2_MH_p05 Teff(hydro) 5780 K, log g 4.438, [M/H] 0.5; K0_MH_00 Teff(hydro) 5240 K, log g 4.609, [M/H] 0.0; K4_MH_00 Teff(hydro) 4543 K, log g 4.609, [M/H] 0.0; M0_MH_00 Teff(hydro) 3707 K, log g 4.826, [M/H] 0.0; M2_MH_00 Teff(hydro) 3518 K, log g 5.0, [M/H] 0.0; M4_MH_00 Teff(hydro) 3196 K, log g 5.0, [M/H] 0.0
Magnetic Δw are interpolated in Teff between the [M/H] = 0 stars (3196, 3518, 3707, 4543, 5240, 5778, 6807 K), clamped at the ends, and applied at every log g and [M/H] (predict.magnetic).

## C14 — grazing transits relative to depth (T3, median over 15 cases)

| p | b | depth (ppm) | hybrid5 err/depth | hybrid4 err/depth | claret4 err/depth | power2 err/depth | hybrid2 err/depth | quadratic err/depth |
|---|---|---|---|---|---|---|---|---|
| 0.1 | 0.00 | 11199 | 0.007% | 0.060% | 0.029% | 0.039% | 0.603% | 0.703% |
| 0.1 | 0.95 | 6725 | 0.006% | 0.047% | 0.017% | 0.053% | 0.572% | 0.574% |
| 0.2 | 0.00 | 44717 | 0.006% | 0.030% | 0.017% | 0.018% | 0.225% | 0.336% |
| 0.2 | 0.90 | 28883 | 0.002% | 0.015% | 0.011% | 0.022% | 0.251% | 0.381% |

Absolute worst errors (ppm, median): p=0.1 b=0.00: hybrid5 0.84, hybrid4 6.5, claret4 3.2, power2 4.5, hybrid2 65, quadratic 76; p=0.1 b=0.95: hybrid5 0.4, hybrid4 3.4, claret4 1.2, power2 3.4, hybrid2 41, quadratic 41; p=0.2 b=0.00: hybrid5 2.6, hybrid4 13, claret4 7.5, power2 8.4, hybrid2 98, quadratic 1.5e+02; p=0.2 b=0.90: hybrid5 0.58, hybrid4 4.6, claret4 3.2, power2 6.2, hybrid2 74, quadratic 1.1e+02

## C15 — Fisher analysis with the exact-quadrature fits

Set-up: P = 3 d, a/R★ = 10, one transit sampled over ±0.75 T14; JWST-like 100 ppm per 1-min point; TESS-like 1000 ppm per 2-min point; p ∈ {0.1, 0.15}, b ∈ {0.3, 0.7}; geometry (p, b, a/R★, t0) free; 4 star-bands (5800 K G141 and TESS, 3400 K PRISM, 4600 K SOSS o1).

| law | case | combinations with s.d. < 0.05, free geometry (min–max over cases) | median s.d. of 1st / 2nd / 3rd direction |
|---|---|---|---|
| hybrid4 | JWST-like | 1–2 | 0.003 / 0.032 / 0.244 |
| hybrid4 | TESS-like | 0–1 | 0.048 / 0.395 / 0.909 |
| hybrid5 | JWST-like | 1–2 | 0.003 / 0.032 / 0.189 |
| hybrid5 | TESS-like | 0–1 | 0.047 / 0.394 / 0.907 |

hybrid5 best-measured direction (free geometry), components (w1…w5), median over cases:
* JWST-like: 1st [0.7, 0.6, -0.0, 0.02, 0.35]; 2nd [-0.57, 0.36, -0.0, 0.06, 0.51]
* TESS-like: 1st [0.7, 0.6, -0.0, 0.02, 0.35]; 2nd [-0.57, 0.37, -0.0, 0.06, 0.5]

## C16 — rms light-curve deviation (T3, b = 0)

From studies/rev/results.md D4 (t3_study.py, exact-quadrature fits): median over
the 15 Table-3 cases of the rms of law − model over the transit chord (ppm),
alongside the worst |ΔF|.

| law | p = 0.1 rms | p = 0.1 worst | p = 0.03 rms | p = 0.2 rms | p = 0.3 rms |
|---|---|---|---|---|---|
| hybrid5 | 0.50 | 0.84 | 0.08 | 1.4 | 2.6 |
| hybrid4 | 2.8 | 6.5 | 0.43 | 6.7 | 11 |
| four-term | 1.5 | 3.2 | 0.19 | 4.3 | 7.3 |
| power-2 | 1.7 | 4.5 | 0.23 | 4.7 | 7.8 |
| hybrid2 | 22 | 65 | 3.6 | 52 | 76 |
| quadratic | 30 | 76 | 3.8 | 85 | 144 |

## E1b — parameter bias, equal footing (noiseless; see script docstring for the set-up)

### phoenix (15 cases)

| law | mode | |Δ(p²)| med / max (ppm) | |Δb| med / max | max |Δv/v| | rms med (ppm) | chromatic spread max (ppm) | binds |
|---|---|---|---|---|---|---|---|
| hybrid5 | fixed | 0.08 / 0.3 | 0.0001 / 0.0004 | 0.0001 | 0.49 | 0.5 | 0/15 |
| hybrid5 | free | 1.40 / 7.2 | 0.0000 / 0.0003 | 0.0003 | 0.12 | 2.9 | 0/15 |
| hybrid5 | prior | 0.34 / 1.7 | 0.0002 / 0.0003 | 0.0001 | 0.19 | 1.2 | 0/15 |
| hybrid4 | fixed | 0.49 / 1.4 | 0.0012 / 0.0018 | 0.0006 | 2.60 | 0.9 | 0/15 |
| hybrid4 | free | 5.95 / 33.8 | 0.0003 / 0.0018 | 0.0016 | 0.36 | 32.7 | 0/15 |
| hybrid4 | prior | 1.58 / 10.6 | 0.0016 / 0.0056 | 0.0017 | 1.36 | 9.6 | 0/15 |
| claret4 | fixed | 0.54 / 1.1 | 0.0001 / 0.0003 | 0.0001 | 1.47 | 1.3 | 0/15 |
| claret4 | free | 1.31 / 2.2 | 0.0000 / 0.0001 | 0.0001 | 0.03 | 1.2 | 1/15 |
| claret4 | prior | 1.45 / 3.5 | 0.0009 / 0.0016 | 0.0004 | 0.67 | 2.9 | 1/15 |
| power2 | fixed | 2.07 / 8.9 | 0.0010 / 0.0051 | 0.0015 | 1.16 | 10.2 | 0/15 |
| power2 | free | 8.26 / 58.0 | 0.0019 / 0.0140 | 0.0027 | 1.26 | 49.3 | 6/15 |
| power2 | prior | 2.72 / 61.6 | 0.0022 / 0.0151 | 0.0030 | 1.75 | 51.9 | 6/15 |
| hybrid2 | fixed | 21.72 / 60.2 | 0.0226 / 0.0368 | 0.0105 | 14.31 | 28.8 | 0/15 |
| hybrid2 | free | 29.48 / 50.0 | 0.0025 / 0.0104 | 0.0041 | 3.63 | 32.2 | 0/15 |
| hybrid2 | prior | 17.22 / 27.7 | 0.0064 / 0.0181 | 0.0072 | 5.02 | 29.7 | 0/15 |
| quadratic | fixed | 30.79 / 96.7 | 0.0268 / 0.0534 | 0.0147 | 21.79 | 52.6 | 0/15 |
| quadratic | free | 40.93 / 112.8 | 0.0079 / 0.0161 | 0.0007 | 4.78 | 57.3 | 0/15 |
| quadratic | prior | 23.79 / 57.9 | 0.0009 / 0.0018 | 0.0027 | 7.95 | 26.1 | 0/15 |

### mps2 (9 cases)

| law | mode | |Δ(p²)| med / max (ppm) | |Δb| med / max | max |Δv/v| | rms med (ppm) | chromatic spread max (ppm) | binds |
|---|---|---|---|---|---|---|---|
| hybrid5 | fixed | 0.08 / 0.3 | 0.0001 / 0.0005 | 0.0002 | 0.49 | 0.5 | 0/9 |
| hybrid5 | free | 1.33 / 6.5 | 0.0001 / 0.0003 | 0.0002 | 0.14 | 5.2 | 0/9 |
| hybrid5 | prior | 0.34 / 2.0 | 0.0001 / 0.0003 | 0.0001 | 0.19 | 1.5 | 0/9 |
| hybrid4 | fixed | 0.73 / 1.3 | 0.0012 / 0.0019 | 0.0006 | 2.47 | 1.0 | 0/9 |
| hybrid4 | free | 7.59 / 25.8 | 0.0003 / 0.0006 | 0.0012 | 0.31 | 14.6 | 0/9 |
| hybrid4 | prior | 1.85 / 4.7 | 0.0018 / 0.0033 | 0.0010 | 1.35 | 3.2 | 0/9 |
| claret4 | fixed | 0.08 / 0.2 | 0.0000 / 0.0002 | 0.0001 | 0.24 | 0.3 | 0/9 |
| claret4 | free | 0.53 / 1.6 | 0.0000 / 0.0000 | 0.0001 | 0.03 | 1.2 | 0/9 |
| claret4 | prior | 0.31 / 1.4 | 0.0002 / 0.0006 | 0.0001 | 0.15 | 1.6 | 0/9 |
| power2 | fixed | 4.28 / 11.2 | 0.0019 / 0.0055 | 0.0016 | 1.97 | 7.0 | 0/9 |
| power2 | free | 5.64 / 14.7 | 0.0005 / 0.0022 | 0.0003 | 0.41 | 9.0 | 0/9 |
| power2 | prior | 3.95 / 8.4 | 0.0018 / 0.0047 | 0.0014 | 1.55 | 9.4 | 0/9 |
| hybrid2 | fixed | 16.60 / 37.0 | 0.0141 / 0.0274 | 0.0079 | 10.70 | 16.1 | 0/9 |
| hybrid2 | free | 20.90 / 35.0 | 0.0016 / 0.0113 | 0.0030 | 3.78 | 38.5 | 0/9 |
| hybrid2 | prior | 10.73 / 30.7 | 0.0061 / 0.0127 | 0.0041 | 4.40 | 33.2 | 0/9 |
| quadratic | fixed | 29.20 / 55.8 | 0.0185 / 0.0357 | 0.0102 | 15.12 | 25.8 | 0/9 |
| quadratic | free | 36.46 / 65.8 | 0.0071 / 0.0112 | 0.0008 | 3.91 | 25.3 | 0/9 |
| quadratic | prior | 22.02 / 37.4 | 0.0019 / 0.0027 | 0.0014 | 5.57 | 14.5 | 0/9 |


