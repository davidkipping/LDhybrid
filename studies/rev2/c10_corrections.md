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
