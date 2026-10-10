# Gate G3: P3 products

Status: **passed** (2026-10-10). `g3_checks.py` → `g3_report.json`.
Products are described in `docs/products.md`.

## Production summary (headline method, r-weighted)

Light-curve error (ppm) of each row against the exact light curve of its own
tabulated profile, production grid (p = 0.1, b ∈ {0, 0.3, 0.6, 0.9}):

| grid (rows) | law | median | 99th pct | max |
|---|---|---|---|---|
| PHOENIX v3.0 (63,916) | hybrid5 | 2.7 | 6.3 | 49.5 |
| | hybrid4 | 3.8 | 16.5 | 60.8 |
| | Claret 4-term | 8.0 | 19.4 | 74.9 |
| | power-2 | 16.3 | 61.7 | 188 |
| | hybrid2 | 44.4 | 114 | 344 |
| | quadratic | 44.1 | 134 | 508 |
| MPS-ATLAS Set 2 (645,012) | hybrid5 | 0.95 | 4.3 | 7.2 |
| | hybrid4 | 6.8 | 13.6 | 17.6 |
| | Claret 4-term | 1.4 | 5.8 | 10.2 |
| | power-2 | 7.2 | 34.2 | 122 |
| | hybrid2 | 43.1 | 83.0 | 110 |
| | quadratic | 51.4 | 99.6 | 131 |
| Kurucz, MIRI LRS (507) | hybrid5 | 1.4 | 2.1 | 2.2 |

## Checks

1. **Physicality.** Every hybrid row is physical in every grid and method
   (PHOENIX 63,916; MPS-ATLAS Set 1 and Set 2 645,012 each; Kurucz 507).
   Constraints bind in ≤ 1% of PHOENIX hybrid5 rows and never for MPS-ATLAS.
   Power-2 (unconstrained, as published elsewhere) is unphysical in 8.6%
   (PHOENIX) and 3.4% (MPS-ATLAS Set 2) of rows, mostly c > 1 in cool stars.
2. **Interpolation** (leave-one-out, linear across a doubled gap; PHOENIX,
   hybrids, ~4,000 trials per axis): median / 90th percentile light-curve
   error 4.7 / 30 ppm across 200 K in Teff, 15 / 67 ppm across 1.0 dex in
   log g, 13 / 61 ppm across 1.0 dex in [M/H]. At the real grid spacing the
   errors are smaller (≈ ¼ for curvature-dominated error).
3. **Coefficient-spectra closure.** Band weights from the spectra equal the
   direct same-settings band fit to ≤ 8 × 10⁻¹⁰ (float32 storage), 25 models
   × 6 bands × 3 laws; 448/450 band results physical (the per-bin fits are
   unconstrained).
4. **Systematic spread** (`systematics.csv.gz`; 1,020 grid points shared by
   PHOENIX, MPS-ATLAS Set 1 and Set 2, 116,280 rows): the largest pairwise
   light-curve difference between the three grids' coefficients is a median
   93 ppm, 90th percentile 285 ppm, identical for all three hybrids.

**Gate:** interpolation error (median 5–15 ppm, 90th 30–67) is well below the
systematic spread between grids (median 93, 90th 285), and every hybrid row
is physical. Passed.

## Notes

* The worst PHOENIX hybrid5 rows (56 of 63,916 above 15 ppm) are 2300–2500 K,
  metal-poor, blue-optical bands, where the PHOENIX profile itself brightens
  toward the limb at μ ≈ 0.01–0.02; the monotonicity constraint binds there by
  design.
* On MPS-ATLAS, hybrid4 is weaker than Claret 4-term (median 6.8 vs 1.4 ppm);
  hybrid5 is better than Claret 4-term on both grids.
* hybrid2 (ε = 0.36) matches quadratic in median and beats it in the tails.
