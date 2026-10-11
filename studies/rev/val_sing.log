# D11: JWST ℓ comparison with Sing et al. (2026), μ_min = 0.2

## 1. Sing et al.'s per-star Δℓ (vector-extracted from their offset-vs-μ_min figures)

| model | HAT-P-11 | WASP-107 | WASP-12 | WASP-121 | WASP-6 | WASP-80 | WASP-94 | mean 7 | mean w/o WASP-12 | Table 3 |
|---|---|---|---|---|---|---|---|---|---|---|
| phoenix | +0.006 | +0.006 | -0.002 | +0.051 | +0.059 | -0.064 | +0.005 | +0.009 | +0.010 | +0.000 |
| mps1 | +0.040 | +0.046 | +0.013 | +0.061 | +0.090 | -0.009 | +0.026 | +0.038 | +0.042 | +0.034 |
| stagger | +0.014 | +0.029 | -0.028 | +0.050 | +0.062 | +0.004 | -0.001 | +0.019 | +0.026 | +0.020 |

(Table 3 means differ from the plain mean of the plotted per-star values; their averaging is not stated.)

## 2. PHOENIX ℓ under three limb rescalings (μ_min = 0.2)

| star | μ_cri M1 (ours) | μ_cri max dI/dμ (Sing) | Δℓ geom | Δℓ lin | Δℓ lin4 | lin − geom |
|---|---|---|---|---|---|---|
| HAT-P-11 | 0.0472 | 0.0472 | +0.037 | +0.005 | +0.005 | -0.032 |
| WASP-107 | 0.0475 | 0.0475 | +0.036 | +0.004 | +0.004 | -0.032 |
| WASP-121 | 0.0733 | 0.0733 | +0.078 | +0.044 | +0.044 | -0.034 |
| WASP-6 | 0.0477 | 0.0478 | +0.087 | +0.058 | +0.058 | -0.029 |
| WASP-80 | 0.0483 | 0.0484 | -0.036 | -0.067 | -0.067 | -0.031 |
| WASP-94 | 0.0684 | 0.0685 | +0.041 | +0.005 | +0.005 | -0.036 |
* geom: mean Δℓ +0.041 ± 0.044 (s.d.), SEM 0.018
* lin: mean Δℓ +0.008 ± 0.043 (s.d.), SEM 0.018
* lin4: mean Δℓ +0.008 ± 0.043 (s.d.), SEM 0.018
* mean shift lin − geom: -0.032; lin4 − lin: +0.0000

## 3. Per-star offset statistics (observed − model)

| case | n | mean | s.d. | SEM | mean/SEM | median |
|---|---|---|---|---|---|---|
| mps1|hybrid5 | 6 | +0.033 | 0.033 | 0.014 | 2.5 | +0.036 |
| mps1|recipe | 6 | +0.041 | 0.033 | 0.013 | 3.1 | +0.043 |
| mps2|hybrid5 | 6 | +0.035 | 0.031 | 0.013 | 2.8 | +0.037 |
| mps2|hybrid5+B100G | 6 | +0.024 | 0.030 | 0.012 | 2.0 | +0.026 |
| mps2|recipe | 6 | +0.043 | 0.031 | 0.013 | 3.4 | +0.044 |
| phoenix|hybrid5 | 6 | +0.041 | 0.044 | 0.018 | 2.3 | +0.040 |
| phoenix|recipe | 6 | +0.041 | 0.044 | 0.018 | 2.3 | +0.039 |
| phoenix recipe, geom rescaling | 6 | +0.041 | 0.044 | 0.018 | 2.3 | +0.039 |
| phoenix recipe, lin rescaling | 6 | +0.008 | 0.043 | 0.018 | 0.5 | +0.005 |
| phoenix recipe, lin4 rescaling | 6 | +0.008 | 0.043 | 0.018 | 0.5 | +0.005 |

dl_err in val_sing_points.csv is the statistical error of the per-star weighted mean (1/√Σσ⁻²), from the extracted per-bin ℓ errors; it ignores correlations between bins and is far smaller than the star-to-star scatter.

Our recipe vs their plotted per-star Δℓ (μ_min = 0.2):

* mps1|recipe − Sing mps1: HAT-P-11 +0.000, WASP-107 -0.000, WASP-121 -0.008, WASP-6 +0.000, WASP-80 -0.000, WASP-94 -0.000; mean -0.001
* phoenix|recipe − Sing phoenix: HAT-P-11 +0.031, WASP-107 +0.030, WASP-121 +0.027, WASP-6 +0.028, WASP-80 +0.028, WASP-94 +0.036; mean +0.030
* phoenix recipe with Sing's lin4 rescaling − Sing phoenix: mean -0.002, max |Δ| 0.007

## 4. Kepler/TESS (Maxted 2023), b < 0.65

* hosts with b < 0.65: Kepler 24, TESS 10 (of 34 stars; 43 in Maxted's Table 2)
* Kepler mps2 hybrid5 B100G: Δh'1 -0.0023 ± 0.0015 (1.5σ), Δh'2 -0.0039 ± 0.0033 (1.2σ); σ_ext 0.0052, 0.0126
* Kepler mps2 hybrid5 B200G: Δh'1 -0.0130 ± 0.0016 (8.3σ), Δh'2 +0.0065 ± 0.0034 (1.9σ); σ_ext 0.0054, 0.0133
* Kepler mps2 hybrid5 none: Δh'1 +0.0080 ± 0.0017 (4.9σ), Δh'2 -0.0132 ± 0.0033 (4.0σ); σ_ext 0.0058, 0.0127
* Kepler mps2 hybrid5 ssd: Δh'1 +0.0066 ± 0.0017 (3.9σ), Δh'2 -0.0123 ± 0.0033 (3.7σ); σ_ext 0.0060, 0.0127
* TESS mps2 hybrid5 B100G: Δh'1 -0.0039 ± 0.0035 (1.1σ), Δh'2 -0.0003 ± 0.0043 (0.1σ); σ_ext 0.0086, 0.0058
* TESS mps2 hybrid5 B200G: Δh'1 -0.0130 ± 0.0037 (3.5σ), Δh'2 +0.0100 ± 0.0041 (2.5σ); σ_ext 0.0094, 0.0046
* TESS mps2 hybrid5 none: Δh'1 +0.0050 ± 0.0033 (1.5σ), Δh'2 -0.0089 ± 0.0045 (2.0σ); σ_ext 0.0079, 0.0069
* TESS mps2 hybrid5 ssd: Δh'1 +0.0041 ± 0.0033 (1.2σ), Δh'2 -0.0079 ± 0.0045 (1.8σ); σ_ext 0.0080, 0.0067
