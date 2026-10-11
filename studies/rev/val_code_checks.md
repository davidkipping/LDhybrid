# E2: user-code checks

1. prior(σ → 0) − coefficients(): max |Δmean| 1.2e-09, max sd 1.1e-08 (6 cases, mps2 CDS tables)
   example hybrid5 TESS (5600±100, 4.4±0.1, 0.1±0.1), + systematic over ['mps2', 'phoenix', 'mps1']: mean [0.0841, 0.1418, 0.0896, 0.1079, 0.2228], sd [0.007, 0.0143, 0.0109, 0.0042, 0.0178]
2. condition(h1: δ = +0.008, σ = 0.006): sd(h1) 0.0076 → 0.0047 (expected 0.0047); h1 shift +0.0049 (expected +0.0049)
   observables() vs direct evaluation (h1, h2; ℓ, δ as p4d_sing): max |Δ| 8.9e-16
3. refit at PHOENIX models vs tabulated (40 random model-bands): hybrid5 max |Δ| 6.3e-07, hybrid4 max |Δ| 5.6e-07, hybrid2 max |Δ| 1.1e-07, claret4 max |Δ| 1.5e-05, quadratic max |Δ| 2.0e-07, power2 max |Δ| 6.2e-06
   rows with an active constraint (3): max |Δ| 7.5e-02
4. leave-one-out across 200 K in Teff (PHOENIX, 200 trials), light-curve error vs the left-out profile, median / 90th pct (ppm):
   hybrid5: interpolated coefficients 4.8 / 27.8; refit of interpolated profiles 4.8 / 27.8; max |Δw| between the two 6.1e-07
   claret4: interpolated coefficients 7.1 / 27.5; refit of interpolated profiles 7.1 / 27.5; max |Δw| between the two 3.1e-01
   power2: interpolated coefficients 16.5 / 47.4; refit of interpolated profiles 14.3 / 43.6; max |Δw| between the two 1.5e-01

For the linear laws an unconstrained r-weighted fit is a fixed linear map of the profile, so refitting an interpolated profile equals interpolating the coefficients wherever no constraint binds; the interpolation error is that of the profiles themselves. Refitting differs only for power-2 (nonlinear) and where a constraint is active.
