# Round-2 code and release checks (studies/rev2/code_checks.py)

1. light_curve() vs ldhybrid.lightcurve (5 laws × p = 0.03, 0.1, 0.3, 400 points): max |ΔF| = 0.0e+00
2. new hybrid columns: expanded form reproduces I(μ) of the weights to 3.4e-06 (6-decimal weights, 8-digit form); q within [-0.000, 1.000]
3. MPS-ATLAS custom throughput at grid models (3 models × TESS/PRISM/MIRI × hybrid5/four-term/power-2): max |Δw| vs table 0.0e+00; off-node (5650 K, 4.42, +0.12, TESS) hybrid5 vs interpolated table: 0.00 ppm in light-curve terms
4. corrections from the release (standalone) vs API: max |Δw| = 4.9e-07 (48 cases); warning outside the MURaM Teff range: True
5. joint condition() = sequential updates for rho = 0 (|Δmean| 1.4e-17, |Δcov| 1.1e-19); sd(h'1) prior / rho=0 / rho=0.5: 0.0153 / 0.0056 / 0.0052
6. prior(): 3000 K G395H PHOENIX, systematic=True: only PHOENIX covers the star → warning; eigen-s.d. without floor [0.0, 0.00014, 0.00112, 0.00405, 0.00768], with the floor at (3500.0, 5.0, 0.0) [0.0, 0.0006, 0.0016, 0.0046, 0.0511]
   cov_floor = 0.01²: eigen-s.d. [7e-05, 0.00034, 0.0008, 0.00194, 0.00823] → [0.0086, 0.0096, 0.0098, 0.0099, 0.0117]
   grid edge (3550 ± 150 K, MPS): 38% of draws dropped, warning True; σ → 0 with B100G reproduces coefficients(): 1.4e-09
7. apply_offsets(): Kepler, MPS2 + B100G prior: h'1 shifted by -0.0016 (offset −0.0023, joint with h'2, rho 0.04); a prior in another state raises: True; JWST ℓ shift +0.0015 (offset +0.024)
8. four-term refit at 30 PHOENIX models: light-curve error median / 90th 3.5 / 6.2 ppm (all μ) → 1.3 / 2.8 ppm (μ ≥ 0.02)
