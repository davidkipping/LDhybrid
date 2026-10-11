# Revision: validation (D8, D11) and user code (E2)

Scripts: `val_claret2025.py` (D8), `val_sing.py` (D11), `val_code_checks.py` (E2).
Detailed outputs: `val_claret2025.md`, `val_sing.md`, `val_code_checks.md`;
figure inputs `val_sing_points.csv`, `val_maxted_points.csv`.

## Appendix draft: validation of the pipeline

**A.1 Exact reproduction, MPS-ATLAS (Kostogryz et al. 2022).** Their published
band CLVs (J/A+A/666/A60; 24 μ, I/I(1) to six decimals) were refitted with our
engine and their recipe (equal weights on the native μ nodes, I(1) = 1). Over
all 273 280 rows (2 sets × 34 160 models × 4 passbands), the median |Δ| is
0.5–2.2 × 10⁻⁵ in the four-term coefficients and 3–4 × 10⁻⁷ in the power-2
coefficients, with maxima of 1.3 × 10⁻⁴ and 7.6 × 10⁻⁶. The four-term floor
is set by the six-decimal rounding of the published CLVs.

Our own integration of the MPS-ATLAS archive spectra matches their CLVs to
1–7 × 10⁻⁴ (Kepler, TESS, CHEOPS).

One exception: their Set 1 CLVs at [M/H] = −1.7 … −0.75 differ by
0.5–2 × 10⁻² from the archive spectra, which are smooth in [M/H]. We use the
archive spectra.

**A.2 ExoTiC-LD (Grant & Wakeford 2024).** On the same profile, the fitting
agrees to ≤ 1.2 × 10⁻⁵ in the four-term coefficients. Band integration differs
by up to 10⁻² in I/I(1), for two causes:

- different TESS and Spitzer curves;
- ExoTiC-LD's 500-node Gauss–Legendre quadrature over the whole mode, which is
  unconverged on high-resolution spectra (GL500 vs GL20000 differ by up to
  1.7 × 10⁻³).

Our exact bin integrals agree with the exact integral of the linear
interpolants to 4 × 10⁻⁶.

**A.3 Claret & Southworth (2023), power-2 from spherical PHOENIX.** Their
PHOENIX-COND intensities (ξ = 2 km s⁻¹) are not public. We therefore applied
their recipe to PHOENIX-ACES v3.0: M1 and M2 edges, 900 uniform μ points,
unweighted Levenberg–Marquardt. The comparison covers 13 dwarfs in Kepler,
TESS and CHEOPS.

- μ_cri(M1) agrees to ~10⁻³ (e.g. 0.0290 vs 0.0286; 0.0532 vs 0.0532).
- Their M2 edge lies inward of M1, as ours does.
- The power-2 coefficients differ by a median |Δc|, |Δα| = 0.016 (M1) and
  0.014 (M2), with a maximum of 0.04. This is the size of the COND-vs-ACES
  model difference.

**A.4 Claret, Hauschildt & Torres (2025), JWST four-term and power-2.**
Their 2025 grid uses new spherical PHOENIX models: high resolution, solar
composition, ξ = 1 km s⁻¹. These are not public, so the comparison measures
the difference between model versions and recipes, not pipeline error. We
compared 241 models at [M/H] = 0 shared with PHOENIX-ACES v3.0 (2400–7400 K,
log g 3.5–5.5) in the nine passbands common to both: NIRCam F322W2/F444W,
NIRISS SOSS 1/2, and NIRSpec G235H/M, G395H/M and PRISM. That makes 2169
model-passbands. All figures below are median / 90th percentile.

- **μ_cri.** Theirs agrees with our M1 edge to a median |Δ| of 0.0012 (90th
  percentile 0.0054; median signed difference 0.0000), so the limb is placed
  in the same way.
- **Coefficients.** Applying their stated recipe (60 000 points uniform in μ,
  unweighted, unconstrained) to our models gives very different four-term
  coefficients (|Δa_k| = 1–3). This is the known degeneracy of the four-term
  basis. The power-2 coefficients differ by |Δc| = 0.17 and |Δα| = 0.09,
  roughly ten times the PHOENIX-COND vs ACES difference of A.3.
- **Light curves** (p = 0.1 production grid, ppm):

| their coefficients against | four-term | power-2 |
|---|---|---|
| our PHOENIX-ACES model profile | 26 / 55 | 24 / 55 |
| our fit with their stated recipe (uniform μ) | 41 / 66 | 39 / 76 |
| our fit uniform in r | 25 / 54 | 22 / 53 |
| our tabulated r-weighted coefficients | 25 / 53 | 27 / 56 |
| our tabulated uniform-μ coefficients | 40 / 65 | 35 / 74 |

- **Size of the difference.** Their coefficients reproduce our model light
  curves to a median of 24–26 ppm. That is comparable to the difference
  between MPS-ATLAS Set 1 and Set 2 (32 ppm median, P2d) and much smaller than
  PHOENIX vs MPS-ATLAS (172 ppm). We read it as the floor set by the
  model-version difference.
- **Outliers.** The maxima (≈ 370 ppm; 2.8% of rows above 100 ppm) are all
  NIRISS SOSS order 2 at 2400–2600 K, where TiO/VO opacity makes the model
  versions differ most.
- **By passband.** Medians range from 12–16 ppm (G395H/M) to 32–41 ppm
  (F322W2, G235).
- **Which recipe matches.** Their coefficients are closer to our uniform-r or
  area fits (22–25 ppm) than to our uniform-μ fit (39–41 ppm). So their fitting
  in effect down-weights the extreme limb more than a plain uniform-μ fit would.
  Their ReadMe describes a μ_crit strategy "to improve numerical accuracy",
  whose details we cannot reproduce.

**Side finding (not part of D8, flagged for production).** The production
resample of 2000 points uniform in r is not converged for steep PHOENIX limbs.
The limb drop at r ≳ 0.998 is covered by only ~4 resample points. Refitting
with denser resamples lowers the light-curve error against the same profile
(`val_resample_check.md`, ppm at n = 2000 → 30 000):

| law | change |
|---|---|
| hybrid5 | 2.8 → 0.8 (solar analogue, PRISM); 3.0 → 0.8 (K dwarf, G141); 2.6 → 0.6 (2600 K, F444W) |
| four-term | 8.9 → 4.1; 9.7 → 5.5; 4.4 → 0.8 (mid M, TESS) |
| power-2 | 34.9 → 12.7 (2800 K, G141); 18.6 → 4.6 (mid M, TESS) |

- In some cases (e.g. 5800 K TESS, 6500 K Kepler) the error is unchanged or
  rises slightly, by 0.2–0.7 ppm.
- At n = 2000 the r-weighted power-2 fit lands in a different optimum, because
  the discretised objective changes with n; on the 2000-point grid it is the
  global optimum.
- The D8 rows show the same effect. Our tabulated r-weighted fits score
  9.8 / 12.7 ppm (four-term / power-2) against the profile. Unconstrained fits
  on 60 000 r-points score 4.0 / 4.7 ppm. The four-term difference mixes this
  effect with the cost of the four-term monotonicity constraint.
- **Implications.** The headline PHOENIX light-curve errors (Tables 1–3) are
  probably upper limits set by the resample. A production rerun with
  n ≈ 30 000 (or exact trapezoid weights in r on a refined grid) would lower
  them. MPS-ATLAS, whose 24-node profiles have no steep spherical limb, is
  presumably unaffected; this has not been checked.

**A.5 Sing et al. (2026), JWST ℓ.** Their own per-star Δℓ (observed − model,
μ_min = 0.2) was recovered by exact vector extraction from their
offset-vs-μ_min figures. These figures include WASP-12, whose ℓ spectrum is
not published.

- **MPS-ATLAS Set 1.** With their recipe (a quadratic fitted uniform in μ for
  μ ≥ 0.2), our per-star Δℓ agree with theirs to ≤ 0.001 for five stars and
  −0.008 for WASP-121. The mean difference is −0.001.
- **PHOENIX.** Rescaled as they do, μ = (μ_ph − μ_cri)/(1 − μ_cri) with μ_cri
  at max dI/dμ, our per-star Δℓ agree with theirs to a mean of −0.002 (max
  |Δ| 0.007). Our geometric rescaling μ′ = √(1 − (r/r_c)²) at the same μ_cri
  instead gives ℓ lower by 0.032 ± 0.002 (range 0.029–0.036 over the six
  stars). The μ_cri of the two edge definitions coincide (same node; e.g.
  0.0472, 0.0733).
- **Removing the four nodes.** Removing the four nodes inward of μ_cri, as
  they do, changes ℓ by < 10⁻⁴ at μ_min = 0.2.

The PHOENIX difference between the two papers is therefore entirely the
rescaling convention.

## D11: JWST ℓ, definitions and statistics

**(i) Definition.** In all our JWST comparisons, ℓ = 1 − (c₁ + c₂) of the
quadratic law fitted by least squares uniform in μ on 161 points over
0.2 ≤ μ ≤ 1 (p4d_sing.py, `quad_ell_delta`), and δ = c₂/4. This is Sing et
al.'s definition with their recommended μ_min = 0.2.

- *recipe* rows: the quadratic is fitted to the grid's band profile in each
  wavelength bin. For MPS-ATLAS it is fitted on the native nodes with μ ≥ 0.2;
  for PHOENIX, uniform in μ after the M1 cut.
- *hybrid5* rows (Fig. 5): the quadratic is fitted to our hybrid5 law (itself
  the v1 r-weighted fit of the band profile), sampled uniformly in μ ≥ 0.2.

hybrid5 and recipe differ by 0.006–0.008 in Δℓ (MPS-ATLAS) and by < 0.001 for
PHOENIX.

**(ii) WASP-12.** It is absent because Sing et al. publish no ℓ(λ) figure or
table for it (only Table 1 and the aggregate Table 3 include it). It cannot be
added. Its per-star Δℓ in their own figures (μ_min = 0.2) is −0.002 (PHOENIX),
+0.013 (MPS1) and −0.028 (Stagger), all below their means. Omitting it raises
the plain mean of their seven per-star values by +0.001 (PHOENIX), +0.004
(MPS1) and +0.007 (Stagger).

Their Table 3 means (0.000, 0.034, 0.020) are not plain means of the plotted
per-star values (0.009, 0.038, 0.019); their averaging is not stated.

**(iii) Rescaling.** See A.5. Per star (lin − geom): HAT-P-11 −0.032,
WASP-107 −0.032, WASP-121 −0.034, WASP-6 −0.029, WASP-80 −0.031,
WASP-94 −0.036. The mean Δℓ for PHOENIX is +0.041 (geom) against +0.008 (lin;
their per-star mean without WASP-12 is +0.010).

**(iv) Statistics** (6 stars; mean of per-star weighted offsets; SEM = s.d./√6):

| case | mean Δℓ | s.d. | SEM | mean/SEM | median |
|---|---|---|---|---|---|
| MPS-ATLAS Set 2, hybrid5 | +0.035 | 0.031 | 0.013 | 2.8 | +0.037 |
| MPS-ATLAS Set 2 + 100 G, hybrid5 | +0.024 | 0.030 | 0.012 | 2.0 | +0.026 |
| MPS-ATLAS Set 1, hybrid5 | +0.033 | 0.033 | 0.014 | 2.5 | +0.036 |
| PHOENIX, hybrid5 (geometric) | +0.041 | 0.044 | 0.018 | 2.3 | +0.040 |
| PHOENIX, recipe, linear rescaling | +0.008 | 0.043 | 0.018 | 0.5 | +0.005 |

- The 100 G correction moves the mean by −0.011, against a SEM of 0.012–0.013.
- WASP-6 is the largest offset in every case: +0.083 (MPS2) and +0.072
  (MPS2 + 100 G).
- Per-star statistical errors from the extracted bin errors
  (`val_sing_points.csv`, column `dl_err`, 1/√Σσ⁻²) are 0.001–0.004. They
  ignore correlations between bins and are much smaller than the 0.03
  star-to-star scatter, so the error bars in Fig. 5 (right) should show them
  and the shading the scatter.

**(v) Kepler/TESS (Maxted 2023).**

- **Sample.** 34 hosts have b < 0.65: 24 Kepler and 10 TESS, out of 43 in his
  Table 2.
- **Errors.** Weighted means with a fitted external scatter (Maxted's
  approach); the ± is the error of the mean.
- **MPS-ATLAS Set 2 + hybrid5:**
  - Kepler: Δh′₁ = +0.0080 ± 0.0017 (4.9σ) and Δh′₂ = −0.0132 ± 0.0033 (4.0σ).
    With 100 G these become −0.0023 ± 0.0015 (1.5σ) and −0.0039 ± 0.0033
    (1.2σ).
  - TESS: +0.0050 ± 0.0033 (1.5σ) and −0.0089 ± 0.0045 (2.0σ). With 100 G,
    −0.0039 ± 0.0035 (1.1σ) and −0.0003 ± 0.0043 (0.1σ).
  - σ_ext is 0.005–0.009 (h′₁) and 0.005–0.013 (h′₂).
- **Figure input.** All cases (grid × law × correction) are in
  `val_maxted_points.csv`.

## Sampling coordinates q (for Sect. 2)

hybrid2 (ε = 0.36): the exact physical region is the triangle with vertices O = (0, 0) and

    V_c = (−0.124615, 1.124615),   V_ℓ = (1.201039, −0.201039),

with V_c = (−g₁/(1 − g₁), 1/(1 − g₁)) and V_ℓ = (g₀/(g₀ − 1), −1/(g₀ − 1)). Here
g₀ = 2/(N ε³) and g₁ = 2/(N (1 + ε)³) are the slopes of Π_ε in ξ = μ² at the
limb and at the centre, and N = ε⁻² − (1 + ε)⁻².

The three edges are:

| edge | condition that holds on it |
|---|---|
| O–V_c | dI/dξ = 0 at the centre (ξ = 1) |
| O–V_ℓ | dI/dξ = 0 at the limb (ξ = 0) |
| V_c–V_ℓ | I(0) = 0 |

O is the uniform disc. A uniform point of the triangle is

    w = √q₁ [(1 − q₂) V_c + q₂ V_ℓ],   q₁, q₂ ~ U(0, 1).

hybrid4/5 (n = 4, 5 weights): a uniform point of the simplex
{w_j ≥ 0, Σ w_j ≤ 1} (a flat Dirichlet with n + 1 components), by stick-breaking:

    w_j = (1 − q_j^{1/(n+1−j)}) ∏_{i<j} q_i^{1/(n+1−i)},   j = 1 … n,   q_j ~ U(0, 1).

Implied priors (200 000 draws), mean ± s.d. and [5th, 95th] percentiles:

| law | h′₁ | h′₂ | ℓ | I(0) |
|---|---|---|---|---|
| hybrid2 | 0.757 ± 0.141 [0.49, 0.95] | 0.247 ± 0.090 [0.08, 0.38] | 0.26 ± 0.27 [−0.12, 0.76] | 0.33 ± 0.24 |
| hybrid4 | 0.754 ± 0.123 [0.52, 0.92] | 0.213 ± 0.080 [0.09, 0.35] | 0.35 ± 0.23 [−0.03, 0.73] | 0.20 ± 0.16 [0.01, 0.53] |
| hybrid5 | 0.789 ± 0.108 [0.58, 0.93] | 0.192 ± 0.073 [0.08, 0.32] | 0.41 ± 0.21 [0.05, 0.74] | 0.17 ± 0.14 [0.01, 0.45] |

So "uninformative" q-sampling is informative in profile space. For hybrid5 it
favours a dark limb (I(0) ≲ 0.45 at 95%) and h′₁ ≈ 0.8. Observed FGK values
(h′₁ ≈ 0.7–0.9, h′₂ ≈ 0.15–0.25) lie inside its bulk.

## E2: user code

New in `ldhybrid/interpolate.py`, which is standalone (numpy and scipy; h5py
only for refit). Thin wrappers are in `ldhybrid/api.py`: `ld.observables`,
`ld.condition`, `ld.prior`, and `ld.coefficients(..., refit=True)`.

- **`observables(law, coeffs=None)`**
  - Returns h1 = I(2/3), h2 = I(2/3) − I(1/3), ell = 1 − (c₁ + c₂) and
    delta = c₂/4. The quadratic is fitted uniform in μ over 0.2 ≤ μ ≤ 1 (161
    points, the D11 definition).
  - Each comes as (const, grad) with obs = const + grad · c, which is exact
    for the linear laws.
- **`LDTables.prior(law, teff, σ, logg, σ, mh, σ, band, grid, n, systematic, seed)`**
  - Monte-Carlo mean and covariance of the interpolated coefficients.
  - `systematic=True` adds the covariance of the phoenix/mps1/mps2 coefficients
    at the nominal point.
- **`condition(mean, cov, law, observable, delta, sigma)`**
  - Exact Gaussian update by an empirical offset in one linear observable:
    K = Cg/(gᵀCg + σ²), m′ = m + Kδ, C′ = C − K gᵀC.
  - This is how the tabulated empirical offsets map onto w: a 1-D (ℓ) or 2-D
    (h′₁, h′₂) constraint on a 5-D weight vector, which leaves the other
    directions at their prior.
- **`LDTables.refit` / `coefficients(..., refit=True)`**
  - Barycentric (Delaunay) interpolation of the post-cut band profiles on one
    r-uniform 2000-point grid (PCHIP per corner), then a refit with the
    production weighting and constraint.
  - The constraint is solved with lstsq, then BVLS for the simplex, then SLSQP;
    power-2 uses LM.
  - The profiles are read from `data/products/profiles/<grid>.h5`, written by
    FORK 1 (phoenix only at the time of testing).

Checks (`val_code_checks.md`):
1. prior with σ → 0 reproduces `coefficients()` to 1 × 10⁻⁹, with s.d.
   ≤ 1 × 10⁻⁸.
2. condition (Δh′₁ = +0.008, σ = 0.006): s.d.(h′₁) falls from 0.0076 to 0.0047
   and h′₁ moves by +0.0049, both exactly as expected. observables() matches
   direct evaluation to 9 × 10⁻¹⁶.
3. refit at the models reproduces the tabulated coefficients to ≤ 7 × 10⁻⁷
   (hybrids, quadratic, power-2) and ≤ 8 × 10⁻⁶ (four-term). This is the
   float32 profile storage. Where the four-term monotonicity constraint binds,
   SLSQP on finite differences differs from the production QP by up to
   1.3 × 10⁻³ in one coefficient.
4. Leave-one-out across 200 K in Teff (PHOENIX, 200 random model-bands); light
   curve vs the left-out profile, median / 90th percentile (ppm):

| law | interpolated coefficients | refit of interpolated profiles |
|---|---|---|
| hybrid5 | 6.0 / 34.5 | 6.0 / 34.5 (identical; max \|Δw\| 7 × 10⁻⁷) |
| four-term | 11.7 / 35.2 | 11.7 / 35.2 |
| power-2 | 20.3 / 56.2 | 19.3 / 52.2 |

**Conclusion for the referees' suggestion (R1 4iii, R2 m15).** For the linear
laws, an r-weighted least-squares fit is a fixed linear map of the profile.
Refitting an interpolated profile is therefore identical to interpolating the
coefficients wherever no constraint binds. The interpolation error is that of
the model intensities themselves, not of the coefficients. Refitting helps only
power-2 (nonlinear), by ~1 ppm in the median and ~4 ppm at the 90th
percentile. The option is provided for completeness.
