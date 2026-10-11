# Response to the referees, second round (draft)

We thank the referees again. All three accept the reframing and recommend
minor revision, conditional on the companion paper being public. It will be
on arXiv before acceptance, and the reference will be completed then. We have
made the changes below.

## Common points

1. **"Elementary" made checkable.** Sect. 2 now carries the argument for a
   circular star. The blocked flux is ∮ G(r) dθ by Green's theorem, with G
   polynomial (even powers) or rational (double poles) in r². A circular or
   elliptical boundary arc has a rational parametrisation, so the integrand is
   rational. Odd powers of μ or simple poles give elliptic or non-elementary
   functions instead. A basis-shape figure was added (new Fig. 1). Until the
   companion code is in use, light curves of any tabulated law can be computed
   with `light_curve()` in our interpolation code, or with any code that
   accepts a radial profile.

2. **Advantages for circular planets.** These are now separated into two
   groups:
   - shared with the quadratic and four-term laws: convexity, analytic
     marginalisation (conditional on the geometry), and exact spectral
     summation;
   - specific to the hybrid laws: the monotonic basis with a simplex
     sufficient region that the fits lie inside, and the elementary light
     curve.

   The Introduction also says that, for circular-planet users, the main value
   may be the uniform tables, the error budget and the conventions analysis.

3. **Limb rescaling (R1 new 2a, R2 new 1).** We have corrected the claim.
   - Linear and geometric rescaling change the light curves by 63–151 ppm
     (median 99) on the Table 3 cases, and by a median of 100 ppm (90th
     percentile 212) on the method study.
   - We state that the geometric rescaling is the physically correct one, and
     that the PHOENIX agreement of Sing et al. rests on the linear remapping.
   - The rescaling now appears in Table 7 as a convention, in the abstract and
     in the Conclusions.

4. **Release (R1 new 9, R2 new 3).**
   - The CDS hybrid tables now carry the q coordinates and the expanded form.
   - corrections.dat, offsets.dat and systematics.dat are added with
     byte-by-byte descriptions.
   - The ReadMe states:
     - that the r-weighted coefficients are for transits and the uniform-μ
       fits for limb-sensitive uses;
     - the bound and unphysical fractions per law, grid and method;
     - how the Δw corrections are constructed and which extrapolations they
       involve;
     - the rule that an offset must match the correction state.
   - The standalone code now applies the Δw corrections, so the recommended
     MPS-ATLAS Set 2 + 100 G prior can be built from the CDS files.
   - It also computes MPS-ATLAS coefficients for any throughput from the model
     archive (about 9.5 GB per set).
   - The kurucz grid is removed.

## Referee 1

- **New 1 (advantages):** see common point 2.
- **New 2b (ℓ route):** the headline now uses the model-profile route, as Sing
  et al. do.
  - MPS-ATLAS Set 2 gives +0.043 (3.4σ), and +0.033 (2.7σ) with 100 G applied
    to the profiles.
  - The hybrid5 route is 0.008 lower, and we explain why.
  - The abstract and Fig. 5 are updated. We state the per-star weighting and
    note that Sing et al.'s tabulated means differ from the plain means.
- **New 3 (four-term recipe):** stated in the abstract, Sect. 5, Sect. 8 and
  the Conclusions. With μ < 0.02 excluded, the four-term law matches hybrid5
  on PHOENIX. The code offers this refit (`refit(..., mu_min=0.02)`: median
  3.5 → 1.3 ppm). We keep the tabulated fits on the stated weighting, so that
  all laws are treated alike.
- **Response 2(d) (light-curve-optimal fits):** not added. The μ_min option
  gives the practical gain. Fitting to a fixed set of light curves would make
  the coefficients depend on that set.
- **New 4 (Claret et al. 2025):** Appendix A now reports the coefficient-level
  differences.
  - Power-2 differs by |Δc| = 0.17 and |Δα| = 0.09; the four-term coefficients
    differ by 1–3.
  - We state that neither spherical comparison is like-for-like, that their
    recipe could not be reproduced, and that we cannot separate model
    differences from recipe differences.
  - The "model versions" and "COND–ACES" attributions are removed.
- **New 5 (Sun):** the recommended combination (MPS-ATLAS Set 2 + 100 G) is
  now confronted with the Sun. It gives 136/55 ppm in the optical/near-IR,
  against 43/53 ppm without correction and 24/40 ppm with the small-scale
  dynamo. Sect. 8 now says the quiet Sun prefers a weaker field and suggests
  the small-scale-dynamo correction for quiet Sun-like stars.
- **New 6 (Δw extrapolations):** both are stated in Sect. 9 and the ReadMe.
  - Holding the ratio constant below μ = 0.1, rather than extrapolating it
    linearly, changes the light curves by ≤ 0.7 ppm.
  - The size of the 100 G correction (median 34 ppm, 90th percentile 133) is
    in Table 7.
  - The MURaM nodes are listed: seven solar-metallicity stars from 3196 to
    6807 K, plus two metallicity variants.
- **New 7 (budget):** Table 7 is restructured into uncertainties and
  conventions.
  - Interpolation is measured at two gap sizes. The error scales as the gap to
    the power 1.7; the native-spacing estimate is 1.3 ppm in the median and
    12 ppm at the 90th percentile.
  - "Smaller than interpolation" is now "comparable to".
  - Photon vs energy weighting is presented as a convention, with its full
    range.
  - The hybrid2/quadratic row now reads 41–51 ppm.
- **New 8 (free fits):** recommendations now distinguish fixed and free
  coefficients. Hybrid4's 34 ppm free bias is stated. The note on power-2 in
  the bias table says it is constrained to c ≤ 1. Howarth (2011) and Espinoza
  & Jordán (2016) are credited.
- **New 10 (power-2):** the text now gives:
  - the unphysical fractions by Teff;
  - the distribution of c;
  - the constrained cost;
  - a warning for codes that assume 0 < c ≤ 1.

  Unphysical entries in Table 3 are marked.

  On whether Claret & Southworth's PHOENIX tables show comparable fractions:
  their μ-uniform, unweighted fits on PHOENIX-COND are not public in a form
  we could check.
- **Minor comments:**
  - quadrature wording: "converged", with the PCHIP nonlinearity noted;
  - reference light-curve convergence: ≤ 0.012 ppm, and ≤ 0.25 ppm in the
    error under chord sampling;
  - lightning approximation citations added (Gopal & Trefethen 2019;
    Herremans et al. 2023), with α = 1 stated;
  - ladder sensitivity: hybrid5 is near-optimal at α = 1. For hybrid4 the
    test showed that α = 1 was far from optimal. We therefore adopted
    α = 25/16 for hybrid4 (ln ε scaled by 0.8; ε = 0.0286, 0.386), tested out
    of sample (held out by grid, random half-splits, and 8100 further
    model–passband combinations), and rebuilt all hybrid4 products. Grid
    medians fall from 4.6/6.8 to 2.1/1.5 ppm (PHOENIX/MPS-ATLAS), so hybrid4
    now matches the four-term law. Its free-coefficient depth bias falls from
    34 to 12 ppm. Sect. 2 explains the choice;
  - the list of four-term wins in TESS is corrected;
  - the MPS-ATLAS < 4600 K remark is added;
  - the Fig. 2 range is limited to the grid;
  - Fig. 4 labels and rounding are fixed, and it now shows two stars;
  - Fig. 5 now has error-bar notes and a mean ± SEM row, and WASP-94 is
    relabelled WASP-94A;
  - the statistics are described;
  - the metric labels for the medians are stated;
  - the chromatic-edge maximum (33 ppm) is given;
  - the metal-poor brightening survives the M2 edge and edges moved further
    in, so it is a model feature;
  - the grazing transit errors are given relative to the depth;
  - Δb is added to Table 5;
  - the Fisher set-up and the best-measured eigenvector are given;
  - the uniform-μ bound fractions are given;
  - ATLAS9 is cited (Castelli & Kurucz 2003), and the Stagger baseline mixing
    is stated;
  - the TESS and IRAC handbooks are cited;
  - the TESS offsets are quoted;
  - the ExoTiC-LD causes are separated, with numbers;
  - I(1) = 1 and the flux-conservation Δw are given;
  - interpolation is now described as Delaunay barycentric;
  - LDTk (Parviainen & Aigrain 2015) is cited;
  - law names at sentence starts are fixed.

## Referee 2

- **New 1:** see common point 3.
- **New 2 (prior safeguards):** done in the code and stated in Sect. 8.
  - When only PHOENIX covers the star, the code warns and substitutes the
    grid spread at 3500 K. Since that floor is rank-deficient, we also
    recommend a covariance floor (`cov_floor`, e.g. 0.01²).
  - A warning reports how many draws fall outside the grid.
  - `condition()` and `apply_offsets()` apply h′₁ and h′₂ jointly with their
    correlation (ρ = 0.02–0.09 for MPS-ATLAS, −0.37 for PHOENIX Kepler) and
    refuse an offset in a different correction state. The docstring example
    is fixed.
  - The passband applicability of each offset is in offsets.dat.
  - The flat-q statement is corrected for the infrared, where we recommend
    table-centred priors.
  - The Fisher details are given.
- **New 3:** see common point 4. The MURaM construction and its extrapolation
  are stated; the custom MPS-ATLAS route is implemented.
- **New 4 (bias test):** redone on an equal footing.
  - Every linear law is free within its exact physical region, and power-2
    within c ≤ 1. The hybrid constraints never bind, so hybrid4's 34 ppm comes
    from the law.
  - Δb is added, and the set-up is stated.
  - A prior-constrained mode (0.05 per coefficient) is added.
  - MPS-ATLAS is repeated for nine FGK and early-M cases: the four-term law is
    slightly better there.
  - Full noisy Fisher depth widths per law are not added. The Fisher
    combinations are given, and the prior mode illustrates prior-constrained
    fitting.
- **New 5 (power-2 in M dwarfs):** see R1 New 10. Recommendations are
  qualified.
- **New 6 (interpolation):** see R1 New 7. 90th percentiles are quoted, and the
  scaling is measured rather than assumed.
- **New 7:** see common point 2.
- **Minor comments:** addressed:
  - rms column added to Table 4;
  - grazing errors given relative to the depth;
  - Table 7 sources and labels fixed;
  - method-study metric labels fixed;
  - the ℓ sign is corrected;
  - the route effect is stated;
  - the interpolation error is now quoted as 4 ppm (the 4.8 ppm was a
    different sample);
  - Fig. 2 range fixed;
  - labels fixed;
  - Fig. 5 mean ± SEM added;
  - ReadMe σ symbols and the transit caveat fixed;
  - release sizes stated (0.9 GB CDS, 1.1 GB Zenodo);
  - uniform-μ bound fractions given;
  - MIRI sentence added;
  - narrow-bin remedy stated;
  - spectra quadrature stated;
  - conditional marginalisation clarified;
  - activity references (Ballerini et al. 2012; Csizmadia et al. 2013) added;
  - chromatic-edge maximum given.

## Referee 3

- **New 1:** see common point 1.
- **New 2:** see common point 2. The K = 3 justification now rests on the
  ~1 ppm floor of the profile ambiguity. "As accurate as the four-term law" is
  qualified to tabulated coefficients.
- **New 3 (abstract):**
  - law accuracy is stated in bands (1–3, 5–8 and 40–50 ppm);
  - the limb-convention result is added;
  - photon vs energy weighting is presented as a convention;
  - "smaller than interpolation" is now "comparable to".
- **New 4:** free-fit recommendations added (Sect. 8).
- **Minor comments:** addressed:
  - μ defined;
  - acronyms expanded;
  - chord velocity defined;
  - ATLAS9 restored and named consistently;
  - "r-weighted" made explicit in the fraction statement;
  - 7524 = 198 × 38 passbands;
  - null-law metric and the 5–30 split given;
  - the hybrid2 rescan sentence rewritten;
  - the four-term wins list fixed;
  - the bias set-up stated;
  - solar sources given, including the recommended combination;
  - Table 7 labels fixed;
  - MURaM wording fixed;
  - Fig. 4 (ΔF, solar analogue, rounding; reduced to two stars);
  - Fig. 5 (error bars, WASP-94A);
  - Fig. 3 caption warning;
  - Fig. 2 now one column;
  - the closed-form check count given (2×10⁻⁸ ppm over 81 000 points);
  - the provenance of the prior intervals given;
  - the cadence qualification added;
  - Barros et al. (2022) and Akinsanmi et al. (2024) cited, and the oblateness
    statement made precise;
  - the Sing et al. table moved to Appendix A.

## Still open (authors)

- Author list, affiliations, acknowledgements.
- Companion arXiv identifier.
- Zenodo DOI and licence.
- Propagate the new hybrid4 poles (α = 25/16) to the companion paper
  (SquishierPlanet and Elómë v0.15.0 already carry them).
