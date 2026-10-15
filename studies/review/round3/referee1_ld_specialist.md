```
Manuscript: Limb-darkening coefficients for rational laws with elementary transit light curves
            (second revision, third round)
Referee: Referee 1 (theoretical limb-darkening specialist; model-atmosphere LDC tabulations)
Recommendation (A&A grade): B   [= Minor revision]
Style and language: (a) well written, concise and self-contained. The main body now runs to
                    about 12.3 pages (Conclusions end on p. 13), slightly over the 12-page
                    limit. Moving Table 1 or Fig. 4 to an appendix would fix this with no loss.
Identify referee: no
```

---

## Confidential comments to the Editor

Most of my second-round report has been dealt with, and dealt with well.

- The limb-rescaling error is corrected and turned into a finding: about 100 ppm, which accounts for Sing et al.'s PHOENIX agreement.
- The JWST ℓ is now computed as described (+0.043, 3.4σ).
- The Claret et al. (2025) appendix reports the coefficient-level disagreement candidly.
- The recommended MURaM 100 G correction is now set against the Sun, which prefers a weaker field.
- The error budget now separates uncertainties from conventions.
- The release matches the paper.

The authors also acted on my minor comment about the pole ladder. They found that the hybrid4 poles were far from optimal, re-tuned them with proper out-of-sample tests, and rebuilt every hybrid4 product. Hybrid4 is now as accurate as the four-term law: 1.8 vs 5.8 ppm out of sample, and the free-fit depth bias falls from 34 to 12 ppm. This is a genuine improvement, and it came from the review.

After their response they changed hybrid4 again, to a "curvature-simplex" basis. The span is the same, so profiles and light curves are unchanged. The constraint region and the sampling map are new. I checked the mathematics independently:
- the five-inequality reduction;
- the vertices S3 and S4 to six digits;
- the volume fractions: 8.8% of the positive–monotonic region in my Monte Carlo, against 9% claimed;
- that 96% of the excluded volume is convex in μ².

All of it is correct. My concerns are with how the new basis is framed and documented, not with the algebra.

1. **The curvature conditions are oversold.** "Concave in μ² with a convex slope" is a regularity of these two grids in these passbands, not a physical requirement. Kipping's conditions are physical. The paper nonetheless calls the region "exact", "complete" and the hybrid4 "physical region", and draws an analogy with Kipping's q1, q2 that holds for the construction but not the scope: the region excludes 91% of positive, monotonic hybrid4 profiles. A simple grey LTE calculation shows that the conditions fail in the Wien regime for cool stars. That is relevant because the coefficient spectra extend to 0.1 μm.

2. **One statement is contradicted by the authors' own files.** The paper says "all 2.9 million tabulated hybrid4 fits satisfy the curvature conditions with a margin", and that only 123 fits would gain from dropping them. According to studies/rev4/rv_bound.md, the curvature constraint is active in 18.6% of the uniform-μ PHOENIX rows (12 541 fits, at all T_eff, mostly infrared bands). Those fits lie on the faces of the region with zero margin. The Conclusions' "which every tabulated fit is" is circular, because the fits are constrained to the region.

3. **Several numbers went stale with the hybrid4 changes:**
   - the weighting factors "1.7–6" and "up to 23" (now 1.9 and 32);
   - the number of MURaM stars;
   - a cross-reference to a number that is not in the paper;
   - hybrid4 on the Sun, which is omitted and was computed with the old poles.

4. **The response again claims changes that are not in the manuscript.** This is the third round in which this has happened:
   - Gopal & Trefethen and Herremans et al. are cited but missing from the bibliography (they print as "?").
   - LDTk is in the bibliography but never cited.
   - The reference-light-curve convergence (≤ 0.012 ppm) is absent.
   - The ≤ 0.7 ppm size of the Δw extrapolation is absent.

   These are small items, but a final check of the response against the manuscript is needed.

5. **Hybrid4 now differs from the companion paper as submitted, in both its poles and its basis.** Acceptance was already to be conditional on the companion being on arXiv. It should also be conditional on the companion carrying the same hybrid4 definition. Otherwise users will meet two different "hybrid4" laws.

None of this needs substantial new computation. I recommend minor revision. I do not need to see the whole manuscript again, but I would be glad to check a revised Sect. 2 and the responses to New comments 1–2 below.

---

## Report to the authors

### Summary of the paper

The authors tabulate coefficients of three rational limb-darkening laws (hybrid2, hybrid4, hybrid5) and of the quadratic, power-2 and four-term laws, computed in the same way.

- **Models and passbands.** PHOENIX-ACES v3.0 (spherical; M1 edge per passband) and MPS-ATLAS Sets 1 and 2 (plane-parallel), in 40–41 passbands.
- **Fitting.** Least squares uniform in r, with exact quadrature, plus fits uniform in μ.
- **Fit quality.** Scored by the largest light-curve error for p = 0.1.

Changes since the second round:
- hybrid4 has re-tuned poles, from α = 1 to α = 25/16 (ε = 0.0286, 0.3858), and a new basis. Its shapes are the vertices of the region where the profile is positive, monotonic, concave in μ², and has a convex slope in μ². That region is the simplex w ≥ 0, Σw ≤ 1, and it maps to [0, 1]⁴ with a Kipping-type q1 = s⁴.
- The linear vs geometric limb-rescaling convention is quantified (about 100 ppm) and shown to account for Sing et al.'s PHOENIX agreement.
- JWST ℓ is computed from the model profiles: +0.043 (3.4σ), and +0.033 (2.7σ) with 100 G.
- The MURaM correction is confronted with the Sun.
- The error budget is restructured.
- Parameter biases now include Δb and a prior mode.
- Appendix A reports coefficient-level comparisons.
- The CDS release now includes q, the expanded form and three auxiliary tables.

### General assessment

**Q1 (why publish).** The case is now clear and well supported. The tables are of lasting value. Several features go beyond current practice and are, in my view, the paper's most durable contribution:
- a light-curve error for every row;
- constraint and physicality flags;
- the grid-to-grid systematic;
- the conventions analysis (limb rescaling; photon vs energy weighting).

**Q2 (assumptions).** Mostly clear. The exception is new: the curvature conditions imposed on hybrid4 are an assumption about stellar profiles. They are presented as if they were physical conditions (New comment 1).

**Q3 (methods).** Adequate. The derivation of the hybrid4 region is correct and clearly written; I verified it independently. The reference-light-curve convergence promised in the response is not in the text (New comment 3).

**Q4 (emphasis).** Good. The abstract is proportionate and does not mention the new basis, which is appropriate. The Conclusions overstate it (New comment 1).

**Q5 (figures and tables).**
- Fig. 1(b) is useful, but its caption does not describe it.
- The S3/S4 equations run beyond the column in the compiled PDF.
- Table 1 and Fig. 4 could go to an appendix to meet the page limit.

**Q6 (electronic material).** Appropriate. The ReadMe needs small corrections for hybrid4 (New comment 4).

**Q7 (IAU).** No problems.

---

### Response to my second-round comments

**New 1 (advantages shared with the standard laws). Addressed, with one new inaccuracy.**

Sects. 1 and 2 now separate the consequences of linearity, which are shared with the quadratic and four-term laws, from what is specific to the hybrid laws. That is what I asked for. The four-term spectra are mentioned as possible ("the same construction would work"), but not provided; that is acceptable.

The new sentence "their physical regions are simplices" (Sect. 1) is, however, not accurate:
- hybrid5's physical region is not a simplex; it contains one (about 8% of it, Sect. 2);
- hybrid4's simplex is not its physical region but a curvature-restricted part of it (New comment 1).

**New 2a (limb rescaling). Addressed.**

The light-curve effect is given:
- 63–151 ppm on the Table 3 cases;
- a median of 100 ppm (90th percentile 212 ppm) on the method study.

It appears in Table 6 and in the abstract, and the geometric convention is identified as the correct one. This turns a convention into a finding, as I hoped.

**New 2b (ℓ route). Addressed.**

The model-profile route is used throughout:
- +0.043 ± 0.013 (3.4σ);
- +0.033 ± 0.012 (2.7σ) with 100 G.

The 0.008 difference of the hybrid5 route is explained, and the averaging is stated. I recomputed Table A.1:
- geometric − linear = 0.029–0.036;
- linear vs Sing et al.: mean |Δ| = 0.002;
- means +0.041 (SEM 0.018), +0.008 and +0.010.

All agree.

**New 3 (four-term recipe). Addressed.**

The μ < 0.02 result is in the abstract, Sects. 5 and 8 and the Conclusions, and the code offers the refit. Keeping the tabulated fits on a common weighting is a defensible choice once it is stated this clearly.

**Response 2(d) (light-curve-optimal fits). Declined, with an acceptable reason.**

**New 4 (Claret et al. 2025). Addressed.**

- The coefficient-level differences are given: |Δc| = 0.17, |Δα| = 0.09, and 1–3 for the four-term law.
- The comparison is stated not to be like-for-like, and the recipe is stated not to be reproducible.
- The "model versions" and "COND–ACES" attributions are removed.

A request to Claret for a handful of intensity files would still be the cleanest test, but that is a suggestion.

**New 5 (Sun). Addressed.**

The recommended Set 2 + 100 G combination gives 136/55 ppm, against 43/53 ppm uncorrected and 24/40 ppm with the small-scale dynamo. Sect. 8 now suggests the small-scale-dynamo correction for quiet Sun-like stars. See minor comment 6 for the wording.

**New 6 (Δw extrapolations). Partly addressed.**

- Both extrapolations are stated in Sect. 9.2 and the ReadMe.
- The size of the 100 G correction is in Table 6.

Two problems remain:
- Sect. 9.2 says that "Sect. 8 gives the size of this assumption" (the constant ratio below μ = 0.1). Neither Sect. 8 nor Table 6 does. The ≤ 0.7 ppm in your response is not in the manuscript.
- The paper says the MURaM simulations exist for "six solar-metallicity stars between 3196 and 6807 K". The ReadMe (Note 10) and your response list seven: 3196, 3518, 3707, 4543, 5240, 5778 and 6807 K.

**New 7 (error budget). Mostly addressed.**

Done:
- uncertainties and conventions are separated;
- photon vs energy weighting is given with its full range;
- the hybrid2/quadratic row reads 41–51 ppm;
- "comparable to" replaces "smaller than".

The native-spacing interpolation estimate covers T_eff only (minor comment 1).

**New 8 (free fits). Addressed.**

- Fixed and free coefficients are distinguished.
- Δb is given.
- The c ≤ 1 bound for power-2 is noted.
- Howarth (2011) and Espinoza & Jordán (2016) are credited.

The re-tuned hybrid4 now has a free bias of 12 ppm instead of 34 ppm. On the ranking and the basis dependence of the prior, see minor comment 3.

**New 9 (release vs manuscript). Addressed for the release; not fully for the response.**

The ReadMe now carries:
- the transit-only caveat;
- the bound and unphysical fractions;
- the auxiliary files;
- q and the expanded form.

The rms row is in Table 4. However, the response again claims several changes that are not in the manuscript (New comment 3).

**New 10 (power-2 with c > 1). Mostly addressed.**

Given now:
- the fractions by T_eff;
- the c distribution (median 1.18, > 2 in 3%, up to 76);
- the constrained cost (4.0 → 17 ppm);
- the daggers in Table 3;
- the warning to codes.

Item (ii), a note in the abstract or in Table 2 that the power-2 statistics include the unphysical rows, is not done (minor comment 2).

**Minor comments of the second round.**

Addressed: 1, 4–18, 20–22 and 25.

Partly addressed:
- 3: the citations are in the text but the references are missing; see New comment 3.
- 19: the small-scale-dynamo numbers are given, but not the field-free MURaM ones.
- 24: lower-case law names still begin Fig. 4's caption and the first recommendation item.

Not addressed in the manuscript, despite the response:
- 2: reference-light-curve convergence and the time sampling of the chord;
- 23: LDTk is in the reference list but not cited; Gopal & Trefethen (2019) and Herremans et al. (2023) are cited but absent from the list.

I acknowledge that my comment 3 on the ladder led to the hybrid4 re-tuning. The test design is sound:
- held out by grid;
- random half-splits;
- 8100 independent model–passband combinations.

The result is a clear improvement.

---

### New major comments

**1. The hybrid4 curvature conditions are a modelling assumption, not physical conditions. State their domain of validity and claim "complete"/"exact" proportionately.**

**The derivation is correct and clearly written.** I verified:
- the reduction of I″ ≤ 0 and I‴ ≥ 0 (in x = μ²) to five linear inequalities. The third derivative involves only the two pole terms, whose ratio is monotonic, so endpoint checks suffice. I″ and I′ then follow at x = 1;
- the vertices S3 and S4 to six digits, in both the paper's and the ReadMe's representations;
- by Monte Carlo, that the region holds 8.8% of the positive–monotonic region (you quote 9%), and that 96.3% of the excluded volume is convex in μ² somewhere (96%);
- the q-map (q1 = s⁴ and Dirichlet stick-breaking), which is uniform and bijective, as stated.

**My concern is the status given to the result.** Kipping's (2013) conditions, positivity and monotonicity, are physical: no real star violates them, chromospheric limb brightening apart. Concavity in μ² and a convex slope in μ² are empirical regularities of LTE photospheric profiles in the near-peak and Rayleigh–Jeans regimes. They are useful regularities, and your evidence that the r-weighted fits satisfy them is good, but they are not physical laws.

*(a) Where they hold.* Your case would be stronger if you said so in a sentence or two:
- Every positive, monotonic quadratic-law profile satisfies both conditions. With I = c₀ + c₁μ + c₂μ², both reduce to c₁ = u₁ + 2u₂ ≥ 0, which is monotonicity at the limb.
- A power-2 profile satisfies them if and only if α ≤ 2. I checked the CDS power-2 tables: no PHOENIX or MPS-ATLAS Set 2 row, r- or μ-weighted, has α > 2.
- More generally, writing I(μ) = f(μ), the third derivative is d³I/dx³ = (μ²f‴ − 3μf″ + 3f′)/(8μ⁵). So a profile that is increasing and concave in μ with f‴ ≥ 0 satisfies both conditions. You state the concavity half of this; the I‴ half would complete it.
- The excluded profiles are those locally steeper than μ² (I ∝ μ^β with β > 2 over part of the disc), such as I = μ⁴.
- The observed Sun satisfies both conditions. I evaluated the Neckel & Labs (1994) fifth-order polynomials, as in your repository: at all 30 wavelengths from 303 to 1099 nm, I″ ≤ 0 and I‴ ≥ 0 hold on 0 < μ ≤ 1. The only exception is a marginal I‴ < 0 at μ > 0.98 at 366 nm. This is direct observational support, and I suggest you include it.

*(b) Where they fail.* The domain should be stated.

- **Wien regime.** For a grey LTE (Eddington) atmosphere with T⁴ = ¾T_eff⁴(τ + ⅔), the monochromatic I(μ) = ∫B[T(t)]e^{−t/μ}dt/μ violates the conditions:
  - T_eff = 3000 K: I″ > 0 for μ ≳ 0.13 at 300 nm, and for μ ≳ 0.06 at 200 nm;
  - 4000 K at 200 nm: I″ > 0 for μ ≳ 0.1;
  - 5800 K at 200 nm: I‴ < 0 near the disc centre.

  This is the expected behaviour. B_ν(T(τ)) becomes strongly convex in τ when hν ≫ kT, and the profile then steepens faster than μ². Your 41 passbands are band-integrated and red-weighted for cool stars, so this does not affect the tables. It does affect:
  - the PHOENIX coefficient spectra, which extend to 0.1 μm;
  - any near-UV band built from them for an M or K dwarf (WFC3/UVIS G280 starts at 0.2 μm; there are dedicated NUV transit programmes).
- **Spherical limbs.** Your own node analysis (studies/rev3/h4_scheme.md) finds the following fractions of PHOENIX band profiles:

  | Condition violated somewhere | all μ | μ ≥ 0.05 | μ ≥ 0.1 |
  |---|---|---|---|
  | I″ > 0 | 15.7% | 4.2% | 0.07% |
  | I‴ < 0 | 99.2% | 13.1% | 0.19% |

  The paper quotes only the last column ("above μ = 0.1, 0.2%"). The 0.05 < μ < 0.1 annulus is where hybrid4's inner pole begins to act (√ε₀ = 0.17). Please give all three columns, and say how the derivatives were estimated: third finite differences on 78 non-uniform nodes are noisy.
- **Chromospheres and NLTE.** Emission cores, such as Ca II H&K, Mg II and Hα in active M dwarfs, are limb-brightened. They violate monotonicity itself, so no constrained law here applies. A sentence would prevent misuse of narrow-bin spectra in those lines.

*(c) Wording.* Please make the wording consistent with the above:
- In Sect. 1 and the Conclusions, replace "physical regions are simplices" and "hybrid4's simplex is complete" with something like:
  > "hybrid4's simplex is exactly the set of hybrid4 profiles that are positive, monotonic, concave in μ² and of convex slope; it holds 9% of the positive–monotonic hybrid4 profiles."
- Say plainly that the analogy with Kipping's q1, q2 is one of construction. Kipping's map covers every physical quadratic profile; yours deliberately excludes 91% of the physical hybrid4 profiles.
- "Which no model produces" (Sect. 2) should read "which none of our band-integrated profiles produces".

*(d) Consistency with hybrid5.* hybrid5 keeps a sufficient simplex that contains profiles convex in μ², such as the vertex I = μ⁴. The two laws therefore now carry different implicit assumptions. A user who changes from one law to the other under uniform-q priors changes the prior on the star.
- Please explain why hybrid5 is not treated the same way. With three poles, I‴ is a sum of three terms, the endpoint argument fails, and the region is no longer a simplex.
- Please also say whether the hybrid5 fits satisfy the curvature conditions anyway.

*(e) Free fits.* In Sect. 8 the tables are recommended mainly as priors for free fits. With q ∈ [0, 1]⁴, the curvature conditions become an assumption about the observed star, not only about the models. Please add this to the "Active stars" item. Unocculted spots and faculae change the effective profile, and it is not obvious that the effective profile stays concave in μ².

**2. Some statements about the curvature simplex are contradicted by your own study files.**

**(a) "With a margin" and the "123 fits".** Sect. 2 states:
> "All 2.9 million tabulated hybrid4 fits satisfy the curvature conditions with a margin. The only fits that would gain from dropping them are 123 PHOENIX fits…"

studies/rev4/rv_bound.md shows a different picture:

| PHOENIX fits | Bound | Fit without curvature conditions leaves the region |
|---|---|---|
| r-weighted | 128 (0.19%) | 123 |
| uniform in μ | 12 541 (18.6%) | all 12 541 |

So in 12 664 tabulated fits the curvature conditions are active. These fits lie on a face of the region, with zero margin. studies/rev3/h4_kipping_check.md likewise counts 13 281 PHOENIX rows with a q coordinate exactly 0 or 1.

The 12 541 uniform-μ rows span 2300–7400 K (median 4900 K) and lie mostly in the infrared: IRAC2, G395H/NRS2, F444W and G395M. This is the spherical-limb inflection, which the μ-weighting resolves.

**(b) Which metric "gain" refers to.** "Would gain" holds for those rows only in light-curve terms: the unconstrained-curvature fits have a worse light curve (median 75.9 against 47.3 ppm). The uniform-μ table, however, is recommended precisely for limb-sensitive uses, where the profile residual σ_μ is the relevant metric, and in that metric every one of these fits would gain.

**(c) What I ask:**
- Correct the statement.
- Attribute the 19% uniform-μ binding fraction in Sect. 4.2 to the curvature conditions, not only to "weighting the limb more".
- Give the cost in profile terms. Your earlier analysis (h4_completeness.md, old basis) found an rms change of median 5 × 10⁻⁵ and maximum 3 × 10⁻³ in I/I(1). If the new basis confirms this, it is small and easily stated.

**(d) Circularity.** The Conclusions' "which every tabulated fit is" is circular, since the fits are constrained to the region. The non-trivial statement is that the fits without the curvature constraint satisfy it: all MPS-ATLAS fits, and all but 123 r-weighted PHOENIX fits. Please say that instead.

**(e) The 123 rows.** These rows (2300–2600 K, B, g, G430L, V, r) look like the metal-poor, limb-brightening rows of Sect. 5. Please say whether they are the same rows. If they are, "give up concavity to follow the spherical limb" should be reconciled with the diagnosis in Sect. 5 that these profiles brighten towards the limb.

**(f) Flags and checks:**
- Table 5's note says the free fits are "within each law's exact physical region". For hybrid4 they are within the curvature simplex (studies/rev4/rv_bias.py). It never binds there, so the numbers stand, but the note is wrong for hybrid4.
- The ReadMe's `Phys` flag ("positive and monotonic") and `Bound` flag ("a physical constraint binds") should say that for hybrid4 the binding constraint is in almost all cases a curvature condition.
- The narrow-bin check of Sect. 9.1 (studies/rev3/prod_checks.py) tests positivity and monotonicity only. A summed hybrid4 coefficient spectrum can lie outside the curvature simplex without being flagged, and then has no q. Please report the fraction of such bins, including the UV, and say what the code does with them.

**3. Numbers made stale by the hybrid4 changes, and response items missing from the manuscript.**

*(a) Weighting factors.*
- In the new Table 1 the smallest uniform-μ → uniform-r factor on MPS-ATLAS is 1.9, for the quadratic law (140/74). The 1.7 came from the old hybrid4. "1.7–6" (abstract, Sect. 4.2, ReadMe abstract) and "1.7–15" (Conclusions) should therefore read 1.9–6 and 1.9–15.
- "Against equal weights on the native nodes the factor is up to 23" is now up to 28 on PHOENIX (108/3.9) and 32 on MPS-ATLAS (74/2.3), both for hybrid4.
- The upper ends of both ranges (15 and 6, and 28 and 32) are now set by hybrid4, whose poles were tuned for the r-weighted fit, as Sect. 4.2 itself notes. The weighting comparison is therefore confounded at its upper end. I suggest quoting the ranges without hybrid4 in the abstract (3–9 and 1.9–3.5), or saying there that the upper value reflects the tuning.

*(b) hybrid4 on the Sun.* Sect. 7.1 quotes 3–5 ppm for hybrid5 and the four-term law fitted to the Neckel & Labs polynomials, but nothing for hybrid4. Your studies/p4/p4a_sun.md gives 8–10 ppm for hybrid4, twice the four-term value. That file was last regenerated before the α = 25/16 change. Since the paper now says hybrid4 is as accurate as the four-term law, please recompute it with the adopted poles and basis, and quote it.

*(c) Missing numbers.* See Response 6 above: the MURaM star count (six vs seven), and the Sect. 9.2 → Sect. 8 cross-reference to a number that is not in the paper.

*(d) Response vs manuscript.* As in the second round, the response claims changes that I cannot find:
- **References.** `gopal2019` and `herremans2023` are cited in Sect. 2 but are not in the bibliography; the compiled PDF prints "of ? and ?". `parviainen2015` (LDTk) is in the bibliography but never cited.
- **Light-curve convergence.** The reference-light-curve convergence ("≤ 0.012 ppm, and ≤ 0.25 ppm in the error under chord sampling") and the chord sampling are not in Sect. 4.3.
- **Δw extrapolation.** The ≤ 0.7 ppm size of the constant-ratio extrapolation is not in the manuscript.

Please check the final response item by item against the compiled PDF.

**4. Document the redefinition of hybrid4 for table users and for the companion paper.**

hybrid4 has changed twice since the version I reviewed in round 2: new poles (α = 1 → 25/16) and a new basis. Profiles, light curves and light-curve errors do not depend on the basis. Everything a user reads per coefficient does:
- w₃ and w₄;
- q;
- the Δw corrections;
- the s_w columns of systematics.dat;
- what "a prior of 0.05 per coefficient" means.

**(a) The companion paper.** The relation to the companion paper's hybrid4 appears only in the ReadMe (Note 15, "weights on T … are w B"). The paper describes the companion as deriving "the placement of the poles". Please:
- give the matrix B, or the expanded form, in the paper;
- state that the companion paper uses the same poles and basis.

As in my confidential comments, I recommend that acceptance be conditional on the two papers defining hybrid4 identically.

**(b) The ladder rationale.** With α = 25/16 the ladder discretises the double-pole representation of μ^{2α} = μ^{25/8}, a profile that is itself convex in μ² and is excluded by the new region. That is not an error, because the light curve is elementary for any ε > 0. But "the ladder follows from discretising an exact representation of μ^{2α}" is now a parametrisation of convenience. Please say that, for hybrid4, α is a tuning parameter.

**(c) Grid dependence of the poles.** Your scan (studies/rev3/h4_analysis.md) shows that the two grids prefer different scales:

| | Best s | Error at best s |
|---|---|---|
| PHOENIX | ≈ 0.87 | 1.62 ppm |
| MPS-ATLAS | ≈ 0.77 | 1.09 ppm |

Poles tuned on PHOENIX give 3.4 ppm on MPS-ATLAS, three times that grid's optimum. The adopted s = 0.80 is a reasonable compromise, but hybrid4's accuracy depends on the grid more than hybrid5's does. A sentence would warn users who apply hybrid4 to future grids (3D, or a new PHOENIX release).

**(d) Smaller documentation items:**
- The letter s is used both for the ladder scale (s = 0.8) and for the depth s = Σw (q₁ = s⁴) within one page. Please rename one of them.
- The paper writes S3 and S4 on (1 − μ², (1 − μ²)², Π₀, Π₁). The ReadMe writes them on (1 − μ⁴, (1 − μ²)², Π₀, Π₁). Both are correct (I checked), but one representation, or an explicit statement of both, would avoid confusion.
- Fig. 1(b) uses S₁ and S₂, which the text does not define.
- Fig. 1's caption describes only panel (a).
- The S3/S4 display overflows the column in the PDF.
- The ReadMe Description says that "the hybrid laws are constrained to their physical regions". For hybrid4 it should say "to its curvature simplex (Note 15)".
- Please state in the ReadMe that the hybrid4 Δw (corrections.dat) and s_w (systematics.dat) refer to the curvature-simplex basis.
- h4_kipping_check.md verifies the magnetic corrected sets for hybrid4, but not the 3D sets (MURaM hydro, Stagger) that Sect. 9.2 says remain physical. Please confirm that they lie in the curvature simplex.

---

### New minor comments

1. **Sect. 5 and Table 6, interpolation.**
   - The native-spacing estimate, "1.3 ppm (12)", is for T_eff only. In log g and [M/H] the PHOENIX native step is 0.5 dex, and the doubled-gap values (15 and 13 ppm) scaled in the same way give roughly 4–5 ppm. The budget row should say "T_eff" or combine all three.
   - The median scales as gap^1.7 but the 90th percentile as about gap^1.2 (studies/rev2/c6_interp.md). The text implies one exponent.
   - The text gives 26 ppm for the 90th percentile in T_eff; Table 6's doubled-gap row gives "(28–65)".
   - The MPS-ATLAS [M/H] test uses a native step of 0.05 dex (c6_interp.md), whereas the text says "metallicity steps are 0.1 dex".
2. **Table 2 / abstract, power-2.** Please note in the Table 2 notes that the power-2 statistics include the unphysical (c > 1) rows: 22% of PHOENIX rows. This was second-round comment 10(ii).
3. **Sect. 6 and Sect. 8, free fits.**
   - The ranking "four-term law least biased, then hybrid5 (7 ppm) and hybrid4 (12 ppm)" is by the maximum. By the median, hybrid4 (0.61) is the least biased. On MPS-ATLAS, hybrid4 free (1.7 / 3.9 ppm, rv_bias.log) is worse than the four-term law (0.53 / 1.6), which is not reported. Please give the basis of the ranking.
   - A Gaussian prior of "0.05 per coefficient" depends on the basis. The hybrid4 prior entry of Table 5 changed (0.52 → 0.33 ppm median) when only the basis changed. For the four-term law, whose coefficients are of order unity and strongly degenerate, 0.05 means something quite different. Please say so, or express the prior on basis-free quantities such as h′₁, h′₂ and I(0).
4. **Sect. 2, the implied prior for hybrid4.** You give the implied 90% intervals for hybrid5 but not for hybrid4, the law whose map is new. For hybrid4, q₁ = s⁴ implies:
   - I(0) between 0.013 and 0.53 (90%);
   - P(I(0) > 0.5) = 6%.

   Please give these and the h′₁ and h′₂ intervals. The G395H remark applies to both laws.
5. **Sect. 4.2, uniform-μ hybrid4.** After the re-tuning, the uniform-μ and equal-per-node hybrid4 fits became much worse in light-curve terms:
   - PHOENIX, uniform in μ: 25 → 58 ppm;
   - MPS-ATLAS, equal per node: 15 → 74 ppm.

   You already say that hybrid5 or the four-term law is better for limb-sensitive uses. Please add this to ReadMe Note 2, and give σ_μ for hybrid4 vs hybrid5 in the μ table, since σ_μ is the relevant metric there.
6. **Sect. 7.1.** Two different "small-scale dynamo" numbers appear for the near-infrared:
   - the MURaM small-scale-dynamo profile itself (23 ppm);
   - Set 2 corrected by the small-scale-dynamo ratio (40 ppm).

   The sentence "MPS-ATLAS Set 1 … and MURaM with a small-scale dynamo … are closest" makes them easy to confuse. Please label them distinctly, and give the field-free MURaM values (72/33 ppm) as second-round comment 19 asked. In Table 6, "Models vs Sun (MPS-ATLAS, MURaM) 23–136" mixes direct MURaM profiles with corrected MPS-ATLAS; please say which.
7. **Sect. 2, constraint comparison across laws.** hybrid4 is now constrained by the curvature conditions, whereas the quadratic and four-term laws are constrained only by positivity and monotonicity. The comparisons in Tables 1 and 2 are therefore slightly conservative for hybrid4. Please say so in one clause.
8. **Sect. 9.1.** The narrow-bin "unphysical" fractions (0.25–0.33% for hybrid4) refer to positivity and monotonicity. Please say so, and give the fraction outside the curvature simplex (New comment 2).
9. **Sect. 2, sampling-check counts.** "400 000 uniform draws from the simplex; 82 000 profiles drawn uniformly from the region" do not match the counts in h4_scheme.md (200 000) or h4_kipping_check.md (100 000 per test). Please make the counts in the text traceable.
10. **Style.**
    - Fig. 4's caption and the first recommendation item ("Law, fixed coefficients. hybrid5, hybrid4, …") begin a sentence with a lower-case law name.
    - The main body slightly exceeds 12 pages.
11. **Placeholders.** These are known: the Zenodo DOI, the licence, the author fields, and the companion arXiv identifier. The companion reference must be completed before acceptance, and it must define hybrid4 as here (New comment 4).

---

*Checks performed by the referee in this round:*
- **Hybrid4 region.**
  - Rebuilt the basis from `ldhybrid/interpolate.py`.
  - Confirmed the S3 and S4 coefficients of Sect. 2 and Note 15.
  - Monte Carlo over a box containing the positive–monotonic region (4 × 10⁶ draws): region volume 0.321 (0.322 quoted); curvature simplex 8.8% of it (9%); 96.3% of the excluded volume convex in μ² somewhere (96%).
  - Verified the endpoint reduction analytically and the uniformity of the q-maps for hybrid2, hybrid4 (q₁ = s⁴; Beta(1,3) and Beta(1,2) stick-breaking) and hybrid5.
  - Checked the hybrid5 prior intervals: I(0) 0.010–0.451; P(I(0) > 0.5) = 3.1%.
- **Physical domain.**
  - Evaluated the conditions on all 30 Neckel & Labs (1994) polynomials (satisfied).
  - Evaluated them on grey-LTE monochromatic profiles: they fail at ≤ 300 nm for 3000 K and at 200 nm for 4000–5800 K.
  - Confirmed that no tabulated power-2 row has α > 2.
- **Against the study files.**
  - Compared the 9%, 96%, 13%, 123, 26 → 15 ppm, 0.2% and 19% figures with h4_scheme.md, h4_completeness.md, h4_kipping_check.md and rev4/rv_bound.md. All agree, except the "margin" and "only 123" statements (New comment 2).
  - Compared the pole tuning with h4_analysis.md, h4_validate.md and rev3/results.md: 5.8 → 1.8 ppm out of sample; held-out factors 1.5–2.1.
  - Compared the bias set-up with rev4/rv_bias.py and its log.
- **Manuscript arithmetic.**
  - Recomputed the Table 1 factors (New comment 3).
  - Counted the Table 3 wins (11/15 and 9/15) and the Table 4 medians at p = 0.1 against Table 3.
  - Recomputed the abstract ranges against Table 2, the JWST statistics and Table A.1.
  - Checked the interpolation extrapolation against rev2/c6_interp.md and the solar numbers against p4/p4a_sun.md.
- **Release and build.**
  - Read the CDS ReadMe in full.
  - Cross-checked the citation keys against the bibliography (two missing, one uncited).
  - Viewed all 14 pages of the compiled PDF.
