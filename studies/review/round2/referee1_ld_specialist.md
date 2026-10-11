```
Manuscript: Limb-darkening coefficients for rational laws with elementary transit light curves
            (revised version, second round)
Referee: Referee 1 (theoretical limb-darkening specialist; model-atmosphere LDC tabulations)
Recommendation (A&A grade): B   [= Minor revision; in substance a moderate revision, see below]
Style and language: (a) well written, concise and now self-contained. No reduction needed; the
                    main body is close to, but within, the 12-page limit. Table 3 and Fig. 4
                    still overlap, and one of them could go to an appendix if space is needed.
Identify referee: no
```

---

## Confidential comments to the Editor

The revision is substantial and made in good faith. The authors took my first-round report seriously.

- They recomputed every table with an exact least-squares quadrature, after finding that their own 2000-point resample did not resolve the spherical limb. They report the consequences candidly, including those that weaken their case: hybrid4 is now worse than the four-term law on both grids, and power-2 is unphysical in 22% of PHOENIX rows.
- They added the null laws, the out-of-sample test of the hybrid2 pole, the parameter-bias table, the planet-size and grazing tests, an error budget and a validation appendix.
- They now credit Morello et al. for the weighted-r fit and discuss Sing et al. (2026) directly.

The reframing deals with my main first-round objection. The paper now claims nothing about non-circular planets that it tests itself, and every result in it can be checked without the companion paper. I recommend that acceptance be conditional on the companion paper being publicly available (arXiv). If possible, the Editor could also send the companion manuscript to the referees now.

Several problems remain. None is fatal, but some affect statements in the abstract and the recommendations, so I would like to see the next version.

1. **A demonstrable error in Sect. 7.3.** The paper asserts that the limb conventions which move the JWST limb intensity ℓ by ~0.03 change a transit light curve "by only a few ppm". The convention in question (linear vs geometric μ rescaling) was never tested on light curves. My own estimate with the authors' light-curve code is ~130–190 ppm for p = 0.1, comparable to the grid-to-grid difference.
2. **The JWST offsets are computed one way and described another.** The quoted numbers (+0.035, 2.8σ → 2.0σ) come from hybrid5 fits, whereas the text describes a quadratic fitted to the model profile. Done as described, the MPS-ATLAS offset is +0.043 (3.4σ), per the authors' own validation notes.
3. **The four-term law is not tabulated at its best.** The PHOENIX "advantage" of hybrid5 over the four-term law is a property of the fitting recipe. A four-term fit that ignores μ < 0.02 (r > 0.9998, about one pressure scale height on the Sun) matches hybrid5. A table that also aims to be "a uniform reference for the standard laws" should not ship four-term coefficients that are three times worse than necessary without saying so prominently.
4. **The Claret, Hauschildt & Torres (2025) comparison omits the coefficient level.** The authors' own notes show power-2 differences of |Δc| = 0.17 and four-term differences of order unity, and they could not reproduce that recipe. The appendix gives only light-curve differences and attributes them to model versions.
5. **The recommended grid correction is not tested against the Sun.** In the authors' own solar comparison, the MURaM 100 G profiles are further from Neckel & Labs than either the field-free MURaM or uncorrected MPS-ATLAS Set 2. This is not mentioned.
6. **Release and response do not match the manuscript in places.** The CDS ReadMe lacks the stated transit-only caveat and the three extra files. The hybrid tables lack the q and expanded-form columns described in Sect. 9. The response says rms values were added, but they were not.

All of these can be fixed with modest computation and rewriting. Points 2 and 3 change numbers in the abstract. Hence grade B, but I ask to see the revised version.

---

## Report to the authors

### Summary of the revised paper

The authors tabulate coefficients of three rational limb-darkening laws: hybrid2, hybrid4 and hybrid5, which combine even polynomials in μ² with fixed double poles (μ²+ε)⁻². They also tabulate the quadratic, power-2 and four-term laws, computed the same way.

- **Models and passbands.** PHOENIX-ACES v3.0 (spherical; M1 edge per passband) and MPS-ATLAS Sets 1 and 2 (plane-parallel), in 40–41 passbands.
- **Fitting.** Constrained least squares uniform in r, now evaluated by exact Gauss–Legendre quadrature over the PCHIP profile. Fits uniform in μ are tabulated as well.
- **Fit quality.** Scored by the largest light-curve deviation for p = 0.1. Full-grid medians for hybrid5 are 1.3 ppm (PHOENIX) and 0.85 ppm (MPS-ATLAS); for the four-term law 3.4 and 1.1 ppm.
- **New material:**
  - null laws without poles;
  - out-of-sample tests of the hybrid2 pole;
  - p- and b-dependence;
  - noiseless parameter biases with fixed and free coefficients;
  - an error budget with recommendations;
  - a discussion of Sing et al. (2026), including a reconstruction of their PHOENIX offsets;
  - an appendix validating the pipeline against Kostogryz et al. (2022), ExoTiC-LD, Claret & Southworth (2023) and Claret et al. (2025).
- **Scope.** The main use is now stated to be circular planets; the elliptical-outline algorithm is deferred to a jointly submitted companion paper.

### General assessment

**Q1 (why publish).** The tables are of lasting value, and are now honestly framed and much better documented. Several features go beyond the current standard of LDC tabulation:
- a light-curve error is tabulated for every row;
- the constraint and physicality flags;
- the explicit grid-to-grid systematic;
- the coefficient spectra;
- the photon/energy warning.

The exact reproduction of Kostogryz et al. (2022) is convincing for plane-parallel models. The case for publication no longer rests on the unpublished companion paper.

**Q2 (assumptions).** Now mostly clear. The remaining exceptions:
- the parameter range over which the magnetic and 3D corrections are applied (New comment 6);
- the extrapolation of the MURaM ratios below μ = 0.1 (New comment 6);
- which route is used to compute ℓ (New comment 2).

**Q3 (methods).** Adequate for reproduction, with the exceptions in the minor comments: light-curve sampling and the convergence of the reference light curve.

**Q4 (emphasis).** Much improved. The abstract still states that law errors are "smaller than the errors from interpolation". That compares native-spacing law errors with doubled-spacing interpolation errors (New comment 7).

**Q5 (figures and tables).** Fig. 1 is now explained and Fig. 5 is legible. Fig. 2 shows something the text does not say: the four-term law beats hybrid5 on MPS-ATLAS below ~4700 K (minor comment 5). Table 3 and Fig. 4 still duplicate each other.

**Q6 (electronic material).** CDS: coefficients, corrections, offsets, systematics. Zenodo: coefficient spectra, profiles, code. This is appropriate, but the release must match the paper (New comment 9).

**Q7 (IAU).** No problems. "WASP-94" in Fig. 5 should read WASP-94A, as in Table 6.

---

### Response to my first-round comments

**Major comment 1 (unpublished foundation; practical motivation). Largely addressed, with conditions.**

The reframing is the right decision.
- Every test in the paper uses a circular planet and a numerical quadrature that does not depend on the companion paper's algorithm.
- The elliptical-outline claim is now explicitly attributed to the companion paper.
- The closed-form statements for the standard laws are qualified correctly: Mandel & Agol series, Maxted & Gill, Chishtie et al. (2025, arXiv:2506.18860, which I confirmed exists and treats arbitrary power laws in hypergeometric form).
- The necessity claim is softened.
- The oblateness literature and numerical tools are cited.

Three things remain.

(i) **Companion availability.** Acceptance should be conditional on the companion paper being on arXiv, with the reference completed. The 10⁻⁸ ppm closed-form check (Sect. 4.3, App. A) and the derivation of the ε ladder still rest on it.

(ii) **The "elementary" claim can be made self-contained at almost no cost.** I suggest a few sentences in Sect. 2 or a short appendix:
- The blocked flux is ∮ G(r) dθ along the boundary of the occulted region, with G(r) = ∫₀ʳ I(r′) r′ dr′.
- For I ∝ (1 − r²)ᵏ, G is a polynomial in r². For a double pole (c − r²)⁻², G is rational in r². For μ itself, G contains (1 − r²)^{3/2}. For a simple pole, G contains log(c − r²).
- On a circle or an ellipse, a rational parametrisation makes r² and dθ rational in the parameter. Even polynomials and double poles therefore give integrands rational in that parameter, hence elementary functions. Odd powers give elliptic integrals.

This explains in one line why the poles must be double. It would let readers verify the central property without the companion paper.

(iii) **The case for the hybrid laws for circular planets is overstated (New comment 1).** The quantitative speed and accuracy comparison I requested is deferred to the companion paper. That is acceptable under the new framing, but then the circular-planet advantages listed in Sect. 1 need to be stated correctly.

**Major comment 2 (r-weighting). (a), (b), (c), (e) addressed; (d) not addressed, and the authors' own new results now make it more pressing.**

- (a) Morello et al. (2017, 2020) and Claret et al. (2012) are credited plainly. The factors are presented as quantification, not discovery.
- (b) Table 1 now covers all six laws on both grids. I checked every ratio: uniform-μ → uniform-r medians give 3.0–8.7 on PHOENIX and 1.7–3.5 on MPS-ATLAS, as stated in the abstract, Sect. 4.2 and the Conclusions. The description of Claret & Southworth's 900-point practice is now correct.
- (c) The circularity is acknowledged, and the cost at the limb is quantified (0.11 vs 0.07; σ_μ larger by 1.2 and 1.6). Uniform-μ fits are offered for limb-sensitive uses. The CDS ReadMe I could inspect does not, however, carry this caveat (New comment 9).
- (e) Settled by the exact quadrature.
- (d) The authors decline the SPAM-type comparison because a light-curve fit depends on geometry. But the light-curve metric already fixes a geometry set (p = 0.1, b = 0–0.9). For the linear laws, fitting to that set is itself a linear least-squares problem and costs almost nothing. The new limb-cut test (Sect. 5) shows why it matters: the r-weighted four-term fit is a factor ~3 from what the law can achieve (New comment 3).

**Major comment 3 (equal footing). Addressed.**

- (a) The parameter count is explicit throughout. Hybrid4 is stated to be worse than the four-term law, and the hybrid5 advantage on MPS-ATLAS is stated to be negligible.
- (b) Bound fractions and the costs of the constraints are given for every law. I verified them against the production statistics.
- (c) The null laws are added. The factor "5–30" is 4.7–20 on PHOENIX and 4.5–30 on MPS-ATLAS. "≈5–30" is acceptable, but see minor comment 9 on which metric is used.
- (d) The out-of-sample test of the hybrid2 pole is convincing.
- (e) "At most 3 ppm" has been removed.
- (f) Physicality is now stated as a property of the method.

**Major comment 4 (systematics, biases, geometry, Sing et al.). Mostly addressed; the parameter-bias request is partly addressed.**

- The error budget (Sect. 8) is a welcome addition, though its construction needs work (New comment 7).
- The p- and b-dependence and grazing transits are now covered. I checked Table 4 against Table 3 (medians 0.84, 6.5, 3.2, 4.5 and 76 ppm) and the quoted log-slopes.
- Profile interpolation with refitting is offered and correctly analysed for the linear laws.
- I asked for the biases in Rp/R*, b and a/R*; only the depth bias is reported. The authors' bias study also computed Δb, which reaches 0.05 for the quadratic law with fixed coefficients. That matters for stellar densities and should be in Table 5.
- The discussion of Sing et al. is now direct and fair (New comment 8 asks for some context).
- My point (ii), the native-spacing interpolation error, is acknowledged in the text ("about four times smaller") but not carried into the budget or the abstract (New comment 7).

**Major comment 5 (spherical models). (a), (b), (d) addressed; (c) addressed for plane-parallel models but not for spherical ones; (e), (f) answered but left open.**

- (a), (b) The model identity, M2 and the μ_cri comparison are now specified.
- (d) The edge per passband vs per model is stated, and the chromatic r_c is quantified. Please also give the maximum (New comment 10).
- (c) The Kostogryz et al. (2022) reproduction (to ~10⁻⁵ over 273 280 rows) is exemplary. For spherical models, however:
  - Neither comparison is like-for-like, because the PHOENIX-COND intensities and those of Claret et al. (2025) are not public.
  - The Claret et al. (2025) comparison is reported only in light-curve terms. The coefficient-level disagreement, and the failure to reproduce that recipe, are omitted (New comment 4).
  - In A.3, the 0.015 power-2 difference is called "the size of the COND–ACES difference" with no independent estimate of that difference, which makes the argument circular.
- (e) The τ_slant = 1 edge is declined for lack of structure files. That is acceptable, but the authors' new discussion of rescaling conventions contains an error (New comment 2).
- (f) The brightening of metal-poor cool models is flagged but not investigated. A cheap diagnostic would be to check whether it survives the M2 edge, i.e. whether it lies in the nodes adjacent to the edge. That would tell users whether these rows are an artefact of where the edge falls.

**Major comment 6 (observations). Mostly addressed.**

- (a) The Sun is now described correctly (polynomial representation; flexibility, not the steep limb). However, the MURaM 100 G comparison with the Sun is not reported (New comment 5).
- (b) Kepler and TESS are now presented as a consistency check, with significances, the sample size and the 200 G overshoot. The TESS offsets are not quoted (minor comment 18).
- (c) JWST:
  - WASP-12 is explained.
  - ℓ is defined, but inconsistently with the numbers quoted (New comment 2).
  - The rescaling claim is now demonstrated. Table 6 reproduces Sing et al.'s PHOENIX offsets to 0.002, which is excellent work.
  - Significances are quoted.
- (d) M-dwarf limitation stated. Done.
- (e) The construction of Δw and the linear-Gaussian use of the empirical offsets are now explained clearly, and the `condition()` update is the right tool. The extrapolations behind Δw are not stated (New comment 6).
- (f) Eclipsing binaries declined. Acceptable; it was a suggestion.

**Major comment 7 (methods). Addressed, with two exceptions.**

- (a) Constraint regions and how they are imposed. Done. The coefficient spectra are unconstrained fits, and their unphysical narrow-bin fractions are given.
- (b) Interpolant dependence quantified: Akima and spline ≤ 1.4 ppm, linear-in-μ 2.5 ppm. The convergence of the reference light-curve quadrature for the C¹ PCHIP profile is still not shown (minor comment 2).
- (c) Grid ranges, Set 1/Set 2 differences and the removal of ATLAS9. Done. ATLAS9 is still used as the 1D reference for the Stagger Δw and needs a citation (minor comment 16).
- (d) Bessell curves and the photon/energy effect. Done. A reference is still needed for overriding the SVO flags for TESS and IRAC (minor comment 17).
- (e) Method study defined; LM set-up given. Done.
- (f) The CDS/Zenodo split is appropriate, but the release does not yet match the paper (New comment 9).

**Minor comments.** Most are addressed. Still open:
- 10: the consequence of fixing I(1) = 1 under r-weighting is not discussed;
- 11: the change in the coefficients under flux conservation is not given;
- 12: light-curve rms values are not given in the paper, although the response says they are;
- 16: how the interpolation handles the gaps in the PHOENIX grid (the code uses Delaunay barycentric interpolation for refits; the text says "multilinear");
- 22: the Fig. 2 x-range;
- 24: the duplication of Table 3 and Fig. 4;
- 27: Parviainen & Aigrain (2015, LDTk) is still not cited;
- 29: law names in lower case at the start of sentences.

---

### New major comments

**1. The advantages claimed for circular planets are mostly shared by the standard laws.**

Sect. 1 (para. 4) and Sect. 2 list the advantages of the hybrid laws for circular planets:
- linearity in the coefficients;
- convex constrained fits;
- analytic marginalisation under Gaussian priors;
- exact summation of coefficient spectra.

The quadratic and four-term laws are also linear in their coefficients, and so are Agol et al.'s polynomial laws, for which analytic marginalisation was introduced. The four-term light curve of Mandel & Agol is linear in its c_n. All four properties therefore hold for them too; your own constrained four-term fits are convex QPs.

What is specific to the hybrid laws is narrower, though still useful:
- every basis shape is monotonic, so a *simplex* is a sufficient physical region;
- the model fits almost always lie inside it (binding in < 1% of rows), so interpolation and uniform q-priors stay physical. Four-term fits, by contrast, routinely have mixed-sign coefficients;
- the light curve is elementary.

Please rewrite these passages to make the distinction. As a corollary, consider providing four-term coefficient spectra as well, since the linear-map argument applies equally.

**2. Sect. 7.3: the light-curve effect of the limb-rescaling convention is not "a few ppm", and ℓ is computed inconsistently.**

*(a) Rescaling.*

The paragraph argues that the rescaling convention moves ℓ by as much as the model–data offset, "whereas the same choices change a transit light curve by only a few ppm (Sect. 4)". Sect. 4.4 tests M1 vs M2, edge per model vs per passband, and coarser μ grids. It does not test the linear rescaling μ′ = (μ − μ_cri)/(1 − μ_cri) against the geometric μ′ = √[1 − (r/r_c)²].

The two are very different over the whole disc. For μ_cri = 0.05 they differ by Δμ′ ≈ 0.036 at μ = 0.2 and 0.024 at μ = 0.5. I estimated the light-curve effect with your own `ldhybrid/lightcurve.py` on a smooth quadratic-like profile:
- p = 0.1, μ_cri = 0.047: max |ΔF| ≈ 137 ppm;
- μ_cri = 0.07: ≈ 195 ppm.

This is comparable to the PHOENIX vs MPS-ATLAS difference. The paragraph's conclusion, that "which atmosphere fits JWST" depends on a convention, therefore holds for light curves as much as for ℓ. Please:
- (i) compute the light-curve effect of the linear rescaling for the Table 3 cases and add it to Table 7;
- (ii) remove or correct the "few ppm" statement;
- (iii) state which convention you consider physically correct.

In my view the geometric one is. It is the projection of a sphere of radius r_c, which is how Claret's r-method and Espinoza & Jordán define it. The linear map has no geometric basis away from the disc centre. Taking a position here would turn a "convention" into a finding: Sing et al.'s PHOENIX agreement is an artefact of their rescaling.

*(b) How ℓ is computed.* Sect. 7.3 says the model ℓ is obtained "in the same way" as Sing et al.: a quadratic fitted uniformly in μ to the model profile at μ ≥ 0.2. The caption of Fig. 5, however, says the offsets are "computed through the hybrid5 fits". According to studies/rev/val_sing.md, the quoted statistics (+0.035, 2.8σ; +0.024, 2.0σ with 100 G) come from the hybrid5 route, which fits a quadratic to the r-weighted hybrid5 law. The route actually described gives:
- +0.043 (3.4σ) for MPS-ATLAS Set 2;
- +0.041 for Set 1, which matches Sing et al.'s own +0.042 for their six stars.

The difference, 0.006–0.008, is half a standard error. It arises because the r-weighted hybrid5 fit departs from the MPS-ATLAS profile by a few 10⁻³ at μ ≥ 0.2, and ℓ amplifies this by extrapolation. Please use the model-profile route throughout, since that is the like-for-like comparison with Sing et al. Then recompute the 100 G figure and update the abstract ("from 2.8σ to 2.0σ") and Fig. 5 accordingly. If you keep the hybrid5 route, say so in the text and explain why.

**3. The PHOENIX advantage of hybrid5 over the four-term law comes from the fitting recipe, not the law.**

Sect. 5 now notes that a four-term fit restricted to μ ≥ 0.02 matches hybrid5 (1.3 vs 1.2 ppm). Your limb-cut study shows more:
- The tabulated r-weighted four-term fit (3.7 ppm on the method study) improves by a factor ~3 when μ < 0.02 is excluded.
- Power-2 and the quadratic law also improve slightly.
- Only hybrid5 gains from including the extreme limb.

The reason is clear. With the exact quadrature, the steep drop at r > 0.9998 carries a large share of the r-weighted cost. The fit then trades interior accuracy for it, and the light-curve error is set by the interior (Sect. 5, D12). Hybrid5's innermost pole (√ε₀ ≈ 0.04) absorbs that region; the four-term law cannot.

Physically, μ < 0.02 corresponds to the outer 2 × 10⁻⁴ R*, about 140 km on the Sun, or roughly one pressure scale height. There the profile of a static 1D spherical model depends on the outer boundary, on where the M1 edge falls, and on the neglect of 3D corrugation. Precision there is not physically meaningful.

I ask that:
- (i) the abstract and Conclusions state that the PHOENIX difference between hybrid5 and the four-term law vanishes when the fit ignores μ < 0.02;
- (ii) you consider tabulating the four-term (and power-2) PHOENIX fits with that limb cut, or with the light-curve-optimal coefficients of Response 2(d), or at least provide them as an alternative column. A table that presents itself as "a uniform reference for the standard laws" should not carry four-term coefficients a factor three worse than the law permits;
- (iii) the error budget lists "Law, hybrid5 / four-term" as equivalent, which in practice it already does.

**4. Appendix A: report the coefficient-level comparison with Claret et al. (2025) and avoid the circular "floor" argument.**

The appendix reports that Claret et al.'s coefficients reproduce your model light curves to 24–26 ppm, and calls this the floor set by model versions. According to the validation notes in the repository (studies/rev/val_results.md, A.4):
- Applying their stated recipe (uniform μ, unweighted, unconstrained) to your models gives four-term coefficients different by |Δa_k| = 1–3, and power-2 coefficients different by |Δc| = 0.17 and |Δα| = 0.09. That is ten times your own COND–ACES estimate from A.3.
- Their coefficients sit closer to your r-weighted fits than to your uniform-μ fits, contrary to their stated recipe.
- You could not reproduce their μ_crit strategy.

These facts are the substance of a validation and belong in the appendix. A 25 ppm agreement that coexists with |Δc| = 0.17 does not show that the pipelines agree; it shows that the light-curve metric tolerates large coefficient degeneracies. Please:
- (i) give the coefficient differences (at least for power-2, which is not degenerate);
- (ii) say plainly that the spherical-model comparison is not like-for-like and that Claret et al.'s recipe could not be reproduced;
- (iii) drop the attribution to "model versions" unless you can support it independently.

In A.3, "the size of the difference between the COND and ACES models" needs a number or a reference. Otherwise it is circular.

In my experience, Claret has shared intensity files for specific comparisons. A request for a handful of PHOENIX-COND or 2025 models, giving a true like-for-like test, would settle this.

**5. The recommended FGK grid correction should be confronted with the Sun.**

Sect. 8 recommends MPS-ATLAS Set 2 with the MURaM 100 G correction. Your solar comparison (studies/p4/p4a_sun.md) shows that the MURaM 100 G profiles are further from Neckel & Labs than the field-free and small-scale-dynamo MURaM profiles, and further than uncorrected MPS-ATLAS Set 2:

| Range | MURaM 100 G | MURaM SSD | MPS-ATLAS Set 2 |
|---|---|---|---|
| 450–700 nm | 83 ppm | 52 ppm | 43 ppm |
| 700–1100 nm | 75 ppm | 23 ppm | 53 ppm |

This is not necessarily a contradiction: the quiet Sun's mean field is closer to the small-scale-dynamo case, and Kepler hosts may be more active. But it should be reported in Sect. 7.1 and discussed in Sect. 8. The recommendation should then read as the correction that best matches the Kepler hosts, with the Sun preferring a weaker field. Users of the tables for quiet, Sun-like stars need to know this.

**6. State the extrapolations behind the Δw corrections.**

Sect. 9.2 and the production notes show two extrapolations:
- **Below μ = 0.1.** The MURaM ratios are held constant there. In that region the hybrid5 pole weights w₃–w₅ act, so their Δw are set by assumption rather than by MURaM.
- **Across log g and [M/H].** The magnetic corrections are applied at all log g and [M/H] of the MPS-ATLAS grid within the MURaM T_eff range. MURaM simulations exist only at solar metallicity and specific gravities.

Please:
- state both extrapolations in the paper and in the CDS ReadMe;
- quantify the first, e.g. by comparing a constant ratio with a linear extrapolation in μ of the ratio, in light-curve terms;
- give the light-curve size of the 100 G correction itself, as a line in Table 7.

**7. The error budget (Table 7) mixes uncertainties with conventions and compares unlike spacings.**

- **Interpolation.** It is listed at doubled gaps (4–15 ppm), although the text says it is about four times smaller at native spacing. At native spacing it is ~1–4 ppm, comparable to the law errors of hybrid5 and the four-term law. It is also negligible for MPS-ATLAS in [M/H] (steps of 0.05–0.1 dex). The abstract statement "these law errors are smaller than the errors from interpolation" should then become "comparable to".
- **Photon vs energy weighting.** This is a convention that is either right or wrong, not an uncertainty. Its range is also given only for the broad bands (44–258 ppm), whereas G141 and 2MASS J give 1–22 ppm. List it separately, as a convention that must be right, with the full range.
- **Missing lines.** The table should include the rescaling convention (New comment 2) and the size of the magnetic correction (New comment 6).
- **Hybrid2/quadratic.** The range should be 41–51 ppm (MPS-ATLAS Set 1, quadratic: 51).

**8. Free-coefficient biases: hybrid4 and the recommendations.**

Table 5 shows that hybrid4 with free coefficients biases the depth by up to 34 ppm, with a chromatic spread of 24 ppm. That is twenty times the four-term law (1.8, 3.2) and worse than power-2 in the median. Since the main use of the tables is now as priors for fits with free coefficients (Sect. 8), the recommendation "hybrid4 where four parameters and an elementary light curve are wanted, at a cost of a few ppm" should add that with free coefficients the cost is tens of ppm.

The free power-2 fits are bounded by c ≤ 1, whereas the tabulated power-2 is unconstrained; please say so in the table note.

The argument that "statistical adequacy is not unbiasedness" is correct. It was made earlier, e.g. by Espinoza & Jordán (2016) and Howarth (2011), and should be credited here.

**9. The release and the response do not yet match the manuscript.**

In the release I could inspect (data/release/cds):
- The ReadMe does not carry the caveat that the r-weighted coefficients are optimised for transits and that uniform-μ fits should be used for limb-sensitive applications. The response says it does.
- The ReadMe does not list or describe `corrections.dat`, `offsets.dat` or `systematics.dat`, and these files are absent from the directory.
- The hybrid tables do not contain the q coordinates or the expanded form c₀ + c₂μ² + c₄μ⁴ + Σ d_k(μ²+ε_k)⁻² that Sect. 9 says each record holds.
- Note (2) should mention that the uniform-μ hybrid5 fits bind their constraint in 10% of PHOENIX rows, and that the uniform-μ power-2 fits are unphysical in 29%. That table is the one recommended for eclipsing-binary and interferometric users.

Please bring the release, the paper and the response into agreement. The response also states that rms light-curve values are "reported alongside the maxima". They are not in the manuscript, although your geometry study computed them. Please add them to Table 4 or correct the response.

**10. Power-2 with c > 1.**

That 22% of PHOENIX r-weighted power-2 fits (29% of uniform-μ fits) are unphysical is a striking regression relative to the tables users are familiar with. According to your notes, the LM solutions reach c = 76 with α → 0, i.e. a logarithmic law. Many light-curve codes assume 0 < c ≤ 1, for instance the (h₁, h₂) parametrisation of Maxted (2018), and will fail silently or reject such rows. Please:
- (i) give the distribution of c among the flagged rows in the text;
- (ii) state in the abstract or Table 2 that the power-2 statistics include these rows, and give the constrained values (median 11 ppm on PHOENIX in your full-grid rerun);
- (iii) consider tabulating the constrained fit alongside, or instead of, the unconstrained one.

It would also help to say whether Claret & Southworth's PHOENIX power-2 tables show comparable fractions of c > 1. If not, the difference is a consequence of the edge treatment and the exact quadrature, which users should know.

---

### Minor comments

1. **Sect. 4.2.** "12-point Gauss–Legendre … exact for this profile to rounding". The residual involves rational pole terms and the Jacobian μ/r, so the quadrature is converged, not exact. Please say "converged (unchanged from 6 to 12 points)". The PCHIP slopes depend nonlinearly on the data, so the exact-quadrature fit is not strictly a linear map of the nodes. Qualify the statements in Sects. 2 and 5 (refitting equals interpolating coefficients) accordingly; your measured difference of ~10⁻⁶ shows the effect is negligible.
2. **Sect. 4.3, reference light curve** (first-round 7b). `lightcurve.py` integrates the C¹ PCHIP profile with 48-point Gauss–Legendre on 4 sub-intervals per segment. These are not aligned with the profile nodes, which crowd at r > 0.998. Please show convergence, e.g. n_sub = 4 vs 16 and 48 vs 96 points, for the model light curves themselves, not only for the smooth hybrid5 closed form. Please also state the time sampling of the chord (120 points per half chord in the code), because the narrow features at second and third contact in Fig. 4 may be under-sampled in the maximum.
3. **Sect. 2.** The pole ladder resembles the tapered exponential clustering of "lightning" rational approximation (Gopal & Trefethen 2019; Herremans, Huybrechs & Trefethen 2023, SIAM J. Numer. Anal.). Their construction, clustered poles augmented with a low-degree polynomial, is precisely the hybrid form. Please cite them. Also state for which exponent α the ladder was derived. Since the hybrid2 pole moved from 0.208 to 0.36 when tuned, a one-line test of the "untuned" ladder would reassure readers: scale ln ε_k by ±20% and report the change in the median light-curve error.
4. **Sect. 6, para. 2.** "The four-term law wins in TESS for the M and K dwarfs and the solar analogue." Per Table 3, hybrid5 wins for the late M dwarf in TESS (1.3 vs 1.6 ppm). The four cases are the mid and early M dwarfs, the K dwarf and the solar analogue.
5. **Sect. 5 and Fig. 2.** "hybrid5 is the most accurate law on both grids" holds for the grid medians. Fig. 2 (right), however, shows the four-term law below hybrid5 on MPS-ATLAS for T_eff ≲ 4700 K. Please say so. The MPS-ATLAS panel extends to 7600 K, beyond the grid (7500 K); if these are bin centres, say so in the caption.
6. **Fig. 4 vs Table 3.** The annotated quadratic maxima for the early M dwarf (76 and 69 ppm in G141 and PRISM) differ from Table 3 (77 and 70). The row label reads "Sun analogue"; the text and table use "solar analogue".
7. **Fig. 5.** The JWST error bars are invisible except for HAT-P-11. Say in the caption that the per-star statistical errors (0.001–0.004) are smaller than the symbols. "WASP-94" should read WASP-94A.
8. **Sect. 7.3, statistics.** Besides New comment 2, please say whether the per-star offsets are weighted means over wavelength bins, and with which weights. Your notes also find that the Table 3 means of Sing et al. are not plain means of their per-star values. Please say so, since the "agreement with their result" depends on which of their means is meant.
9. **Sect. 5, null laws.** The medians quoted (41, 4.3, 1.2 ppm for the hybrids) use the p = 0.1 production metric. Table 1 uses p ∈ {0.03, 0.1, 0.15} on the same method study (1.9, 6.5, 56 ppm). Please state the metric in both places, or readers will see two "method study" medians for the same law.
10. **Sect. 4.1.** Give the maximum chromatic r_c variation across the JWST passbands (1.6 × 10⁻³, i.e. 33 ppm for a 1% transit, per your production statistics). Users need to know that some models have a strongly chromatic edge.
11. **Sect. 5, the last paragraph before Table 2.** For the flagged brightening rows, see Response 5(f): a check with the M2 edge would be informative.
12. **Table 4.** For grazing transits the depth is smaller. Please also give the errors relative to the transit depth, or note that the relative error is not smaller.
13. **Table 5.** Add Δb (and Δv/v or Δa/R*) with fixed coefficients; see Response 4. Fixed-coefficient quadratic Δb reaches 0.05.
14. **Sect. 8, Fisher analysis.** Specify which noise level corresponds to TESS and which to JWST, the cadence, and the planet parameters assumed. "To better than 0.05": in which normalisation of the combinations?
15. **Sect. 4.2, constraints.** State that the uniform-μ hybrid5 fits bind in 10% of PHOENIX rows (0.5% for r-weighted). This bears on the limb-sensitive uses recommended for that table.
16. **Sect. 9.2.** The Stagger Δw uses ATLAS9 profiles from ExoTiC-LD as the 1D reference. Please cite the ATLAS9 grid (Castelli & Kurucz 2003, or the version in ExoTiC-LD) and state that the Stagger–ATLAS9 difference is assumed to transfer to MPS-ATLAS Set 2.
17. **Sect. 3.** Please cite the instrument documentation behind treating the SVO-flagged TESS and IRAC curves as photon-counting responses (e.g. the TESS Instrument Handbook; the IRAC Instrument Handbook).
18. **Sect. 7.2.** Quote the TESS offsets as well as Kepler (MPS-ATLAS Set 2: +0.0050 ± 0.0033 and −0.0089 ± 0.0045; with 100 G, −0.0039 and −0.0003). The abstract mentions both missions, but the "known offset" is significant only for Kepler.
19. **Sect. 7.1.** Give the field-free and small-scale-dynamo MURaM numbers alongside the 100 G ones (New comment 5), and the PHOENIX numbers per range.
20. **Appendix A, ExoTiC-LD.** Separate the two causes of the 10⁻² profile difference: the curve differences, and the quadrature, which your notes give as ≤ 1.7 × 10⁻³ (500 vs 20 000 nodes). State the modes concerned. A statement about another public code's convergence should carry its numbers.
21. **First-round minor 10 and 11, still open.** Please state that I(1) = 1 is an interpolation constraint at the disc centre, which carries negligible weight under r-weighting. Please also give the change in the coefficients when flux conservation is imposed, for eclipsing-binary users.
22. **Sect. 6, interpolation.** Describe how the gaps in the PHOENIX grid and the irregular MPS-ATLAS spacing are handled. The code's refit path uses Delaunay barycentric interpolation, the text says "multilinear". Both preserve physicality, but the text should match the code.
23. **References.**
    - Complete the companion-paper reference (arXiv identifier) before acceptance.
    - Add Parviainen & Aigrain (2015, LDTk), Gopal & Trefethen (2019) and Herremans et al. (2023).
    - Chishtie et al. (2025): please check the journal status before publication.
24. **Style.** Several sentences still begin with a lower-case law name ("hybrid5 is the better…", "hybrid4 (3.8–15 ppm)…", "hybrid2 beats…", "hybrid4 and power-2 come next…"). A&A style capitalises or rephrases.
25. **Placeholders.** The red \tbd{} entries (Zenodo DOI, licence), the author and institution fields and the companion arXiv identifier must be completed. These are known and need no further comment.

---

*Checks performed by the referee in this round:*
- **Abstract, text and tables against each other.** I recomputed the Table 1 improvement factors, Table 3 → Table 4 medians, Table 6 differences, the JWST mean, scatter and significance, the Kepler significances, the chromatic-edge depth equivalents and the record counts. Everything agrees except the items in minor comments 4, 6 and 9.
- **Against the authors' study files.** I compared the manuscript with studies/rev/results.md, val_results.md, val_sing.md, production_stats.md and studies/p4/p4a_sun.md. Five items are in the files but missing from the manuscript:
  - the ℓ route (hybrid5 vs recipe: +0.035 vs +0.043);
  - the Claret et al. (2025) coefficient-level differences;
  - the MURaM 100 G solar comparison;
  - the free-coefficient Δb values;
  - the power-2 c distribution.
- **Release and code.** I inspected data/release/cds/ReadMe and the hybrid5 byte-by-byte description, and ldhybrid/lightcurve.py and profile.py (quadrature).
- **Rescaling estimate.** I estimated the light-curve effect of linear vs geometric limb rescaling with ldhybrid/lightcurve.py: 137 ppm (μ_cri = 0.047) and 195 ppm (μ_cri = 0.07) for p = 0.1, b = 0.
- **Literature.** I confirmed Chishtie, Saeed & Goderya (arXiv:2506.18860), exact power-law transits in hypergeometric form, and Herremans, Huybrechs & Trefethen (arXiv:2302.02743), tapered exponential pole clustering with polynomial augmentation.
