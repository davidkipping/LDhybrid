```
Manuscript: Limb-darkening coefficients for rational laws with elementary transit light curves (revised version)
Referee: Referee 3 (independent generalist; not a limb-darkening or transit specialist), second round
Recommendation (A&A grade): B   [= Minor revision], conditional on the companion paper being publicly available (arXiv) before acceptance
Style and language: (a) Well written and concise, and now largely self-contained. About 5% could be saved without losing content: merge Fig. 4 and Table 3 (or move one to an appendix), and move Table 6 with the second paragraph of Sect. 7.3 to Appendix A. No language editing is needed.
Identify referee: no
```

---

## Confidential comments to the Editor

The authors have taken the first round seriously. The revision is substantial and, in my view, honest.

- They found and corrected a real error in their own pipeline: the least-squares integral under-resolved the steep limb of spherical models. They recomputed everything and report the less favourable comparisons that followed. For example, hybrid4 is now worse than the four-term law on both grids, and the hybrid5 advantage on MPS-ATLAS is below their own interpolant ambiguity.
- Almost all the numerical inconsistencies I listed are resolved. I re-derived the main ratios and counts from the tables (weighting factors, passband and record counts, Table 3 ranges, JWST means) and found them consistent, apart from a few small slips listed below.
- The new error budget (Sect. 8, Table 7) and the bias table (Table 5) do what I asked. They turn ppm residuals into consequences and state frankly that the choice of atmosphere matters more than the choice of law when coefficients are fixed.

On the central issue of my first report, dependence on the companion paper, the authors chose to narrow the scope rather than import the companion's content. The paper is now explicitly a coefficient-table paper whose main use is circular planets, and every test uses a circular planet. I consider this an acceptable editorial choice. On that footing, the paper largely stands on its own: the laws are fully defined, the fits and their quality are verifiable from the paper and the CDS tables, and the external validation (Appendix A, against Kostogryz et al. 2022, Claret & Southworth 2023, Claret et al. 2025 and Sing et al. 2026) is good.

Three things still depend on the companion:
- the claim in the title and abstract that the light curves are elementary;
- the light-curve code a user needs to apply the hybrid tables;
- the ladder derivation.

The companion reference still has placeholder authors and no arXiv number. I recommend that the Editor make its public availability a condition of acceptance. I also ask the authors (Major comment 1) for one paragraph that makes the "elementary" claim checkable here, at least for a circular planet. This is a small addition.

What remains is mostly a matter of proportion and wording:
- The abstract still ranks laws at the sub-ppm level and says that law errors are "smaller than the errors from interpolation". That holds only for the doubled-gap interpolation test; at native spacing the authors' own estimate makes the two comparable.
- The Introduction lists, as advantages of the hybrid laws for circular planets, linearity properties that the quadratic and four-term laws share.
- The free-coefficient results (Table 5) favour the four-term law over hybrid5 and strongly over hybrid4. This is not reflected in the recommendations.
- ATLAS9 is still used but is no longer cited.

None of these needs new large computations, so I move from C to B.

The product (uniform, light-curve-scored tables for six laws on two modern grids and 41 passbands, magnetic and 3D corrections, coefficient spectra, and an interpolation and prior code) is of lasting value. The demonstration in Sect. 7.3 that the JWST "which atmosphere fits" conclusion depends on a limb-rescaling convention is a useful result in its own right.

---

## Report to the authors

### Summary of the revised paper

The paper tabulates coefficients for three rational limb-darkening laws (hybrid2, hybrid4, hybrid5). Each law is an even polynomial in μ plus one to three normalised double-pole terms (μ²+ε)⁻² with fixed ε > 0.

- **Grids and passbands.** Coefficients are given for PHOENIX-ACES (spherical) and MPS-ATLAS Sets 1 and 2 (plane-parallel), in 41 passbands.
- **Fitting.** The fits are uniform in projected radius, evaluated exactly between profile nodes, and constrained to physical profiles. Quadratic, power-2 and four-term coefficients are computed the same way.
- **Figure of merit.** Fit quality is measured by the light-curve error. This is supplemented by:
  - parameter biases for fixed and free coefficients (Table 5);
  - the dependence on planet size and impact parameter (Table 4);
  - null laws without poles;
  - out-of-sample tests of the hybrid2 pole;
  - an error budget (Table 7).
- **Comparison with observations.** The predictions are compared with the Sun, with Kepler/TESS (Maxted 2023), and with JWST (Sing et al. 2026). The comparison shows that the PHOENIX JWST discrepancy between studies is entirely due to the limb-rescaling convention.
- **Products.** Magnetic and 3D corrections, coefficient spectra and an interpolation and prior code are provided.
- **Scope.** The paper states that the light-curve algorithm and the elliptical-planet application are left to a jointly submitted companion paper, and that the main use of these tables is for circular planets.

### Response to my first-round comments

**Major comment 1 (self-containment; dependence on unpublished work): partly addressed.**

What has been done:
- The scope is now explicit (Sect. 1, third paragraph), and the companion is described as submitted jointly.
- Sect. 2 now gives:
  - the properties of Π_ε;
  - the monotonicity of every shape;
  - the reason for the polynomial basis;
  - a short explanation of the ladder, with what K controls;
  - the hybrid2 triangle and its edges;
  - the explicit q-maps;
  - the profile prior that a flat q prior implies.
- "Elementary" is defined, and the closed-form statements for the standard laws are now correct.
- The new null-law test (Sect. 5) shows internally what the poles contribute (factors of 5–30). This is the most useful addition for a reader who cannot consult the companion.

What remains:
- (i) The paper still gives no argument, even a schematic one, for why even powers of μ plus double poles in μ² give elementary light curves. Even for the circular planet, which is now the main use, it is asserted with a forward reference. See new Major comment 1.
- (ii) The basis-shape figure was not added. With the explicit definitions this is no longer essential, but a one-column panel of T_j(μ) for hybrid5 would still help a reader see where each pole acts.
- (iii) The companion reference has no authors or arXiv number. Acceptance should wait until it is public.

**Major comment 2 (the elliptical-planet application is never shown): addressed by re-scoping.**

The authors decline to add the oblate-planet demonstration, the timings against squishyplanet/JoJo, and the elliptical-outline accuracy test. They move these to the companion paper and instead restrict the claims of this paper to circular planets. I accept this. The statement in the old Sect. 6 that hybrid5 is "accurate to a few ppm for any elliptical planet outline" has gone, the Conclusions state "We tested only circular planets", and the oblateness literature and the existing tools are now cited.

Two small residual points:
- The statement "Since no oblate exoplanet has yet been detected" (Sect. 1) should be made precise. Rotational oblateness has so far yielded only upper limits (e.g. Zhu et al. 2014; Liu et al. 2025). There are tentative detections of tidal deformation (Barros et al. 2022; Akinsanmi et al. 2024), which I suggested in round 1 and which are still not cited.
- The re-scoping raises the question of what the hybrid laws offer for circular planets. See new Major comment 2.

**Major comment 3 (inconsistent headline numbers): addressed.** In detail:
- (a) "3–10" is replaced by 3–9 (PHOENIX) and 1.7–3.5 (MPS-ATLAS). Both agree with Table 1, now given for all six laws and both grids.
- (b) The 93 versus 170 ppm discrepancy is reconciled. They are different statistics (p = 0.1 over shared grid points versus the worst case over p ≤ 0.15 in the method study), and the text now says so.
- (c) The method study is defined once: (34 + 16) × 8 = 400.
- (d) The passbands are listed, and I count 41. The per-grid coverage (40 for PHOENIX, 41 for MPS-ATLAS) and the record counts in Table 8 are consistent.
- (e) The T_eff ranges are now given per grid.

A few new small inconsistencies are listed under the minor comments.

**Major comment 4 (proportionality of the abstract): largely addressed.** In detail:
- (a) The abstract now quotes grid-wide medians and 99th percentiles, and Sect. 6 says that the 15 cases are examples.
- (b) Physicality is described as holding by construction. Power-2 is now treated fairly: the cost of constraining it is given (4.0 → 17 ppm in the affected rows), as are the bound fractions for every law.
- (c) The comparison is now made at equal parameter count. The four-term law beats hybrid4 on both grids, and hybrid5 beats it on MPS-ATLAS by only 0.25 ppm, below the interpolant ambiguity. The unclear word "bracket" is gone.
- (d) The 100 G result is called a consistency check.

The abstract still has two proportionality problems, which I take up in new Major comment 3.

**Major comment 5 (figure of merit and its practical meaning): addressed.** In detail:
- (a) Morello et al. (2020) is credited directly. The metric is acknowledged to favour r-weighting partly by construction, and the cost at the limb is quantified (limb residual 0.11 versus 0.07; σ_μ 1.2–1.6 times larger).
- (b) The comparison of r and area weighting is settled with the exact quadrature.
- (c) Table 5 now gives the depth biases, including the chromatic spread, and Table 4 the p dependence. The p² scaling was corrected to p^1.1–1.4, which I confirm from Table 4.
- (d) The grid-density weighting is stated.

**Major comment 6 (accuracy below the systematic floors): largely addressed.** The new Sect. 8 and Table 7 are what I asked for: the error sources are ordered, and the fixed versus fitted cases are distinguished. The Fisher-matrix statement that one transit constrains one or two combinations of the weights is useful guidance.

What remains:
- The interpolation error at native spacing is stated as an expectation ("should be about four times smaller"), not measured, and the abstract and Conclusions use the doubled-gap figure.
- There is no explicit advice on which law to fit with free coefficients, although Table 5 answers it.

See new Major comments 3 and 4.

**Major comment 7 (robustness of the reference light curve and of the limb): addressed.** In detail:
- (a) The interpolant sensitivity is quantified (≤ 1.4 ppm for smooth interpolants, 2.5 ppm for linear), and "exact" is qualified.
- (b) The limb weight under r-weighting (0.08% at μ < 0.04) is stated. The text explains that the hybrid5 advantage on PHOENIX comes from the extreme limb (four-term fitted to μ ≥ 0.02 matches hybrid5), which reconciles Fig. 1 with the narrative.
- (c) The sensitivity of ℓ versus the light curve is reconciled in Sect. 7.3 and now forms a clear, demonstrated result (Table 6).

**Major comment 8 (circularity in tuning and testing): addressed.** The hybrid2 pole is tested out of sample by grid and by random half-splits, and the analytic poles of hybrid4 and hybrid5 are noted as untuned. One sentence in this part is hard to follow (minor comment 6).

**Major comment 9 (observational comparisons): largely addressed.** In detail:
- (a) The Sun: the Neckel & Labs representation is now described correctly, and the solar claim is gone from the abstract. Newer solar CLV data reaching closer to the limb were not considered. That is acceptable given the qualified wording.
- (b) Kepler/TESS is presented as a pipeline and atmosphere check, with significances and sample sizes.
- (c) JWST: ℓ is defined, significances are given (2.8σ → 2.0σ), WASP-12 is explained, and error bars have been added.
- (d) The M-dwarf limitation is stated plainly. The preference for Set 2 over Set 1 rests on the Kepler match after the 100 G correction. Since that field strength was itself inferred from the same Kepler data, the preference shares the circularity that the authors acknowledge in Sect. 7.2, and one clause saying so would be enough. It would also help to give, in Sect. 7.1, the solar comparison for the recommended combination (MPS-ATLAS Set 2 + 100 G), which the text does not report.

**Major comment 10 (corrections and physicality): addressed.** The Δw construction, the base models and the physicality check (none of 2.1 million corrected sets leaves the physical region) are now given.

**Major comment 11 (external validation): addressed, and well.** Appendix A reproduces Kostogryz et al. (2022) to rounding. It also compares with ExoTiC-LD, Claret & Southworth (2023), Claret et al. (2025) and Sing et al. (2026). Claret's procedure is described correctly (900 points uniform in μ).

**Minor comments, first round.** Most are addressed. Still open:
- No. 9: several acronyms are still not expanded (see minor comment 1 below).
- No. 28: Figure 2 could be one column.
- No. 30: Figure 4 still duplicates Table 3.
- No. 35: the "Sun analogue"/"solar analogue" naming is still inconsistent in Fig. 4.
- No. 26 and 36: placeholders (author list, Zenodo DOI, licence, companion reference). These are known to the authors, so I note them only for completeness.

### General assessment (Q1–Q7)

**Q1. Why should this paper be published?** As a tables paper it now earns its place. Its strengths are:
- uniform, physically constrained, light-curve-scored coefficients for six laws on two modern grids and 41 passbands, including the JWST modes;
- coefficient spectra that integrate exactly to any bandpass;
- magnetic and 3D corrections;
- a code that turns the tables into priors.

Two findings have value beyond the hybrid laws: that photon versus energy weighting of the passband matters more than any fitting choice, and that the JWST model–data comparison depends on the limb-rescaling convention. The hybrid laws themselves are presented proportionately. hybrid5 is about as accurate as the four-term law, hybrid4 costs a few ppm, and hybrid2 is comparable to the quadratic law. Their distinctive value lies in the elementary light curve, whose proof is in the companion, and in the simplex physical region.

**Q2. Are the assumptions spelled out clearly?** Yes, now. The limb definition, the interpolant, the weighting and its cost, the constraint regions, the tuning set, the grid-density weighting of statistics, and the circular-planet scope are all stated.

**Q3. Are the methods fully described?** Yes for the coefficients. The light-curve method is that of the companion, apart from the reference quadrature, which is described. ATLAS9 is still used (μ resampling test in Sect. 4.4, and the Stagger Δw baseline) but is no longer described or cited.

**Q4. Are the new results adequately emphasised?** Mostly. The abstract over-weights sub-ppm law rankings and does not mention the limb-convention result of Sect. 7.3, which is the most transferable new finding.

**Q5. Are all figures and tables necessary and properly laid out?** Figures 1, 2, 3 and 5 are now properly labelled, and every line style is explained. Table 3 and Fig. 4 still present the same 15 cases. Table 6 and the Sing et al. comparison detail could go to Appendix A. Tables 2 and 3 are set as `table*` but would fit in one column, and Fig. 2 could be stacked into one column. The paper has 8 tables and 5 figures in about 10 pages of main text, so some consolidation would help readability.

**Q6. Electronic-only material.** The CDS tables (18 plus 3 auxiliary files) and the Zenodo archive (coefficient spectra, profiles, code) are appropriate. The DOI and licence remain to be filled in. The CDS ReadMe should repeat the warning that the r-weighted coefficients are optimised for transits.

**Q7. IAU designations.** These are fine. Use WASP-94A consistently: the table says WASP-94A, but Fig. 5 says WASP-94.

**Length.** The main text ends on page 10 of 12, within the 12-page limit. The growth from about 6 to about 10 pages is mostly justified (Sects. 2, 4.2, 6 and 8 and Appendix A).

### Major comments (second round)

**1. Make "elementary" checkable in this paper, at least for the circular case.**

The title, the abstract's first sentence and Sect. 1 rest on the claim that even polynomials in μ plus double poles in μ² give elementary light curves. In this paper the only evidence is the 10⁻⁸ ppm agreement of the quadrature with a closed form that is itself in the companion. One short paragraph in Sect. 2 would make the claim plausible without importing the algorithm. It could say:

- On the sky, μ² = 1 − x² − y², so these intensities are rational functions of the sky coordinates.
- Green's theorem reduces the occulted flux to integrals along the boundary arcs, which are circular or elliptical.
- For rational integrands these integrals reduce to algebraic, logarithmic and inverse-trigonometric functions.
- An odd power of μ brings in √(1 − x² − y²), and with it elliptic integrals.

The authors should check this wording against the companion; the point is to give the reader the mechanism.

Since the light-curve code is that of the companion (Sect. 9.3), please also say how a user who wants these tables for a circular planet computes a hybrid light curve today. Options are the archived reference quadrature or any code that accepts a numerical radial profile. Ideally the interpolation package would include such a function, so that the tables are usable before the companion code is released.

**2. Distinguish the advantages unique to the hybrid laws from those shared with the standard linear laws.**

Sect. 1 (fourth paragraph) and Sect. 2 ("Linearity has four consequences") list, as reasons to use the hybrid laws for circular planets:
- convex constrained fits;
- interpolated coefficients that stay physical;
- analytic marginalisation under Gaussian priors;
- coefficient spectra that integrate exactly to any passband.

The quadratic and four-term laws are also linear in their coefficients. Three of these four properties therefore hold for them too; only power-2 is non-linear. What is genuinely distinctive is narrower:
- (i) the light curve is elementary;
- (ii) all shapes are monotonic, so the simplex is a sufficient physical region, and in practice the fits fall inside it. This is not true of four-term fits, whose coefficients usually alternate in sign.

Please say "as for the quadratic and four-term laws" where it applies, and state the distinctive points explicitly.

The sentence "Their five-parameter member is as accurate as the four-term law or more so" (Sect. 1) needs a qualifier: it holds with tabulated coefficients, but not with free coefficients (Table 5; see Major comment 4). The rationale for stopping at K = 3 (Sect. 2) refers to free fits, whereas the tables are mainly for fixed coefficients. A cleaner justification is that hybrid5 already reaches the ~1 ppm floor set by the interpolant ambiguity and by interpolation.

**3. Proportionality of the abstract and Conclusions.**

(a) *"These law errors are smaller than the errors from interpolation"* (abstract, Conclusions, Sect. 8). The interpolation figures (4–15 ppm) come from the leave-one-out test across doubled gaps. The authors themselves expect the error at native spacing to be about four times smaller, i.e. 1–4 ppm. That is comparable to, not larger than, the hybrid5 and four-term law errors (1–3 ppm). Please either:
- measure the scaling, for example by comparing leave-one-out (gap 2h) with leave-two-out (gap 3h) and extrapolating to h; or
- rephrase as "comparable to the interpolation error".

Table 7 should then list the native-spacing estimate, labelled as an estimate.

(b) *Sub-ppm rankings in the abstract.* The abstract gives eight medians to two significant figures, including 0.85 versus 1.1 ppm on MPS-ATLAS. Sect. 4.3 shows that this difference is below the interpolant ambiguity. For a general reader, a simpler statement would be both clearer and more accurate. For example: hybrid5 and the four-term law reproduce the model light curves to about 1–3 ppm (99th percentile below 6 ppm), hybrid4 and power-2 to 5–8 ppm, and hybrid2 and the quadratic law to 40–50 ppm.

(c) *The limb-convention result.* The finding that the PHOENIX JWST discrepancy between studies is entirely due to the μ-rescaling convention, and that "which atmosphere fits JWST depends on a convention", is well demonstrated (Table 6). It deserves a clause in the abstract. It interests more readers than the law ranking does.

(d) *Photon versus energy weighting as an "error".* In the abstract, Sect. 8 and Table 7 this is listed alongside model-grid differences as a source of error. It is not an uncertainty: for photon-counting detectors there is a right answer, and getting it wrong is a mistake. Please present it as such, e.g. "an error 5–100 times larger than any fitting choice arises if an energy response is used for a photon-counting detector". Separate it in Table 7 from the genuine systematics.

**4. Recommendations for free-coefficient fits.**

Table 5 shows that with free coefficients the four-term law is the least biased (1.3 / 1.8 ppm), ahead of hybrid5 (1.9 / 6.0 ppm), while hybrid4 (10 / 34 ppm) is little better than power-2. Sect. 8 ("Law") recommends hybrid5 or the four-term law "where the light curve must be reproduced to a few ppm" and hybrid4 "where four parameters and an elementary light curve are wanted". These recommendations do not distinguish fixed from free use. Please add one sentence. For example: when coefficients are fitted freely, use the tables as priors (as Sect. 8 already advises), and note that hybrid4 is a poor choice for free fits.

### Minor comments

1. **Definitions (Sects. 1–3).** μ is never defined. Please give μ = cos θ at first use, and its relation μ² = 1 − r² (now only in Sect. 4.2). Still not expanded: ACES, COND, SpecInt, NGTS, IRAC, SOSS, LRS, MURaM (only cited) and Stagger (not described). "Chord velocity" (Sect. 6) needs a short definition.
2. **ATLAS9.** It is used in Sect. 4.4 ("ATLAS (17 points)") and as the baseline of the Stagger Δw (Sect. 9.3), but it is no longer cited or described since it was removed from the release. Please restore the reference (Castelli & Kurucz 2003, or the ExoTiC-LD bundle source) and use one name consistently ("ATLAS" versus "ATLAS9"). Please also comment briefly on adding a Stagger − ATLAS9 shift to MPS-ATLAS Set 2 weights, i.e. mixing baselines.
3. **Sect. 2.** "It binds in fewer than 1% of them" holds for the r-weighted fits only. For the tabulated μ-uniform fits the hybrid5 constraint binds in about 10% of PHOENIX rows (quadratic about 13%, according to the authors' own production statistics). Please say "r-weighted" here, and give the μ-uniform fractions in Sect. 4.2, since those rows are also tabulated.
4. **Sect. 2, hybrid2 tuning.** "7524 model–passband combinations (198 … models in all their passbands)": 7524 = 198 × 38, i.e. the earlier passband set. Please say "in 38 passbands".
5. **Sect. 5, null laws.** The hybrid medians quoted here (41, 4.3, 1.2 ppm "on the method study") differ from the r-uniform medians for the same method study in Table 1 (56, 6.5, 1.9 ppm). Presumably the former are for p = 0.1 and the latter are the worst case over p ≤ 0.15. Please say so. Also, "PHOENIX; MPS-ATLAS similar" hides the fact that the upper end of "5–30" comes from MPS-ATLAS (24 → 0.79); on PHOENIX the factors are 4.7–20.
6. **Sect. 5, hybrid2 pole.** The sentence beginning "Rescanned with the exact quadrature…" is hard to follow. It needs to say:
   - which statistic is being minimised;
   - which sample is used (the 99th percentile of 115 ppm does not match Table 2's 112 ppm);
   - what "102 and 99 ppm at the 90th percentile for ε = 0.32 and 0.40" is to be compared with.
7. **Sect. 6.** "The four-term law wins in TESS for the M and K dwarfs and the solar analogue": in Table 3, hybrid5 wins for the late M dwarf in TESS (1.3 versus 1.6 ppm). The count of 11 of 15 is right, but the description should read "for the mid and early M dwarfs, the K dwarf and the solar analogue".
8. **Sect. 6 / Table 5.** Please state the noise model used in the fits (unweighted least squares on noiseless curves?), the cadence, and which parameters were free besides p, b and the chord velocity (e.g. was a/R⋆ or the mid-time fixed?).
9. **Sect. 7.1.** Table 7 lists "Models vs Sun (MPS-ATLAS, MURaM) 23–150 ppm", but the text gives only 25 and 23 ppm for the best models and 156–307 ppm for PHOENIX. Please give the source of the upper value, and add the result for the recommended MPS-ATLAS Set 2 + 100 G combination.
10. **Table 7.** The caption says "medians", but the photon/energy and Sun rows are ranges over cases. "Model grid (PHOENIX vs MPS-ATLAS) 87" is the median of the largest pairwise difference among all three grids, including Set 1 versus Set 2 (Sect. 5). Please label it accordingly. The "hybrid2 / quadratic 41–49" row should be 41–51 if MPS-ATLAS Set 1 is included (Table 2).
11. **Sect. 9.3.** "The ratio of the hydrodynamic MURaM profile to the MPS-ATLAS Set 2 profile (the code that synthesised the MURaM spectra)" reads as if MPS-ATLAS were a hydrodynamics code. Please rephrase, e.g. "MPS-ATLAS, which also computed the MURaM spectra". Please also state the T_eff matching tolerance ("nearest grid point within 150 K", if that is what was done), and that the MURaM ratios held constant below μ = 0.1 make the corrected inner-pole weights unconstrained by MURaM.
12. **Fig. 4.** Several points:
    - the row label is "Sun analogue" while the text and Table 3 say "solar analogue";
    - the panel annotations "max |ΔI|" should be "max |ΔF|" (flux, not intensity);
    - the annotations differ from Table 3 by rounding in two places (early M: quadratic 69 versus 70 in PRISM, 76 versus 77 in G141).

    With 15 dense panels the figure largely repeats Table 3. Showing one star (e.g. the mid M dwarf) in the three passbands, with Table 3 kept for the numbers, would be clearer and would save about half a page.
13. **Fig. 5.** In the JWST panel the per-star error bars are visible only for HAT-P-11. If the others are smaller than the symbols, please say so in the caption. Please also use "WASP-94A", as in Table 6.
14. **Fig. 3.** The PHOENIX w₃ panel shows the jumps described in the text. A sentence in the caption warning against physical interpretation of the individual weights would help readers who look only at figures.
15. **Fig. 2 and Tables 2–3.** These could be set at one-column width (Fig. 2 stacked vertically), in line with A&A guidance.
16. **Sect. 4.3.** Please state how many cases the 10⁻⁸ ppm quadrature check covers (first-round minor comment 15).
17. **Sect. 2, priors.** "Cover the values of FGKM dwarfs": please say which grid and passbands these values are drawn from.
18. **Sect. 8, last paragraph.** The cadence remark is welcome. "Upper limits for long cadences" is true only for the residual, not necessarily for the bias in fitted parameters; please qualify.
19. **References.** Please complete the companion reference before acceptance. Barros et al. (2022) and Akinsanmi et al. (2024) are relevant to the statement on deformed planets in Sect. 1.
20. **Placeholders.** The author list, affiliations, acknowledgements, Zenodo DOI and licence (red \tbd{} marks) remain to be completed.
