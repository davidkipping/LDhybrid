```
Manuscript: Limb-darkening coefficients for rational laws with elementary transit light curves (second revision)
Referee: Referee 3 (independent generalist; not a limb-darkening or transit specialist), third round
Recommendation (A&A grade): B   [= Minor revision], conditional on (i) the companion paper being public (arXiv) and defining hybrid4 exactly as this paper does, and (ii) the main body being brought within the 12-page limit
Style and language: (b) Well written, and now self-contained on the points I raised. The main body runs to about 12.4 pages in the A&A layout, so it must be cut by at least 0.5 page, and I would aim for 8-10% (about 1 page). Fig. 4 (hybrid5 weights), the hybrid4 ladder tuning, the derivation and numerical checks of the hybrid4 curvature simplex, and the hybrid2 out-of-sample test can all move to appendices without loss (details in Major comment 1). No language editing is needed.
Identify referee: no
```

---

## Confidential comments to the Editor

The second revision answers my second-round comments well. Four substantive changes stand out:

- **"Elementary" is now checkable.** Sect. 2 gives a short and correct argument (Green's theorem, a rational parametrisation of the boundary arc) for why even powers of μ and double poles give elementary light curves, and why odd powers and simple poles do not. A basis-shape figure has been added.
- **Advantages are separated.** Those that the hybrid laws share with the quadratic and four-term laws are now distinguished from those specific to them.
- **The abstract is more proportionate.** Law accuracy is given in bands. The interpolation comparison now reads "comparable to", and the scaling with gap size was measured rather than assumed. The limb-rescaling result is in the abstract.
- **Free-coefficient recommendations have been added** (Sect. 8). The four-term law is named as the least biased.

Since my last report the authors have also re-tuned the hybrid4 poles (α = 25/16) and put hybrid4 on a new basis, the "curvature simplex". The pole re-tuning is a real improvement: hybrid4 now matches the four-term law with fixed coefficients and was tested out of sample. The curvature simplex is mathematically neat. It is, however, the part of the revision I am least comfortable with, for four reasons:

1. **It is over-sold.** The text calls the region "exact" and "complete". It is exact only for a set of conditions the authors chose so that the region would be a simplex. That set contains two curvature conditions beyond positivity and monotonicity, and it keeps 9% of the positive, monotonic hybrid4 profiles.
2. **It is not explained for a general reader.** The third-derivative condition has no stated physical motivation.
3. **It costs space the paper does not have.** Counted from the PDF in the official A&A layout, the main body now ends about 40% of the way down page 13.
4. **It creates a coherence risk between the two jointly submitted papers.** The authors' own response lists the hybrid4 change as still to be propagated to the companion paper. Two papers that define "hybrid4" differently (poles or basis) would confuse users permanently. I recommend that the Editor check this before acceptance.

The hybrid4 change has also left stale numbers. Some are in the abstract: the weighting factors "1.7-6" and "up to 23" no longer match Table 1. There are also small production problems: two citations print as "?", and an equation runs off the column. None of this is serious, and I see no reason to change my recommendation. The paper remains a useful and honest tables paper, and its finding on the limb-rescaling convention (Sect. 7.3) has value beyond the hybrid laws. All the remaining points are presentation, proportion and consistency, and none needs new large computations.

---

## Report to the authors

### Summary of the revised paper

The paper tabulates coefficients for three rational limb-darkening laws: hybrid2, hybrid4 and hybrid5. Each is an even polynomial in μ plus one to three normalised double-pole terms (μ²+ε)⁻² with fixed ε > 0.

- **Grids and passbands.** The models are PHOENIX-ACES (spherical) and MPS-ATLAS Sets 1 and 2 (plane-parallel), in 41 passbands.
- **Fits.** The fits are least squares uniform in projected radius, constrained to a physical region. Quadratic, power-2 and four-term coefficients are computed in the same way.
- **Scoring.** Every fit is scored by its transit light-curve error.
- **New in this revision:**
  - a Green's-theorem argument for the elementary light curve;
  - a basis-shape figure;
  - re-tuned hybrid4 poles (α = 25/16, tested out of sample);
  - a new hybrid4 basis whose non-negative simplex is exactly the set of hybrid4 profiles that are positive, monotonic, concave in μ² and have a convex slope, with a unit-hypercube map (q₁-q₄) in the manner of Kipping (2013);
  - a restructured error budget that separates uncertainties from conventions;
  - a measured scaling of the interpolation error;
  - a corrected and quantified treatment of the limb-rescaling convention in the JWST comparison;
  - free-coefficient recommendations.
- **Main results.**
  - hybrid5, hybrid4 and the four-term law reproduce model light curves to 1-3 ppm (median, p = 0.1) with tabulated coefficients.
  - With free coefficients the four-term law is the least biased.
  - Model-grid differences (median 87 ppm) and two conventions (limb rescaling, about 100 ppm; photon versus energy weighting, up to 258 ppm) dominate the error budget.

### Response to my second-round comments

**Major 1 (make "elementary" checkable; how users compute light curves now): addressed.**
- The new paragraph in Sect. 2 is what I asked for. I checked it:
  - G(r) is polynomial in r² for (1−r²)^k and rational for (1+ε−r²)⁻².
  - A circle or ellipse has a rational parametrisation in which r² and dθ are rational.
  - A simple pole gives a logarithm in G, which leads to non-elementary integrals along the arc.
- It also explains, in passing, why the poles are double rather than simple, which a non-specialist would otherwise wonder about.
- Fig. 1 shows the shapes.
- The interpolation code now computes the light curve of any tabulated law for a circular planet (Sect. 9.2), which answers the practical question.
- Two small problems with the new material are listed under Minor 1-3.

**Major 2 (advantages unique to the hybrid laws): largely addressed, but the new hybrid4 material re-introduces an over-statement.**
- The split into shared properties (convexity, exact spectral summation, conditional analytic marginalisation) and specific properties is now clear.
- The four-parameter and five-parameter members are said to be "as accurate as the four-term law" only "with tabulated coefficients".
- K = 3 is now justified by the ~1 ppm floor.
- However, the Introduction now states that the hybrid laws' "physical regions are simplices". This is not accurate as written:
  - For hybrid5 the simplex is a sufficient subset filling about 8% of the physical region (Sect. 2 says so).
  - For hybrid4 it is the exact region of a stricter set of conditions, filling 9% of the positive-monotonic region.
  - For hybrid2 the region is an exact triangle, but so is that of the quadratic law (Kipping 2013, whose q₁, q₂ the paper adopts). A triangular region is therefore not specific to the hybrid laws.
- What is genuinely specific is still what I wrote last time: all basis shapes are monotonic, so a simple sufficient region exists, and the fits actually fall inside it. The four-term law could be given the same monotonic-shape simplex (each 1 − μ^{k/2} is monotonic), but its fits have mixed-sign coefficients and fall outside it. See Major comment 2.

**Major 3 (proportionality of the abstract and Conclusions): largely addressed.**
- (a) "Comparable to the interpolation error" is now used, and the scaling was measured. The native-spacing estimate, however, covers T_eff only and does not follow from the stated scaling at the 90th percentile (Minor 9).
- (b) The bands (1-3, 7-8, 40-50 ppm) are clearer. The Results paragraph of the abstract has nonetheless grown to about 17 lines and now carries specialist detail (the μ < 0.02 four-term refit). See Major comment 4.
- (c) The limb-convention result is now in the abstract, and correctly qualified ("accounts for the reported agreement of PHOENIX with JWST").
- (d) Photon versus energy weighting is now called "an error, not an uncertainty" (Sect. 3) and listed under "Conventions" in Table 6. Thank you.

**Major 4 (free-coefficient recommendations): addressed in Sect. 8, not in the Conclusions.**
- The recommendation list now separates fixed from free coefficients and names the four-term law as least biased with free coefficients (2.2 ppm, against 7 for hybrid5 and 12 for hybrid4).
- The Conclusions still say only that "freeing the coefficients of a two-parameter law does not remove its depth bias". The second Conclusions bullet ("hybrid4 ... is as accurate as the four-term law") is true only with fixed coefficients. Please add the free-fit ranking there (one clause).

**Minor comments (second round).** In brief:

| No. | Topic | Status |
|---|---|---|
| 1 | Definitions | μ, COND, SpecInt, NGTS, IRAC, SOSS, LRS and the chord velocity are now defined. Partly: ACES is explained but not expanded, MURaM is not expanded, and Stagger is still not described. |
| 2 | ATLAS9 | Addressed: restored, cited, and the baseline mixing stated. The naming is still mixed: "ATLAS (17 points)" in Sect. 4.4 and "Kurucz" in Appendix A. |
| 3-8 | | Addressed. |
| 9 | Sun | Addressed: sources and the recommended combination are given. This raises a new question (Major 5). |
| 10 | Table 6 labels | Addressed. |
| 11 | MURaM wording | Partly. The 150 K nearest-point tolerance is still not stated. The sentence "Sect. 8 gives the size of this assumption" (Sect. 9.2) points to nothing: Sect. 8 and Table 6 give the size of the 100 G correction, not the effect of holding the MURaM ratio constant below μ = 0.1. The ≤0.7 ppm figure quoted in your response is not in the paper. |
| 12-14, 16, 18 | Figures, closed-form check, cadence | Addressed. |
| 15 | One-column figures and tables | Partly. Fig. 3 is now one column. Tables 1-3 are still full-width; Table 2 would fit in one column. |
| 17 | Prior provenance | Addressed (NIRSpec G395H example). |
| 19 | References | Barros and Akinsanmi are cited. The companion reference is still a placeholder (known). |
| 20 | Placeholders | Still present (known). |

### General assessment (Q1-Q7)

**Q1. Why should this paper be published?** For the reasons I gave in round 2, which still hold:
- uniform, physically constrained, light-curve-scored tables for six laws on two modern grids;
- coefficient spectra;
- corrections;
- a prior code;
- two transferable findings: passband photon weighting, and the limb-rescaling convention that decides whether PHOENIX "fits" JWST.

The re-tuned hybrid4 is a useful four-parameter law with an elementary light curve. Its new basis is a convenience for sampling, not a scientific result, and should be presented as such.

**Q2. Are the assumptions spelled out clearly?** Mostly. The new hybrid4 conditions (concavity in μ² and a convex slope) are assumptions. They restrict free fits and uniform priors, but the physical reason for the third-derivative condition is not given (Major 2).

**Q3. Are the methods fully described?** Yes. One exception is the ladder parameter α, which is used with two meanings (Minor 6).

**Q4. Are the new results adequately emphasised?** Yes, with one caveat. The Conclusions give the hybrid4 basis a full bullet while the abstract (rightly) does not mention it. The Conclusions omit the free-coefficient ranking, which a user needs more (Major 4).

**Q5. Are all figures and tables necessary and properly laid out?**
- Fig. 1 is useful, but its caption describes only panel (a) (Minor 3).
- Fig. 4 (hybrid5 weights, full page width) is the least necessary figure in the main text. Its caption says the weights "should not be interpreted physically", so it belongs in an appendix.
- Table 4 (geometry) could also go to an appendix, with one sentence kept in Sect. 6.
- The other figures and tables are necessary.

**Q6. Electronic-only material.** The CDS and Zenodo plan is appropriate. Because the hybrid4 shapes S₃ and S₄ are fixed numerical combinations, the full-precision constants must be in the ReadMe (as stated) and identical in the companion paper and its code.

**Q7. IAU designations.** These are fine. WASP-94A is now used consistently.

**Length.** I count the main body (Sects. 1-10, with Tables 1-7 and Figs. 1-6) as running from page 1 to about 40% of the way down page 13, i.e. about 12.4 pages in the official A&A layout. My round-2 estimate (page 10) was made on the previous, non-A&A layout. This version is over the 12-page limit, and the new author list and affiliations will add a few lines. See Major comment 1.

### Major comments (third round)

**1. Bring the main body within 12 pages; move derivations to appendices.**

A&A's limit is firm, and this is a tables paper whose core results fit comfortably in 11 pages. I suggest moving the following to a new Appendix B ("Construction of the hybrid laws"), keeping one or two summary sentences of each in Sect. 2:
- (a) the hybrid4 ladder tuning (the paragraph beginning "The value for hybrid4 is a tuned constant");
- (b) the derivation of the five inequalities, the vertex list, the numerical checks (400 000 and 82 000 draws), and the 9%/96%/123-fit statistics of the curvature simplex;
- (c) the explicit q-maps of hybrid4 and hybrid5;
- (d) the hybrid2 out-of-sample paragraph of Sect. 5.

To Appendix C or the CDS ReadMe:
- (e) Fig. 4 (the weights), since a transit constrains only one or two combinations of them;
- (f) the itemised record list of Sect. 9, which belongs in the ReadMe.

Items (a)-(c) and (e) alone save about one page.

This would also improve the order of Sect. 2, which now runs as follows: definitions, ladder, hybrid4 tuning, K, linearity, the hybrid5 simplex, the hybrid4 curvature derivation, interpolation, the elementary argument, the q-maps. The paper's title claim (the elementary light curve) comes after a page of hybrid4 region geometry. A general reader would be better served by this order: definitions, elementary light curve, linearity, physical regions (briefly), sampling.

**2. The hybrid4 curvature simplex: explain why it matters, state what it assumes, and do not over-claim "complete".**

For a non-specialist the motivation can be given in two sentences:
- if the physical region is exactly {w ≥ 0, Σw ≤ 1}, then a box prior on q gives a uniform prior with no rejected draws;
- linear interpolation of tabulated weights cannot leave the region.

That is a real convenience, and the paper should say so plainly. As written, however:

- (a) **"Exact" and "complete" are relative to conditions chosen for the purpose.** The Introduction ("the simplex is exact ... and nothing else") and the Conclusions ("hybrid4's simplex is complete") will be read by most readers as "covers all physical hybrid4 profiles". It covers 9% of them. Please say each time that the region is exact *for the four stated conditions*, two of which (concavity in μ², convex slope) are additional assumptions. Say also that the region excludes 91% of positive, monotonic hybrid4 profiles.

- (b) **Motivate the third-derivative condition.** Concavity in μ² is motivated: it follows from concave and increasing in μ. "Convex slope" (d³I/dx³ ≥ 0) is justified only by "supported by the models" and by the fact that it makes the region a simplex.
  - What does it mean physically for a stellar profile?
  - How strongly do the models support it? The authors' own study (studies/rev3/h4_scheme.md) finds I‴ < 0 somewhere in 99% of PHOENIX band profiles, and I″ > 0 in 16%. These violations are almost all at μ < 0.1, which the text summarises only as "violate them mainly at μ < 0.1".
  - The paper should give these fractions. They are the same limb region in which the paper credits hybrid5 with its PHOENIX advantage.

- (c) **"All 2.9 million tabulated hybrid4 fits satisfy the curvature conditions with a margin."** There are two problems with this sentence.
  - It is true by construction, since the fits are restricted to the curvature simplex (Sect. 4.2), as I noted for "every hybrid coefficient set is physical" in round 1.
  - It appears inconsistent with two other statements. The hybrid4 constraint binds in 19% of the uniform-μ PHOENIX fits (Sect. 4.2). And the next sentence says that 123 PHOENIX fits "would gain from dropping" the curvature conditions (26 → 15 ppm). In a convex fit, if no curvature constraint is active at the optimum, removing it cannot improve the fit.
  - Please report instead which faces bind, how often (per grid and weighting), and at what cost. This is the same information the paper gives for the other laws.

- (d) **Priors.** A uniform prior on q₁-q₄ now imposes concavity in μ² and a convex slope in free fits. For hybrid5 the paper gives the profile intervals that a uniform prior implies (h′₁, h′₂, I(0)). Please give them for hybrid4 as well, and note in the Sect. 8 "Priors" item that the hybrid4 uniform prior carries these shape assumptions.

- (e) **Proportion.** The curvature simplex now occupies about half a page, a figure panel, two numerical equations, five lines of the Introduction and a Conclusions bullet. That is more than the paper gives the error budget's main finding. I suggest:
  - one sentence in the Introduction;
  - one short paragraph in Sect. 2 (motivation, the conditions, the 9% and the cost);
  - the rest in an appendix (Major 1);
  - removing the Conclusions bullet or shortening it to a clause.

**3. Consistency of hybrid4 between the two papers.**

hybrid4 is now defined by new poles (α = 25/16, ε = 0.028601, 0.38582) and a new basis (S₃, S₄ given to six decimals, full precision in the ReadMe). The response to round 2 lists "propagate the new hybrid4 poles to the companion paper" as still open.

The companion is cited in this paper for the ladder and for the closed-form light curves. That includes the 2 × 10⁻⁸ ppm check of the hybrid4 light curves against "the closed-form light curves of SquishierPlanet (2026)". Please confirm in the text that the companion uses the same poles, and either the same basis or the basis-free expanded form, so that tabulated hybrid4 coefficients can be used directly with the companion code. If the companion keeps α = 1, the name "hybrid4" will refer to two different laws.

**4. Coherence of the abstract, Introduction and Conclusions.**

The story is now coherent:
- the laws and why they are elementary;
- tables in the Claret manner, for circular planets mainly;
- r-weighting;
- the law accuracies;
- the dominance of model grids and conventions;
- comparisons with observations;
- the products.

Three adjustments would make it proportionate:

- (a) **Abstract.** The Results paragraph is now the longest part of the abstract (about 17 lines). I suggest dropping the μ < 0.02 four-term sentence (a recipe detail) and the "even polynomials" clause, and adding the free-coefficient ranking in half a sentence.
  - The JWST statement "reduces the offset from 3.4σ to 2.7σ" suggests a meaningful reduction. Sect. 7.3 says the shift is "less than one standard error", so "does not significantly reduce" (or giving both offsets with their errors) would be more proportionate.
  - Likewise, the circularity of the Kepler/TESS 100 G result (the field strength was inferred from the same data) was signalled in the previous abstract ("as found before") but has been lost. Please restore "as found by Kostogryz et al. (2024)".
- (b) **Introduction**, fourth paragraph: see Major 2(e).
- (c) **Conclusions**: add the free-coefficient ranking (see Major 4 of round 2), and say "not significantly" for JWST.

**5. Which magnetic correction should a user apply to an FGK star?**

Sect. 8 recommends MPS-ATLAS Set 2 + 100 G for FGK stars. Sect. 7.1 now shows that this combination is the *worst* MPS-ATLAS option for the Sun in the optical (136 ppm, against 43 uncorrected and 24 with the small-scale dynamo). The text then suggests the small-scale dynamo "for quiet Sun-like stars". Most Kepler and TESS hosts in Maxted's sample are also old, quiet FGK stars, so a user cannot tell from the paper which recommendation applies to their star.

Please either:
- give a criterion (e.g. an activity indicator, or T_eff/age); or
- say that the choice is unresolved and that the difference between the two corrections (up to about 100 ppm in the optical) should go into the prior width.

The current wording leaves the main user-facing recommendation ambiguous.

### Minor comments

1. **Missing references.** "\citet{gopal2019}" and "\citet{herremans2023}" (Sect. 2, the lightning approximations) have no bibliography entries and print as "?" in the PDF. Conversely, Parviainen & Aigrain (2015) is in the reference list but is never cited (the response says LDTk is now cited).

2. **Eq. for S₃, S₄ (Sect. 2).** The two lines overflow the column and are cut off in the PDF ("+0.000778(1 − μ…"). Please break them, or move them to the appendix (Major 1). With the printed six-digit constants, both shapes are very slightly non-monotonic near the disc centre (|dI/dx| ~ 10⁻⁷-10⁻⁶ at μ > 0.985, by my check). Please say that implementations must use the full-precision ReadMe values.

3. **Fig. 1.**
   - The caption describes only panel (a). Panel (b) (the hybrid4 vertex profiles) needs a sentence.
   - The legend uses S₁ and S₂, which are not defined in the text (T = {1−μ², (1−μ²)², S₃, S₄}).
   - Please also correct "falls most steeply near μ ≈ √ε" and "the three poles of hybrid5 act at μ ≈ 0.04, 0.2, and 0.6". For Π_ε ∝ (μ²+ε)⁻², the steepest descent in μ is at μ = √(ε/5) ≈ 0.45√ε, and Π_ε = 1/2 at μ ≈ 0.64√ε. For the hybrid5 poles these are μ ≈ 0.018, 0.09, 0.28 (steepest) and 0.026, 0.13, 0.37 (half drop). Panel (a) itself shows Π_0.3874 already near 0.2 at μ = 0.6. The round-2 wording "at μ ≲ √ε" was correct.

4. **"Conditions I-II" (Sect. 2)** are used three times but never defined in this paper. Please define them once: (I) I(μ) ≥ 0, (II) I non-decreasing from limb to centre, following Kipping (2013).

5. **Stale numbers after the hybrid4 change** (abstract, Sect. 4.2, Conclusions). From Table 1, the uniform-μ to uniform-r reduction on MPS-ATLAS ranges from 140/74 = 1.9 (quadratic) to 14/2.3 = 6.1 (hybrid4). The lower bound "1.7" came from the earlier hybrid4 (20/12). The abstract's "1.7-6" and the Conclusions' "1.7-15" should read 1.9-6 and 1.9-15. Against equal weights on the native nodes, the factor reaches 108/3.9 = 28 (PHOENIX) and 74/2.3 = 32 (MPS-ATLAS) for hybrid4, not "up to 23".

6. **The ladder parameter α (Sect. 2).**
   - α also denotes the power-2 exponent (Sects. 4.2 and 7). Please use another symbol for one of them.
   - The text says that the ladder discretises "an exact representation of μ^{2α}" and also that "α sets how strongly [the poles] cluster". A reader cannot tell whether α = 25/16 means that hybrid4 now approximates μ^{25/8}, or only that the spacing changed. One clarifying clause would help.

7. **hybrid4 medians "on the method study".** Sect. 2 gives 5.0 → 1.9 ppm, Sect. 5 gives 2.1 ppm, and Table 1 gives 3.9/2.3 ppm. These are presumably different metrics or grid combinations. Please state the metric and grid for the Sect. 2 figures.

8. **The 0.05 prior (Sects. 6 and 8, Table 5).** A prior of fixed width "on each coefficient" depends on the basis. It means something different for hybrid4 weights on the new S₃, S₄ shapes than on the Π_ε shapes they replace (same span). It means something quite different again for four-term coefficients, whose values differ by 1-3 between recipes (Appendix A). The "Prior" column of Table 5 is therefore not a like-for-like comparison across laws. Please say so. Better, express the recommended prior in a basis-free quantity (e.g. widths in h′₁, h′₂, which the code already handles).

9. **Interpolation (Sect. 5, Table 6).**
   - The text gives a 90th percentile of 26 ppm for the 200 K T_eff gap, but Table 6 gives "28-65".
   - Under the stated scaling (gap^1.7), 26 ppm extrapolates to 26 × 2^−1.7 ≈ 8 ppm at the native 100 K spacing, not the 12 ppm quoted. Please explain.
   - The native-spacing estimate covers T_eff only. By the same scaling, the doubled-gap log g (15 ppm) and [M/H] (13 ppm) errors for PHOENIX would be about 4-5 ppm at the native 0.5 dex spacing, which is larger than the law errors. Table 6 should say "T_eff" or give all three.

10. **Table 6.**
    - "Magnetic correction, 100 G (size)" is listed under "Uncertainties", but it is a correction, not an uncertainty. Its uncertainty is better represented by the spread between the 100 G and small-scale-dynamo corrections (Major 5).
    - Sect. 8 says that the power-2 law error (7-8 ppm) is "an order of magnitude below ... the size of the magnetic correction" (34 ppm). That is a factor of about 4-5.

11. **Sect. 9.2.** "Sect. 8 gives the size of this assumption": see the round-2 minor 11 entry in the table above. Please add the ≤0.7 ppm figure from your response, or remove the cross-reference. Please also reconcile the number of MURaM stars: the paper says six solar-metallicity stars, the response says seven.

12. **Hybrid4 and hybrid5 polynomial bases differ.** hybrid4 uses {1−μ², (1−μ²)²} and hybrid5 uses {1−μ⁴, (1−μ²)²}, so w₁ means different things in the two laws. A clause saying so, and pointing to the basis-free expanded form in the tables, would prevent mistakes.

13. **Sect. 2, "specific to the hybrid laws".** For hybrid2, the exact triangular region is shared with the quadratic law (Kipping 2013). See the Major 2 entry in "Response to my second-round comments".

14. **Naming.** Please use ATLAS9 consistently: Sect. 4.4 has "ATLAS (17 points)", and Appendix A has "four Kurucz ... models". Please expand MURaM once, and give one descriptive clause for Stagger (a 3D radiation-hydrodynamics grid, Magic et al. 2015).

15. **Table 2** could be set at one-column width. **Table 7** floats onto page 13 after the Conclusions have begun; the typesetter will handle this, but it counts towards the page limit.

16. **Placeholders and companion reference** (author list, affiliations, acknowledgements, Zenodo DOI, licence, companion arXiv number). These are known; I list them only for completeness.
