```
Manuscript: Limb-darkening coefficients for rational laws with elementary transit light curves
Referee: Referee 3 (independent generalist; not a limb-darkening or transit specialist)
Recommendation (A&A grade): C   [= Major revision]
Style and language: (a) well written and concise, but NOT yet self-contained. The text is compressed to the point of opacity. It needs about 1.5-2 pages more (definitions, motivation, limitations), not fewer, and still stays well inside the 12-page limit.
Identify referee: no
```

---

## Confidential comments to the Editor

This is a competent and carefully executed tabulation paper in the tradition of Claret's coefficient tables. It covers three model grids, 38 passbands, six laws, two fitting weightings, magnetic and 3D corrections, and coefficient spectra. A product like this would be useful to the transit community and probably to the CDS. The numerical work looks careful: limb treatment, photon weighting, constrained convex fits, interpolation tests, sensitivity to limb choices. I found no sign of a fundamental error in the method.

My concerns are about **self-containment, motivation, and proportionality of the claims**, which the Editor asked a general referee to judge.

1. **The paper rests on an unpublished paper.** The paper exists to tabulate coefficients for laws that are defined, motivated, and validated (closed-form light curves, the pole "ladder", the claim that standard laws are not elementary for ellipses, the q-sampling) only in a reference "in prep." with no authors. As submitted, a reader cannot check the paper's raison d'être. Either the companion paper must be publicly available (arXiv at least) before acceptance, or this paper must carry enough of the argument to stand on its own. I recommend both.
2. **The headline application is never shown.** No light curve of an elliptical planet appears anywhere. All tests use circular planets. The "why should one care" case (how large oblateness signals are compared with these law errors, and why an elementary solution beats existing numerical tools such as squishyplanet) is not made.
3. **Several headline numbers are not supported by, or are inconsistent with, the paper's own tables.** The abstract's "factors of 3-10" from the weighting does not match Table 1, which gives factors of about 2-4. There are also inconsistent sample sizes and model-to-model differences, and a passband count that does not add up. The "at most 3 ppm" abstract claim comes from 15 hand-picked cases, while the full-grid maximum is 50 ppm.
4. **The ppm-level ranking of the laws is at or below the paper's own systematics:** limb definition (about 5 ppm), interpolation (5-15 ppm), model-to-model differences (about 100-170 ppm), and observed offsets (hundreds of ppm). The paper does not discuss what the few-ppm advantage of hybrid5 over the four-term law means in practice.
5. **The "observed solar limb darkening" test is weaker than presented.** Neckel & Labs (1994) give fifth-order polynomial fits that exclude μ < 0.12, which is exactly where the hybrid poles act. The Kepler/TESS and JWST comparisons test the model atmospheres, not the laws, and largely reproduce published results.

None of this is fatal, and all of it is fixable with analysis the authors can clearly do. I therefore recommend **major revision (C)** rather than rejection. The acceptability of the paper depends strongly on the status of the companion paper, and the Editor may wish to make public availability of that work a condition.

---

## Report to the authors

### Summary of the paper

The authors give theoretical limb-darkening coefficients for three "hybrid" rational laws (hybrid2, hybrid4, hybrid5). Each law is an even polynomial in μ plus one to three normalised double-pole terms (μ²+ε)⁻² with fixed ε. The light curves of these laws are said to be elementary for a planet with an elliptical outline. The coefficients come from PHOENIX-ACES (spherical), two MPS-ATLAS sets and ATLAS9 (plane-parallel) intensities, integrated over 38 passbands. The fits are least squares uniform in projected radius, constrained to positive, monotonic profiles. Quadratic, power-2, and four-term coefficients are computed in the same way for comparison.

The main results are:

- Weighting uniformly in r rather than in μ reduces the light-curve error of all laws.
- hybrid5 reproduces model light curves to a few ppm, better than the four-term law.
- All hybrid coefficients are physical.
- The model-to-observation offsets (Kepler, TESS, JWST, Sun) belong to the atmospheres, not to the laws, and a 100 G MURaM correction removes the Kepler/TESS offset.

Magnetic and 3D corrections, coefficient spectra at R = 3000, and an interpolation code are also provided.

### General assessment

**Q1. Why should this paper be published?** Coefficient tables of lasting value are explicitly welcome in A&A. If the hybrid laws become the standard for non-circular planets, these tables (and the coefficient spectra) will be used. Two parts also have value independent of the hybrid laws: the systematic study of fit weighting measured by light-curve error, and the uniform recomputation of the standard laws with physical constraints. Validity appears high at the level of the numerical method.

However, the case for *why* the community needs these particular laws is asserted, not demonstrated (Major comments 1-3). The practical significance of the accuracy gains is not put in context (Major comments 4-6). As written, a non-specialist A&A reader cannot tell whether this is a niche product for a handful of oblateness measurements or a general replacement for Claret's tables. The paper must say which, and argue for it.

**Q2. Are the assumptions spelled out clearly?** Partly. Several assumptions that matter at the claimed ppm level are not stated or tested:

- the monotone-cubic interpolant used to define the "exact" model profile and to resample before fitting;
- the choice of a single limb definition;
- that accuracy for circular planets transfers to elliptical outlines;
- that a maximum light-curve deviation for p = 0.1 is the right figure of merit;
- the tuning of the hybrid2 pole on the same grids used to evaluate it.

See Major comments 3, 6, 7 and 8.

**Q3. Are the methods fully described?** The fitting and limb procedures are described adequately for a specialist. The laws themselves are not: the origin of the pole ladder, the "exact double-pole representation of μ^{2α}", and the light-curve solution are deferred to the in-prep paper. Undefined items include method M2, ℓ, SPAM, ODF, "relative flux error", the "method study" subset, and the grid on which hybrid2 was tuned.

**Q4. Are the new results adequately emphasised?** Over-emphasised in places (abstract numbers; "every hybrid coefficient set is physical", which is true by construction) and under-emphasised in others. In particular, the finding that limb-cut and rescaling conventions move JWST-relevant quantities by as much as the model-data offsets is potentially important and is buried in one sentence (Sect. 7.3).

**Q5. Are all figures and tables necessary and properly laid out?** Mostly necessary. Several have legend or caption problems (minor comments). Table 3 and Fig. 4 largely duplicate each other. All five figures are set at full page width; Figs. 2 and 5 would fit in one column. A figure of the basis shapes T_j, and one of an elliptical-planet light curve, are missing and would do more for the reader than Fig. 3.

**Q6. Electronic-only material.** The coefficient tables and coefficient spectra belong at the CDS, as proposed. The interpolation code should be archived with a persistent identifier (e.g. Zenodo DOI), not only "made available". Table 3 could move to an appendix if Fig. 4 is kept, or vice versa.

**Q7. IAU designations.** Object names (HAT-P-11, WASP-6, -80, -94, -107, -121, Kepler, TESS hosts) are standard. No issues.

### Major comments

**1. Self-containment and dependence on unpublished work.** The paper cites the in-prep work (SquishierPlanet 2026, no authors) for every load-bearing statement about the laws:

- that odd powers of μ give elliptic integrals and fractional powers give no closed form for an elliptical occultor;
- that double poles "restore the steep limb";
- that the hybrid laws "rival the four-term law";
- the pole ladder ε_k = exp[2π(√(k+1/2) − √K)] and its derivation from "an exact double-pole representation of μ^{2α}";
- the q-sampling coordinates;
- the closed-form solution used to validate the quadrature.

A referee cannot assess these, and A&A policy is that in-preparation references are not evidence. I ask the authors to:

(a) make the companion paper publicly available (arXiv) and cite it properly before this paper is accepted; and

(b) independently of (a), add a short section (about one page) that lets this paper stand alone. That section should:
- state the light-curve integral for an elliptical occultor and indicate, at least schematically, why even powers of μ and double poles in μ² yield elementary functions, while μ^{odd} gives elliptic integrals;
- explain in two or three sentences where the ladder comes from (it looks like an exponentially spaced quadrature of a Stieltjes-type integral representation of μ^{2α}; if so, say so, and say what K controls);
- show the basis shapes T_j(μ) in a small figure;
- explain what "poles beyond the stellar limb" means (poles at μ² = −ε, outside the physical domain).

**2. Motivation: the elliptical-planet application is asserted but never shown.** The title, abstract and introduction justify the laws by elementary light curves for elliptical planets. Yet every test in the paper uses a circular planet of p = 0.1, and the hybrid closed form is used only once, as a check on the quadrature. The sentence in Sect. 6 that hybrid5 is "accurate to a few ppm for any elliptical planet outline" is therefore an extrapolation. Please:

(a) show at least one worked example: an oblate (f ≈ 0.1) planet with an inclined spin, for a solar analogue and an M dwarf. Give the oblateness signal (oblate minus equal-area circle) next to the light-curve error of each law for the same geometry. That is the comparison a reader needs to judge whether a 40-100 ppm (quadratic) or 5-10 ppm (four-term) law error matters for the science case. It also addresses whether law errors, which also peak at ingress and egress, are degenerate with oblateness;

(b) verify that the law-versus-model error for an elliptical outline is comparable to that for a circle. Ingress of an elongated outline samples the limb differently, and the profile at μ ≲ 0.04, where hybrid5 departs most from the model (Fig. 1 inset), may be weighted differently;

(c) explain the advantage of an *elementary* solution over existing practice. The authors themselves compute light curves "by Gauss-Legendre quadrature with endpoint substitutions, accurate to 10⁻⁸ ppm", which is apparently fast and precise for any profile. Existing codes already model non-spherical planets with polynomial limb darkening, e.g. squishyplanet (Cassese et al. 2024, JOSS 9, 6972) and the algorithm of Liu et al. (2025, AJ; arXiv:2406.11644). If the gain is speed, differentiability, numerical stability, or the absence of model-dependent quadrature error, please quantify it (timings and precision against these tools);

(d) cite the relevant oblateness literature, e.g. Seager & Hui (2002), Barnes & Fortney (2003), Carter & Winn (2010), Zhu et al. (2014), Akinsanmi et al. (2020), Cassese et al. (2024), Liu et al. (2025), and tidally deformed ultra-hot Jupiters (e.g. Barros et al. 2022; Akinsanmi et al. 2024). The introduction should also state how many systems this is expected to matter for.

**3. Inconsistencies between headline numbers and the paper's own tables.**

(a) *"Factors of 3-10" (abstract, Sect. 4.2, Conclusions).* From Table 1, going from "uniform in μ" to "uniform in r" reduces the median error by 13/4.9 = 2.7 (hybrid5), 55/13 = 4.2 (four-term) and 57/22 = 2.6 (power-2). The worst-case reduction is 3.0, 3.7 and 2.1. The factor-of-10 range is reached only against "equal per native node" (e.g. 247/22 = 11), which is not "uniform in μ". Either the claim refers to laws or samples not shown (quadratic, hybrid2/4, the full grid), in which case show them, or it must be corrected to about 2-4.

(b) *Model-to-model difference.* Sect. 4.4 states that PHOENIX and MPS-ATLAS "for the same star differ by a median of 170 ppm". Sect. 5 states that PHOENIX and the two MPS-ATLAS sets differ by "a median of 93 ppm (90th percentile 285 ppm)". Please reconcile, or state precisely what is compared in each case (which pairs, which subset, which law).

(c) *Sample sizes.* Table 1 is based on "34 models × 8 passbands" (272 combinations), while Sect. 4.2 cites "400 model-passband combinations" in "our method study". The hybrid2 pole was tuned on "7524 model-passband combinations" from an unspecified grid. Please define the method-study sample once and say which grid(s) and passbands were used for tuning.

(d) *Passbands and grid coverage.* I count 37 passbands in the list of Sect. 3, not 38. Yet Table 4 and the record counts (1682 × 38 × 2 = 127 832; 16 974 × 38 × 2 = 1 290 024) imply 38 for both PHOENIX and MPS-ATLAS. This conflicts with the statement that MIRI LRS "lies beyond the PHOENIX data" and is covered only by ATLAS9. Please say which passbands each grid covers.

(e) *Parameter ranges.* Sect. 3 restricts all grids to 2300 ≤ T_eff ≤ 7400 K, but Table 4 and Fig. 2 give 3500-7500 K for MPS-ATLAS.

**4. Proportionality of the abstract claims.**

(a) "For a central transit of M to G dwarfs its light curves deviate from those of the model profiles by at most 3 ppm, against up to 115 ppm for the quadratic law." This rests on 15 hand-picked cases (five solar-metallicity PHOENIX stars, three passbands). Over the full PHOENIX grid (Table 2) the hybrid5 maximum is 50 ppm and the quadratic maximum is 508 ppm. Please quote grid-wide statistics in the abstract (median / 99th percentile), or state clearly that these are selected examples. Please also explain how the five stars were chosen.

(b) "Every hybrid coefficient set is physical." This holds by construction, because the constraint is imposed. It is not a result. The meaningful statement is how often the constraint is active and what it costs (Sect. 4.2 says at most 1%). Please rephrase. The contrast drawn with power-2 ("unphysical in 8.6% of PHOENIX rows") is fair only if power-2 is given the same opportunity. It was fitted *unconstrained*, while every other law was constrained. Either fit power-2 with its (simple) physical constraint, or report the light-curve error of the unconstrained power-2 fits separately for the physical and unphysical rows. Also report whether the 8.6% depends on the r-weighting (i.e. the fraction for μ-uniform fits).

(c) "Better than the four-term law" with medians of 2.7 vs 8.0 ppm (PHOENIX) and 0.95 vs 1.4 ppm (MPS-ATLAS). On the recommended grid the difference is 0.45 ppm, which is negligible in practice. hybrid5 also has one more free parameter than the four-term law. The like-for-like comparison is hybrid4 against four-term, and there the four-term law wins on MPS-ATLAS (1.4 vs 6.8 ppm median) and for all M dwarfs in Table 3. The abstract should state this balance plainly. The phrase in Sect. 6 that the two- and four-parameter hybrids "bracket the standard laws they replace" is unclear and should be replaced by a direct statement.

(d) The abstract states that the 100 G correction "removes" the Kepler/TESS offset. This reproduces Kostogryz et al. (2024), who chose the field strength using essentially the same Kepler measurements (Maxted 2023). It is a consistency check, not an independent validation, and the abstract should say so.

**5. Choice of figure of merit and its practical meaning.** The "light-curve error" (maximum |Δflux| for p = 0.1, four impact parameters) is a sensible and welcome choice, but some context is needed.

(a) The weighting "uniform in r" (or area) is effectively selected to minimise this metric. Its superiority over μ-uniform weighting *on this metric* is therefore partly by construction, and largely anticipates the "weighted-r" recipe of Morello et al. (2020), which should be credited more prominently than "consistent with". Please state whether r-weighting worsens any other metric (e.g. σ or the intensity residual near the limb, the secondary-eclipse or phase-curve brightness distribution, or the h′₁, h′₂ that observations constrain).

(b) Table 1 shows area-uniform weighting to be better than r-uniform (3.9 vs 4.9 ppm median for hybrid5, about 20%). Yet r-uniform is adopted on the grounds that the two are "equivalent within the noise". What noise? Please justify the choice or adopt area weighting.

(c) What users ultimately need is the bias in fitted parameters (R_p/R_⋆, b, a/R_⋆, and the *chromatic* bias in transmission spectra), not a peak residual. A short demonstration of bias for one representative case (in the spirit of Espinoza & Jordán 2016) would turn ppm numbers into astrophysical consequences. Please also note explicitly that errors scale as p² (Sect. 6), so that for a typical sub-Neptune (p ≈ 0.03) all laws except quadratic and hybrid2 are at the sub-ppm level.

(d) Grid-wide medians weight models by the grid density (e.g. 46 metallicities for MPS-ATLAS), not by astrophysical relevance. A sentence acknowledging this would suffice.

**6. The claimed accuracy sits below the paper's own systematic floors.** The paper reports several systematic errors larger than the law-to-law differences it emphasises:

- limb-definition choices change light curves by a median of ≲5 ppm (90th percentile ≲8 ppm) (Sect. 4.4);
- leave-one-out interpolation errors are 5, 15 and 13 ppm in T_eff, log g and [M/H] (Sect. 5);
- PHOENIX and MPS-ATLAS differ by about 100-170 ppm;
- the observed offsets correspond to hundreds of ppm.

hybrid5 and the four-term law differ by a median of about 5 ppm. A reader is left asking what accuracy a user actually gets from these tables, and whether the law choice matters at all once model coefficients are adopted. Please add a short Discussion or "Limitations and recommendations" section that orders these error sources and states the practical conclusion. My reading is: when coefficients are *fixed* from models, the atmosphere dominates and the law is secondary. When coefficients are *fitted*, the law's flexibility matters, but then the weak determination of w₃ (Sect. 4.4, Fig. 3) and the five-dimensional degeneracy are a real issue. Guidance on which hybrid law to fit freely, and with what priors, is needed. Please also note that the leave-one-out test doubles the grid spacing, and estimate the interpolation error at native spacing, which should be roughly four times smaller for linear interpolation of smooth functions.

**7. Robustness of the "exact model light curve" and of the limb region.**

(a) Both the reference light curve ("exact light curve of the tabulated model profile") and the 2000-point resample used for fitting rely on monotone cubic interpolation between the μ nodes. For MPS-ATLAS (24 nodes, plane-parallel) the differences between laws at the 0.5-1 ppm level (Table 2) could depend on this interpolant. Please repeat a subset with a different interpolant (e.g. linear, or Akima/PCHIP variants) and report the change. "Exact" should then be qualified as "exact for the interpolated profile".

(b) Fig. 1 seems to contradict the narrative. The text says double poles "restore the steep limb", but in both insets hybrid5 departs from the PHOENIX nodes at μ ≈ 0.005-0.03 by up to ~0.08 in I/I(1), visibly more than the four-term law. Under r-weighting the region μ < 0.04 (r > 0.9992) carries <0.1% of the weight, which explains both this and the weakly determined w₃. Please say so explicitly. Note also that this region is exactly where spherical 1D models are least reliable and where the M1/M2 choice operates, so the hybrid5 advantage cannot be credited to "following the limb".

(c) Sect. 7.3 states that a linear versus geometric rescaling of μ moves ℓ by ≈0.03, as large as the model-data offset itself. Sect. 4.4 states that limb choices change light curves by only ≲5 ppm. Please reconcile these for the reader (presumably ℓ is far more limb-sensitive than the light curve, but ℓ is not defined). This point deserves more emphasis: it implies that conclusions about which atmosphere "fits JWST" depend on a convention.

**8. Circularity in tuning and testing.** The hybrid2 pole (ε = 0.36) was tuned on 7524 model-passband combinations and then evaluated on grids that presumably include them. The weighting was selected on a "method-study" subset that is also part of the evaluation grid. The structure of the laws was apparently developed on PHOENIX J-band profiles. None of this invalidates the results, but the paper should (i) state which models and passbands were used for each tuning step, and (ii) report the performance on a held-out set: e.g. tune ε on PHOENIX and report the hybrid2 error on MPS-ATLAS (or the reverse), or use a cross-validation split. The analytically placed poles of hybrid4/5 are not affected, which is a point in their favour worth stating.

**9. The observational comparisons test the atmospheres, not the laws, and the solar test is weaker than presented.**

(a) *Sun.* Neckel & Labs (1994) give fifth-order polynomial fits in μ to their CLVs, and to my understanding the points at μ < 0.12 were excluded from those fits. "Fitting hybrid5 to the observed solar limb darkening" therefore means fitting a smooth, odd-and-even polynomial that is an *extrapolation* in exactly the region where the hybrid poles act. The 3-5 ppm agreement measures how well each law represents that polynomial, not the Sun's limb. Please state what was fitted and over what μ range, and remove or qualify "the observed solar limb darkening to 3-5 ppm" in the abstract. The model-versus-Sun differences (PHOENIX 156-307 ppm) depend on the same extrapolated region and on the limb cut, so they also need qualifying. If limb data at μ < 0.1 are needed, more recent solar CLV observations extending close to the limb exist and should be considered.

(b) *Kepler/TESS.* The quantities h′₁ = I(2/3) and h′₂ = I(2/3) − I(1/3) probe the disc interior, where all reasonable laws agree. The authors themselves conclude that the offsets "belong to the models, not to the law". This section mainly reproduces Maxted (2023) and Kostogryz et al. (2024). It is useful as a pipeline check and should be presented as such, perhaps more briefly.

(c) *JWST.* Six stars with a star-to-star scatter (0.031) comparable to the mean offset (0.035) support only a modest statement. The standard error of the mean is about 0.013, so the corrected offset of +0.024 is about 2σ. "Removes about a third" should carry this uncertainty. Please also define ℓ (at which μ, for which law), and give per-star error bars in Fig. 5.

(d) *Recommendation for cool stars.* The recommendation "MPS-ATLAS Set 2 + 100 G for FGK stars" is reasonable but rests on the Kepler calibration that fixed the field strength (see 4d). Why Set 2 rather than Set 1, which Sect. 7.1 finds closest to the Sun in the optical? There is no recommendation for M dwarfs (T_eff < 3500 K), where only PHOENIX is available. Yet PHOENIX is "furthest" from the Sun and shows the largest Kepler/TESS offsets, and M dwarfs are a main JWST target class and four of the five showcase stars of Sect. 6. Please give guidance, and state this limitation explicitly.

**10. Corrections and physicality.** The weight shifts Δw for magnetic and 3D corrections are added linearly. Is the corrected set guaranteed to remain in the physical region? Adding Δw can leave the simplex even if both underlying fits are inside it. Please state how often this happens, and what the interpolation code does in that case. Please also specify the "matching 1D models" for the Stagger and MURaM-hydro corrections, and their T_eff/log g coverage.

**11. External validation of the pipeline.** Apart from μ_cri (agreement to ~10⁻³), the paper does not compare its coefficients directly with published tables. A short comparison of the μ-uniform four-term (and power-2) coefficients with Claret & Southworth (2022, 2023) for a few passbands, in light-curve terms, would give users confidence that the "uniform in μ" rows can be used as a bridge to earlier work. Please also check that Claret's procedure is correctly described as "uniform in μ": if it fits the native nodes, the relevant row of Table 1 is "equal per native node".

### Minor comments

1. **Abstract, Context.** "Polynomials in μ² plus double poles (μ²+ε)⁻² beyond the stellar limb" is opaque to a general reader. Suggest: "laws whose intensity is an even polynomial in μ plus terms of the form (μ²+ε)⁻², with fixed ε > 0".
2. **Abstract.** "Elementary" is used throughout without definition. Please define it (no elliptic or special functions) and distinguish it from "closed form". In Sect. 1, "For a circular planet every standard law has a closed-form light curve" is inaccurate for the power-2 law. The analytic qpower2 (Maxted & Gill 2019) is an approximation, and an exact form in terms of the Appell F₁ function appeared only recently (e.g. Chishtie et al. 2025, arXiv:2506.18860; also Short et al. 2019, RNAAS). Mandel & Agol's results involve elliptic integrals, so they are "closed form" but not "elementary".
3. **Sect. 1.** Please add one or two sentences on why limb-darkening accuracy matters (transmission spectroscopy biases, R_p/R_⋆, stellar density), with references beyond Claret's tables (e.g. Csizmadia et al. 2013; Espinoza & Jordán 2016; Patel & Espinoza 2022).
4. **Sect. 2, Eq. (2) onward.** State explicitly that Π_ε(0) = 1 and Π_ε(1) = 0, and that Π_ε is monotonic on [0, 1].
5. **Sect. 2.** "For hybrid2 this region is exactly a triangle". Please give the vertices or the inequalities.
6. **Sect. 2.** "The simplex … is a sufficient region that contains every model fit in this work". Was this checked on the coefficient spectra and the corrected (Δw) sets as well?
7. **Sect. 2.** The note that "linear interpolation … stays inside the physical region" holds for multilinear interpolation with non-negative weights. Note that it fails for higher-order (e.g. cubic) interpolation, which users may apply.
8. **Sect. 2.** Please say why the polynomial part stops at μ⁴, and why the basis is {1−μ⁴, (1−μ²)²} rather than {1−μ², 1−μ⁴}. They span the same space, so presumably the choice is for the physical region.
9. **Sect. 3.** Define or expand: ACES, COND, SpecInt, ODF (opacity distribution functions), MURaM, SVO, NGTS, IRAC, SOSS, LRS, and "set 1/set 2" (how do abundances and mixing length differ?).
10. **Sect. 3.** The reduction error "below 1.2 × 10⁻⁵ in I/I(1)": what is the corresponding light-curve error?
11. **Sect. 3.** The TESS energy-counter remark is useful. Please say what "shift TESS profiles by ~1%" means (in I/I(1)? in coefficients?).
12. **Sect. 4.1.** Method M2 is mentioned in Sect. 4.4 but never defined. Please define both M1 and M2.
13. **Sect. 4.2.** "Gives the outer 0.5% of the radius (μ < 0.1) a tenth of the weight". Please add "about ten times its share of the projected area (~1%)". As written, "a tenth" reads as an under-weighting.
14. **Sect. 4.2.** "No constraint was active for 400 …; over the full grids they are active in at most 1%". Please give the fraction per law and per grid.
15. **Sect. 4.3.** The quadrature accuracy "10⁻⁸ ppm" (10⁻¹⁴ relative) is at the level of double-precision round-off. Please state what was compared and over how many cases. Also say whether the "exact model light curve" is computed with the same quadrature.
16. **Sect. 4.3.** Are the four impact parameters of the merit function sufficient for grazing geometries (b + p > 1)? Sect. 6 discusses b ≤ 0.85 only.
17. **Sect. 4.4.** "The light curves of hybrid5 are the least sensitive of all laws". Please quantify against the others.
18. **Sect. 5.** "56 of 63 916 above 15 ppm": 63 916 = 1682 × 38, which is the PHOENIX record count per method. Please say so, to help the reader.
19. **Sect. 5 / Table 2.** Why is MPS-ATLAS Set 1 not shown? Add a column, or state that it is similar.
20. **Sect. 6.** The "worst error changes by a median factor of 1.01 … for 0 ≤ b ≤ 0.85" statement is welcome. Please give the corresponding statement for p (you say errors scale ≈ p²; over what range?).
21. **Sect. 6 / Table 3.** "hybrid2 is comparable to the quadratic law". In fact it is better for all M dwarfs and worse for the K dwarf and solar analogue in TESS (39 vs 27; 34 vs 29 ppm). The pattern deserves a sentence. Note also that Fig. 1 does not show hybrid2 or power-2 at all.
22. **Sect. 7.2.** Define SPAM ("synthetic-photometry/atmosphere-model", Howarth 2011) at first use. "Too noisy and degenerate to discriminate further": please quantify (e.g. the typical uncertainty of the Patel & Espinoza coefficients).
23. **Sect. 7.3.** Define ℓ and the "r-method" and "linear rescaling of μ".
24. **Sect. 8.** Define "relative flux error" and specify the "light-curve error" column (p = 0.1, which b?). Give the CDS ReadMe column list in an appendix or table.
25. **Sect. 8.1.** "Fitted to the unnormalised intensity with one stellar edge per model". Please explain why a per-model edge is used here but per-passband choices were tested earlier. Also explain the reason for the 2-10 ppm difference from the passband tables.
26. **Sect. 8.2.** Code and data availability: please give a repository URL and an archived DOI, plus version and licence.
27. **Fig. 1.** The grey band around the profile is not explained. Hybrid2 and power-2 are absent from the legend and the plot. The caption should give log g and [M/H] for both models. The residual axis is in units of 10⁻³ of I/I(1) on a symmetric-log scale; please state the linear threshold. A&A requires every line style to be explained.
28. **Fig. 2.** The MPS-ATLAS panel is Set 2 only; please say so in the caption. The text under the panels duplicates the caption and can be removed. A single-column version would suffice.
29. **Fig. 3.** The PHOENIX w₃ curves are visibly noisy with T_eff. Please comment on whether this noise propagates into interpolated light curves (presumably not, given the degeneracy) and whether users should be warned against interpreting the individual weights physically. The large differences in w₃-w₅ between PHOENIX and MPS-ATLAS reinforce that point.
30. **Fig. 4.** With six laws on a symmetric-log axis the panels are dense. Since Table 3 gives the numbers, consider showing only hybrid5, four-term and quadratic, or move the full figure to an appendix.
31. **Fig. 5.** The legend uses the same dark filled circle for "Kepler" and "MPS-ATLAS 2". The left-panel shading is centred on zero, whereas the text quotes a mean offset. Say what the shading represents (± scatter about zero?). The JWST panel has no error bars. The caption says "hybrid5 predictions", but the point of the figure is the atmosphere models.
32. **Table 1.** Please define "equal per native node", and state p and b for this table explicitly (p = 0.03-0.15 here versus p = 0.1 elsewhere, which explains why its medians exceed those of Table 2).
33. **Table 4.** The number of MPS-ATLAS metallicities (46) and models (16 974) is unusual. A note on the grid steps would help users.
34. **Conclusions.** These are a short repetition of the abstract. Please add the limitations: model systematics dominate, w₃ is poorly constrained, no M-dwarf recommendation, and the elliptical-planet accuracy is untested (if it remains so).
35. **References.** Complete the in-prep reference (authors, arXiv). Check the entries flagged in the source (Hestroffer 1997; Howarth 2011; Kostogryz et al. 2024; Sing et al. 2026). Note the inconsistent naming "Sun analogue"/"solar analogue" in Table 3 and the text.
36. **General.** The manuscript still contains placeholders (author list, institution, acknowledgements). These must be completed before the paper can be processed further.
