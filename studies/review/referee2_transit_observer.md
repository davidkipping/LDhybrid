```
Manuscript: Limb-darkening coefficients for rational laws with elementary transit light curves
Referee: Referee 2 — observational transit / transmission-spectroscopy analyst (TESS, Kepler, HST, JWST light-curve fitting)
Recommendation (A&A grade): C   [= Major revision]
Style and language: (a) well written and concise; however, it is not fully self-contained (see Major comment 1), and the added material requested below should fit within the 12-page limit
Identify referee: no
```

---

## Confidential comments to the Editor

This paper tabulates theoretical coefficients for three new "hybrid" limb-darkening laws. The laws are even polynomials in μ plus fixed double poles beyond the limb. They are presented as giving elementary closed-form transit light curves for planets with elliptical outlines (oblate or rotating planets). The paper also re-tabulates the quadratic, power-2 and four-term laws with a different weighting, uniform in projected radius. It covers PHOENIX, two MPS-ATLAS sets and ATLAS9 over 38 (in fact, 37 are listed) passbands, and it adds 3D and magnetic corrections from MURaM and Stagger. Technically the work looks careful:

- the limb treatment follows Claret & Southworth (2023);
- the photon-weighting and passband conventions are handled with care;
- the fits are constrained to physical profiles;
- I verified the quoted pole positions from the ladder formula;
- the numbers in Table 3 agree with the statements in Sect. 6.

The argument for radius-uniform weighting is sound and useful. It applies to every law, and it may turn out to be the paper's most widely used result.

My concerns, as someone who would have to use these tables in a fitting pipeline:

1. **The whole motivation rests on an unpublished paper.** The paper's reason to exist is elementary light curves for elliptical planets, and this depends entirely on "SquishierPlanet (2026), in prep." That paper is cited for the closed form, the pole ladder, the claim about which laws are elementary, and the q-sampling. None of it can be checked. The authors' own light curves are computed by Gauss–Legendre quadrature accurate to 10⁻⁸ ppm. So the paper itself shows that accurate non-elementary evaluation is available, and it never demonstrates the practical advantage it claims: no oblate-planet light curve, no speed comparison, no numerical-stability test. squishyplanet (Cassese et al. 2024) already models non-spherical planets with arbitrary-order polynomial limb darkening, and it is not cited.

2. **The headline metric is not what observers face.** It is the largest deviation, with fixed tabulated coefficients, from the light curve of the model profile. Observers usually fit limb darkening freely or with priors. In that setting the relevant quantity is the bias in Rp/R*, b, a/R* (and oblateness) at realistic noise and cadence. Sing et al. (2026), which the manuscript cites, conclude from JWST data that the quadratic law is statistically preferred for most FGK stars, with negligible depth bias. The manuscript does not engage with this.

3. **Law errors are tiny next to model errors.** The gains are a few ppm against tens of ppm for the quadratic law. The model-to-model differences are 93–170 ppm, and PHOENIX departs from the Sun by 156–307 ppm, both by the paper's own figures. The paper says this honestly in one sentence, but it never draws the practical conclusion for a user.

4. **There is no usable guidance for priors.** Nothing covers coefficient uncertainties, covariances, the degeneracy of w3 (visibly jagged in Fig. 3), or how to apply the scalar empirical offsets (Δh′₁, Δh′₂, Δℓ) to a five-weight law.

5. **Several internal inconsistencies need fixing:**
   - the passband count (37 listed, 38 claimed);
   - the PHOENIX record counts, which imply MIRI LRS rows beyond the PHOENIX wavelength limit;
   - "3–10×", which Table 1 supports only against per-node weighting (it is about 2–4× against uniform-in-μ);
   - a "tabulated" model systematic that does not appear in the table schema;
   - six of the seven Sing et al. hosts are used, and the omission is not explained.

None of this is fatal. The tables are of lasting value, which is the kind of material A&A welcomes, and the radius-weighting result is publishable on its own. But in its present form the paper does not show a transit observer why or how to switch from quadratic or power-2 to these laws. I recommend major revision. I would want to see the revised version.

The manuscript still carries placeholder author and affiliation fields and an "in prep." self-citation. The Editor may want the companion paper made available to referees, at least as an arXiv preprint, before the revision is assessed.

---

## Report to the authors

### Summary of the paper

The authors tabulate limb-darkening coefficients for three "rational" laws (hybrid2, hybrid4, hybrid5). Each combines terms in even powers of μ with normalised double-pole terms Π_ε at fixed pole positions. The companion paper states that these laws give elementary transit light curves for planets of elliptical outline.

The models are:

- PHOENIX-ACES v3.0 (spherical; the limb is placed at the maximum intensity gradient);
- MPS-ATLAS Sets 1 and 2 (plane-parallel);
- ATLAS9, used for MIRI LRS only.

The coefficients cover 38 passbands. The fits are constrained least squares, uniform in projected radius, with uniform-in-μ fits also tabulated. Quadratic, power-2 and four-term coefficients are computed the same way.

The authors report the following:

- radius-uniform weighting lowers the light-curve error of every law;
- hybrid5 is the most accurate law (median 2.7 ppm for PHOENIX and 0.95 ppm for MPS-ATLAS, for p = 0.1);
- the hybrid tables are fully physical;
- interpolation errors are 5–15 ppm;
- the PHOENIX–MPS-ATLAS model systematic is about 100 ppm.

Against observations, the MPS-ATLAS offsets found by Maxted (2023) for Kepler and TESS are removed by a 100 G MURaM correction, and about a third of the JWST offset found by Sing et al. (2026) is removed. The release consists of CDS tables (24 files), PHOENIX coefficient spectra at R = 3000, 3D and magnetic Δw corrections, and an interpolation code.

### General assessment (Q1–Q7)

**Q1. Why should this paper be published?**
Tabulated coefficients from several model grids, with documented fitting choices, are useful to the community, and A&A has a long tradition of such tables. Two contributions are of value beyond the hybrid laws:

- a quantified case for radius- or area-uniform weighting;
- the explicit model-to-model systematic.

The case for the hybrid laws themselves is less well made. The paper shows they are accurate representations of model profiles. It does not show the following:

- that the closed form gives a practical advantage (the derivation is unpublished, and no oblate application is shown);
- that the accuracy gain matters for real data analysed in the usual ways;
- how an observer should use them as priors or free parameters.

These points are needed for the tables to be adopted rather than merely archived.

**Q2. Are the assumptions spelled out clearly?**
Mostly yes for the atmosphere-side choices: limb placement, photon weighting, passband conventions and constraints. No for the light-curve side. The paper does not state, or does not justify:

- the choice of max |Δ| over p = 0.1 and b ≤ 0.9;
- that coefficients are held fixed;
- that cadence and noise are ignored;
- that the reference "exact" light curve rests on monotone cubic interpolation between μ nodes.

The paper also assumes, without stating it, that the steep spherical-model drop at μ ≲ 0.05 is real. The Sun comparison (Sect. 7.1) suggests otherwise.

**Q3. Are the methods fully described?**
Not fully. Missing are:

- the closed-form light curve;
- the origin of the pole ladder;
- the q-coordinates;
- the exact constraint region used for hybrid4/5 (the exact physical region, or the sufficient simplex?);
- how ℓ, h′₁ and h′₂ are computed from hybrid5 for the comparisons with observations;
- how the tabulated "prior widths" are meant to be applied.

See Major comments 1, 5 and 7.

**Q4. Are the new results adequately emphasised?**
The radius-weighting result is emphasised, though its size is overstated (Major comment 8). The model-systematic result, which matters most for users, appears in one sentence at the end of Sect. 5.

**Q5. Are the figures and tables necessary and properly laid out?**
Mostly. All five figures are full width. Figs. 2 and 3 would work at one column. Several legends are incomplete or ambiguous (Minor comments 20–24).

**Q6. What should be electronic only?**
The coefficient tables and coefficient spectra should go to the CDS, as planned. The 3D and magnetic Δw tables and the empirical-offset table should also be at the CDS, each with a ReadMe that defines exactly how they are applied. The code should be archived with a DOI (e.g. Zenodo) and cited in the paper.

**Q7. IAU designations.**
No problems: WASP-, HAT-P- names and passband names are standard.

---

### Major comments

**1. The motivation and the closed form need to be in this paper, or in a citable preprint, and the practical advantage must be shown.**

The paper rests on the claim that only laws polynomial in μ² plus poles of order ≥ 2 give *elementary* light curves for elliptical planets. It also rests on the pole ladder and the q-sampling. All of these are attributed to "SquishierPlanet (2026), in prep." A referee or reader cannot evaluate work in preparation. At a minimum I ask for:

- (a) a short appendix giving the form of the hybrid light curve, or the key integrals, for a circular and an elliptical occultor, plus the origin of the ε_k ladder (one paragraph would do);
- (b) the companion paper on arXiv before acceptance.

Beyond that, the practical advantage of "elementary" needs to be demonstrated, not asserted.

- The authors compute their own light curves by Gauss–Legendre quadrature to 10⁻⁸ ppm. So an accurate non-elementary route already exists.
- squishyplanet (Cassese et al. 2024, JOSS / arXiv:2409.00167) already models transits of triaxial planets with polynomial limb darkening of nearly arbitrary order. This is directly relevant prior work and should be cited and compared with.
- Please give at least: run-time per light-curve point for the closed form against quadrature against squishyplanet; numerical precision near the contact points and near the innermost pole (ε = 0.0016); and differentiability for gradient-based samplers.
- Please state whether the light-curve code itself is released. Only an interpolation code is mentioned. An observer cannot use hybrid coefficients in batman, exoplanet/starry, juliet or Eureka! without a light-curve implementation. A circular-planet implementation, or at least a mapping onto existing polynomial LD solvers, would make the tables usable tomorrow.

**2. Demonstrate the oblate-planet use case.**

Oblateness signals are typically tens of ppm, concentrated at ingress and egress (Seager & Hui 2002; Barnes & Fortney 2003; Carter & Winn 2010; Zhu et al. 2014; Akinsanmi et al. 2020; recent JWST work by Cassese et al.). They are known to be degenerate with limb darkening. The quadratic law's ~100 ppm ingress/egress residuals for M dwarfs (Fig. 4) are the same size and sit in the same part of the light curve. That makes a compelling case, but the paper never makes it.

I ask for one injection–recovery experiment. Inject an oblate planet using the PHOENIX (or MPS-ATLAS) model profile, then recover it with:

- (i) quadratic or power-2 with fixed coefficients;
- (ii) the same laws with free coefficients;
- (iii) hybrid2/4/5 with fixed coefficients;
- (iv) hybrid2/4/5 with coefficients free or under priors.

Report the bias in oblateness, obliquity, Rp/R* and b at JWST-like precision. This would turn the paper's motivation into a result.

**3. The light-curve-error metric does not represent how coefficients are used, and the paper should report parameter biases under free and prior-constrained fitting.**

The metric is max |Δ| between the fixed-coefficient law light curve and the model-profile light curve, for p = 0.1 and b ∈ {0, 0.3, 0.6, 0.9}. This measures how well a law *represents* a model. It does not measure the error an observer incurs. In practice:

- (a) For JWST and HST spectroscopy, limb darkening is usually fitted freely or with priors. The comparison that matters is then a fit of each law with free coefficients to the model light curve, together with Rp/R*, b, a/R* and t₀. The quantities to report are the residual structure and the parameter biases (cf. Howarth 2011; Espinoza & Jordán 2016; Morello et al. 2017; Maxted 2018). A quadratic law with free coefficients absorbs most of the 30–115 ppm in Table 3. Sing et al. (2026) find empirically that the quadratic law is statistically preferred for most FGK stars and gives minimal depth bias. Please engage with that conclusion directly. If the hybrid laws are no better under free fitting, say so. If they are better, for example by reducing depth bias in M dwarfs, show it.
- (b) A maximum deviation is a poor proxy for χ². Please also report rms deviations over the transit, and the implied depth bias.
- (c) Please put the numbers in context with noise and cadence. Examples: TESS 120 s and 200 s, Kepler 30 min (Kipping 2010 integration smoothing), JWST white-light and R ~ 100 bins. Integration smoothing reduces sharp ingress and egress features, and supersampling requirements for the steep pole terms should be noted.
- (d) The scaling "∝ p²" should be shown, not stated, and extended to large p (0.15–0.3, giant planets around M dwarfs, e.g. WASP-80b, TOI-5205b) and to grazing geometries (b > 1 − p). In grazing geometries only μ ≲ 0.4 is sampled and the steep-limb terms matter most. The current b ≤ 0.9 and p = 0.1 excludes exactly the cases where law choice is most consequential.

**4. Put the law accuracy in proportion to the model systematics, and say what a user should conclude.**

By the paper's own numbers, the systematics dwarf the law errors:

| Quantity | Light-curve error (p = 0.1) | Where in the paper |
|---|---|---|
| hybrid5 law error | ~1–3 ppm | Tables 2–3 |
| Quadratic law error | 30–115 ppm | Table 3 |
| PHOENIX vs MPS-ATLAS (same star) | 93 ppm median, 285 ppm 90th percentile | Sect. 5 |
| PHOENIX vs MPS-ATLAS (same star) | 170 ppm median | Sect. 4.4 |
| PHOENIX vs the Sun | 156–307 ppm | Sect. 7.1 |
| Interpolation over grid gaps | 5–15 ppm | Sect. 5 |
| MPS-ATLAS vs JWST hosts | Δℓ ≈ 0.035 | Sect. 7.3 |

Please also reconcile the two PHOENIX–MPS-ATLAS figures (93 and 170 ppm), or define them more precisely if they are different quantities.

The practical message for a user with fixed coefficients is that grid choice matters about 30–100 times more than law choice. That should be stated in the abstract and conclusions, together with what follows from it: when to fix, when to fit, and what prior widths to use.

Relatedly, consider the M dwarfs, where the quadratic law is worst against PHOENIX (72–115 ppm). Much of that deficit probably comes from μ ≲ 0.05, where the spherical-model limb (located by method M1) drops steeply. The Sun comparison (PHOENIX furthest from observation) and Fig. 1 (residuals of 10⁻¹ in I at μ ≈ 0.01) both suggest this. Please quantify which part of the quadratic error comes from μ < 0.05, and discuss whether that steep drop is observable or a consequence of how the limb is defined. Otherwise the main accuracy advantage over quadratic is an advantage at fitting a feature that may not be real.

**5. Make the tables usable as priors: uncertainties, correlations, degeneracies and offsets.**

As a user, I need to know how to turn a table row into a prior. At present:

- (a) No uncertainty accompanies any coefficient. Please provide a documented recipe, ideally implemented in the code, that propagates (T_eff, log g, [M/H]) uncertainties into a coefficient mean and covariance, as LDTk (Parviainen & Aigrain 2015) does. The recipe should optionally add the model-to-model systematic.
- (b) Fig. 3 shows that the hybrid5 weights, especially w3 and w4 for PHOENIX, are jagged functions of T_eff, with jumps of about 0.1 between neighbouring models. The text concedes that w3 is "weakly determined". So individual weights are not meaningful as Gaussian priors. Please give the correlation structure, for example the principal directions of the light-curve Fisher matrix in w-space for typical p and b. Please also recommend a reduced parametrisation for free fitting. One option: fix the pole weights (w3–w5) at their tabulated values and fit w1 and w2, which are the smooth, polynomial, large-scale terms, with priors. Another option: two "observable" combinations such as h′₁ and h′₂.
- (c) Sect. 8.2 tabulates "the empirical offsets of Sect. 7 with their star-to-star scatter as prior widths". But the offsets are defined in (h′₁, h′₂) for Kepler/TESS and in ℓ for JWST, while the law has five weights. The mapping from a two-dimensional (or one-dimensional) offset to Δw is not unique. Please state exactly how the offsets are applied. Please also show that applying them, for example as a prior in a fit, reproduces the observed h′₁ and h′₂ or ℓ.
- (d) The q-coordinates "provided with the code" are the obvious prior for free fits, but they are not defined in the paper. A uniform prior over the 5-simplex is informative in profile space, because a Dirichlet(1,…,1) distribution concentrates away from the vertices. Please define the mapping and show the prior it implies on h′₁ and h′₂ (or on I(0)), so users know what an "uninformative" hybrid5 fit assumes.
- (e) Linearity in w is a real practical advantage, and the paper understates it. For fixed geometry the light curve is linear in the weights, so with Gaussian priors the limb darkening can be marginalised analytically, as for the polynomial laws in starry. This is worth a paragraph. The positivity and monotonicity constraints break pure Gaussianity, and how to handle that should be discussed.
- (f) Sect. 5 says the model systematic is something "we tabulate", but the record description in Sect. 8 has no such column. Please add it, or say where it lives.

**6. Coverage is thinnest where the hybrid laws are claimed to matter most.**

- (a) **M dwarfs.** The largest gains are for M dwarfs (Table 3), but these are covered only by PHOENIX v3.0 below 3500 K. That is the grid furthest from the Sun (Sect. 7.1) and has no empirical validation in the paper: the Maxted (2023) and Sing et al. (2026) samples are FGK. No recommendation is given for M dwarfs. Please say what an M-dwarf observer should use, and what prior widths are justified. Kostogryz et al. (2026) state that the magnetic effect is negligible in M dwarfs, which should be mentioned.
- (b) **Coefficient spectra.** These are given only for PHOENIX, but the recommended FGK grid is MPS-ATLAS Set 2 + 100 G. A JWST user who wants R ~ 100 spectroscopic-bin coefficients from the recommended grid therefore has no route. Please either provide MPS-ATLAS coefficient spectra (on its native ODF sampling, with a warning about how coarse that sampling is for narrow bins) or explain why not.
- (c) **Narrow bins.** In narrow bins across strong molecular bands (H₂O, CO, CH₄ in M dwarfs; Na and K cores), limb brightening or non-monotonic profiles can occur. Please state whether the per-bin coefficients are constrained, what happens when the summed band coefficients fall outside the physical region, and how often this occurs.
- (d) **MIRI LRS.** This is covered only by ATLAS9 at 3500–6500 K. MPS-ATLAS intensities are computed on the standard ATLAS ODF grid, which extends to the far-IR; the ExoTiC-LD mps1/mps2 bundles use the same 1221-point wavelength grid as its ATLAS9 bundle. So why is MPS-ATLAS not used for MIRI? More importantly, MIRI LRS transit and eclipse work is dominated by M-dwarf hosts below 3500 K, which are not covered at all. Limb darkening is weak at 5–12 μm, so a remark that a simple law suffices there would itself be useful guidance.
- (e) **Passbands.** Please consider adding PLATO, Ariel (AIRS channels and FGS), Pandora (VIS and NIR) and NIRSpec G395H split into NRS1/NRS2, as JWST observers use them. Note that whole-band G395H coefficients are rarely what is needed for spectroscopy.

**7. The observational comparisons do not test the hybrid laws' distinguishing feature, and the JWST analysis needs clarification.**

- (a) Maxted's h′₁ = I(2/3) and h′₂ = I(2/3) − I(1/3) probe the disc interior. Agreement between hybrid5 and four-term offsets ("they belong to the models, not to the law") is therefore expected by construction. No comparison probes μ ≲ 0.1, where the laws differ. The Neckel & Labs (1994) profiles are themselves fifth-order polynomials in μ, so a fit to them tests flexibility on a smooth profile, not the steep limb. Please acknowledge this limitation. Grazing or high-b transits would be the natural observational test.
- (b) Sing et al. (2026) define ℓ = 1 − (c₁ + c₂) from *quadratic* fits restricted to μ ≥ 0.2, and they analyse **seven** hosts (WASP-121, WASP-12, WASP-94A, WASP-6, HAT-P-11, WASP-107, WASP-80). The manuscript uses six "FGK hosts" and calls ℓ a "limb intensity". Please:
  - state which star was dropped and why (WASP-12, presumably);
  - state exactly how ℓ was computed from your models (quadratic fit to the hybrid5 profile over μ ≥ 0.2, or hybrid5 I(0)?). These are very different quantities;
  - check that your rescaling argument for the PHOENIX discrepancy (+0.041 against their 0.000 ± 0.041) also accounts for the μ_min = 0.2 truncation.
- (c) Statistical significance:
  - the reduction from +0.035 ± 0.031 to +0.024 (scatter, n = 6) is a ~0.01 shift against a standard error of the mean of ~0.013;
  - the residual offset after the 100 G correction is ~2σ from zero;
  - "removes about a third" therefore overstates what six stars show, with WASP-6 as an obvious outlier;
  - please quote the standard error and use a robust statistic.
- (d) The 100 G level was in effect calibrated on these same Kepler/TESS data (Maxted 2023; Kostogryz et al. 2024). So the left panel of Fig. 5 is a consistency check, not an independent validation. Please also credit Maxted (2023), who already attributed the Δh′₁ ≈ 0.006 offset to ~100 G fields.
- (e) Stellar activity (unocculted and occulted spots and faculae) alters the effective limb darkening measured in transits and is relevant to several of these hosts (e.g. WASP-6, HAT-P-11, WASP-107) and to M dwarfs in general. One or two sentences on how the tabulated priors should be widened for active stars would help users.

**8. Correct the size of the weighting improvement.**

The abstract, Sect. 4.2 and the conclusions state that radius-uniform weighting reduces the light-curve error "by factors of 3–10 for every law". Table 1 shows that this range holds only against "equal per native node". Against uniform-in-μ, which is the relevant comparison with Claret's tables, the median improvement is:

| Law | Median error, μ-uniform → r-uniform | Factor |
|---|---|---|
| hybrid5 | 13 → 4.9 ppm | 2.7× |
| four-term | 55 → 13 ppm | 4.2× |
| power-2 | 57 → 22 ppm | 2.6× |

For the worst cases the factors are 3.0×, 3.7× and 2.1×. Please quote the improvement relative to uniform-in-μ, about 2–4×.

Please also state how the r-weighted quadratic, power-2 and four-term coefficients differ, in coefficient units, from Claret's and from ExoTiC-LD/LDTk outputs for the same models. Users switching tables need to know the size of the shift. Claret, Hauschildt & Torres (2025) provide JWST four-term and power-2 coefficients from the same PHOENIX models, so a direct comparison is possible and would validate the pipeline.

**9. Treat power-2 on the same footing as the other laws.**

Power-2 is fitted without constraints, and the paper then reports that it is "unphysical (c > 1) in 8.6% of PHOENIX rows". Every other law is fitted within its physical region. This penalises power-2 for a choice the authors made, and it inflates the "every hybrid coefficient set is physical" contrast. Please fit power-2 with c ≤ 1 and α > 0 (the region is simple), or else state clearly that the unphysical rows are an artefact of the fitting choice. Also note that "every hybrid set is physical" holds by construction, because the fits are constrained. It is not an empirical result and should not be listed among the Results in the abstract.

**10. Internal consistency of the grids, passbands and tables.**

- (a) Sect. 3 lists 37 passbands, but the text, the abstract and Table 4 say 38. Please list all 38.
- (b) Table 4 gives 127,832 PHOENIX records. That equals 1682 models × 38 passbands × 2 methods, and Sect. 5 uses 63,916 = 1682 × 38 PHOENIX rows. This implies the PHOENIX table contains MIRI LRS, which lies beyond the PHOENIX wavelength coverage (5.5 μm), as Sect. 3 itself states. Please clarify which passbands each grid covers.
- (c) Sect. 3 restricts the grids to T_eff ≤ 7400 K, but Table 4 gives MPS-ATLAS up to 7500 K, and the MPS-ATLAS curves in Fig. 2 appear to extend to about 7600 K.
- (d) hybrid4 is markedly worse than four-term on MPS-ATLAS (median 6.8 against 1.4 ppm). That is the recommended grid. Please explain this, for example the K = 2 pole ladder being suited to the steep spherical limb but not to the plane-parallel one. Please also give an explicit law recommendation by grid: "for MPS-ATLAS, use hybrid5 or four-term, not hybrid4".
- (e) The hybrid2 pole (ε = 0.36) was tuned on 7524 model–passband combinations and then evaluated on the same grids. Please identify the tuning set, and show that hybrid2's performance is not sensitive to ε near 0.36, for example the error at ε = 0.3 and at 0.42.

---

### Minor comments

1. **Abstract.** "for a central transit of M to G dwarfs … at most 3 ppm" rests on five PHOENIX stars at solar metallicity (Table 3). Please qualify it, or support it with grid-wide statistics restricted to b = 0.
2. **Abstract / Sect. 5.** "against JWST it removes about a third". See Major comment 7c.
3. **Sect. 1.** "the four-term law represents model intensities best" is not quite true for spherical models without a limb cut. Please qualify.
4. **Sect. 1.** Please cite the oblateness literature (Seager & Hui 2002; Barnes & Fortney 2003; Carter & Winn 2010; Zhu et al. 2014; Akinsanmi et al. 2020), squishyplanet (Cassese et al. 2024), and the main empirical and LD-bias works (Sing 2010; Csizmadia et al. 2013; Parviainen & Aigrain 2015; Espinoza & Jordán 2016; Morello et al. 2017). For power-2 sampling, cite Short et al. (2019) and Maxted (2018) for the h₁/h₂ parametrisation.
5. **Sect. 2, Eqs. (3)–(5).** Please state explicitly that every T_j is monotonically decreasing in μ, which is why w ≥ 0 implies monotonicity. Please show the hybrid2 "exact triangle" constraints, i.e. the two vertices that are not at the origin.
6. **Sect. 2.** "the simplex … is a sufficient region that contains every model fit in this work". Was the simplex itself the constraint used in the QP, or was it the exact physical region? This matters for users setting priors (Major comment 5d).
7. **Sect. 2.** The claim that "a least-squares fit is a fixed linear map" holds only for unconstrained fits. Please note that it fails when a constraint is active. This is relevant to Sect. 8.1.
8. **Sect. 3.** "would otherwise shift TESS profiles by ~1%". One per cent of what: I(μ) at a given μ, or a coefficient? Please be specific.
9. **Sect. 3.** Please state which wavelength ranges of the NIRSpec and NIRISS modes are used: whole grating range or detector-limited, NRS1 + NRS2 combined, SOSS order 2 truncation.
10. **Sect. 4.2.** "gives the outer 0.5% of the radius (μ < 0.1) a tenth of the weight" reads as under-weighting. In fact this region holds about 1% of the area but receives 10% of the weight under μ-uniform sampling. Please rephrase.
11. **Sect. 4.2 / Table 1.** Area-uniform weighting is slightly better than r-uniform for hybrid5 (3.9 against 4.9 ppm). Why tabulate r-uniform? Please also say whether the μ-uniform fits use Claret's 100-point grid, or the same 2000-point resample.
12. **Sect. 4.2.** "Constraining … costs nothing: … no constraint was active for 400 model–passband combinations". Please give the grid-wide fraction per law and grid. Sect. 5 mentions blue passbands at 2300–2500 K.
13. **Sect. 4.3.** The reference light curve depends on monotone cubic interpolation between μ nodes. For MPS-ATLAS (24 nodes) the choice of interpolant could itself account for a large part of the 0.95 ppm median. Please quantify the reference's own uncertainty, for example cubic against linear interpolation in μ and in r.
14. **Sect. 4.4.** Please give the corresponding numbers for the other laws, and define "90th percentile ≲ 8 ppm" by law.
15. **Sect. 5.** The interpolation errors (5, 15 and 13 ppm) exceed the hybrid5 law error by a factor of 2–5. Consider offering an option to interpolate the intensity profiles and refit, which is cheap for a linear QP, rather than interpolating coefficients. Also say whether the leave-one-out tests were done at grid edges.
16. **Sect. 6.** "Errors scale approximately as p²". See Major comment 3d. Please also report b > 1 − p.
17. **Sect. 7.2.** The Patel & Espinoza (2022) comparison uses SPAM coefficients, which are not what is tabulated. As written it validates the PHOENIX intensities, not the tables. Please say so.
18. **Sect. 8.** Please give the total size of the CDS release. Twenty-four files, two of them with 1.29 × 10⁶ records each, is large. Consider also releasing the band-integrated I(μ) profiles on the resampled r grid, so users can refit with any law or weighting.
19. **Sect. 8.2.** Please specify the sign convention and the reference 1D model for each Δw. Please also say whether Δw are additive in w (the hybrid laws are linear, so presumably yes) and whether adding them can leave the physical region.
20. **Fig. 1.** The thick grey band is not explained. hybrid2 and power-2 are omitted without comment. The panel shows hybrid5 residuals of about −10⁻¹ in I at μ ≈ 0.01–0.02, larger than four-term there. This sits uneasily with Sect. 1's "restores the steep limb". Please comment.
21. **Fig. 3.** Five shades of blue are hard to tell apart. The jaggedness of w3 and w4 for PHOENIX deserves discussion in the text, not only the remark in the caption (Major comment 5b). The non-monotonic T_eff dependence (minimum near 3500 K, maximum near 4300 K in w1 and w2) should be explained physically (onset of molecular opacity?).
22. **Fig. 4.** Please give exposure-time and cadence assumptions. Consider adding a curve, or a column in Table 3, for the quadratic law with *free* coefficients (Major comment 3a).
23. **Fig. 5.** The legend uses the same dark filled circle for "Kepler" and for "MPS-ATLAS 2", and "PHOENIX" grey appears only in the right panel. Error bars are not defined (1σ of the mean?). The right panel omits WASP-12 (Major comment 7b). The caption should say how ℓ is computed.
24. **Figs. 2–3.** These would read well at one-column width, per A&A guidance.
25. **Table 2.** Is "MPS-ATLAS Set 2" representative of Set 1? Please add Set 1 or a remark.
26. **Table 4.** Please state the [M/H] values for MPS-ATLAS (46 values between −2 and +1 implies non-uniform spacing) and the log g values.
27. **Sect. 9.** "hybrid5 is the most accurate law tested". Please add "against 1D model profiles". Against the Sun, four-term is within 1 ppm of it.
28. **References.** "SquishierPlanet 2026, in prep." is not acceptable as the basis of the method (Major comment 1). The author and affiliation placeholders need completing.
29. **Data and code availability.** Please give the repository URL and archive DOI for the interpolation code, and for the light-curve code if released. Please state the licence and the dependencies' versions.

---

*Sources consulted to check claims and citations:* Sing et al. 2026, arXiv:2609.00263 (sample of seven hosts; definition ℓ = 1 − (c₁ + c₂) from quadratic fits at μ_min = 0.2; MPS-ATLAS Set 2 Δℓ = +0.036 ± 0.031; PHOENIX 0.000 ± 0.041); Kostogryz et al. 2026, arXiv:2606.21912 (magnetic effect negligible in M dwarfs); Kostogryz et al. 2024, Nat. Astron. 8, 929; Maxted 2023, MNRAS 519, 3723 (43 FGK hosts; offset Δh′₁ ≈ 0.006 attributed to ~100 G); Claret, Hauschildt & Torres 2025 (JWST four-term and power-2 coefficients from PHOENIX); Cassese et al. 2024, squishyplanet (arXiv:2409.00167, JOSS 10.21105/joss.06972); ExoTiC-LD stellar-grid documentation (mps1/mps2 and kurucz share a 1221-point wavelength grid).
