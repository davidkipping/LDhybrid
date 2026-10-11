```
Manuscript: Limb-darkening coefficients for rational laws with elementary transit light curves
Referee: Referee 1 (theoretical limb-darkening specialist; model-atmosphere LDC tabulations)
Recommendation (A&A grade): C   [= Major revision]
Style and language: (a) well written and concise; if anything over-condensed. No reduction is
                    needed. Sections 3 (Method) and 6 (Observations) need expanding, which fits
                    easily within the 12-page limit (the paper is currently about 8 pages).
Identify referee: no
```

---

## Confidential comments to the Editor

The paper provides tables of limb-darkening coefficients (LDCs) for three new "hybrid" laws: even polynomials in μ² plus fixed double poles outside the limb. It covers PHOENIX-ACES, MPS-ATLAS (Sets 1 and 2) and ATLAS9 models in 38 passbands used for exoplanets. For comparison it also provides quadratic, power-2 and four-term coefficients, together with coefficient spectra and 3D/magnetic corrections. The work is competent and carefully done. The tables would be useful if the laws are taken up, and the CDS product described is broadly what A&A welcomes. The light-curve error is a sensible metric to tabulate. The authors' figures are mostly consistent with their tables: I checked the claims in Sect. 5 entry by entry against Table 3 and they hold.

I nevertheless recommend major revision, for these reasons.

1. **Unverifiable foundation.** The reason the laws exist is that they give elementary light curves for elliptical planet outlines. That claim, the derivation of the pole positions, and the closed-form solution used to validate the quadrature all rest on a single reference marked "in prep." (SquishierPlanet 2026). A referee cannot evaluate it. This needs to be fixed before acceptance: either the companion paper is available, or the essentials go into an appendix.
2. **Headline claims are overstated.**
   - **Weighting.** Fitting uniformly in r is presented as a new finding that cuts light-curve errors by "factors of 3–10". It is the weighted-r quasi-spherical fit of Morello et al. (2017, 2020), and the authors' own Table 1 supports only a factor of about 2–4 relative to uniform-μ fitting.
   - **hybrid5 vs four-term.** The claim that hybrid5 beats the four-term law compares a five-parameter law with a four-parameter one. The like-for-like four-parameter hybrid is five times worse than the four-term law on MPS-ATLAS.
   - **"At most 3 ppm".** This abstract claim comes from 15 hand-picked cases, while the grid maximum is 50 ppm.
3. **Accuracy gains are within the systematics.** The gains are at the 2–5 ppm level, which equals the authors' own limb-placement systematic and is smaller than their interpolation error. The manuscript does not address the recent empirical result (Sing et al. 2026, which it cites) that quadratic laws are statistically preferred at JWST precision.
4. **Methods and observations need work.** Several methodological points need clarification or extra tests:
   - whether the stellar edge is placed per passband or per model, which implies a chromatic stellar radius;
   - no demonstration that the pipeline reproduces existing published coefficients;
   - no comparison with Claret, Hauschildt & Torres (2025) for the overlapping JWST passbands.

   The observational section largely repeats Maxted (2023) and Kostogryz et al. (2024). It uses 6 of the 7 JWST targets of Sing et al. without explanation, and it criticises Sing et al.'s rescaling without showing the test.

None of these is fatal. With the companion paper available, honest framing, and the extra validation, I expect the paper to be publishable in A&A. I would be glad to see the revised version.

---

## Report to the authors

### Summary of the paper

The authors define three limb-darkening laws that are linear in their weights: hybrid2 (μ² term plus one pole), hybrid4 (two even-polynomial terms plus two poles) and hybrid5 (two even-polynomial terms plus three poles). The poles Π_ε are placed at fixed values of ε, the hybrid4 and hybrid5 positions from an "analytic ladder". The authors state that these laws give closed-form, elementary transit light curves for planets of elliptical outline.

They fit these laws, and the quadratic, power-2 and four-term laws, to band-integrated intensities of:
- spherical PHOENIX-ACES v3.0 models, with the edge placed at the maximum of |dI/dr| (method M1);
- plane-parallel MPS-ATLAS Sets 1 and 2;
- ATLAS9, for MIRI LRS only.

Fits use least squares uniform in projected radius r and are constrained to positive, monotonic profiles. Fit quality is measured mainly by the largest deviation of a p = 0.1 transit light curve from that of the model profile.

The main results:
- r-weighting reduces this error by "factors of 3–10".
- hybrid5 reaches medians of 2.7 ppm (PHOENIX) and 0.95 ppm (MPS-ATLAS), better than the four-term law.
- hybrid4 is comparable to the four-term law on PHOENIX; hybrid2 is comparable to the quadratic law.

The authors compare model predictions with the solar CLV of Neckel & Labs, with Kepler/TESS measurements (Maxted 2023) and with JWST measurements (Sing et al. 2026). They recommend MPS-ATLAS Set 2 with a 100 G MURaM correction for FGK stars. The products for the CDS are: 24 tables, coefficient spectra at R = 3000 for PHOENIX, Δw corrections and an interpolation code.

### General assessment

**Q1. Why should this paper be published?**

If the hybrid laws are adopted for modelling non-circular planets, tables of their coefficients over standard grids and passbands will be needed, and this paper supplies them in a careful and well-documented way. Several features are welcome and above the usual standard:
- every record carries a light-curve error;
- unphysical coefficients are explicitly flagged;
- the gap between models (PHOENIX vs MPS-ATLAS) is quantified;
- the coefficient spectra exploit the linearity of the laws.

However, the paper's value depends entirely on the laws being useful. That case is made only by reference to an unpublished paper, and is not weighed against existing numerical tools or against recent empirical evidence that higher-order laws are not needed at present precision (Major comments 1 and 4). The validity of the tables themselves is broadly convincing. The comparative claims are over-stated (Major comments 2 and 3).

**Q2. Are the assumptions spelled out clearly?**

Partly. Several key choices are stated without justification or are ambiguous:
- whether the edge is set per model or per passband;
- whether hybrid4/5 are constrained to the simplex or to the exact physical region;
- why ATLAS9 is used when MPS-ATLAS covers MIRI;
- how the Bessell curves are integrated;
- the microturbulence and convection parameters of each grid;
- the reference 1D models for the 3D/MURaM corrections.

**Q3. Are the methods fully described?**

No (Major comment 7). Several elements are not described at all: method M2, the constraint implementation, the Levenberg–Marquardt setup for power-2, the accuracy of the reference light curves, and the construction of the 3D/magnetic corrections. The paper is too terse to be reproduced from the text alone.

**Q4. Are the new results adequately emphasised?**

The new results (coefficients for the hybrid laws and the light-curve error of every fit) are emphasised. The emphasis is distorted, though, by presenting r-weighting as a new finding and by over-claiming in the abstract (Major comments 2 and 3).

**Q5. Are all figures and tables necessary and properly laid out?**

Mostly, yes. All five figures are full-page width; Figs. 1 and 3 could be one column. Several legend and line-style problems remain (Minor comments 21–25). Table 3 could be shortened, or Fig. 4 moved to an appendix, because they carry the same information.

**Q6. What should be published in electronic form only?**

The coefficient tables, coefficient spectra, Δw corrections and empirical offsets belong at the CDS, each with a proper ReadMe. Table 4 currently lists only the 24 coefficient files. The code needs a persistent archive (Zenodo or ASCL with a DOI). Fig. 3 and possibly Fig. 4 could go to an appendix.

**Q7. IAU designations.**

No problems: WASP-, HAT-P- and the Sun follow standard usage. The Kepler/TESS hosts are not named, which is fine if the source catalogue is cited.

---

### Major comments

**1. The foundation of the paper is an unpublished companion paper.**

Every central premise is attributed to "SquishierPlanet (2026), in prep.":
- that the hybrid laws give elementary light curves for elliptical outlines;
- that *only* laws built from even polynomials plus poles of order ≥ 2 do so (a necessity claim, Sect. 1);
- that odd powers lead to elliptic integrals and fractional powers to "no closed form at all";
- the derivation of the pole ladder ε_k = exp[2π(√(k+1/2) − √K)];
- the earlier tuning of the hybrid2 pole (0.208);
- the "closed-form hybrid solution" used to validate the quadrature to 10⁻⁸ ppm (Sect. 3.3).

A&A policy is clear that work "in preparation" cannot support results. I checked the pole values and they follow from the formula given. However, I cannot judge whether the formula is the right one, whether the light curves are indeed elementary, or how stable they are numerically.

I ask the authors to do one of the following:
- (a) make the companion paper available (at least on arXiv) and cite it properly; or
- (b) add a short appendix that states the closed-form light-curve result for one Π_ε term and one even-polynomial term for an elliptical occultor. It should also give the origin of the ε ladder and how the stated accuracy of the reference solution was checked.

The practical motivation also needs to be argued, not assumed:
- Elliptical-outline transits can already be computed to high accuracy with arbitrary polynomial limb darkening by numerical one-dimensional boundary integrals. Examples are Cassese et al. (2024, *squishyplanet*), the earlier oblate-planet work of Seager & Hui (2002), Barnes & Fortney (2003) and Carter & Winn (2010), and general eclipse codes.
- What is gained in speed, differentiability or accuracy by an elementary solution? A short quantitative comparison (runtime and accuracy against such a code with the four-term or power-2 law) would make the case.
- The authors should also state the amplitude of the oblateness signals the laws are meant to model. The required accuracy of the LD law (≲ 1–3 ppm vs ~10 ppm) only matters if the signal is of that order.

The statement in Sect. 1 that "for a circular planet every standard law has a closed-form light curve" needs qualifying. The nonlinear (four-term) solution of Mandel & Agol (2002) is a series in hypergeometric-type functions. The power-2 law is handled by approximations or series (e.g. Maxted & Gill 2019; Short et al. 2019). Neither is "closed form" in the sense used here for the hybrid laws.

**2. The r-weighting result is not new, and Table 1 does not support "factors of 3–10".**

(a) *Priority.* Fitting the intensity profile with weights proportional to the sampling interval in r, i.e. uniform in r, is exactly the "weighted-r" (quasi-spherical) fit of Morello et al. (2020, ExoTETHyS, Sect. 2.1.4). It in turn generalises the quasi-spherical fits of Claret et al. (2012). Morello et al. (2017, AJ 154, 111) had already shown that sampling uniformly in r rather than μ is what reconciles model coefficients with empirical ones, and that the choice affects light-curve residuals at the tens-of-ppm level. The manuscript says only that its finding is "consistent with" Morello et al. (2020). It should say plainly that its adopted method is that one. The abstract and conclusions should then present the factor as confirmation and quantification, not as a result of this paper.

(b) *The number.* Table 1 shows that moving from uniform in μ to uniform in r changes the median/worst light-curve error as follows:

| Law | Median | Worst |
|---|---|---|
| hybrid5 | 13 → 4.9 (×2.7) | 33 → 11 (×3.0) |
| Four-term | 55 → 13 (×4.2) | 123 → 33 (×3.7) |
| Power-2 | 57 → 22 (×2.6) | 152 → 72 (×2.1) |

Factors of 3–10 appear only against "equal weight per native node". That is not "Claret's method": Claret & Southworth (2023), for example, resample to 900 equally spaced μ points before a Levenberg–Marquardt fit, so their weighting is the "uniform in μ" row. The abstract, Sect. 3.2 and the Conclusions must quote the factor actually supported (about 2–4 against uniform-μ fitting), or show the full-grid statistics that support 3–10.

(c) *Circularity of the metric.* The figure of merit is the light-curve deviation of a transit, which samples the disc by area. A fit weighted in r (or area) is by construction closer to optimising that figure of merit. That r-weighting wins on it is therefore expected, not a demonstration that the established procedure is "suboptimal".

The cost elsewhere is visible in the authors' own Fig. 1. With r-weighting, hybrid5 misses the PHOENIX limb nodes by up to ~0.07–0.1 in I/I(1) for μ < 0.03, where the four-term fit follows them. The tables are used well beyond transits: eclipsing binaries (especially partial and grazing eclipses), interferometric visibilities, microlensing caustic crossings, and Rossiter–McLaughlin modelling. For those uses the limb intensities matter and the μ-weighted fits are the appropriate product. Please:
- report Claret's merit function σ (uniform in μ) for the r-weighted fits, alongside the light-curve error;
- state explicitly in Sect. 7 and in the CDS ReadMe that the r-weighted coefficients are optimised for transit light curves.

(d) If the goal is the light curve, the logical end point is to fit the coefficients to synthetic light curves directly (the SPAM approach of Howarth 2011), which the authors use only for the Patel & Espinoza comparison. Please show how close r-weighting gets to the SPAM optimum for at least the Table 3 cases. This would also show whether the ~2–3 ppm floor of hybrid5/hybrid4 on PHOENIX is a property of the law or of the weighting.

(e) Table 1 shows uniform-in-area fitting is better than uniform-in-r for every law (3.9 vs 4.9, 12 vs 13, 19 vs 22 ppm median). The text calls them "equivalent within the noise", but these are deterministic computations with no noise. Why was r adopted rather than area?

**3. The comparison between laws is not on an equal footing, and the abstract over-claims.**

(a) *Number of parameters.* hybrid5 has five free weights, the four-term law four. The fair comparison with the four-term law is hybrid4. By the authors' Table 2:
- On MPS-ATLAS, hybrid4 is worse than the four-term law by a factor of ~5 in the median (6.8 vs 1.4 ppm) and ~2 at the 99th percentile.
- On PHOENIX, Fig. 2 shows hybrid4 is worse than the four-term law below ~4000 K, i.e. for all M dwarfs.

The abstract states that hybrid5 is "better than the four-term law" but does not mention that the like-for-like hybrid is worse. On MPS-ATLAS the hybrid5 advantage (0.95 vs 1.4 ppm median) is negligible in practice. Please rewrite the abstract and Sect. 4 to make the parameter count explicit, and consider an information-criterion-style or per-parameter comparison.

(b) *Constraints applied unequally.*
- Power-2 is fitted *without* constraints, then described as "unphysical (c > 1) in 8.6% of PHOENIX rows", while all other laws are constrained to their physical regions. For c > 1 the intensity becomes negative only for μ < (1 − 1/c)^{1/α}, an annulus that may be vanishingly thin. Please either fit power-2 with c ≤ 1 (and α > 0) and report the cost in light-curve error, or quantify the width and flux of the negative-intensity annulus.
- Conversely, the four-term law is constrained to its "exact positive–monotonic region", whereas published four-term tables are unconstrained. Please report how often the four-term constraint is active and what it costs. The four-term figures in Tables 1–3 may be penalised by a constraint that users of existing tables do not apply.

(c) *Missing null model.* The hybrid laws' raison d'être is that even polynomials in μ² are elementary but cannot follow the steep limb, and that the poles fix this. The natural baseline is therefore an even polynomial of the same order (e.g. 1 − a(1−μ²) − b(1−μ²)², and the analogous 4- and 5-term versions), which are also elementary for an ellipse. Please add these to Tables 2 and 3. Without them the reader cannot tell how much the poles contribute.

(d) *hybrid2.* hybrid2 is no more accurate than the quadratic law: worse for the K dwarf and solar analogue in TESS in Table 3, and similar in median in Table 2. Its single pole was tuned on the same 7524 model–passband combinations on which it is then evaluated. Please give a cross-validated estimate (tune on half the grid, evaluate on the other half) or acknowledge that the hybrid2 numbers are in-sample.

(e) *"At most 3 ppm".* The abstract claims hybrid5 light curves "deviate from those of the model profiles by at most 3 ppm" for M to G dwarfs. This comes from 15 hand-picked cases (Table 3) at solar metallicity on grid nodes. Over the full PHOENIX grid, Table 2 gives a 99th percentile of 6.3 ppm and a maximum of 50 ppm. Please state the result in a way that reflects the full grid.

(f) "Every hybrid coefficient set is physical" (abstract) holds by construction, because the fits are constrained. It is a property of the method, not a result, and should be phrased as such.

**4. The accuracy gains must be set against the systematics and against the empirical evidence.**

The manuscript's own numbers show that the differences between the best laws lie inside its error budget:
- Changing the limb treatment (M1 vs M2, edge per model vs per passband, coarser μ grids) moves light curves by a median of ≲ 5 ppm, 90th percentile ≲ 8 ppm (Sect. 3.4).
- Linear interpolation (leave-one-out) costs 5, 15 and 13 ppm in Teff, log g and [M/H] (Sect. 4).
- PHOENIX and MPS-ATLAS differ by a median of 170 ppm (Sect. 3.4), or 93 ppm (Sect. 4). Please reconcile these two numbers.

Against these, hybrid5 vs hybrid4 vs four-term on PHOENIX is 2.7 vs 3.8 vs 8.0 ppm median. In Table 3, hybrid5 is in fact slightly *worse* than hybrid4 for four of the K/G cases. The ranking of the three best laws is therefore not robust to the authors' own systematics, and after interpolation it is likely to vanish for any star off a grid node.

Please:
- (i) state this clearly;
- (ii) note that a leave-one-out test with doubled spacing overestimates linear-interpolation error by roughly a factor of 4, and give the expected error at the true grid spacing;
- (iii) consider recommending interpolation of the band-integrated *intensities* followed by refitting where ppm accuracy matters.

The maximum light-curve deviation at p = 0.1 is a reasonable summary, but it is not what an observer cares about. What matters is the bias it induces in fitted parameters. Please add, for the Table 3 cases at least, the bias in Rp/R*, b and a/R* when the model light curve is fitted with each law using the tabulated coefficients (in the spirit of Espinoza & Jordán 2015). For JWST work, the chromatic bias in (Rp/R*)² across a spectrum is the key quantity.

The geometry should also be extended:
- The innermost hybrid5 pole acts only at μ ≲ 0.04 (r ≳ 0.9992), and its weight w₃ is admitted to be weakly determined (Fig. 3 shows jumps of ~0.1 between adjacent 100 K PHOENIX models).
- The test set stops at b = 0.9 for p = 0.1, i.e. it excludes grazing geometries (b ≳ 1 − p), where the light curve is set entirely by the limb.
- The claim that errors scale as p² should be shown, and tested for p = 0.2–0.3, which is relevant for giant planets around M dwarfs.

Finally, Sing et al. (2026), cited in Sect. 6.3, conclude from JWST white-light curves of seven FGK hosts that the quadratic law is statistically preferred over higher-order laws. They also find that the extra parameters are not constrained by the data and that limb structure at μ ≲ 0.2 does not measurably affect transit depths at JWST precision. This bears directly on whether a five-parameter law tuned to the extreme limb of 1D models is useful. It should be discussed openly rather than citing that paper only for an offset.

**5. Spherical models: definition and consistency of the limb, and validation against existing tables.**

(a) *Model identity.* "PHOENIX-ACES-AGSS-COND SpecInt v3.0" is not the PHOENIX-COND grid used in Claret et al. (2012, 2013), Claret (2017, 2018) and Claret & Southworth (2023). I checked the Göttingen server: the v3.0 README states that these are the SpecInt intensities of Husser et al. (2013) extended to 5.5 μm. Please say so, and state the microturbulence and other fixed parameters.

(b) *Agreement of μ_cri.* The claim that the authors' μ_cri "agree with those of Claret & Southworth (2023) to ~10⁻³" needs specifying, because Claret & Southworth (2023) used PHOENIX-COND at solar metallicity with Vξ = 2 km s⁻¹, resampled to 900 μ points. Which models, passbands and quantities were compared? Method M2 is invoked in Sect. 3.4 but never defined; please define it (the midpoint between μ_cri and the adjacent node, in Claret & Southworth's notation).

(c) *Pipeline validation.* The paper never shows that the pipeline reproduces existing tabulated coefficients when the same choices are made. Please show, for a few representative models, that with uniform-μ weighting, the same integration (photon or energy) and the same edge, the quadratic, power-2 and four-term coefficients agree with published values. Suitable comparisons are Claret & Southworth (2022, 2023) for Kepler/TESS, and Claret, Hauschildt & Torres (2025, A&A 699, A97) for the JWST NIRCam/NIRISS/NIRSpec passbands. The last of these is directly overlapping work: it tabulates four-term and power-2 coefficients for 11 JWST passbands from high-resolution spherical PHOENIX models. It is listed in the introduction but not compared with. A short table of differences, in coefficients and in light-curve terms, would establish the pipeline's credibility far better than the μ_cri statement.

(d) *Edge per model or per passband: a chromatic stellar radius.* Sect. 3.1 does not say which edge the tables use. Sect. 3.4 tests both. The coefficient spectra use "one stellar edge per model" and differ from the passband tables by 2–10 ppm, which suggests the tables use a per-passband edge. A per-passband edge means a passband-dependent definition of R*. For μ_cri ≈ 0.05, a change of 0.01 in μ_cri changes r_c by ~5 × 10⁻⁴. For a 1% transit that corresponds to ~10 ppm in depth, comparable to the features sought in transmission spectra. Please state which edge each product uses and quantify the chromatic r_c variation across the JWST passbands. Please also recommend, for spectroscopy, a single edge per model (as in the coefficient spectra) or explain why a per-passband edge is preferable.

(e) *Alternative edge definitions.* Sing et al. (2026) introduce a physically motivated edge at τ_slant = 1, reconstructed from the hydrostatic structure, as an alternative to the maximum-gradient criterion. It should be discussed alongside M1/M2. If feasible, include it in the sensitivity tests of Sect. 3.4.

(f) *Limb brightening in metal-poor cool models.* For metal-poor PHOENIX models at 2300–2500 K in blue passbands, the authors find that the profile brightens towards the limb at μ ≈ 0.01–0.02 and the monotonicity constraint becomes active. Is this brightening physical (e.g. a temperature inversion in the outer layers) or an artefact of where the edge falls? If physical, a monotonicity constraint misrepresents the model. The same question applies more strongly to the R = 3000 coefficient spectra, where limb brightening in strong lines or molecular band heads is plausible.

**6. The comparisons with observations.**

(a) *Sun.*
- The Neckel & Labs (1994) profiles are fifth-order polynomial fits in μ, smooth and with no steep limb, over a limited μ range. That hybrid5 reproduces them to 3–5 ppm (four-term 4–6 ppm) therefore says little about the laws' ability to follow a real limb. Please state the μ range of the observations and that they are polynomial representations.
- For the model–Sun comparison (PHOENIX "furthest", 156–307 ppm), please state the solar parameters and interpolation used, and how the model edge is reconciled with the observational definition of the solar limb. A large part of the PHOENIX discrepancy may come from the edge definition rather than the atmosphere.

(b) *Kepler/TESS.* Reproducing Maxted's (2023) offsets and their removal by a 100 G MURaM correction confirms Kostogryz et al. (2024). Maxted (2023) had already ascribed the h′₁ ≈ 0.006 offset to ~100 G fields. That hybrid5 gives the same h′ as the four-term law is expected, since h′₁ and h′₂ sample μ = 2/3 and 1/3, where all good fits agree. The section can be shortened and framed as a consistency check. The Patel & Espinoza paragraph (agreement "to 0.01" in quadratic coefficients that are "too noisy and degenerate to discriminate") adds little and could be removed.

(c) *JWST.*
- Sing et al. (2026) analyse seven targets. Fig. 5 shows six; WASP-12 is missing. Please explain.
- Sing et al. define ℓ = 1 − (c₁ + c₂) from a *quadratic* fit to the model intensities, and recommend μ_min = 0.2. Fig. 5 is labelled "offsets of hybrid5 predictions". How was ℓ obtained from a hybrid5 profile, and with which μ_min? Δℓ depends strongly on μ_min (Sing et al., Sect. 4), so the comparison is meaningless unless this is specified.
- The claim that Sing et al.'s smaller PHOENIX offset "follows from their linear rescaling of μ" should be demonstrated, with both rescalings applied to the same models, not asserted. Note that their numerical method also removes the four nodes adjacent to μ_cri, which the manuscript does not mention.
- Please show the measurement uncertainties in Fig. 5 (right).
- The mean offset is +0.035 ± 0.031 (scatter), i.e. ±0.013 in the mean, about 2.7σ. After the 100 G correction it is +0.024, about 1.9σ. Please quote significances rather than "removes about a third".

(d) *Recommendation for FGK stars.* The recommendation of MPS-ATLAS Set 2 plus a 100 G correction rests on ~40 Kepler/TESS hosts and six JWST hosts. That is reasonable as a default, but the paper should state:
- the parameter range over which the MURaM correction exists (3200–6800 K per Kostogryz et al. 2026);
- that it is negligible for M dwarfs, so M-dwarf users have no observationally validated choice of grid. This is precisely where the hybrid laws gain most over the quadratic law (72–115 ppm), and the advantage is unvalidated there.

(e) *Δw corrections and empirical offsets.*
- Adding a weight shift Δw, computed for MURaM or Stagger "relative to the matching 1D models", to coefficients from a *different* 1D grid (MPS-ATLAS Set 2) assumes the correction transfers between grids. Which 1D models are "matching" for Stagger (Magic et al. 2015 used their own 1D references)? What is the μ sampling of the 3D intensities, and does w + Δw remain in the physical region?
- The "empirical offsets ... as prior widths" are measured in h′₁, h′₂ (Kepler/TESS) or ℓ (JWST). Mapping two numbers onto five hybrid5 weights is under-determined. Please describe exactly how users should apply them.

(f) Persona remark rather than a requirement: tests against well-studied detached eclipsing binaries with measured LD (e.g. from the DEBCat systems) are a traditional and limb-sensitive check of LDC tables. They would strengthen a paper that presents itself "in the manner of the tables of Claret and collaborators".

**7. Methods that need a fuller description (reproducibility).**

(a) *Constraints.*
- For hybrid4/5, was the fit restricted to the simplex (w_j ≥ 0, Σw_j ≤ 1), which is sufficient but not necessary, or to the exact positive–monotonic region? If the simplex, negative weights that still give a physical profile are excluded and the constraint may bias some fits.
- How are infinitely many monotonicity inequalities enforced for the four-term law: on a μ grid, or analytically?
- The "linear map" argument (Sect. 2 and Sect. 7.1) holds only for unconstrained fits. Are the coefficient spectra fitted with or without constraints? If with constraints, the 10⁻⁹ summation identity cannot hold where a constraint is active. If without, the summed band coefficients may be unphysical in those cases.

(b) *Reference light curves.*
- The "exact" model light curve is computed from a monotone cubic interpolant of the tabulated nodes. Its derivative discontinuities limit the quadrature, and the 10⁻⁸ ppm accuracy is quoted only for the hybrid closed form. Please give the accuracy for the reference model light curve itself.
- For MPS-ATLAS with 24 μ nodes, the "truth" between nodes near the limb is the interpolant. The sub-ppm hybrid5 figures partly measure agreement with a particular interpolant. Please show how they change with a different interpolant, e.g. linear in I(r).

(c) *Grids.*
- Table 4 gives MPS-ATLAS as 3500–7500 K and log g 3.5–5.0, against "2300 ≤ Teff ≤ 7400 K, 3.5 ≤ log g ≤ 5.5" in Sect. 3. Fig. 2 extends to 7600 K.
- State the MPS-ATLAS Set 1/Set 2 differences precisely (abundance sets, mixing-length/overshoot prescriptions), the microturbulence of each grid, and the number and placement of μ points in ATLAS9.
- Give the provenance of the ATLAS9 grid (Castelli & Kurucz 2003? which version?) beyond "the ExoTiC-LD subset".
- MPS-ATLAS spectra extend into the mid-infrared, and Table 4 lists 38 passbands for the mps grids. Why is ATLAS9 needed for MIRI LRS at all?

(d) *Passband integration.*
- State explicitly how the Bessell Johnson–Cousins curves were treated, photon or energy, and why.
- Justify, with reference to the TESS instrument documentation, overriding the SVO energy-counter flag for TESS.
- Previous tables differ on photon vs energy integration, so give the size of the effect for a few passbands.

(e) *Tests and fits.*
- The "method study" is described as 400 model–passband combinations in Sect. 3.2 but as 34 × 8 = 272 in the Table 1 caption. Define it.
- Give the Levenberg–Marquardt starting values and convergence criteria for power-2.

(f) *Availability.*
- No URL or DOI is given for the interpolation code or the coefficient spectra.
- Table 4 lists only the 24 coefficient files. The coefficient spectra, Δw corrections, empirical offsets and the "model systematic, which we tabulate" (Sect. 4) need their own CDS entries and ReadMe descriptions.
- The coefficient spectra (1682 models × ~10⁴ bins × 3 laws) may exceed what CDS normally hosts; please confirm with CDS or provide a Zenodo archive.

---

### Minor comments

1. **Abstract.** Correct "factors of 3–10" (Major comment 2b) and "at most 3 ppm" (3e). Qualify "better than the four-term law" with the parameter count (3a). Rephrase "every hybrid coefficient set is physical" (3f). The Methods paragraph should mention that M1 is one of two edge definitions tested.
2. **Sect. 1, para. 1.** The power-2 law is attributed to Hestroffer (1997) (verified: A&A 327, 199). Its adoption for transits and the name "power-2" are due to Morello et al. (2017), which should be cited here and in Sect. 3.2.
3. **Sect. 1, para. 2.** "Only laws whose intensity is a polynomial in μ² plus terms with poles of order two or more ... give elementary light curves" is a necessity claim. Either prove or reference it, or soften it to "laws of this form give ...".
4. **Sect. 2.** Give one sentence of motivation for the ε ladder and the meaning of K. State why K = 2 and 3 were chosen and not K = 4.
5. **Sect. 2.** "For hybrid2 this region is exactly a triangle": please give its vertices.
6. **Sect. 2.** The q-coordinates "in the manner of Kipping (2013)" for a 5-simplex should be given explicitly or referenced to a verifiable source.
7. **Sect. 3.** The R = 3000 reduction error (1.2 × 10⁻⁵): over which wavelength range and passbands was it measured? Note that the NIRSpec G395H resolution (R ≈ 2700) is comparable to R = 3000, so coefficient spectra cannot be integrated over bins much narrower than a few resolution elements.
8. **Sect. 3.1.** Is μ_cri determined from the band-integrated I(r), from a reference wavelength, or from the bolometric intensity (Major comment 5d)?
9. **Sect. 3.2.** "Claret's least-squares method weights points uniformly in μ": please make this specific to the papers emulated. Claret et al. (2012) quasi-spherical fits and Morello et al. (2020, with Claret as co-author) already use a different weighting.
10. **Sect. 3.2.** With r-weighting the disc centre carries negligible weight, so fixing I(1) = 1 is an interpolation constraint. Please confirm that the central intensity of the fit is therefore exact, and comment on any consequence for the flux normalisation.
11. **Sect. 3.2.** "Imposing flux conservation changed the light-curve error by less than 0.5 ppm". Please also give the change in the coefficients and in the flux error, since the eclipsing-binary community relies on flux conservation (Claret & Bloemen 2011).
12. **Sect. 3.3.** Normalise "ppm" explicitly (relative to out-of-transit flux). Consider tabulating the rms deviation as well as the maximum.
13. **Sect. 3.4.** Reconcile 170 ppm with 93 ppm (Sect. 4) for the PHOENIX vs MPS-ATLAS difference, and define the subsets.
14. **Sect. 4.** "56 of 63 916 above 15 ppm": this is clearly the r-weighted PHOENIX rows only (1682 × 38). Say so.
15. **Table 2.** Add MPS-ATLAS Set 1 and a column for Claret's σ. Consider adding the flux error.
16. **Sect. 4, interpolation.** How are missing PHOENIX nodes (1682 of the full rectangular grid) and the irregular MPS-ATLAS log g spacing handled by the interpolation code? Does linear interpolation in a non-rectangular grid still guarantee physicality (it does for convex combinations, but say how the simplices are built)?
17. **Sect. 5.** "Errors scale approximately as p²": show this (one line in Table 3 or a sentence with numbers).
18. **Sect. 5.** "A central transit ... sweeps the full range of μ": true only to within p of the limb. The grazing case is untested (Major comment 4).
19. **Sect. 6.1.** List the 30 wavelengths or cite the table used. State that model intensities were evaluated at the same monochromatic wavelengths rather than band-integrated.
20. **Sect. 7.** Specify units and formats for every column in the ReadMe: σ in which weighting, the definition of "largest intensity residual", and the definition of "relative flux error" (relative to what?). Explain the two flags.
21. **Fig. 1.**
    - The thick grey band behind the PHOENIX nodes is not explained in the legend or caption.
    - Give the linear threshold of the symmetric-log residual axis.
    - The text should say that hybrid5 and hybrid4 deviate from the limb nodes by up to ~0.1 in I/I(1) for μ ≲ 0.03 with r-weighting; this is visible in both insets.
    - Power-2 and hybrid2 are omitted from this figure but appear in all others.
22. **Fig. 2.** Bands are shown only for the hybrid laws; say so in the legend (currently only in a footnote-like annotation). The right panel's x-axis extends beyond the MPS-ATLAS grid range quoted in Sect. 3.
23. **Fig. 3.**
    - The five shades of blue are hard to tell apart; use distinct line styles or colours.
    - The jagged PHOENIX w₃ (and to a lesser extent w₄) shows what users fitting with priors will face. Please discuss the correlations between weights and recommend which linear combinations are well determined, e.g. via principal components of the coefficient covariance across the grid.
24. **Fig. 4.** Give the linear threshold of the symmetric-log axis. The figure duplicates Table 3; one of the two could move to an appendix.
25. **Fig. 5.**
    - The legend uses the same filled dark circle for "Kepler" and for "MPS-ATLAS 2".
    - The right panel has no error bars on the measurements.
    - The row label "+ small-scale dynamo" is truncated in the rendered figure.
    - The 200 G row is shown but not discussed.
26. **Title and subtitle.** The subtitle says "ATLAS" and the text "ATLAS9"; be consistent. As the ATLAS9 grid covers a single passband, it hardly merits a place in the subtitle.
27. **References.**
    - "SquishierPlanet (2026), in prep." is not acceptable as a citation (Major comment 1).
    - Sing et al. (2026): the first authors are Sing, D. K., Lothringer, J. D., Valenti, J. A., et al. (arXiv:2609.00263).
    - I checked Hestroffer (1997) A&A 327, 199; Howarth (2011) MNRAS 418, 1165; Kostogryz et al. (2024) Nat. Astron. 8, 929; Claret, Hauschildt & Torres (2025) A&A 699, A97; Maxted (2023) MNRAS 519, 3723; Morello et al. (2020) AJ 159, 75. These appear correct; the "% VERIFY" comments in the source should be removed.
    - Please add: Morello et al. (2017, AJ 154, 111); Claret et al. (2012) as the origin of the quasi-spherical fit (already cited, but not for this); Cassese et al. (2024, *squishyplanet*); Carter & Winn (2010); Seager & Hui (2002); Barnes & Fortney (2003); Parviainen & Aigrain (2015, LDTk); Vögler et al. (2005) for MURaM; Castelli & Kurucz (2003) for ATLAS9; Bessell & Murphy (2012) for the passband curves.
28. **Sect. 6.2.** Give the size of the Maxted sample used (number of Kepler and TESS hosts after the b < 0.65 cut).
29. **General.** "hybrid5" etc. begin several sentences in lower case. A&A style would capitalise at the start of a sentence, or rephrase.

---

*Literature checks performed by the referee:*
- Sing et al. (2026, arXiv:2609.00263): 7 targets, the ℓ = 1 − (c₁+c₂) definition, μ_min = 0.2, the numerical rescaling μ = (μ_ph − μ_cri)/(1 − μ_cri) after removing four nodes, the τ_slant = 1 method, and the conclusion preferring the quadratic law.
- Kostogryz et al. (2026, arXiv:2606.21912): MURaM coverage 3200–6800 K, effect negligible for M dwarfs.
- Claret & Southworth (2023): PHOENIX-COND, solar metallicity, Vξ = 2 km s⁻¹, 900 μ points, Levenberg–Marquardt, definitions of M1 and M2.
- Claret, Hauschildt & Torres (2025): 11 JWST passbands, four-term and power-2.
- Morello et al. (2020): "weighted-r" fit with weights proportional to the sampling interval in r.
- Göttingen PHOENIX v3.0 SpecInt README: identical to earlier SpecInt, extended to 5.5 μm.
- Cassese et al. (2024): *squishyplanet*, elliptical/triaxial planets with polynomial limb darkening.
