```
Manuscript: Limb-darkening coefficients for rational laws with elementary transit light curves (second revision)
Referee: Referee 2 — observational transit / transmission-spectroscopy analyst (TESS, Kepler, HST, JWST light-curve fitting)
Recommendation (A&A grade): B   [= Minor revision. Not yet acceptable as it stands: the new hybrid4 material contains one incorrect statement and leaves out the prior that its sampling scheme implies. Both can be fixed with short calculations, and I would not need to see the full manuscript again]
Style and language: (b) well written and self-contained. The main body runs to about 12.4 pages, slightly over the 12-page limit, so it could be cut by ~5% without losing content. Candidates: shorten the hybrid4 ladder-tuning paragraph of Sect. 2 and move its out-of-sample details to the appendix; remove the repetition between Sect. 7.3 and Appendix A
Identify referee: no
```

---

## Confidential comments to the Editor

The authors have dealt with my second report thoroughly and honestly. Points now settled:

- **Limb-rescaling convention.** It is now quantified (median 99–100 ppm), and the authors state which convention is geometrically correct and what this means for the PHOENIX–JWST agreement reported by Sing et al. It also appears in the abstract, the error budget and the Conclusions.
- **Release.** The CDS release now carries the auxiliary tables, the q columns and the expanded form. I checked the expanded form against the weights: they agree to 2×10⁻⁶ in I/I(1). The recommended MPS-ATLAS Set 2 + 100 G prior can now be built from the CDS files with the standalone code.
- **Prior safeguards.** The code now warns below 3500 K, warns about draws that fall outside the grid, applies offsets jointly with their correlation, and refuses an offset in the wrong correction state. I tested all four.
- **Bias test, interpolation tails and power-2 in M dwarfs.** All three are now treated properly.

Since the second round the authors have also changed hybrid4 substantially. They have moved its poles (α = 25/16) and given it a new basis, the "curvature simplex", together with q1–q4 sampling coordinates presented as the analogue of Kipping's q1, q2. The pole change is well tested out of sample, and it makes hybrid4 as accurate as the four-term law, a real improvement. The mathematics of the curvature simplex is correct; I verified the reduction to five linear inequalities and the q-map. The conditions are also well supported: by MPS-ATLAS, by the smooth solar profiles (I checked 30 Neckel & Labs wavelengths), and by every physical quadratic profile and every power-2 profile with α ≤ 2.

From the viewpoint of someone who would use these tables in a fitting pipeline, three problems remain:

1. **One statement is wrong.** The paper says that all 2.9 million tabulated hybrid4 fits satisfy the curvature conditions "with a margin". This number comes from a study of the earlier tables. In the released tables, about 12,700 PHOENIX fits lie exactly on a third-derivative face, including 62–74% of the uniform-μ fits in the 3–5 μm bands. This contradicts the paper's next sentence, which says that some fits would gain from dropping the conditions.
2. **The promoted prior is strongly informative.** The new q1–q4 are offered as the Kipping analogue, i.e. as an uninformative default for free fits. Under that uniform prior, however, I(0) > 0.7 has a probability of 0.8%. The tabulated hybrid4 weights in NIRSpec G395H and MIRI LRS sit at q1 ≈ 10⁻³–10⁻². The paper gives the implied prior for hybrid5 but not for hybrid4.
   - I made a linearised (Fisher-matrix) estimate for one JWST-like transit. A free hybrid4 fit under the uniform-q prior shifts the depth by about −40 ppm, about 0.8 of the posterior width. A prior flat in I(0), with the same Dirichlet direction, reduces this to about +5 ppm.
   - This is an estimate, not a full injection–recovery test, but it shows that the issue is not academic.
3. **hybrid4 has now had three definitions.**
   - The three are: α = 1 poles in the monomial basis; α = 25/16 poles in the monomial basis; α = 25/16 poles in the curvature basis.
   - The ReadMe and the documentation handle the conversion, but the paper does not mention it.
   - The two papers are submitted jointly, so the companion paper must define hybrid4 identically.

These issues do not affect the main conclusions: the tables, the error budget and the conventions analysis. I recommend minor revision. The Editor may wish to make acceptance conditional on:

- the corrected statement;
- the implied hybrid4 prior, with a recommended alternative;
- an explicit basis/version note that agrees with the companion paper.

As before, acceptance should also wait until the companion paper is public. The placeholders (authors, DOI, licence, companion arXiv identifier) remain.

---

## Report to the authors

### Summary of the paper

The authors tabulate coefficients for three linear "rational" limb-darkening laws (hybrid2, hybrid4 and hybrid5). Each law is an even polynomial in μ plus fixed double poles in μ², and its transit light curve is elementary. The coefficients cover spherical PHOENIX-ACES v3.0 models (40 passbands) and plane-parallel MPS-ATLAS models, Sets 1 and 2 (41 passbands). The fits are constrained least squares, uniform in projected radius, with a converged quadrature. The standard laws are computed in the same way, and every fit is scored by the light-curve error it produces.

New in this version:

- **hybrid4.** Its two poles move to an α = 25/16 ladder (ε = 0.0286 and 0.3858), tuned on the light-curve error and tested out of sample. It is tabulated in a new basis, (1 − μ², (1 − μ²)², S3, S4), in which the region w ≥ 0, Σw ≤ 1 is exactly the set of hybrid4 profiles that are positive and monotonic, concave in μ², and have a convex slope. The coordinates q1–q4 sample this region uniformly.
- **Corrected JWST statement.** The limb-rescaling convention is now quantified (about 100 ppm in light-curve terms).
- **Restructured error budget.** It now separates uncertainties from conventions.
- **Code.** The prior code gains safeguards (an M-dwarf floor, a covariance floor, joint offsets and a correction-state check). The standalone code now applies the Δw corrections, and MPS-ATLAS coefficients can be computed for custom throughputs.
- **Bias test.** It is redone with each law free in a constrained region, with a prior mode and with Δb.

### General assessment (Q1–Q7)

**Q1 (why publish).** The case is now clear and honestly made. These tables of all six laws are uniform and optimised for transits. Together with the error budget and the analysis of the limb and passband conventions, they will be useful to transit observers whether or not they adopt the hybrid laws. The re-tuned hybrid4 gives a four-parameter law with an elementary light curve and the accuracy of the four-term law, a genuine improvement over the second version.

**Q2 (assumptions).** Mostly clear. The new curvature conditions (I″ ≤ 0 and I‴ ≥ 0 in μ²) are an assumption added to Kipping's conditions. Their support and their cost need a more accurate statement (Major comment 1). The prior that the new sampling scheme implies is not stated (Major comment 2).

**Q3 (methods).** Fully described for the tables. The q1–q4 map can be implemented from the paper alone, which I did. The S3 and S4 shapes are given to six decimals, and in full precision in the ReadMe, though in a different basis there (minor comment 2).

**Q4 (emphasis).** Good. The abstract and Conclusions give the main results in proportion.

**Q5 (figures and tables).** All are useful. The S3/S4 equation runs off the column. Two citations print as "?". The caption of Fig. 1 does not describe panel (b) (minor comments 1–3).

**Q6 (electronic material).** The CDS ReadMe is now complete and well documented. The q1 column needs more significant digits (minor comment 8).

**Q7 (IAU designations).** No issues. WASP-94A is now labelled correctly.

---

### Response to my second-round comments

**Major 1 (the limb-rescaling convention changes the light curve by ~100 ppm).** *Addressed.*

- Sect. 7.3 now gives 63–151 ppm (median 99) for the Table 3 cases, and a median of 100 ppm (90th percentile 212) on the method study.
- The geometric rescaling is identified as the projection of a sphere. The PHOENIX "agreement" of Sing et al. is attributed to the linear remapping.
- The convention has its own row in Table 6 and appears in the abstract and the Conclusions.

This is exactly what I asked for.

**Major 2 (prior safeguards).** *Largely addressed.* Two items remain.

- **(a) M dwarfs: partly addressed.** The code now warns and substitutes the spread at 3500 K. I confirmed this at 3000 K in G395H. The paper does not, however, give the recommended M-dwarf prior widths in coefficient units, as I asked. Nor does it say whether the 3500 K spread is plausibly representative of 2800 K. The authors' own Fig. 4 shows that w3 to w5 behave differently below 3500 K.
- **(b) Covariance floor: partly addressed.** A floor is recommended (0.01²) and implemented (`cov_floor`). As implemented, however, it only raises diagonal variances that lie below the floor. It does not bound the eigenvalues. With `systematic=True` the diagonal already exceeds the floor, so the near-singular directions survive: for hybrid5 at 5600 K in TESS, the smallest standard deviation along an eigenvector is 0.0018, not ≥ 0.01. See Major comment 4.
- **(c) Flat-q prior: addressed for hybrid5, but not carried over to the new hybrid4 map.** The new map is more informative still (Major comment 2).
- **(d) Offsets: addressed.**
  - h′₁ and h′₂ are applied jointly, with their correlation.
  - An offset in the wrong correction state is refused. I tested this.
  - The passband applicability is given in `offsets.dat` and Note 11.
- **(e) Fisher analysis: addressed for hybrid5.**
  - The set-up is stated.
  - The analysis was redone with the exact quadrature, and for the new hybrid4 basis (`studies/rev4/rv_fisher.md`).
  - The best-measured hybrid5 combination is given.
  - The corresponding hybrid4 combination, now in a new basis, is not given, and there is still no explicit reduced parametrisation for free fits (minor comment 9).

**Major 3 (the release does not match the paper).** *Addressed*, with two small residual inconsistencies.

- **(a)** The auxiliary tables are at the CDS with byte-by-byte descriptions.
- **(b)** The q and expanded-form columns are present and correct. I verified the expanded form against w (2×10⁻⁶) and the q round trip.
- **(c)** The standalone `interpolate.py` applies Δw.
- **(d)** The MURaM construction and its extrapolation are stated, and the size of the correction is in Table 6.
- **(e)** The custom MPS-ATLAS route is implemented.
- **(f)** The kurucz grid is gone, and the code warns about dropped draws.

Residual inconsistencies:

- Sect. 9.2 says the MURaM simulations exist for "six solar-metallicity stars between 3196 and 6807 K". The ReadMe, the response and `corrections.dat` list seven magnetic stars (3196–6807 K). For the 3D (hydro) correction there are six (3518–6807 K).
- The Δw corrections are absent for WFC3 UVIS G280 (both orders) and MIRI LRS. The recommended FGK prior therefore cannot be built in those modes: `prior(..., "JWST_MIRI_LRS", correction="B100G")` raises an error. Neither the paper nor the ReadMe says so (minor comments 6–7).

**Major 4 (bias test).** *Partly addressed.*

- **(a) Equal footing.** Addressed in substance. All laws are now free within a constrained region, and the constraints do not bind in the noiseless fits. The wording is now inaccurate for hybrid4, however. Its "exact physical region" in the Table 5 footnote is the curvature simplex, not the region of conditions I–II, so hybrid4 and hybrid5 are again constrained differently. This does not matter numerically here, but the footnote should say so.
- **(b) Δb added.** Δv or Δ(a/R★) is still not reported, although it was fitted.
- **(c) MPS-ATLAS.** Repeated for nine cases, in the text. Adequate.
- **(d) Noise.** Not done; the authors declined. The new q1–q4 scheme makes this more pressing (Major comment 2).
- **(e) Set-up.** Stated, and Howarth (2011) and Espinoza & Jordán (2016) are credited.

**Major 5 (power-2 in M dwarfs).** *Addressed.*

- The fractions by T_eff, the distribution of c and the cost of the constrained refit are given.
- The unphysical entries in Table 3 are marked, and the recommendation is qualified.
- Table 4, which takes medians over the same cases, is not marked (minor comment 10).

**Major 6 (interpolation tails).** *Addressed.*

- The 90th percentiles are quoted (26 ppm at 200 K; 28–65 ppm in Table 6).
- The scaling is measured at two gap sizes (exponent 1.7) rather than assumed.
- The [M/H] test on MPS-ATLAS (0.3 ppm) is the direct check I asked for.

**Major 7 (advantages specific to the hybrid laws).** *Addressed*, but one new sentence overstates the case.

- The Introduction now separates the advantages shared with the quadratic and four-term laws from those specific to the hybrid laws.
- It says that the broader value for circular-planet users may lie in the uniform tables, the error budget and the conventions.
- The new sentence "their physical regions are simplices" is not accurate (Major comment 1d).

**Minor comments 1–20.**

Addressed:

- 1: rms row added to Table 4.
- 2: grazing errors given relative to the depth.
- 3: Table 6 sources and labels fixed.
- 4: metric labels stated.
- 5: list of four-term wins fixed. I checked the 11/15 and 9/15 counts against Table 3.
- 6: sign of the ℓ shift corrected.
- 7: ℓ-route effect stated.
- 9: Fig. 3 range now within the grid.
- 10: labels fixed.
- 11: mean ± SEM row added.
- 12: ReadMe σ symbols, transit caveat and release size fixed.
- 13: uniform-μ binding fractions given.
- 14: MIRI sentence added.
- 15: remedy for unphysical narrow bins stated.
- 16: quadrature of the coefficient spectra stated.
- 17: marginalisation now stated as conditional on the geometry.
- 18: activity references added.
- 19: maximum chromatic edge (33 ppm) given.

Not changed, with reason:

- 8: the authors explain that 4.8 ppm came from a different sample. Accepted.
- 20: placeholders remain, as acknowledged.

Two items remain open from my first round:

- **Round-1 comment 1.** Sect. 9.2 should still state the language, gradient support and archive of the companion light-curve code. It should also say how the tables can be used with existing transit codes in the meantime (minor comment 13).
- **Round-2 Major 2(a).** The recommended M-dwarf prior widths in coefficient units are still missing (minor comment 14).

---

### New major comments

**1. The curvature simplex: one statement is incorrect, and its support and cost need to be stated with numbers (Sect. 2, Sect. 4.2, Conclusions).**

The construction is correct and attractive. I verified each step:

- the third derivative in x = μ² lies in a two-function span whose ratio is monotonic;
- the reduction to I‴(0) ≥ 0, I‴(1) ≥ 0, I″(1) ≤ 0, I′(1) ≥ 0 and I(0) ≥ 0;
- the vertices;
- the q-map.

The text around it needs correcting in four places.

- **(a) "All 2.9 million tabulated hybrid4 fits satisfy the curvature conditions with a margin" is not true of the released tables.**
  - In the curvature basis, the face w3 = 0 is I‴(0) = 0 and the face w4 = 0 is I‴(1) = 0.
  - I read `phoenix_hybrid4.dat.gz` from the release. Of the 128 bound r-weighted rows:
    - 123 lie exactly on w3 = 0 (2300–2700 K; Johnson B and V, Sloan g and r, STIS G430L);
    - 6 lie on Σw = 1.
  - Of the 12,541 bound uniform-μ rows, 12,351 lie on w4 = 0 and 200 on w3 = 0. They span 2300–7400 K. The fraction is 74% of all PHOENIX rows in IRAC2, 62% in NIRSpec G395H, 9% in PRISM and 4% in WFC3 G141.
  - The third-derivative condition therefore binds in about 12,700 tabulated fits, with zero margin.
  - The statement appears to come from `studies/rev3/h4_scheme.md`. That study was run on the earlier w-simplex tables, where the minimum margins were 0.029 and 0.025. It is also logically inconsistent with the next sentence of the paper. If the curvature constraints were strictly inactive at the optimum of a convex programme, dropping them could not change the fit, so no fit could "gain from dropping them".
  - Please correct the statement and give the binding statistics by condition and method.
- **(b) The uniform-μ hybrid4 fits are shaped by the curvature condition, not by physicality.** Sect. 4.2 says that in the uniform-μ fits the constraints "bind more often on PHOENIX: … 19% for hybrid4". For hybrid4 the binding constraint is almost always I‴(1) ≥ 0, and in the mid-infrared it binds in most rows. The uniform-μ fits are recommended for limb-sensitive uses (eclipsing binaries, interferometry, microlensing), so users should be told that, for PHOENIX at 3–5 μm, these hybrid4 coefficients are set by an added curvature condition. The advice to prefer hybrid5 or the four-term law for such uses is right; please give this as an additional reason.
- **(c) Support by the models: give the numbers.** "PHOENIX profiles violate them mainly at μ < 0.1 … above μ = 0.1, 0.2% of them do" is correct, but it hides the scale. According to the authors' own node finite differences (`h4_scheme.md`):
  - 99.2% of PHOENIX band profiles have I‴ < 0 somewhere, and 13.1% do so at μ ≥ 0.05;
  - 15.7% have I″ > 0 somewhere.

  "They are supported by the models" holds for MPS-ATLAS and, away from the limb, for PHOENIX. Please give these fractions. I would also suggest adding two pieces of supporting evidence, which strengthen the case:
  - every physical quadratic-law profile, and every power-2 profile with α ≤ 2, satisfies both conditions identically;
  - the Neckel & Labs (1994) solar polynomials satisfy both conditions at 29 of the 30 wavelengths. I checked this; the exception is a marginal I‴ < 0 at μ > 0.98 at 366 nm.
- **(d) "Exact" and "complete" are relative to the added conditions.**
  - The Introduction says the hybrid laws' "physical regions are simplices". The physical region of conditions I–II is a simplex only for hybrid2.
  - For hybrid4 the simplex is the curvature region, which fills 9% of the I–II region. That is less than the 13% of the previous, merely sufficient, simplex.
  - For hybrid5 the simplex is sufficient and fills about 8%.
  - The comparison with Kipping's q1, q2 should make the same distinction. Kipping's coordinates cover the whole region of conditions I–II for the quadratic law; q1–q4 cover only the curvature region.
  - Please reword the Introduction, Sect. 2 and the Conclusions accordingly, and define "conditions I–II" (with attribution to Kipping 2013) at first use. At present they are used before they are stated.

**2. The prior implied by q1–q4 is strongly informative against weak limb darkening, and the paper does not say so (Sect. 2, Sect. 8).**

The paper presents q1–q4 as the analogue of Kipping's q1, q2, and "all three maps give uniform priors over the regions". Many users will read this as an invitation to sample q1–q4 ∈ U(0, 1) in a free fit, as they do with Kipping's coordinates. The uniform measure on a 4-simplex concentrates at large depth s = 1 − I(0): the marginal density is 4s³ (Kipping's 2D case gives 2s). From 2×10⁵ draws of the paper's map I find, for hybrid4:

| Quantity under uniform q1–q4 | Value |
|---|---|
| I(0), 5th / 50th / 95th percentiles | 0.013 / 0.16 / 0.53 |
| P(I(0) > 0.5) | 6.3% |
| P(I(0) > 0.7) | 0.8% |
| h′₁, 90% interval | 0.65–0.94 |
| h′₂, 90% interval | 0.11–0.38 |

The tabulated hybrid4 fits lie in the tails of this distribution in the infrared:

| Grid, star | Band | I(0) | q1 |
|---|---|---|---|
| MPS-ATLAS Set 2, 5600 K | NIRSpec G395H | 0.76 | 0.0035 |
| PHOENIX, 3000 K | NIRSpec G395H | 0.68 | 0.011 |
| MPS-ATLAS Set 2, 5600 K | MIRI LRS | 0.86 | 4×10⁻⁴ |

In G395H, h′₁ = 0.953 lies above the 95th percentile of the prior. The hybrid5 statement (P(I(0) > 0.5) ≈ 3%) is now given; the corresponding hybrid4 statement is missing.

To gauge the size of the effect, I made a linearised estimate:

- **Set-up.** One transit with p = 0.1, b = 0.3, a/R★ = 10, P = 3 d. I used the authors' `light_curve` and Fisher matrix, with p, b, a/R★ and t₀ free, at their JWST-like (100 ppm per minute) and TESS-like (1000 ppm per 2 min) noise. The data were taken at the true (tabulated) profile.
- **Method.** I importance-sampled the prior, weighted by the Gaussian likelihood in w marginalised over the geometry, and propagated the posterior mean of w to the depth through the conditional regression of p on w.

| Case | Depth shift, uniform q1–q4 | Posterior width of depth | Depth shift, prior flat in I(0) |
|---|---|---|---|
| MPS-ATLAS 5600 K, G395H, JWST-like | −42 ppm | 49 ppm | +5 ppm |
| MPS-ATLAS 5600 K, PRISM, JWST-like | −39 ppm | 60 ppm | +5 ppm |
| PHOENIX 3000 K, G395H, JWST-like | −36 ppm | 52 ppm | +17 ppm |
| MPS-ATLAS 5600 K, TESS, TESS-like | −204 ppm | 496 ppm | −1 ppm |

The flat-in-I(0) prior is q1′ = s, with the same Dirichlet direction. This is only a Gaussian approximation, and a full MCMC injection–recovery would be better. It indicates, however, that the uniform-q default shifts the depth by roughly 0.4–0.8 of its own posterior width, coherently across bands and with a band-dependent size. That is the kind of systematic a transmission spectrum is sensitive to.

Please:

- **(a)** give the implied hybrid4 prior on I(0), h′₁ and h′₂ (as for hybrid5), and say plainly that U(0, 1) on q1–q4 is not an uninformative default in the near- and mid-infrared;
- **(b)** offer, and implement in `from_q`/`to_q`, a variant flat in the depth s (or in I(0)). Because the hybrid4 map already separates s from the direction, this costs one line, and it is arguably closer to what users expect from "Kipping's q1";
- **(c)** run a small injection–recovery test with noise. A few star/band cases at JWST-like and TESS-like precision would do. Compare the depth bias and coverage of free hybrid4 under (i) uniform q1–q4, (ii) the depth-flat variant and (iii) the tabulated prior with a floor. Comparing the quadratic law under Kipping's q would put the numbers in context. I asked for noisy tests in round 2 (4d). Now that the paper promotes a sampling scheme for free fits, I think they are needed to support it.

**3. Versioning and interoperability of hybrid4 (Sects. 2, 9; ReadMe; companion paper).**

hybrid4 has now had three definitions:

- α = 1 poles (ε = 0.01176, 0.3041) in the monomial basis;
- α = 25/16 poles in the monomial basis (Elómë v0.15.0, and LDhybrid until 2026-10-14);
- α = 25/16 poles in the curvature basis, as tabulated now.

Weights from one definition give a wrong profile, without any warning, in code that expects another. Old q values do not convert one-to-one: the authors' own upstream note says they must be recomputed from w. The ReadMe (Note 15), `docs/interpolator.md` and the standalone files handle this well, for example by giving `w @ BASIS` for codes that take monomial weights. The paper itself says nothing. Please:

- **(a)** state in Sect. 9 that the tabulated hybrid4 weights refer to the curvature basis with α = 25/16 poles, and give the conversion to the monomial-pole weights (the 4×4 matrix B and its inverse, or a pointer to Note 15). State that q values from any earlier basis must be recomputed from w;
- **(b)** carry a version tag ("LDhybrid v1.0; hybrid4 = curvature basis, α = 25/16") in the ReadMe header and in the code outputs, so that chains and priors from different definitions are not mixed;
- **(c)** ensure that the jointly submitted companion paper defines hybrid4 identically, with the same poles and the same basis or an explicit conversion. If the companion paper keeps the monomial basis, a reader of the two papers will find two different sets of "hybrid4 weights" for the same law. I could confirm that the locally installed SquishierPlanet already uses the curvature basis, but I cannot check the companion manuscript;
- **(d)** explain why hybrid5 does not get the same treatment. Presumably, with three poles, I‴ spans three functions and the ratio argument fails. Also explain why its q-map draws the weights in sequence rather than the depth first. A user moving from hybrid4 to hybrid5 meets two different parametrisations of the "same" idea. Offering the depth-first map for hybrid5 as well, which is the same flat Dirichlet reparametrised, would make the two consistent.

**4. The recommended prior widths depend on the basis; the covariance floor does not do what the text says (Sects. 6, 8, 9.2; `interpolate.prior`).**

- **(a)** Sect. 8 recommends "a prior of 0.05 per coefficient" and "a covariance floor (e.g. 0.01²)". Both are defined per weight, and so depend on the basis. The paper's own revision shows this. The hybrid4 "Prior" entry of Table 5 changed from 0.52 to 0.33 ppm (median) when only the basis changed: same span, same poles, same fits (addendum; `studies/rev3/results.md`). The same nominal width means different things for hybrid4, hybrid5 and the four-term law.
  - Please express the recommendation in basis-independent terms, for example as widths on I(0), h′₁ and h′₂, or as a floor on the light-curve covariance mapped to w.
  - Alternatively, state the basis explicitly and show that the conclusion ("removes most of the bias") is insensitive to it.
  - Note also that for the two-parameter laws the prior does not "remove most of the bias": for power-2 the maximum is 62 ppm with the prior against 58 ppm free, and for the quadratic law 58 against 113 ppm.
- **(b)** `cov_floor` raises only the diagonal variances that lie below the floor (`cov + diag(clip(floor − diag(cov), 0))`). It does not bound the eigenvalues. When the diagonal already exceeds the floor, which is the usual case with `systematic=True`, nearly singular directions remain. The smallest standard deviations along eigenvectors at 5600 K are:
  - hybrid5, TESS: 0.0018;
  - hybrid4, TESS: 0.0049.

  The intended floor was 0.01. Since the paper's own argument is that the narrow directions are the problem, please add the floor isotropically (cov + floor·I) or as an eigenvalue floor, and say which in the text.
- **(c)** Narrow Gaussian priors centred on small infrared weights put mass outside the simplex. The documentation says "truncate there if your sampler needs it". The paper should say this too, and note that truncation shifts the effective prior mean in the infrared. In my tests, 0.2% of the hybrid4 prior mass in G395H falls outside the region, so the effect is small there, but it grows with wider floors.

---

### Minor comments

1. **References.** \citet{gopal2019} and \citet{herremans2023} (Sect. 2) have no entries in the bibliography and print as "?" in the PDF. The response states that they were added.
2. **S3/S4 (Sect. 2, Eq. after 5).**
   - The equation overflows the right column in the compiled PDF.
   - The paper writes S3 and S4 on (Π, Π, 1 − μ², (1 − μ²)²). The ReadMe Note 15 writes them on (1 − μ⁴, (1 − μ²)², Π, Π), so the polynomial coefficients differ: 0.001252/0.000778 in the paper against 0.000626/0.001404 in the ReadMe. Both are correct, but a user checking one against the other will be confused. Please use one basis, or state both.
   - Please also give the defining conditions of S3 and S4 in the paper, as the ReadMe does: S3 has I(0) = 0 and I′ = I″ = I‴ = 0 at x = 1; S4 has I(0) = 0, I′ = I″ = 0 at x = 1, and I‴ = 0 at x = 0. Users can then compute them to machine precision.
3. **Fig. 1.** The caption describes only panel (a). Please describe panel (b), the vertices of the hybrid4 region, and its legend.
4. **"Conditions I–II".** These are used before they are defined (Sect. 2, hybrid2 paragraph). Please define them at first use and attribute them to Kipping (2013).
5. **Weighting factors (abstract, Sect. 4.2, Conclusions).**
   - In Table 1 the smallest plane-parallel factor is now the quadratic law's (140/74 = 1.9). The quoted "1.7–6" is a leftover from the earlier hybrid4 (20/12). Please correct to 1.9–6, and 1.9–15 in the Conclusions.
   - The upper bound of 15 on PHOENIX is set by hybrid4, whose poles were tuned on the r-weighted light-curve error. Its uniform-μ error rose from 25 to 58 ppm with the new poles. The factor is therefore partly by construction, and this deserves a clause in the abstract. Without hybrid4 the factors are 3–9 and 1.9–3.5.
6. **MURaM stars (Sect. 9.2).** "six solar-metallicity stars between 3196 and 6807 K": there are seven for the magnetic corrections (3196, 3518, 3707, 4543, 5240, 5778 and 6807 K) and six for the 3D hydro correction (3518–6807 K). Please distinguish the two.
7. **Missing corrections.** `corrections.dat` covers 38 passbands. WFC3 UVIS G280 (±1) and MIRI LRS have no Δw, and the code raises an error when `correction="B100G"` is requested there. Please state this in Sect. 9.2, in the MIRI bullet of Sect. 8 and in the ReadMe.
8. **Precision of q1 (ReadMe).** q1 = s⁴ is printed as F9.6. In the mid-infrared q1 ≈ 4×10⁻⁴–10⁻³, so only 3 significant figures survive. The profile reconstructed from the printed q differs from that of the printed w by up to 8×10⁻⁵ in I/I(1); I checked this on a sample of the MPS-ATLAS Set 2 table. Please print q in E format, or state that w is the primary quantity.
9. **Fisher (Sect. 8).** Please give the best-measured hybrid4 combination in the new basis, since this is the law most users of the elementary light curve will fit. Please also give an explicit recommendation for a reduced free-fit parametrisation, for example fitting the one or two best-measured combinations and holding the rest at their prior.
10. **Table 4.** The power-2 medians and maxima include the unphysical M-dwarf fits marked in Table 3. Please add a note.
11. **Table 5.** Please add Δv (or Δa/R★), which was fitted (round-2 comment 4b). The footnote should say which region hybrid4 and hybrid5 were free in (Major comment 1d).
12. **Coefficient spectra (Sect. 9.1).** The unphysical-bin fractions for hybrid4 (0.25–0.33%) come from `studies/rev3/prod_checks.md`, which was run after the pole change but before the change of basis. Please confirm whether "unphysical" for hybrid4 bins now means outside the region of conditions I–II or outside the curvature simplex. These are different sets, and an unconstrained summed bin may violate the curvature conditions while remaining physical.
13. **Light-curve code (Sect. 9.2).** "The elementary light-curve code … is that of the companion paper": please give its language, whether it provides gradients (relevant for HMC samplers such as exoplanet or juliet with dynesty), and where it is archived. Please also note that, until then, any code accepting a numerical radial profile can use the basis-free expanded-form columns. This is round-1 comment 1.
14. **M-dwarf prior widths (Sect. 8).** Please give the recommended widths in coefficient units (or in h′₁, h′₂ and I(0)) for T_eff < 3500 K, as asked in round 2.
15. **MIRI LRS (Sect. 8).** "A two-parameter law suffices": the quoted 1.7 ppm is for power-2. The quadratic law's median error at 5–12 μm is about 17 ppm. Please name the law.
16. **"With a margin" (Conclusions, `docs/interpolator.md`).** The Conclusions say that every tabulated fit is a profile of the curvature region. This is true, but some lie on its boundary (Major comment 1a). Please keep the wording consistent with Sect. 2 once corrected.
17. **Placeholders.** Authors, affiliations, the Zenodo DOI, the licence, the companion arXiv identifier and the acknowledgements remain to be completed.

---

*Checks performed for this report.* All were made on the released files, read-only; scratch output only.

- Compared the round-2 and current manuscripts word by word, and read the compiled PDF and all figures.
- Recomputed:
  - the Table 1 weighting factors;
  - the hybrid4 and hybrid5 win counts in Table 3;
  - the medians in Table 4 against Table 3.
- Verified the curvature-simplex reduction and the q-map analytically. Checked that the paper's S3/S4 coefficients are consistent with ReadMe Note 15.
- Read `phoenix_hybrid4.dat.gz` (all rows) and samples of the MPS-ATLAS Set 2 table. From them:
  - classified the bound rows by face;
  - computed the distances of the weights to the faces, and q, by passband;
  - verified the expanded-form and q columns against w.
- Drew 2×10⁵ uniform q for the implied hybrid4 and hybrid5 priors on I(0), h′₁ and h′₂.
- Ran `LDTables.prior` (with and without `systematic`, with `cov_floor`, and at 3000 K), `apply_offsets` and the Δw path, including the MIRI and G280 cases.
- Made a linearised Fisher/importance-sampling estimate of the prior-induced depth shift with `ldhybrid.interpolate.light_curve` for four star/band cases.
- Checked the Neckel & Labs (1994) solar polynomials against the curvature conditions.
- Read:
  - `studies/rev3/{results,h4_scheme,h4_completeness,h4_kipping_check,prod_checks}.md`;
  - `studies/rev4/{rv_bound,rv_fisher,rv_validate}.md` and `rv_bias.log`;
  - `docs/upstream/hybrid4_*.md`, `docs/interpolator.md`, `standalone/README.md` and the CDS ReadMe.
