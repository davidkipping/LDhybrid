```
Manuscript: Limb-darkening coefficients for rational laws with elementary transit light curves (revised version)
Referee: Referee 2 — observational transit / transmission-spectroscopy analyst (TESS, Kepler, HST, JWST light-curve fitting)
Recommendation (A&A grade): B   [= Minor revision, at the upper end ("moderate"): several items need short new calculations, none needs new physics]
Style and language: (a) well written, concise and now self-contained. At 12 pages including references and appendix, it is within the limit; no reduction needed
Identify referee: no
```

---

## Confidential comments to the Editor

The authors have responded seriously and in good faith to my first report. The revision is much better.

- The paper is now self-contained on the law side: basis, ladder, constraint regions and sampling map.
- It reports parameter biases with fixed and free coefficients, an error budget, explicit recommendations, a prior/offset recipe in code, the p-scaling and grazing geometries, and a pipeline validation against Kostogryz et al. (2022), Claret & Southworth (2023) and Claret et al. (2025).
- It corrects most of the internal inconsistencies I listed.
- While checking our points, the authors found and fixed a real convergence error in their own fitting quadrature. They recomputed everything and reported the less flattering comparisons honestly (hybrid4 is worse than four-term on both grids; hybrid5's lead on MPS-ATLAS is negligible). This speaks well for the work.

On the companion paper, I accept the reframing. This paper now presents the hybrid laws as limb-darkening laws for circular planets, every test uses a circular planet, and the oblate-planet application is deferred to the jointly submitted SquishierPlanet paper. That answers my first-round comments 1 and 2, provided the companion is public before acceptance and the light-curve code is archived. The Editor may wish to make acceptance conditional on that.

From a user's perspective, a few problems remain. They are fixable, but they matter, because they affect what a pipeline user would actually do with the release:

1. **A factual claim in the JWST section is wrong in size.** The authors say that the limb-rescaling convention which reconciles their PHOENIX results with Sing et al. (2026) changes a transit light curve "by only a few ppm". I checked this with their own released band profiles. Linear and geometric rescaling at the same edge give model light curves that differ by 60–130 ppm for p = 0.1. That is comparable to the PHOENIX–MPS-ATLAS difference. It belongs in the error budget, and it changes how the JWST comparison should be read.
2. **The prior machinery has pitfalls.**
   - The Monte-Carlo covariance is nearly singular. Some directions have standard deviations of 10⁻⁵–10⁻⁴, so it effectively fixes the combinations the data cannot constrain.
   - The "grid-to-grid" term is silently absent for M dwarfs, exactly where the paper tells users to widen priors.
   - The flat "uninformative" q-prior strongly disfavours the weak limb darkening typical of JWST near-infrared bands.
3. **Power-2 for M dwarfs.** 70–95% of the tabulated PHOENIX power-2 rows below 3500 K are unphysical (c > 1). Yet power-2 is recommended as the accurate two-parameter law.
4. **The release does not match the text.** The CDS directory and ReadMe lack the corrections, offsets and systematics files and the q / expanded-form columns described in Sect. 9. The recommended FGK prior (MPS-ATLAS Set 2 + 100 G) cannot be built from the CDS release with the standalone code.
5. **The bias test.** Table 5 does not treat the laws on an equal footing (quadratic and four-term are free and unconstrained, while the hybrids are confined to the simplex). It reports only the depth, and only for PHOENIX.

None of this changes the main results, and I would not need to see the full manuscript again. A point-by-point response, with the corrected Sect. 7.3 statement, the corrected release, and the added rows/columns, would suffice. Placeholders (authors, Zenodo DOI, licence, companion arXiv identifier) remain, as the authors acknowledge.

---

## Report to the authors

### Summary of the revised paper

The authors tabulate coefficients for three linear "rational" limb-darkening laws (hybrid2/4/5). Each is an even polynomial in μ plus fixed double poles in μ², whose light curves are elementary. The coefficients cover spherical PHOENIX-ACES v3.0 and plane-parallel MPS-ATLAS (Sets 1 and 2) models. PHOENIX covers 40 passbands and MPS-ATLAS 41, now including MIRI LRS and the G395H NRS1/NRS2 ranges.

The method is as follows:

- fits are constrained least squares, uniform in projected radius, with an exact quadrature;
- standard laws are computed the same way;
- each fit is scored by the light-curve error it produces.

New in the revision:

- grid-wide statistics (hybrid5 median 1.3 / 0.85 ppm);
- null-law tests;
- p-scaling and grazing tables;
- a noiseless bias test with fixed and free coefficients;
- an error budget with recommendations;
- a demonstration that the PHOENIX discrepancy with Sing et al. (2026) is a rescaling convention;
- a validation appendix;
- prior/offset code.

The main use is stated to be circular planets. The elliptical-outline application is deferred to the companion paper.

### General assessment (Q1–Q7)

**Q1 (why publish).** The paper is now useful to transit observers beyond the hybrid laws. Four things stand out:

- the uniform, transit-optimised r-weighted tables for all six laws;
- the photon-versus-energy warning (44–258 ppm, larger than any fitting choice);
- the error budget;
- the demonstration that ℓ-based "which atmosphere fits JWST" conclusions depend on a limb convention.

The case for adopting the hybrid laws themselves for circular planets is honest but modest (Major comment 7). That is acceptable for a tables paper.

**Q2 (assumptions).** Much improved. The light-curve metric, the reference profile and its interpolant are now stated, along with the fixed-coefficient assumption and the cadence caveat. The remaining unstated assumptions concern the bias test (Major comment 4) and the prior recipe (Major comment 2).

**Q3 (methods).** Largely complete. The exceptions are the free-fit set-up in Table 5, the Fisher analysis (one sentence; no assumptions or eigenvectors), and the construction and coverage of the MURaM Δw (Major comment 3).

**Q4 (emphasis).** Good. The abstract and Sect. 8 now say that atmosphere choice dominates law choice for fixed coefficients.

**Q5 (figures and tables).** All are useful. Table 5 needs Δb and Δv columns, and Fig. 2 has an axis-range issue (minor comments). Figures 2 and 3 could still be one column wide.

**Q6 (electronic material).** The CDS should carry the three auxiliary tables (corrections, offsets, systematics), each with an exact usage note. At present they are absent from the release directory and the ReadMe (Major comment 3).

**Q7 (IAU designations).** No issues. WASP-94A is labelled "WASP-94" in Fig. 5.

---

### Response to my first-round comments

**1. Motivation and closed form in this paper or a citable preprint; show the practical advantage.**
*Largely addressed, by reframing.*

What the revision now does:

- Sect. 2 gives the basis, the meaning and origin of the ladder (its discretisation of a double-pole representation of μ^{2α}), the hybrid2 triangle with vertices and edges, and the q-maps.
- "Elementary" is defined. The closed-form statements for the standard laws are correctly qualified (Mandel & Agol; Maxted & Gill 2019; Chishtie et al. 2025).
- squishyplanet and Liu et al. (2025) are cited.
- The split of the work between the two papers is explicit, and the joint submission with arXiv posting before acceptance is acceptable to me.

What is still missing:

- The light-curve code for the hybrid laws is "that of the companion paper". The paper should state in Sect. 9.2 that this code evaluates circular-planet light curves, its language, whether it gives gradients, and where it is archived.
- It should also note that, until users adopt that code, the tables can be used with any code that accepts a custom radial profile (e.g. the numerical-integration mode of batman).

I no longer ask for timings against squishyplanet in this paper. They belong to the companion.

**2. Demonstrate the oblate-planet use case.**
*Not done here, and I accept the authors' reasoning.*

No oblate exoplanet has been detected. This paper now claims nothing about oblate planets beyond what the companion shows, and the Conclusions state that only circular planets were tested. The remaining condition is that the companion, including its injection–recovery tests, is public before this paper is accepted.

**3. The metric does not represent how coefficients are used; report parameter biases under free and prior-constrained fits.**
*Partly addressed.*

- (a) Table 5 is a real improvement. In particular, it shows that a free quadratic law fits M-dwarf light curves to ~5 ppm rms yet biases the depth by up to 113 ppm, with 57 ppm chromatic spread. That is the result I hoped the authors would test, and it directly engages Sing et al. (2026).
- However, the free fits are not on an equal footing, only depth is reported, only PHOENIX is used, and prior-constrained fits are not shown (Major comment 4).
- (b) rms deviations. The response states that "rms values are reported alongside the maxima", but the manuscript does not do so; only the 5 ppm rms of the free quadratic appears. The numbers exist in the authors' studies (rms at b = 0 by law and p) and should be added to Table 4 or Table 3.
- (c) Cadence. The Kipping (2010) sentence is adequate.
- (d) p-scaling and grazing geometries. Fully addressed (Table 4; p^1.1–1.4, not p²). See minor comment 2 on expressing the grazing errors relative to the depth.

**4. Put law accuracy in proportion to the model systematics; say what a user should conclude.**
*Addressed.*

- Sect. 8 and Table 7 do this clearly, and the abstract and Conclusions state it.
- The 93/170 ppm discrepancy is reconciled (87 ppm median at p = 0.1 over shared points; 172 ppm is the worst case over p ≤ 0.15).
- The μ < 0.05 decomposition (8–13 of 70–113 ppm) answers my question about whether the quadratic deficit is an artefact of the steep limb. It is not, and that is a useful result.
- Remaining issues with the budget itself are in Major comment 1 (limb rescaling) and minor comment 3.

**5. Make the tables usable as priors: uncertainties, correlations, degeneracies, offsets.**
*Partly addressed.*

- (a) A Monte-Carlo prior over stellar parameters, with an optional grid-to-grid term, is implemented. It is well verified for its intended mathematics, but it has practical pitfalls (Major comment 2).
- (b) The correlation structure is summarised in a single sentence: one (TESS) or two (JWST) combinations constrained to better than 0.05. The eigenvectors, the assumed geometry and noise, and a recommended reduced parametrisation are not given (Major comment 2d).
- (c) Mapping h′₁, h′₂ and ℓ offsets onto w as exact linear-Gaussian constraints is elegant and correct, and it answers my question in principle. Usage caveats are missing (Major comment 2c).
- (d) The q-maps and their implied profile prior are given. The claim that the prior "covers the values of FGKM dwarfs" holds only in the optical (Major comment 2b).
- (e) Analytic marginalisation is noted. Addressed.
- (f) The systematic now has its own file in the text, but the file is not in the CDS release (Major comment 3).

**6. Coverage is thinnest where the hybrid laws matter most.**
*Partly addressed.*

- (a) M dwarfs. The recommendation is now explicit and honest: no grid is validated, and priors should be at least as wide as the grid spread. However, the code cannot supply that spread below 3500 K (Major comment 2a).
- (b) MPS-ATLAS coefficient spectra. Their absence is explained. The text says that "the code computes MPS-ATLAS coefficients for any throughput from the model files", but I find no user-facing function that does so; `coefficients_custom` is PHOENIX-only (Major comment 3).
- (c) Narrow bins. Quantified (3.6% of hybrid5 R = 100–1000 bins in M dwarfs) and flagged. Please also say what the user should do (minor comment 15).
- (d) MIRI LRS now comes from MPS-ATLAS. M dwarfs below 3500 K remain uncovered, which is explained. The remark I suggested on weak limb darkening at 5–12 μm is not made, although the authors' own statistics support it (minor comment 14).
- (e) NRS1/NRS2 are added. That PLATO, Ariel and Pandora curves are not public is an acceptable reason, given the custom-throughput route (but see Major comment 3 on that route being PHOENIX-only).

**7. The observational comparisons do not test the hybrid laws' distinguishing feature; JWST clarifications.**
*Addressed, with one new problem.*

- (a) The interior-probe limitation of h′ and the polynomial nature of the Neckel & Labs profiles are acknowledged.
- (b) The absence of WASP-12 is explained, ℓ is defined exactly, and the rescaling claim is demonstrated per star (Table 6).
- (c) Significances, standard errors and the WASP-6 outlier are reported.
- (d) Maxted (2023) is credited, and the 100 G test is described as a consistency check.
- (e) Activity: one sentence (minor comment 18).
- The new problem is the claim, in the same subsection, that the conventions change the light curve by only a few ppm (Major comment 1).

**8. Correct the size of the weighting improvement.**
*Addressed.*

- 3–9× (spherical) and 1.7–3.5× (plane-parallel) are now quoted, consistently with Table 1, which I recomputed from the table.
- The Morello et al. credit and the limb cost are stated.
- The comparison with Claret et al. (2025) in light-curve terms (24–26 ppm) is what I asked for, and it validates the pipeline.

**9. Treat power-2 on the same footing.**
*Addressed in substance, but it raises a new user issue.*

- The constrained refit (all rows end at c = 1, with a median error of 17 ppm) quantifies the cost.
- Keeping the unconstrained, flagged rows "as in published tables" is defensible.
- "Every hybrid row is physical, by construction" is now phrased as a property of the method.
- However, the unphysical rows are concentrated in exactly the stars where power-2 is popular (Major comment 5).

**10. Internal consistency.**
*Addressed.*

- 41 passbands are listed, with per-grid coverage.
- The record counts check out: 1682 × 40 × 2 = 134,560 and 16,974 × 41 × 2 = 1,391,868.
- The T_eff ranges are corrected in the text. Fig. 2 still extends to 7600 K (minor comment 9).
- hybrid4 on MPS-ATLAS is now stated plainly.
- The hybrid2 tuning set and its out-of-sample sensitivity are given.

**Minor comments of round 1.**
Nearly all are addressed, including the Fig. 1 grey band, the legends, Set 1 in Table 2, the M2 definition, the SPAM caveat, the NIRSpec ranges and the TESS photon/energy wording. Those not addressed are carried over below.

---

### New major comments

**1. The limb-rescaling convention changes the light curve by ~100 ppm, not "a few ppm" (Sect. 7.3, Table 7).**

Sect. 7.3 shows convincingly that the PHOENIX Δℓ difference between this paper and Sing et al. (2026) comes entirely from linear versus geometric rescaling of μ at the same edge. The text then states that the same choices "change a transit light curve by only a few ppm (Sect. 4), because ℓ extrapolates a fit to the limb". Sect. 4.4, however, tested M1 versus M2, per-model versus per-passband edges, and μ-node resampling. It did not test the rescaling convention.

The two conventions differ well inside the disc, not only at the limb. For example, at μ_cri ≈ 0.05, a node at original μ = 0.2 maps to μ′ ≈ 0.19 (geometric) but ≈ 0.16 (linear).

I checked this with the released post-cut PHOENIX profiles (`data/products/profiles/phoenix.h5`):

- I reconstructed the original μ from r_c and remapped the nodes linearly.
- I computed light curves with the authors' `lightcurve.lc_errors` (p = 0.1, b = 0–0.9).
- The model light curves of the two conventions differ by 132 ppm (5800 K, TESS), 91 ppm (5800 K, PRISM), 100 ppm (4600 K, PRISM), 108 ppm (3300 K, TESS) and 63 ppm (3300 K, PRISM), with rms 50–100 ppm.

This is comparable to the PHOENIX–MPS-ATLAS difference (87 ppm). Please:

- (a) correct the statement and quantify the light-curve effect of the rescaling convention over the grid;
- (b) state which convention is geometrically correct for a projected disc. I believe the geometric one, given an edge at r_c. Then say what this implies for reading the PHOENIX "agreement" with JWST in Sing et al.: it is obtained with a profile whose interior has been remapped non-geometrically;
- (c) add the convention to Table 7. "Limb definition and μ sampling, 1–3 ppm" is true only for the choices tested in Sect. 4.4. A user who takes PHOENIX coefficients from a tool using the linear convention (and some public tools do) will see differences of this order, and should be told so;
- (d) adjust the conclusion bullet "where the limb conventions themselves matter" accordingly. The point is stronger than the paper currently makes it: the conventions matter for the light curve, not only for ℓ.

**2. The prior recipe needs safeguards before observers can use it as recommended (Sects. 2, 8, 9.2; `interpolate.prior`, `condition`).**

I ran `LDTables.prior` on the released tables. For hybrid5 at 5600 ± 100 K, log g 4.4 ± 0.1, [M/H] 0.1 ± 0.1 in TESS (MPS-ATLAS Set 2), the square roots of the covariance eigenvalues are 7×10⁻⁵, 3×10⁻⁴, 8×10⁻⁴, 2×10⁻³ and 8×10⁻³. With `systematic=True` they become 3×10⁻⁴ … 2.6×10⁻². In NIRSpec G395H they are as small as 9×10⁻⁶ (4×10⁻⁵ with the systematic term).

The paper says the tables' value "lies in the combinations the data do not constrain". A Gaussian prior this narrow in those directions is equivalent to fixing them at the model value. That is far tighter than the model uncertainty the paper itself documents, for example:

- w3 moves by 0.012–0.017 under limb choices;
- Major comment 1 above;
- the offsets from the Sun and from JWST.

Specific requests:

- **(a) M dwarfs.** Below 3500 K only PHOENIX exists, so `systematic=True` silently adds nothing. I confirmed this at 3000 K, G395H: `grids = ['phoenix']`. This is exactly where Sect. 8 tells users to "adopt prior widths of at least the grid-to-grid spread". Please:
  - have the code warn in this case;
  - provide a documented floor for T_eff < 3500 K, for example the spread at 3500–3800 K or the `systematics` file values at the nearest shared point, expressed in w;
  - state the recommended M-dwarf widths in the paper in coefficient units.
- **(b) A covariance floor.** Recommend, and implement as an option, a floor on the prior covariance that represents the known non-parametric uncertainties (limb convention, μ sampling, model–observation offsets). Alternatively, recommend explicitly that users inflate the prior along the poorly constrained Fisher directions.
- **(c) Flat-q prior.** The flat Dirichlet ("uninformative") prior is informative against weak limb darkening. The authors give 90% intervals of 0.58–0.93 (h′₁), 0.08–0.32 (h′₂) and 0.01–0.45 (I(0)) for hybrid5 and state that they "cover the values of FGKM dwarfs". In the infrared they do not. From the tables:
  - the solar-type MPS-ATLAS hybrid5 fit in G395H has h′₁ = 0.955, h′₂ = 0.066 and I(0) = 0.72;
  - the 3000 K PHOENIX fit in G395H has h′₁ = 0.953 and I(0) = 0.56;
  - in PRISM (solar), I(0) = 0.52.

  Under the flat prior, P(I(0) > 0.5) ≈ 3%. A JWST near-infrared user fitting hybrid5 "freely" in q is therefore pulled towards stronger limb darkening. Please correct the statement, show the implied prior against the tabulated values across wavelength, and suggest a better default for free fits (for example, the tabulated prior with inflated covariance, or a Dirichlet concentration < 1).
- **(d) Offsets.**
  - `condition()` applies one observable at a time, "treating them as independent". Maxted's h′₁ and h′₂ offsets are correlated, so please provide a joint 2-D update (or the correlation) in the offsets file.
  - The offsets file has separate rows for uncorrected and 100 G-corrected predictions. The paper should state that the offset row must match the correction state. Applying the uncorrected Kepler offset (+0.008) on top of the recommended MPS-ATLAS Set 2 + 100 G prior double-counts the magnetic correction. The code's own docstring example does exactly this.
  - Kepler/TESS offsets are passband-specific, and the JWST ℓ offset is a single value averaged over 0.86–5.5 μm. Please say over which passbands each row may be applied.
- **(e) Fisher analysis.** Sect. 8 rests a central recommendation on a Fisher analysis given in one sentence. Please give:
  - the assumed geometry: the study files indicate a/R* = 10, P = 3 d, one transit, with p, b, a/R* and t₀ free;
  - that the "JWST-like" and "TESS-like" cases are 100 ppm per minute and 1000 ppm per 2 min;
  - the eigenvectors of the constrained combinations, in a short table or at the CDS;
  - a recommended reduced parametrisation for free fits, for example fitting the one or two best-measured combinations and holding the rest at their prior.

  Please also confirm that the analysis was repeated with the exact quadrature.

**3. The release does not match the paper (Sect. 9; CDS ReadMe; code).**

As a user I would download the CDS tables and the standalone `interpolate.py`. I find:

- **(a) Missing auxiliary files.** The ReadMe and the release directory contain only the 18 law tables. The `corrections`, `offsets` and `systematics` files described in Sect. 9 and Table 8 exist only as CSV files in the project's product directory. Please add them to the CDS with byte-by-byte descriptions and usage notes.
- **(b) Missing columns.** Sect. 9 states that the hybrid records also contain the sampling coordinates q and the expanded form c₀ + c₂μ² + c₄μ⁴ + Σd_k(μ² + ε_k)⁻². The ReadMe byte descriptions have no such columns. The expanded form is exactly what a user implementing the law in another code needs. Please add the columns, or remove the statement.
- **(c) The recommended prior cannot be built from the release.** The magnetic correction is applied only by the full API (`ldhybrid.api.coefficients(correction=...)`), not by the standalone `interpolate.py`. A CDS user therefore cannot reproduce the paper's recommended FGK prior (MPS-ATLAS Set 2 + 100 G).
- **(d) Construction of the MURaM correction.** The MURaM Δw rest on six stars at [M/H] = 0 (3196–6807 K, each at one log g) plus two metallicity variants near solar T_eff. The code interpolates them linearly in T_eff only, clamps them at the ends, and applies them at all log g and [M/H]. The paper should say this, give the T_eff nodes, and state that applying them away from solar metallicity and the tabulated gravities is an extrapolation. The size of the correction should also appear in Table 7. For example, the 100 G Δw changes the light curve by ~58 ppm in 2MASS H at 6800 K.
- **(e) MPS-ATLAS for arbitrary throughputs.** "The code computes MPS-ATLAS coefficients for any throughput from the model files" (Sect. 9.1) is not supported by the user API. `coefficients_custom` uses the PHOENIX coefficient spectra only. Since MPS-ATLAS is the recommended FGK grid, and since PLATO, Ariel and Pandora users are directed to the custom route, please either expose such a function (documenting the model download it needs and its size) or revise the statement.
- **(f) Smaller code points.** The code still offers a `kurucz` grid, which is removed from the paper. `prior()` silently drops draws outside the grid, which biases the mean near grid edges (e.g. T_eff ≈ 3500 K on MPS-ATLAS); a warning is needed.

**4. Bias test (Table 5): equal footing, more parameters, the recommended grid, and noise.**

Table 5 is the most useful addition for observers. To make it reliable:

- **(a) Equal footing.** The text says each law was fitted "with them free within their physical regions". According to the authors' study notes, however, the free quadratic and four-term fits were unconstrained, the hybrid4/5 fits were confined to the simplex via q (a sufficient, not the exact, physical region), and power-2 was confined to c ≤ 1. Please fit all laws within their exact physical regions, or all without constraint, and report whether the constraint binds in the free hybrid4 and hybrid5 fits. The 34 ppm free-hybrid4 bias may partly be the simplex rather than the law.
- **(b) Other parameters.** Please add Δb and Δ(a/R*) (or the chord velocity), which the authors already computed. Several of these matter to observers:
  - the fixed-coefficient quadratic biases b by up to 0.053 at b = 0.3, and v by 1.5%;
  - stellar density from transits is a standard by-product;
  - b and a/R* biases propagate into the depth chromaticity.
- **(c) MPS-ATLAS.** Please repeat for MPS-ATLAS Set 2, the recommended FGK grid. The test currently uses only PHOENIX, whose steep limb is the most demanding case.
- **(d) Noise.** The "bias" of free hybrid5 (rms residual 0.12 ppm, yet up to 6 ppm in depth) shows a near-degeneracy between depth and limb darkening. In real data this appears as variance, not bias. The quantity an observer needs is the depth uncertainty and bias at realistic noise, for each law: fixed, free, and under the tabulated prior. These could be Fisher widths at the TESS-like and JWST-like noise levels already used in Sect. 8. A prior-constrained row would also directly test the recommendation that "the tables are most useful as priors".
- **(e) Set-up and context.** Please state that t₀ and the baseline were fixed, and that the light curves were sampled at 600 points without integration. Please also cite Espinoza & Jordán (2016) and Howarth (2011) where the free-quadratic result is presented: that freeing the coefficients does not remove the depth bias is consistent with their findings, and the new numbers quantify it for M dwarfs.

**5. Power-2 is unphysical for most M-dwarf rows (Sects. 4.2, 6, 8).**

The paper reports that 22% of PHOENIX power-2 rows have c > 1. The fraction is strongly concentrated. From the released `phoenix_power2` table (r-weighted):

| T_eff (K) | Unphysical rows |
|---|---|
| 2300–3000 | 74% |
| 3000–3500 | 69% |
| 3500–4000 | 30% |
| ≥ 4000 | ≤ 2% |

Below 3500 K, the fraction is 80% in TESS, 95% in WFC3 G141, 81% in NIRSpec PRISM and 94% in NIRISS SOSS order 1.

These rows have I(0) < 0, so a user who samples power-2 in the usual (h₁, h₂) or q parametrisation with a tabulated prior centre will place that centre outside the support. Fig. 1 (right panel) shows the effect: the power-2 profile reaches zero at the limb for the 3000 K star.

The power-2 entries for the M dwarfs in Tables 3 and 4 are therefore, I presume, for unphysical coefficients. The free power-2 fits, which are constrained to c ≤ 1, give depth biases up to 58 ppm (late M, TESS). Please:

- give the fraction by T_eff in Sect. 4.2;
- mark the affected entries in Tables 3 and 4;
- qualify the recommendation in Sect. 8 ("power-2 is more accurate") with "except for M dwarfs with PHOENIX models, where most tabulated power-2 coefficients are unphysical and the constrained law errs by ~17 ppm in median".

**6. Interpolation error, not law error, is the largest pipeline term. Quote its tail and verify the native-spacing claim (Sect. 5, Table 7).**

The leave-one-out medians (4–15 ppm) are given. The tails are not. The authors' own check gives a 90th percentile of 28 ppm for hybrid5 across a 200 K gap (35 ppm in another run), which is 20 times the law error.

The statement that "at the native spacing the error should be about four times smaller" is an assumption (quadratic scaling of linear-interpolation error) and is not verified. MPS-ATLAS, with its 0.05–0.1 dex [M/H] spacing, offers a direct test. Please:

- quote the median and the 90th percentile in Table 7;
- either verify the factor of four or drop it.

**7. Be precise about which circular-planet advantages are specific to the hybrid laws (Sect. 1, Sect. 2).**

The Introduction says that for circular planets the hybrid laws offer the following:

- linearity;
- convex constrained fits;
- analytic marginalisation under Gaussian priors;
- exact summation of coefficient spectra;
- a physical region containing a simplex.

The first four apply equally to the quadratic and four-term laws, which are also linear in their coefficients. analytic marginalisation for polynomial laws is exactly the starry case the authors cite. Only two advantages are specific to the hybrid laws:

- a simple sufficient physical region, which helps interpolation and priors and contrasts with the awkward exact region of the four-term law;
- the elementary light curve.

For the accuracy gain, the authors' own budget shows that for fixed coefficients hybrid5 versus four-term (1.3 vs 3.4 ppm on PHOENIX; 0.85 vs 1.1 on MPS-ATLAS) is negligible against the 87 ppm grid spread. Please reword the paragraph so that a circular-planet user can see what they gain by switching, and what they do not. This would strengthen rather than weaken the paper. Its broadest value to circular-planet users is arguably the uniform, transit-optimised tables for all six laws, together with the budget and the conventions analysis, and the paper could say so.

---

### Minor comments

1. **rms deviations (Sect. 6, Tables 3–4).** The response says these are reported. Please add the rms column (already computed: e.g. hybrid5 0.50, four-term 1.5, quadratic 30 ppm median at p = 0.1, b = 0). See round-1 comment 3b.
2. **Grazing transits (Sect. 6, Table 4).** "Grazing transits give smaller errors than full ones" is true in absolute ppm. Relative to the much shallower grazing depth, however, the errors are larger. For example, the quadratic's 41 ppm at b = 0.95 against a depth of roughly 0.2–0.3% is a larger fractional error than 76 ppm against 1%. Grazing geometries are also where depth, b and limb darkening are most degenerate. Please express these errors relative to the depth, or add a grazing row to the bias test.
3. **Table 7 (error budget).**
   - "Models vs Sun (MPS-ATLAS, MURaM) 23–150": the upper value does not appear in Sect. 7.1, which gives 23–25 ppm for the closest models and 156–307 for PHOENIX. Please give the source.
   - "Photon vs energy weighting" is a user error rather than an uncertainty once the convention is right. Consider listing it separately as a pitfall.
   - Add the size of the magnetic correction and the rescaling convention (Major comments 1 and 3d).
4. **"Method study" medians (Sects. 4.2, 5).** Table 1 gives r-weighted medians of 1.9 (hybrid5) and 56 ppm (hybrid2) on the method study. Sect. 5 gives 1.2 and 41 ppm "on the method study". The difference is the p-range of the metric (0.03–0.15 versus 0.1). Please state the metric in each place.
5. **Sect. 6.** The four-term law wins in TESS for the mid and early M dwarfs, the K dwarf and the solar analogue. For the late M dwarf in TESS, hybrid5 wins (1.3 vs 1.6 ppm, Table 3). The count of 11 hybrid5 wins is correct, but the list "the M and K dwarfs and the solar analogue" names five cases and so is inconsistent with it.
6. **Sect. 7.3, sign.** "Their rescaling lowers ℓ by 0.029–0.036" should be "raises the model ℓ (lowers Δℓ = observed − model) by 0.029–0.036". Table 6 and the authors' study notes (lin − geom in Δℓ) imply this sign.
7. **Sect. 7.3 / Fig. 5.** ℓ is computed through the hybrid5 fits. The authors' notes show that this differs from ℓ of the model profile by 0.006–0.008 for MPS-ATLAS, which is comparable to the 100 G shift (0.011) and to the SEM (0.012–0.013). Please state this. It reinforces the point that ℓ is ill-conditioned.
8. **Sect. 5.** The T_eff interpolation error is 4.8 ppm in the authors' check, not 4 ppm.
9. **Fig. 2.**
   - The MPS-ATLAS curves extend to ~7600 K, but the grid ends at 7500 K. The PHOENIX curves start at 2400 K, but the grid starts at 2300 K.
   - If these are bin centres, say so in the caption and use bins that do not extend beyond the grid. This was round-1 comment 10c.
   - The figure would read well at one-column width.
10. **Labels.** Fig. 4 row label "Sun analogue" should be "solar analogue", as in the text and Table 3. The Fig. 5 label "WASP-94" should be "WASP-94A".
11. **Fig. 5 (right).** Most per-star error bars are smaller than the symbols. Say so, and consider adding a mean ± SEM marker per model to match the numbers quoted in Sect. 7.3.
12. **CDS ReadMe.**
    - The abstract text has lost the σ symbols ("from 2.8to 2.0.").
    - The response says the ReadMe explains that r-weighted coefficients are optimised for transits and that uniform-μ fits are for eclipsing binaries, interferometry and microlensing. Note (2) says only "recommended". Please add the sentence.
    - Please give the total size of the release in the paper (≈ 0.5 GB compressed for the CDS tables; 1.0–1.1 GB on Zenodo; the paper and the response differ).
13. **Uniform-μ tables.** The uniform-μ hybrid5 fits on PHOENIX are constrained in 10% of rows (against 0.5% for r-weighted). Users of the secondary table should be told, since a binding constraint breaks the exact spectral summation and the linear-map arguments.
14. **MIRI LRS (Sects. 3, 8).** On MPS-ATLAS the median light-curve errors in MIRI LRS are 0.33 ppm for hybrid5 and 17 ppm for the quadratic law. A sentence stating that a two-parameter law suffices at 5–12 μm for FGK hosts would be useful guidance. Please also state explicitly in Sect. 3 that no MIRI coefficients exist below 3500 K.
15. **Unphysical narrow bins (Sect. 9.1).** Flagging is good. Please also recommend a remedy, for example a constrained refit of the bin's summed profile, and consider implementing it, since the band-profile machinery exists.
16. **Coefficient spectra (Sect. 9.1).** These are fitted with trapezoidal r-weights on the native nodes, not with the exact quadrature that motivated the recomputation. The 1–2 ppm agreement with the band tables shows this is harmless, but please say explicitly that the spectra use a different quadrature.
17. **Sect. 2, analytic marginalisation.** For the hybrid laws, the light curve is linear in w only at fixed geometry. With geometry free, the marginalisation is conditional, as in starry. A clause would avoid over-reading.
18. **Activity (Sect. 8).** The sentence is welcome but gives no scale. Please cite work that quantifies how unocculted and occulted spots and faculae shift fitted limb darkening (e.g. Ballerini et al. 2012, A&A 539, A140; Csizmadia et al. 2013), or give an indicative widening.
19. **Sect. 4.1.** The chromatic-edge statistic is useful. The maximum across JWST bands (1.6×10⁻³ in r_c, i.e. 33 ppm in depth) should be quoted alongside the median and 90th percentile, since outliers matter for spectra.
20. **Placeholders.** Authors and affiliations, the Zenodo DOI, the licence, the companion arXiv identifier and the acknowledgements remain to be completed.

---

*Checks performed for this report.* All on the authors' released files, read-only; scratch output only:

- recomputed the Table 1 improvement factors and the Table 6 means;
- confirmed the CDS record counts and the band lists from the ReadMe;
- ran `LDTables.prior` and `observables` on the released hybrid5 tables (covariance spectra; flat-Dirichlet prior percentiles from 2×10⁵ draws against tabulated IR fits);
- counted unphysical power-2 rows by T_eff and band in `phoenix_power2.dat.gz`;
- compared linear and geometric limb rescaling on five PHOENIX profiles from `data/products/profiles/phoenix.h5` with `ldhybrid.lightcurve.lc_errors`;
- inspected `interpolate.py`, `api.py`, `predict.py`, `qspace.py`, `corrections.csv`, `empirical_offsets.csv` and `systematics.csv.gz`;
- read `studies/rev/results.md`, `val_results.md`, `val_code_checks.md`, `production_stats.md` and `studies/p2/p2f_fisher.md`.
