# Response to the referees (draft)

We thank the three referees for careful and constructive reports. We revised the
manuscript throughout. Below, "Sect." and "Table" refer to the revised version.

## Changes common to all reports

1. **Companion paper.** The companion paper (SquishierPlanet) is submitted
   jointly with this one and will be on arXiv before acceptance.
   - The split is deliberate.
     - The companion paper derives the light-curve algorithm for circular and
       elliptical outlines and the pole ladder, and it applies the laws to
       oblate planets.
     - This paper treats the hybrid laws as limb-darkening laws and tabulates
       their coefficients.
   - The new Introduction says so explicitly (Sect. 1, third paragraph).
   - No oblate exoplanet has been detected, so the main use of these tables is
     circular planets. Every test in this paper uses a circular planet, and the
     Conclusions say so.
   - We therefore do not include an oblate-planet demonstration, injection–
     recovery of oblateness, timings against squishyplanet, or a stability
     analysis of the closed form. All of these belong to the companion paper.
   - We now cite the oblateness literature and the existing numerical tools
     (Seager & Hui 2002; Barnes & Fortney 2003; Carter & Winn 2010; Zhu et al.
     2014; Akinsanmi et al. 2020; Cassese et al. 2024; Liu et al. 2025) as
     context.
   - The light-curve code is that of the companion paper and is cited in
     Sect. 9.2.
2. **Self-containment.** Sect. 2 now gives:
   - the properties of Π_ε;
   - the monotonicity of every shape;
   - the reason for the polynomial basis;
   - the origin and meaning of the ε ladder and why K ≤ 3;
   - the vertices of the hybrid2 triangle and its edges;
   - the explicit q-coordinates;
   - the profile prior that a flat q prior implies.

   We define "elementary" and qualify the closed-form statements for the standard
   laws (Mandel & Agol 2002; Maxted & Gill 2019; Chishtie et al. 2025). We
   softened the necessity claim.
3. **A correction to our own pipeline.** In checking the referees' points we
   found that the 2000-point resample uniform in r did not resolve the steep
   limb of spherical models. Only about four points fell at r > 0.998.
   - The least-squares integral is now evaluated exactly: Gauss–Legendre on
     every interval between profile nodes (Sect. 4.2). It is converged at six
     points per interval.
   - All tables, statistics, studies and figures were recomputed with it.
   - The effect is largest on PHOENIX. Median light-curve errors fall from 2.7
     to 1.3 ppm for hybrid5, 8.0 to 3.4 for the four-term law, and 16 to 7.9
     for power-2.
   - MPS-ATLAS changes little.
   - Some comparisons changed as a result, and we report them as they are now:
     - hybrid4 is worse than the four-term law on both grids;
     - hybrid5 beats the four-term law on PHOENIX, but on MPS-ATLAS only by an
       amount below the interpolant ambiguity;
     - power-2 is unphysical in 22% of PHOENIX rows.
4. **Coverage.** We made three changes to the passband and grid coverage.
   - MIRI LRS now comes from MPS-ATLAS (both sets), and ATLAS9 is no longer
     part of the release.
   - We added the NIRSpec G395H NRS1 and NRS2 detector ranges, for 41
     passbands in all.
   - PLATO, Ariel and Pandora have no publicly available response curves, so
     we did not add them. The coefficient spectra give coefficients for any
     response.

## Referee 1

**M1 (unpublished foundation).** See common point 1.
- The ladder origin is now stated in Sect. 2. The full derivation is in the
  companion paper.
- The quadrature validation against the closed-form hybrid5 solution is
  described in Appendix A.

**M2 (r-weighting).**
- **(a) Priority.** We now say plainly that we adopt the weighted-r fit of
  Morello et al. (2020), descended from Claret et al. (2012) and motivated by
  Morello et al. (2017). We present our factors as quantification, not
  discovery.
- **(b) Size of the gain.** With the exact quadrature, Table 1 (now for all six
  laws and both grids) shows median gains over uniform-μ fitting of 3–9 on
  PHOENIX and 1.7–3.5 on MPS-ATLAS. We quote these figures everywhere and
  describe Claret & Southworth's 900-point uniform-μ practice correctly.
- **(c) Circularity, and the cost at the limb.** Both are acknowledged.
  - Sect. 4.2 now gives the limb cost: the residual at μ < 0.05 is 0.11 vs
    0.07, and σ_μ is 1.2–1.6 times larger.
  - It states that the r-weighted coefficients are optimised for transits and
    points other uses to the uniform-μ table.
  - The CDS ReadMe says the same.
- **(d) SPAM optimum.** Not added. Fitting light curves depends on the assumed
  geometry, which is why we fit profiles. Table 5 now reports the parameter
  biases that the tabulated coefficients induce.
- **(e) r vs area.** With the exact quadrature, r is better than or equal to
  area for every law except the quadratic (Table 1).

**M3 (equal footing).**
- **(a) Parameter count.** The parameter count is explicit throughout. We now
  state that the four-term law beats hybrid4 on both grids and that the
  hybrid5 advantage on MPS-ATLAS is negligible.
- **(b) Constraints.**
  - We report the bound fractions for every law and grid.
  - We give the cost of the four-term constraint in bound rows: 4.6 ppm
    against 1.8 ppm for the unconstrained fit, which is unphysical.
  - We give the cost of constraining power-2: all unphysical rows end at
    c = 1, and their median error rises from 4.0 to 17 ppm.
  - Power-2 stays tabulated unconstrained and flagged, as in published tables.
- **(c) Null laws.** Added (Sect. 5). Even polynomials in μ² with 2, 4 and 5
  terms give medians of 193, 44 and 24 ppm. The poles improve on them by
  factors of 5–30.
- **(d) hybrid2 out of sample.**
  - Tuned on PHOENIX alone, ε = 0.34; on MPS alone, 0.39.
  - Half-splits give 0.35–0.38, with a held-out excess of 0.5 ppm.
  - With the exact quadrature, 0.36 is still optimal.
- **(e) "At most 3 ppm."** Removed. The abstract now quotes grid-wide medians
  and the 99th percentile.
- **(f) Physicality.** Rephrased as a property of the method.

**M4 (systematics and empirical evidence).**
- New Sect. 8 with an error budget (Table 7). It includes the interpolation
  error at doubled spacing and its expected size at native spacing.
- Interpolating profiles and refitting is offered in the code. For the linear
  laws it equals coefficient interpolation unless a constraint binds, and it
  helps power-2 (16.5 → 14.3 ppm).
- Parameter biases: new Table 5 (noiseless; fixed and free coefficients; with
  the chromatic spread).
- Geometry: new Table 4 covers p = 0.03–0.3 and grazing transits. Errors grow
  as about p^1.1–1.4, not p² (corrected). Grazing transits give smaller errors
  than central ones.
- We discuss Sing et al. (2026) directly. Our noiseless test shows that a free
  quadratic fits the light curve to a 5 ppm rms yet biases M-dwarf depths by
  up to 113 ppm. Statistical preference and unbiasedness are different
  questions.

**M5 (spherical models).**
- **(a)** The model identity is stated: Husser et al. (2013) SpecInt extended
  to 5.5 µm, a different grid from PHOENIX-COND.
- **(b)** We specify the μ_cri comparison (13 dwarfs; Kepler, TESS and CHEOPS)
  and define M2.
- **(c)** Appendix A contains the pipeline validation:
  - exact reproduction of Kostogryz et al. (2022);
  - comparison with ExoTiC-LD;
  - comparison with Claret & Southworth (2023);
  - the new comparison with Claret, Hauschildt & Torres (2025). Their
    coefficients reproduce our model light curves to 24–26 ppm, the size of
    the model-version difference, and their μ_cri agree with ours to 0.0012.
- **(d)** We state per-band and per-model edges. The chromatic r_c across the
  JWST bands is 1.2 × 10⁻⁴ in the median (2.3 ppm for a 1% transit) and
  3.7 × 10⁻⁵ between NRS1 and NRS2. For spectroscopy we recommend the
  coefficient spectra, which use one edge.
- **(e)** We do not adopt the τ_slant = 1 edge. The SpecInt files do not
  include the model structure it needs. Table 6 shows how strongly the
  conventions at the limb affect ℓ.
- **(f)** The metal-poor limb-brightening rows are flagged and quantified. We
  do not claim to know whether the brightening is physical.

**M6 (observations).**
- **(a) Sun.** We state that Neckel & Labs give polynomial representations, so
  the test checks flexibility, not the steep limb. We give the solar
  parameters and band definition, and note that the limb definition
  contributes.
- **(b) Kepler/TESS.** Presented as a consistency check. We give
  significances and the sample size (24 + 10). We note that 200 G
  overshoots.
- **(c) JWST.**
  - WASP-12 is explained: no per-wavelength ℓ is published for it, and
    omitting it changes their means by ≤ 0.007.
  - ℓ is defined.
  - The rescaling claim is now demonstrated (Table 6): it is −0.029 to −0.036
    for every star, and removing the four nodes has a negligible effect.
  - Significances are 2.8σ → 2.0σ.
  - Fig. 5 has error bars.
- **(d) M dwarfs.** The recommendation gives the MURaM coverage and says that
  no grid is validated for M dwarfs.
- **(e) Corrections and offsets.**
  - The Δw construction and base models are stated.
  - None of the 2.1 million corrections leaves the physical region.
  - The offsets map onto w as exact linear-Gaussian constraints, which the
    code implements.
- **(f) Eclipsing binaries.** Not added. Beyond the transit scope, and the
  uniform-μ table serves eclipsing-binary users.

**M7 (methods).**
- **(a)** We describe the constraint regions and how they are imposed. The
  coefficient spectra are unconstrained fits; we give the unphysical
  narrow-bin fractions (Sect. 9.1).
- **(b)** The accuracy of the reference light curve and its dependence on the
  interpolant are quantified (Sect. 4.3).
- **(c)** We give the grid ranges and the Set 1/Set 2 differences. ATLAS9 is
  removed.
- **(d)** Bessell energy curves are stated. The photon vs energy effect is now
  quantified (44–258 ppm) and highlighted.
- **(e)** The method study is defined (400 = 272 + 128). The power-2 LM setup
  is given.
- **(f)** The CDS listing includes corrections, offsets and systematics. The
  spectra and profiles go to Zenodo. DOI pending.

**Minor comments.** Addressed throughout:
- abstract;
- citations (Morello et al. 2017; Castelli & Kurucz removed with ATLAS9;
  Bessell 1990; Vögler et al. 2005; and others);
- figure legends, symlog thresholds, line styles, Fig. 5 markers and labels;
- capitalisation and "solar analogue";
- "56 of …" rows defined;
- MPS Set 1 added to Table 2.

## Referee 2

**M1, M2.** See common point 1. The oblate-planet demonstration is the subject
of the companion paper.

**M3 (metric vs use).**
- New Table 5 gives the parameter biases with fixed and free coefficients,
  including the chromatic spread.
- rms values are reported alongside the maxima.
- The cadence remark is added (Kipping 2010).
- The p-scaling and grazing geometries are in Table 4.

**M4 (proportion).**
- The error budget and the user conclusions are in Sect. 8.
- 93 vs 170 ppm is reconciled: they are different p ranges; the new value
  is 87 ppm.
- The M-dwarf quadratic deficit lies in the disc interior. Only 8–13 of the
  72–114 ppm comes from μ < 0.05, so it is not an artefact of the steep limb.

**M5 (priors).**
- New code: `prior()` (Monte Carlo over stellar parameters, optional
  grid-to-grid covariance), `observables()` (h′₁, h′₂ and ℓ as exact linear
  functionals of w) and `condition()` (Gaussian update by an empirical offset).
- The q-coordinates and the profile prior they imply are given in Sect. 2.
- The analytic marginalisation is noted in Sect. 2.
- The "tabulated systematic" now has its own CDS file.

**M6 (coverage).**
- **(a)** M-dwarf guidance is given, with its limitation.
- **(b)** MPS-ATLAS coefficient spectra are not provided because of the ODF
  sampling. The code computes MPS-ATLAS coefficients for any throughput.
- **(c)** The narrow-bin unphysical fractions are quantified.
- **(d)** MIRI LRS now comes from MPS-ATLAS. M dwarfs below 3500 K have no
  public grid beyond 5.5 µm.
- **(e)** NRS1 and NRS2 are added. PLATO, Ariel and Pandora have no public
  curves.

**M7.**
- h′ probes the disc interior; this is now stated.
- WASP-12, the definition of ℓ, the significances and the rescaling are
  covered as in R1 M6.
- Maxted (2023) is credited for the offset. The 100 G result is a
  consistency check with Kostogryz et al. (2024).
- A sentence on stellar activity is added.

**M8.** See R1 M2. Our r-weighted coefficients differ from Claret,
Hauschildt & Torres (2025) by 24–26 ppm in light-curve terms (Appendix A).

**M9.** See R1 M3(b, f).

**M10.** Passbands, record counts, Teff ranges, hybrid4 on MPS-ATLAS, and the
hybrid2 tuning set and its sensitivity are all corrected or added.

**Minor comments.** Addressed. Notable points:
- photon vs energy (the TESS 1% statement is clarified);
- NIRSpec ranges;
- the μ < 0.1 weighting wording;
- the SPAM caveat for Patel & Espinoza;
- release size: 486 MB of CDS tables plus 1.1 GB at Zenodo.

## Referee 3

**M1, M2.** See common points 1 and 2. We added the basis definitions and
explained the ladder.

**M3.** Every inconsistency is corrected:
- 3–10 → 3–9 and 1.7–3.5;
- 93 vs 170 reconciled;
- 272 vs 400 defined;
- 41 passbands with per-grid coverage;
- Teff ranges.

**M4.** Abstract rephrased: grid-wide statistics, physicality as a property of
the method, parameter count, and the 100 G result as a consistency check.

**M5.**
- The metric and Morello credit are as in R1 M2.
- r vs area is settled with the exact quadrature.
- The bias table is added.
- Grid-density weighting is stated.

**M6.** Sect. 8 (error budget, fixed vs fitted, recommendations). Interpolation
at native spacing is addressed.

**M7.**
- The interpolant sensitivity is quantified.
- Fig. 1 is explained by the limb weight under r-weighting.
- ℓ vs light-curve sensitivity is reconciled in Sect. 7.3.

**M8.** hybrid2 held-out tests (Sect. 5). The ladder poles are analytic.

**M9.** Sun, Kepler/TESS and JWST are as for R1 M6. On Set 2 vs Set 1: the
corrected Set 2 matches the Kepler hosts better (Sect. 8).

**M10.** Δw physicality: none of the corrections leaves the physical region;
the construction is stated.

**M11.** Appendix A.

**Minor comments.** Addressed:
- definitions of acronyms, M1/M2, ℓ, SPAM and the relative flux error;
- Π properties;
- triangle vertices;
- why the basis was chosen;
- the cubic-interpolation caveat;
- figures and captions;
- Set 1 in Table 2;
- the Table 1 caption;
- Conclusions limitations;
- references;
- placeholders (authors and DOIs to be completed by the authors).
