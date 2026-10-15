# Response to the third round of mock review

All three referees recommended minor revision (B). New numbers come from
`studies/rev5/r5_checks.md` (script `studies/rev5/r5_checks.py`); the others
come from the study files named below.

## Points raised by all three referees

### 1. The "with a margin" sentence was wrong

The sentence is removed from Sect. 2 and replaced by the binding statistics,
which come from `rev4/rv_bound.md` and `rev5` items 8–9.

**r-weighted PHOENIX fits.** 128 (0.19%) are bound, and 123 of these bind on
condition IV at the limb (w3 = 0).
- These are the metal-poor rows at 2300–2600 K whose profiles brighten towards
  the limb (the rows of Sect. 5). They also bind hybrid5. The phrase "give up
  concavity" is gone.
- Fitting them under conditions I–II alone would lower their light-curve
  error from 26 to 15 ppm.

**Uniform-μ PHOENIX fits.** 12,541 (18.6%) are bound, almost all on IV at the
disc centre (w4 = 0), mostly at 3–5 µm.
- Relaxing the condition lowers σ_μ by a median of only 0.2% (at most 25%).
- It raises the light-curve error from 47 to 76 ppm.

**Other changes.**
- Sect. 4.2 now attributes the 19% uniform-μ binding fraction to condition IV.
- The circular Conclusions clause is replaced by the non-trivial statement:
  unconstrained, the r-weighted fits meet the conditions in all but 0.2% of
  PHOENIX rows.

### 2. Over-claims ("exact", "complete", "physical", "simplices")

Conditions I–II are now defined at first use (Kipping 2013). Conditions III–IV
are named and called an assumption about the profile, not a physical
requirement.

The Introduction no longer says that the physical regions are simplices.
- What is specific to the hybrid laws is that all their basis shapes are
  monotonic, so a sufficient simplex exists and the fits fall inside it.
- hybrid2's exact triangle is shared with the quadratic law.

Sect. 2 now states where the conditions hold:
- every physical quadratic profile;
- every power-2 profile with α ≤ 2;
- any profile that is increasing and concave in μ with f‴ ≥ 0;
- Neckel & Labs (1994) at 29 of 30 wavelengths (marginal failure at 366 nm, μ > 0.98);
- MPS-ATLAS profiles.

It also states where they fail:
- **PHOENIX limb.** Violation fractions at all μ / μ ≥ 0.05 / μ ≥ 0.1:
  - III: 15.7 / 4.2 / 0.07%;
  - IV: 99.2 / 13.1 / 0.19%.
- **Wien regime.** A grey atmosphere at 3000 K fails at 300 nm for μ > 0.14.
- **Chromospheric lines.**

It also gives the volume: the hybrid4 simplex fills 9% of the region of
conditions I–II. "Which no model produces" is gone. The Kipping analogy is
called one of construction.

**Why hybrid5 is not treated the same way.** With three poles the endpoint
argument fails. On its simplex IV holds automatically, and all tabulated
hybrid5 fits also satisfy III.

### 3. Prior implied by uniform q1–q4 (Referee 2 major 2; Referee 1 minor 4; Referee 3 2d)

Sect. 2 gives the implied hybrid4 prior: h′₁ 0.65–0.94, h′₂ 0.11–0.38 and
I(0) 0.01–0.52, with P(I(0) > 0.7) = 0.8%. The density of the depth s grows as
s^(n−1), so weak limb darkening is disfavoured. The tabulated G395H fits have
I(0) = 0.68–0.76 (q1 ~ 1e-3 to 1e-2).

**What we added.** The paper notes that q1 = s instead gives a prior flat in
I(0). It recommends a table-centred prior in the infrared, and Sect. 8
repeats this.

**Not done in this revision.**
- The depth-flat variant is not implemented in `from_q`/`to_q`.
- The noisy injection–recovery test was not run.

The referee's linearised estimate (about −40 ppm, JWST-like) is not quoted,
because we have not reproduced it.

### 4. The 0.05 prior and the covariance floor depend on the basis

**The 0.05 prior.** Sect. 6 says that a width per coefficient depends on the
basis: on the monomial-pole weights the same width gives a hybrid4 median of
0.52 ppm instead of 0.33. It also says that the Prior column compares laws
only roughly, and that the prior does not remove the bias of the
two-parameter laws. Sect. 8 recommends widths on h′₁, h′₂ and I(0) as the
portable choice.

**The covariance floor.** Sect. 8 now recommends adding the floor to every
direction (0.01² I), and notes that truncating at the region shifts the prior
mean. The code still lifts only the diagonal (`cov_floor`), so the paper and
the code disagree until the code is changed. **Open.**

### 5. Three definitions of hybrid4, and the companion paper

**Version note.** Sect. 9 now says:
- hybrid4 = the curvature basis with α = 25/16 poles (v1.0);
- earlier definitions give wrong profiles in a code that expects another;
- weights on T are w B;
- q from earlier bases must be recomputed;
- Δw corrections and weight spreads are in the curvature basis.

**Basis matrix.** B is given in full in the new Appendix B, with the defining
vertex conditions of S3 and S4, the face-to-condition map, and a warning to
use the full-precision ReadMe constants. S3 and S4 are now written in the
ReadMe basis, so the paper and the ReadMe agree.

**Companion paper.** That it uses the same poles and basis is marked
`\tbd{confirm}`. This is for the authors to confirm.

### 6. Stale numbers and production errors

**Weighting factors.**
- Abstract: 3–9 and 1.9–3.5, with up to 15 and 6 for hybrid4 (tuned).
- Sect. 4.2: 1.9–6, plus the ranges without hybrid4.
- Equal-per-node factor: "up to 23" is now 28 (PHOENIX) and 32 (MPS-ATLAS).
- Conclusions: 1.9–15.

**Citations.**
- gopal2019 and herremans2023 are added. The volume and page numbers need
  checking: SIAM J. Numer. Anal. 57, 2074 and 61, 2580.
- LDTk is now cited in Sect. 3.

**Equations and figures.**
- The S3/S4 equation is moved to Appendix B and broken over lines; there are
  no overfull boxes.
- The Fig. 1 caption now describes panel (b) and defines S1 and S2.
- The steepness statement is corrected: a pole falls to half at
  0.64√ε (μ ≈ 0.03, 0.13, 0.37 for hybrid5).

**MURaM and Δw.**
- Seven MURaM stars underlie the magnetic corrections and six the 3D ones.
- "Sect. 8 gives the size" is replaced by the number: the constant ratio below
  μ = 0.1 changes the light curves by ≤ 0.7 ppm.
- There are no Δw for UVIS G280 and MIRI LRS (Sects. 8 and 9).

**Convergence.** The reference light curves are converged: doubling the
quadrature nodes changes them by ≤ 0.012 ppm, and finer chord sampling changes
the light-curve errors by ≤ 0.25 ppm (Sect. 4.3).

**Interpolation.**
- The median scales as gap^1.7 and the 90th percentile as gap^1.2.
- The native-spacing estimate is labelled T_eff only. log g and [M/H] give
  about 4–5 ppm, and are given their own budget row.
- MPS-ATLAS step: 0.05 dex.
- Table 6: 26–65 ppm.

**Sun.** hybrid4 is recomputed with the current poles: 2.8–3.8 ppm median, so
"hybrid5, hybrid4, and the four-term law … 3–5 ppm". The solar numbers are
relabelled: the small-scale-dynamo MURaM profiles themselves versus Set 2
corrected by them. The field-free MURaM values (72/33 ppm) are added.

**Other text changes.**
- Table 4 and Table 2 notes on unphysical power-2 rows.
- Table 5 footnote: hybrid4 is free in its curvature simplex.
- Δv bias: < 2×10⁻⁴ relative for the 4- and 5-parameter laws (fixed).
- Free-fit ranking stated as by the maximum; by the median hybrid4 is least
  biased; MPS-ATLAS hybrid4 free bias (1.7 ppm) added.
- Fisher: hybrid4 best-measured combination (0.71, 0.65, 0.00, 0.28).
- Narrow bins: the old fractions were I–II at 0.6–5.3 µm. Added: 0.4% of
  M-dwarf hybrid4 bins outside the curvature simplex; below 0.3 µm, up to 18%
  violate I–II and 20–35% lie outside the curvature simplex (Wien regime).
- 3D-corrected hybrid4 sets: none leave the curvature simplex.
- Magnetic correction (Referee 3 major 5): we know of no criterion, so the
  difference between the 100 G and small-scale-dynamo corrections should go
  into the prior width.
- Abstract: JWST "does not significantly reduce"; Kostogryz et al. credited
  for Kepler/TESS.
- Naming: ATLAS9, and the MURaM and Stagger descriptions.
- Active stars: the effective profile need not satisfy III–IV.
- The letter s for the ladder scale is renamed κ.
- The sampling-check counts now match `h4_kipping_check.md`.

### 7. Page limit

The main body now ends on p. 12 (it was about 12.4 pages).

**Moved to the appendices.**
- New Appendix B: the hybrid4 ladder tuning (with the grid dependence: best
  κ 0.87 for PHOENIX and 0.77 for MPS-ATLAS), the hybrid2 out-of-sample test,
  the curvature-simplex derivation and checks, B, and the q-maps.
- New Appendix C: Fig. 4 and Table 4.

**Condensed.** Sect. 2 is reordered (definitions, elementary light curve,
linearity, regions, sampling). The Introduction roadmap, the Sect. 9 record
list and code feature list, the even-polynomial paragraph, the Sing et al.
paragraph and the Conclusions are condensed.

The margin is small; the author list will add a few lines.

## Not done, or for the authors

- **Placeholders:** authors, DOI, licence, the companion arXiv id, and the
  companion code's language, gradients and archive.
- **Release (needs a re-upload):**
  - version tag in the ReadMe header;
  - Phys/Bound flag wording for hybrid4;
  - Δw basis note;
  - missing G280/MIRI corrections noted in the ReadMe;
  - "constrained to their physical regions" in the ReadMe Description;
  - q printed in E format (the paper now says w is primary).
- **Code:** an isotropic `cov_floor`; the depth-flat q variant.
- **Not done:**
  - noisy injection–recovery tests;
  - M-dwarf prior widths in coefficient units;
  - a separate symbol for the ladder α versus the power-2 α (kept, because
    the ReadMe and the companion paper use α).
