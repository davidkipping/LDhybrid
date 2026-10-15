# Addendum to the second-round response: changes made after it

The response to the second round (`studies/review/round2/response_round2.md`) covers:
- the revision for the limb conventions, the completeness of the release and the
  prior safeguards;
- the new hybrid4 poles (α = 25/16).

Since then we have made two further changes:

1. **A&A macro package.** The manuscript now uses the official A&A macro package
   (v9.4). The layout has changed; the content has not.

2. **hybrid4 basis: the curvature simplex.** hybrid4 is now tabulated on a new basis
   with the same span: the same poles, profiles, light curves and light-curve errors.
   - **The new basis.** The shapes are 1 − μ², (1 − μ²)² and two pole-dominated
     shapes S3 and S4 (Sect. 2). They are the vertices of the region where the
     profile is:
     - positive and monotonic (Kipping's conditions);
     - concave in μ² (I″ ≤ 0);
     - of convex slope in μ² (I‴ ≥ 0).
   - **What it gives.** The region is exactly w ≥ 0, Σw ≤ 1. Sampling coordinates
     q1–q4 ∈ [0, 1] cover it uniformly and completely, as an analogue of
     Kipping's q1, q2.
   - **Why we changed it.**
     - The previous region (the simplex on the monomial-pole shapes) was
       sufficient but not complete. It covered 13% of the region allowed by
       Kipping's conditions, and 44% of it was profiles convex in μ².
     - The new region contains 9% of the Kipping region. 96% of what it excludes
       is convex in μ² somewhere.
     - All 2.9 million tabulated fits satisfy the two curvature conditions.
     - Where the constraint binds, 123 PHOENIX fits at 2300–2600 K in blue bands
       would gain from dropping it; their median light-curve error would fall from
       26 to 15 ppm.
   - **What was rerun, and what changed.** The tables, the method-study and bias
     tables and the binding fractions were rerun in the new basis. The numbers
     that changed are:
     - the hybrid4 column of Table 1 (equal-per-node and uniform-in-μ rows only);
     - the hybrid4 prior entry of Table 5;
     - the binding fraction for the uniform-μ PHOENIX fits (20% to 19%);
     - the limb-sensitivity bound (2.7 to 2.4 ppm).
   - **The release.** The CDS hybrid4 tables, the ReadMe (Notes 7 and 15), the
     interpolators and the coefficient spectra now use the new basis. Binding fits
     lie exactly on the faces. hybrid2 and hybrid5 are unchanged.
