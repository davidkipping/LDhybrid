# Revision plan from the three mock referee reports

Framing adopted by the authors (2026-10-10):
- This paper is submitted jointly with the companion SquishierPlanet paper, and the intro says so.
  The problem is split in two: this paper = the laws as LD laws + their coefficients;
  companion = the light-curve algorithm (circular and elliptical outlines), the derivation of
  the pole ladder, and the oblate-planet application.
- The main use of the tables is spherical (circular-outline) planets. No oblate planet has been
  announced. Every test in this paper uses circular planets, and the paper says so.
- The light-curve algorithm stays out of this paper.

Referee keys: R1 = LD specialist, R2 = transit observer, R3 = generalist; M = major, m = minor.

## A. Answered by the framing (intro rewrite + response letter)

A1. Companion paper: cite it as jointly submitted (arXiv ID once posted), and state the division
    of labour in the last paragraph of Sect. 1 and in Data availability. [R1-M1, R2-M1, R3-M1]
A2. Sharpen Sect. 1 for circular planets. What the hybrid laws give a circular-planet user:
    elementary light curves; linearity in w (convex fits, analytic marginalisation, exact
    coefficient spectra); a physical region that is a simplex, so interpolation stays physical;
    accuracy that matches or beats the four-term law. Elliptical outlines are one paragraph of
    motivation, pointing to the companion. [R3-Q1 "niche product or general replacement?"]
A3. Wording:
    - define "elementary" (no elliptic or special functions);
    - soften the necessity claim ("only laws of this form…" → "laws of this form…");
    - qualify "every standard law has a closed form": Mandel & Agol (2002) four-term uses
      special functions; power-2 has qpower2 (Maxted & Gill 2019) and an exact Appell-F1 form
      (Chishtie et al. 2025). [R1-M1, R3-m2, R1-m3]
A4. Context citations, no comparison:
    - oblateness (Seager & Hui 2002; Barnes & Fortney 2003; Carter & Winn 2010; Zhu+ 2014;
      Akinsanmi+ 2020; Liu+ 2025);
    - squishyplanet (Cassese+ 2024).
    [R1-M1, R2-m4, R3-M2d]
A5. Law definition only, nothing on the light-curve algorithm, to keep the paper self-contained:
    - one sentence on the origin of the ε ladder and what K controls; why K = 2, 3;
    - state that Π_ε(0) = 1, Π_ε(1) = 0, and that every T_j decreases monotonically;
    - give the vertices of the hybrid2 triangle;
    - explain why the basis is {1−μ⁴, (1−μ²)²};
    - define the q-coordinates;
    - optional: a small panel showing the basis shapes.
    [R1-m4–6, R2-m5, R3-M1b, R3-m4, m5, m8]
A6. Response letter only:
    - oblate demo, injection–recovery, timings against squishyplanet, stability near the
      innermost pole, differentiability → companion paper; [R1-M1, R2-M1–2, R3-M2a–c]
    - light-curve code release → companion code, cited in Data availability. [R2-M1]

## B. Factual fixes and consistency (text only)

B1. "Factors of 3–10" → ×2.1–4.2 vs uniform μ, and up to ×11 vs per-node weighting. Fix in the
    abstract, Sect. 3.2 and the conclusions. Describe Claret's practice accurately: Claret &
    Southworth resample to 900 μ points, which is the uniform-μ row. [all three]
B2. 93 vs 170 ppm: define both quantities, or report one. [all three]
B3. Method study: 272 vs 400 (PHOENIX 272 + MPS2 128). Define it once, give p and b for Table 1,
    and name the hybrid2 tuning set (7524). [R1-7e, R3-M3c, m32]
B4. Passbands: list all 39 and give per-grid coverage. PHOENIX/MPS have 38 (no MIRI); MIRI LRS
    is ATLAS9 only. Record counts then match. [all three]
B5. Teff and log g ranges per grid (PHOENIX 2300–7400 K; MPS 3500–7500 K); MPS [M/H] steps;
    Fig. 2 axis extent. [all three]
B6. Define M1 and M2. State that the tables use a per-band edge and the spectra a per-model
    edge, and explain the 2–10 ppm difference. [all three]
B7. Define ℓ and how we compute it (quadratic fit, μ ≥ 0.2, from our profile). Explain why
    WASP-12 is missing, or add it. Give significances: +0.035 ± 0.013 (2.7σ) → +0.024 (1.9σ);
    drop "removes a third". [all three]
B8. Define acronyms (ACES, SpecInt, ODF, SPAM, SOSS, LRS, MURaM, SVO, NGTS). Fix wording:
    - "~1%" for TESS (1% of what);
    - "outer 0.5% … a tenth of the weight";
    - "equivalent within the noise" (these are deterministic numbers);
    - "56 of 63,916" (PHOENIX, r-weighted).
    [R1-m14, R2-m8, m10, R3-m9, 11, 13, 18]
B9. Grid provenance:
    - PHOENIX v3.0 = Husser+ 2013 SpecInt extended to 5.5 µm; microturbulence;
    - MPS Set 1 vs Set 2 differences;
    - ATLAS9 version, μ points, why ATLAS9 for MIRI;
    - Bessell photon/energy convention; justification for the TESS flag.
    [R1-M5a, 7c–d]
B10. Methods detail:
    - constraint region per law (hybrid4/5: simplex);
    - how the four-term monotonicity is enforced;
    - the linear-map identity holds only for unconstrained fits (are spectra constrained?);
    - power-2 LM setup;
    - ReadMe column definitions;
    - Δw sign convention and the reference 1D models.
    [R1-M7a, m20, R2-m6–7, m19, R3-M10]
B11. Figures:
    - Fig. 1: grey band, symlog threshold, hybrid2/power-2 absent, log g and [M/H];
    - Fig. 2: Set 2 label; one column;
    - Fig. 3: distinguishable lines;
    - Fig. 4: symlog threshold; possibly move to an appendix;
    - Fig. 5: duplicate marker, error bars, truncated label, 200 G row, shading meaning,
      caption ("hybrid5 predictions").
    [all three]
B12. Housekeeping:
    - subtitle ATLAS → ATLAS9 (or drop it);
    - capitalise "hybrid5" at sentence starts;
    - "Sun analogue" vs "solar analogue";
    - author placeholders; remove "% VERIFY" where checked.
B13. Citations to add (verify each before citing; the mock referees can be wrong):
    - Morello+ 2017;
    - Castelli & Kurucz 2003;
    - Bessell & Murphy 2012;
    - Vögler+ 2005;
    - Parviainen & Aigrain 2015;
    - Espinoza & Jordán 2015/16;
    - Short+ 2019;
    - Maxted 2018;
    - Maxted & Gill 2019;
    - Chishtie+ 2025;
    - Claret+ 2012 (quasi-spherical fit);
    - Sing+ 2026 full author list.
    Also credit Maxted (2023) with the original 100 G attribution.

## C. Claims to temper (text)

C1. Abstract "at most 3 ppm": give grid-wide figures (PHOENIX median 2.7, 99th percentile 6.3,
    max 50) and say the 15 cases are selected examples. [all three]
C2. "Every hybrid set is physical" is true by construction. Report the constraint-active
    fraction instead (needs D6). [all three]
C3. Parameter count:
    - hybrid5 has 5 parameters, four-term 4; the fair comparison is hybrid4 vs four-term;
    - hybrid4 loses on MPS (6.8 vs 1.4) and below 4000 K;
    - give a law recommendation per grid;
    - replace "bracket the standard laws". [all three]
C4. r-weighting: say plainly that this is Morello+ 2017/2020's weighted-r fit, which we
    quantify. Acknowledge that the metric favours it. r-weighted coefficients are for transits;
    point EB, interferometry, microlensing and RM users to the uniform-μ table. [R1-M2, R3-M5a]
C5. Fig. 1: under r-weighting the region μ < 0.04 carries < 0.1% of the weight, so hybrid5
    misses the limb nodes by up to 0.1. Reword "restores the steep limb" (that is what the
    poles allow, not what the r-weighted fit does). [R1-M2c, R2-m20, R3-M7b]
C6. Error budget: a new short "Limitations and recommendations" subsection.
    - Rank: law error < limb convention ≈ interpolation < grid-to-grid < observed offsets.
    - Fixed vs fitted coefficients; what follows for users.
    - Leave-one-out at doubled spacing overestimates interpolation error by about 4×.
    - Grid-density weighting of medians.
    [all three]
C7. Engage with Sing+ 2026, who find the quadratic law preferred:
    - their test is free fits of FGK white-light curves at current precision;
    - our tables serve fixed or prior use, M dwarfs and spectroscopic bins;
    - state where higher-order laws matter and where they do not. [R1-M4, R2-M3a]
C8. Kepler/TESS is a consistency check: 100 G was calibrated on the same data, and h′ probes
    the disc interior. Shorten; drop or recast the Patel & Espinoza paragraph. [all three]
C9. Sun: the Neckel & Labs values are 5th-order polynomials over μ ≥ ~0.12, so this tests
    flexibility on a smooth profile. Qualify the abstract. [all three]
C10. "Exact" → "exact for the PCHIP-interpolated profile". [R1-M7b, R3-M7a]
C11. Recommendations:
    - say why Set 2 rather than Set 1;
    - M dwarfs: PHOENIX only, no validated choice, magnetic effect negligible;
    - MURaM coverage 3200–6800 K;
    - one sentence on stellar activity;
    - one sentence on cadence/supersampling.
    [R1-M6d, R2-M6a, M7e, R3-M9d]
C12. How to use the empirical offsets: map ℓ / h′ to w, or tabulate them in observable space
    and say so. Where the "tabulated model systematic" column lives. [R1-M6e, R2-M5c, f]
C13. Physical notes:
    - the jagged w3/w4 in Fig. 3: weakly constrained, do not read the weights physically;
    - non-monotonic Teff dependence of w1/w2 (molecular opacity onset?);
    - limb brightening in metal-poor 2300–2500 K blue bands.
    [R1-M5f, m23, R2-m21, R3-m29]

## D. Small new analyses (existing pipeline, hours not days)

D1. Even-polynomial null laws (2/4/5 terms in μ², no poles) in Tables 2–3: shows what the poles
    buy. [R1-M3c]
D2. Constrained power-2 (c ≤ 1, α > 0): cost in light-curve error; how often it binds.
    DECISION: also regenerate the power-2 table columns? [R1-M3b, R2-M9, R3-M4b]
D3. hybrid2 pole held-out test:
    - tune on PHOENIX, test on MPS (and the reverse);
    - error at ε = 0.30 / 0.42.
    [R1-M3d, R2-M10e, R3-M8]
D4. Geometry: p = 0.2, 0.3 and grazing b (b > 1 − p) for the Table 3 cases; show the p² scaling
    in numbers; add rms. [all three]
D5. Reference-interpolant sensitivity (PCHIP vs linear in r) on an MPS + PHOENIX subset.
    [R1-M7b, R2-m13, R3-M7a]
D6. Constraint-active fraction per law and grid, including the four-term constraint and its
    cost against an unconstrained four-term fit. [R1-M3b, R2-m12, R3-m14]
D7. Uniform-μ merit σ for the r-weighted fits (the limb cost), as a Table 2 column. Add MPS
    Set 1 to Table 2. [R1-M2c, m15]
D8. Pipeline-validation appendix from G1:
    - Kostogryz exact match;
    - μ_cri vs Claret & Southworth 2023 (specify what was compared);
    - uniform-μ four-term and power-2 vs Claret & Southworth 2022/23 (TESS/Kepler) and
      Claret, Hauschildt & Torres 2025 (JWST), in coefficient and light-curve terms.
    [R1-M5b–c, R2-M8, R3-M11]
D9. Chromatic r_c: spread of the per-band edge across the JWST bands, and its ppm-depth
    equivalent. [R1-M5d]
D10. Δw corrections: how often w + Δw leaves the simplex, and what the code does then.
     [R1-M6e, R2-m19, R3-M10]
D11. Demonstrate the Sing rescaling claim: both rescalings, the same models, μ_min = 0.2,
     including their 4-node removal. [R1-M6c, R2-M7b, R3-M7c]
D12. Fraction of the M-dwarf quadratic error that comes from μ < 0.05 (is the gain about the
     steep limb?). [R2-M4]
D13. Photon vs energy integration: size of the effect for a few bands. [R1-M7d]
D14. Constraint and physicality behaviour in narrow spectral bins (coefficient spectra).
     [R1-M7a, R2-M6c]

## E. Larger items: decide

E1. Parameter bias with free and fixed coefficients (noiseless, Table 3 cases: depth, b, a/R*
    bias per law). Asked by all three; not standard Claret. RECOMMEND a compact version
    (one table), because it also answers Sing+.
E2. Priors (recipe in code):
    - Teff/log g/[M/H] uncertainty → w mean and covariance (LDTk-style);
    - a paragraph on analytic marginalisation;
    - principal directions of the weights;
    - the Dirichlet prior implied on h′.
    RECOMMEND the recipe + paragraph; PCA optional. [R2-M5, R3-M6]
E3. r vs area: area is ~20% better in Table 1.
    - Option 1: keep r and justify it (Morello convention; also report area).
    - Option 2: switch the headline to area and regenerate.
    DECISION. [R1-M2e, R2-m11, R3-M5b]
E4. MIRI LRS from MPS-ATLAS (check that the archive spectra reach 12 µm). Add PLATO, Ariel and
    NRS1/NRS2; Pandora optional. RECOMMEND MPS MIRI if available + NRS1/2 + PLATO + Ariel.
    [R1-M7c, R2-M6d–e]
E5. MPS-ATLAS coefficient spectra: provide them at native sampling, or point to
    api.coefficients_custom. RECOMMEND pointing to the API plus a short ODF caveat. [R2-M6b]
E6. Interpolate intensities and refit, as a code option. RECOMMEND (cheap). [R1-M4iii, R2-m15]
E7. Decline with reasons:
    - SPAM optimum (R1-M2d; optional);
    - DEBCat EB tests (R1-M6f);
    - τ_slant = 1 edge (R1-M5e; discuss only);
    - newer near-limb solar CLVs (R3-M9a).
E8. Release logistics:
    - CDS entries for spectra, corrections and offsets, with a ReadMe;
    - Zenodo DOI, licence, versions;
    - total size;
    - optionally release the band-integrated profiles.
