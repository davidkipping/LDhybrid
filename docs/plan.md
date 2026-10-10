# LDhybrid — plan of work

Status (2026-10-10): **G0–G4 passed** (P4: `studies/p4/results.md`) (`docs/data.md`, `studies/g1/results.md`,
`studies/p2/results.md`, `studies/g3/results.md`); v1.0 frozen in
`ldhybrid/v1.py`; products in `data/products/` (`docs/products.md`);
next P5.
Deferred from P3: Jacobians ∂w/∂(Teff, log g, [M/H]) (users can difference the
tables) and PLATO/Ariel/Pandora passbands (curves not yet sourced). Goal: Claret-style theoretical tables of
limb-darkening coefficients for the rational **hybrid2, hybrid4 and hybrid5**
laws of SquishierPlanet (`../squishierplanet/docs/notes`,
`squishierplanet/laws.py`), regressed against simulated specific-intensity
profiles, with the physical constraints imposed by construction.

## The laws (inherited, frozen at G2)

    I(μ)/I(1) = 1 − Σ_j w_j T_j(μ),   T_j: 1 at the limb, 0 at the centre

| law | shapes T_j | poles ε | free |
|---|---|---|---|
| hybrid2 | 1 − μ², Π_ε | **0.36** (re-tuned in P2c; was 0.208, J band) | 2 |
| hybrid4 | 1 − μ⁴, (1 − μ²)², Π_ε0, Π_ε1 | analytic ladder K=2: 0.01176, 0.3041 | 4 |
| hybrid5 | 1 − μ⁴, (1 − μ²)², Π_ε0..2 | analytic ladder K=3: 0.001597, 0.04128, 0.3874 | 5 |

Poles are universal constants, never fitted per star: this keeps the laws
linear in w, so every fit is a convex problem and rows are comparable.

## Decisions taken

*G2 (2026-10-10), from the P2 evidence: hybrid2 ε = 0.36; headline fits
r-weighted (2000-pt PCHIP resample), simplex / exact triangle, no FCM, with a
uniform-μ secondary table; PHOENIX limb M1 per band for band tables and one
r_c per model for coefficient spectra; no custom radiative transfer; MIRI LRS
for FGK only (Kurucz), M dwarfs dropped from v1.0.*

1. **Weighting.** Fit the *intensity profile*, weighted by radius or disc area
   (ExoTETHyS "weighted-r"; area weighting is the uniform-b, small-planet
   limit of the light-curve error). **No** fitting to synthetic light curves:
   Howarth's SPAM coefficients depend on the adopted impact parameter. A
   uniform-μ table is kept as a secondary product for comparison with Claret.
2. **hybrid2's pole** is re-tuned once over the full grid and all bands, then
   frozen.
3. **Custom radiative transfer** (Korg/Turbospectrum on MARCS) is a P2
   go/no-go, decided by the μ-sampling test (P2d).
4. **Scope.** FGKM dwarfs and subgiants only (no hot stars, no giants), in the
   standard passbands used for transit observations, with particular
   attention to **HST WFC3** and **JWST**.
5. **Both hybrid4 and hybrid5** are tabulated. hybrid5's robustness at the
   extreme limb is measured in P2d; if it fails, its rows carry a
   systematic-error flag, not removal.

Proposed grid scope (confirm): Teff 2300–7500 K, log g 3.5–5.5,
[M/H] −2 to +1.

## Phases and gates

### P0 — Set-up and data

* Package skeleton `ldhybrid/` (pyproject, tests, `.venv`); depend on
  `squishierplanet` (laws, oracle) and `metalplanet ≥ 0.10.7` (fast
  circular-planet hybrid light curves) as editable installs. `data/` is
  git-ignored.
* **Intensity grids** (download scripts, checksums, provenance log):

  | grid | role | source | approx. volume |
  |---|---|---|---|
  | PHOENIX-ACES SpecInt **v3.0** (Husser+2013; 78 μ, 1 Å, ≤ 5.5 µm, spherical) | baseline, all Teff, JWST/HST spectroscopy | ftp.gwdg.de/pub/misc/phoenix/v3.0/SpecIntFITS/ | ~60 GB for the scope above |
  | MPS-ATLAS Set 1 + Set 2 (Kostogryz+2022; 24 μ, plane-parallel) | baseline FGK broadband; convection/abundance systematic | doi:10.17617/3.NJ56TR (manual download), or ExoTiC-LD bundle | ~17 GB |
  | MURaM 3D-MHD (Kostogryz+2026; HD, SSD, 100–300 G; 10 μ ≥ 0.1) | magnetic/3D correction layer | doi:10.17617/3.FBTIYY | small |
  | Stagger 3D (Magic+2015; 10 μ, via ExoTiC-LD) | non-magnetic 3D correction | Zenodo 10.5281/zenodo.7436954 | 0.5 GB |
  | PLATO MARCS profiles (Morello+2022) | independent test bed | Zenodo 10.5281/zenodo.7339706 | small |
  | Kurucz ATLAS9 `ip*.pck` (17 μ, to 160 µm) | MIRI fallback only | kurucz.harvard.edu/grids/ | small |

* Ask the Hamburg group for PHOENIX NewEra centre-to-limb data (not public).
* **Throughputs:** SVO Filter Profile Service (photometric bands), stsynphot
  (HST WFC3 G102/G141/UVIS G280, STIS G430L/G750L), JWST throughputs (pandeia
  or the ExoTiC-LD bundle), Ariel and Pandora where published.
* **Reference coefficient tables** from VizieR: Claret 2017 (600/A30),
  Claret & Southworth 2022/2023 (664/A128, 674/A63), Claret+2025 (699/A97),
  Kostogryz+2022 (666/A60).

**Gate G0:** every source downloaded and loadable, provenance recorded.

### P1 — Pipeline, and reproduce the published tables

* **Loaders** → one `Profile` object: μ nodes, I(λ, μ), stellar parameters,
  geometry flag (spherical/plane-parallel), provenance.
* **Limb treatment** for spherical models, both variants of Claret & Southworth
  2023: M1 (radius at max |dI/dr|, rescale μ) and M2 (midway to the next node
  inward, as Claret & Southworth's tabulated μ_cri show). Record r_c / μ_cri on every row.
  * Subtlety: the cut point depends on λ. Passband integration stays linear
    only if one μ-mapping is used for all λ in a model. Default: one r_c per
    model, from a reference band; compare with per-band r_c in P2d.
* **Band integration:** photon-weighted, I_band(μ) ∝ ∫ I(λ, μ) λ S(λ) dλ.
* **Resampling:** monotone interpolation (PCHIP) in r onto a dense quadrature
  grid, or native nodes with trapezoid weights; both are kept to test in P2d.
* **Fitting engine:** weighted least squares in the shape basis, with I(1) = 1.
  Pluggable weights (uniform-μ, r, area) and constraints, solved as a convex QP
  (cvxpy/OSQP; NNLS where enough):
  * the simplex w ≥ 0, Σw ≤ 1;
  * the exact physical region (I(0) ≥ 0, ∂I/∂ξ ≥ 0 on a dense ξ grid; the
    exact hybrid2 triangle);
  * optional flux conservation, which is a linear equality for these laws.
* **Control laws** through the same pipeline: quadratic, power-2 and Claret
  4-term (Levenberg–Marquardt).
* **Metrics** on every row:
  * weighted intensity rms in a continuous norm (independent of grid N);
  * max |ΔI|;
  * transit light-curve error in ppm, max and rms over
    p ∈ {0.03, 0.1, 0.15} × b ∈ [0, 0.9], from the exact light curve of the
    tabulated profile (a 1D radial integral for a circular planet,
    cross-checked against the squishierplanet oracle);
  * which constraints bind.

**Gate G1:**
* Reproduce the Kostogryz+2022 MPS-ATLAS power-2 and 4-term coefficients
  (same public profiles) to ≲ 1e-3.
* Reproduce ExoTiC-LD and ExoTETHyS weighted-r output on PHOENIX.
* Agree with the Claret PHOENIX tables to the level expected from different
  model versions; his PHOENIX-COND intensities are not public, so exact
  agreement is impossible.

### P2 — Method study and law freeze

Work on a representative subset: ~30 stars × ~8 bands that span Teff, log g,
[M/H] and wavelength (TESS, Kepler, WFC3 G141, NIRISS SOSS, NIRSpec G395H,
MIRI LRS, …).

* **P2a Weighting.** Uniform-μ vs r vs area. Score on light-curve error
  averaged over, and worst over, the (p, b) grid; no single b is privileged.
  Choose the headline weighting.
* **P2b Constraints.** Unconstrained vs simplex vs exact region, ± flux
  conservation. Measure the accuracy cost, how often each constraint binds,
  and the fraction of the grid inside the simplex (so far 92/92 dwarfs).
* **P2c hybrid2 pole.** Re-tune ε over the full grid × all bands, minimising
  the worst and median light-curve error under the headline weighting. Check
  that the hybrid4/5 ladder with α_eff = 1 holds across bands (verify only,
  no re-tuning).
* **P2d Limb sensitivity**, the key study:
  * M1 vs M2;
  * one r_c per model vs per band;
  * spherical vs plane-parallel (PHOENIX vs MPS-ATLAS at matched parameters);
  * μ subsampling of PHOENIX 78 → 24/17/10 nodes;
  * interpolation scheme.

  Report Δw per weight against the spread between grids.
  * hybrid5 verdict: are its weights, especially the innermost pole's, stable
    to well within that spread?
  * Custom-transfer go/no-go: **go** if μ-sampling-induced Δw is comparable to
    the spread between grids, or if MIRI coverage beyond Kurucz is needed.
* **P2e (if go)** Korg on MARCS with 50–100 μ clustered at the limb, at native
  resolution, for the subset first; repeat P2d on it.
* **P2f Conditioning.** Fisher analysis of w under the headline weighting
  (which combinations a profile or transit constrains); define the
  near-degeneracy flags.

**Gate G2:** law definitions frozen as v1.0 (constants in the package, never
changed afterwards); method frozen and written up as an internal note.

### P3 — Production

* **Tables:** for each law (hybrid2/4/5 + controls) × grid (PHOENIX M1 and M2,
  MPS-ATLAS Set 1 and Set 2, custom if go) × passband × method (headline
  weighting; uniform-μ secondary).
* **Passbands:**
  * Kepler/K2, TESS, CHEOPS, PLATO, Spitzer IRAC 1–2;
  * Sloan g′r′i′z′, Johnson-Cousins BVRI, 2MASS JHKs;
  * HST WFC3 G102/G141/UVIS G280, STIS G430L/G750L;
  * JWST NIRISS SOSS (orders 1–2), NIRSpec PRISM/G140H/G235H/G395H, NIRCam
    F322W2/F444W, MIRI LRS;
  * Ariel, Pandora.
* **Coefficient spectra.** Unnormalised c(λ) in fine bins over 0.2–5.5 µm
  (plus MIRI from the fallback), so users can integrate over any spectroscopic
  bin.
  * Closure test: integrated vs direct band fit.
  * Document the gap when inequality constraints bind; the result is still
    physical, because the cone is convex.
* **3D / magnetic correction layer.**
  * Δw = fit(I₁D · C) − fit(I₁D), with C = I₃D/I₁D smoothed by a low-order
    function in μ.
  * MURaM pairs with MPS-ATLAS (shared microphysics), tabulated against Teff,
    ⟨B⟩ ∈ {HD, SSD, 100, 200, 300 G} and band. Stagger is the non-magnetic
    alternative.
  * MURaM has no μ < 0.1 data, so extrapolation below 0.1 is an explicit,
    flagged assumption.
* **Uncertainties:**
  * systematic σ_w from the spread between grids (PHOENIX, MPS1, MPS2);
  * Jacobians ∂w/∂(Teff, log g, [M/H]) for propagating stellar-parameter
    errors.
* **Sampler coordinates:** q-space values on every row (the hybrid2 triangle
  map; stick-breaking for hybrid4/5).

**Gate G3:**
* every row physical, or flagged;
* held-out-node interpolation error below the systematic σ_w;
* coefficient-spectrum closure test passed.

### P4 — Validation against observations

* **Sun:** fit Neckel & Labs (1994) centre-to-limb data and compare with the
  tabulated solar-parameter rows.
* **Kepler/TESS:** the Maxted (2023) sample, compared through law-agnostic
  profile values h₁′ = I(2/3), h₂′ = h₁′ − I(1/3); the Patel & Espinoza (2022)
  TESS sample.
* **JWST:** the Sing et al. (2026) (ℓ, δ) constraints for 7 FGK hosts.
* **Reanalysis:** one or two benchmark planets with HST WFC3 G141 plus JWST
  (e.g. WASP-39b, HD 209458b). Fit hybrid4/5 freely with simplex priors vs
  fixed or table-prior coefficients, and quantify the bias in the
  transmission spectrum.
* Deliver recommended empirical offsets / prior widths per grid, after
  Maxted (2023).

**Gate G4:** validation note written; offsets fixed.

### P5 — Release

* **Paper** (A&A, Claret-style): method, the convexity and linearity results,
  limb sensitivity, 3D/magnetic layer, validation.
* **CDS tables** in long format, one row per (model, band, method, limb
  variant):
  * stellar parameters and band;
  * w_j, plus the expanded c_j;
  * q-space values;
  * r_c / μ_cri;
  * intensity rms and LC ppm;
  * σ_w;
  * flags.

  Publish the band-integrated profiles too, as Kostogryz+2022 and PLATO do.
* **`ldhybrid` package:**
  * on-the-fly coefficients for any user throughput;
  * table interpolation;
  * coefficient-spectrum integration;
  * priors in q-space;
  * hand-off to squishierplanet and MetalPlanet.
* **Upstream:** offer the hybrid laws to ExoTiC-LD.

## Risks

| risk | mitigation |
|---|---|
| hybrid5's inner pole sits at μ ≈ 0.02–0.04, where models disagree and grids are sparse | P2d measures it; custom transfer if needed; flag rather than drop |
| A λ-dependent limb radius breaks linear passband integration | one r_c per model by default; quantify in P2d |
| Weak constraint on some weights under r/area weighting (near-degenerate fits) | P2f Fisher analysis; convex constraints keep fits physical; flag rows |
| PHOENIX grid holes (seen in squishierplanet) | log substitutions; interpolate only within a grid, never across grids |
| 3D grids are sparse in parameters and μ | used only as a smoothed correction layer, with stated assumptions |
| 1D (and non-magnetic 3D) models over-predict limb darkening | P4 empirical offsets; magnetic layer |

## To verify before citing

* Kostogryz+2024 Nat. Astron. volume and page.
* Kipping 2016 (MNRAS 455, 1680?).
* Howarth 2011b, Parviainen & Aigrain 2015, Patel & Espinoza 2022 volume and
  page.
* The PHOENIX SpecInt v3.0 Teff ceiling.
* MPS-ATLAS Set 2 mixing-length calibration.
* Whether ExoTiC-LD reads MURaM natively.
