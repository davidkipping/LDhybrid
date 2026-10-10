# Gate G1: reproducing published limb-darkening products

Status: **passed** (2026-10-09). Scripts in this directory; JSON outputs
alongside.

## G1(a) Kostogryz et al. (2022), MPS-ATLAS — exact reproduction

`g1_kostogryz.py`. Their published CLVs (J/A+A/666/A60, 24 μ, I/I(1) to 6
decimals) refitted with our engine using their recipe (equal weights on the
native μ, I(1) = 1). All 273,280 rows (2 sets × 34,160 models × 4 bands):

| law | median \|Δ\| | 99th pct | max |
|---|---|---|---|
| 4-term (a1…a4) | 0.5–2.2 × 10⁻⁵ | 1.9–8.2 × 10⁻⁵ | 1.3 × 10⁻⁴ |
| power-2 (c, α) | 3–4 × 10⁻⁷ | 1.2–1.5 × 10⁻⁶ | 7.6 × 10⁻⁶ |

The 4-term floor is the 6-decimal rounding of the published CLVs.
**Catalogue note:** the ReadMe labels both coefficient tables "Set1"; by
reproduction, tablee5 belongs to tablee3 (Set 1) and tablee6 to tablee4
(Set 2).

## G1(b) ExoTiC-LD v3.2 on its own Kurucz and Stagger files

`g1_exoticld.py`. 4 Kurucz + 4 Stagger models × 7 modes.

* **Fitting** on the same profile: ≤ 1.2 × 10⁻⁵ in the 4-term coefficients
  (LM tolerance).
* **Band integration** differs by up to 1 × 10⁻² in I/I(1). Two causes, both
  isolated:
  1. *Different curves* for TESS and Spitzer (ours from SVO, theirs bundled).
  2. *ExoTiC-LD's quadrature*: 500-node Gauss–Legendre over the whole mode.
     Our emulation of it reproduces ExoTiC-LD to 10⁻¹⁶, and it is
     unconverged on high-resolution spectra: GL500 vs GL20000 differ by up to
     1.7 × 10⁻³ (Stagger, STIS G430L). Our exact bin integration agrees with
     the exact product integral of the linear interpolants to 4 × 10⁻⁶
     (3 × 10⁻⁸ with 8× bin refinement, now the default for coarse grids).

## G1(c) Claret & Southworth (2023), PHOENIX power-2 — consistency

`g1_claret.py`. Their PHOENIX-COND intensities are not public, so this
compares our PHOENIX-ACES v3.0 run of their recipe (M1/M2 cut, 900 uniform μ,
unweighted LM) with their table, for 13 dwarfs × Kepler/TESS/CHEOPS.

* μ_cri (M1) agrees to ~10⁻³ (e.g. 0.0290/0.0286, 0.0532/0.0532), so the
  limb-cut definition matches.
* **M2 direction:** their tabulated μ_cri(M2) > μ_cri(M1), i.e. the cut moves
  *inward* of the max-|dI/dr| node. Our first implementation moved it
  outward; fixed.
* Power-2: median |Δg|, |Δh| = 0.016 (M1), 0.014 (M2); max 0.04. This is the
  size of the COND-vs-ACES difference, not pipeline error.

## Also verified

* **PHOENIX reduction** (R = 3000 photon bins, 0.1–5.5 µm) against the raw
  1 Å files, real throughputs: ≤ 1.2 × 10⁻⁵ in I/I(1) (≈ 0.1 ppm in light
  curves for p = 0.1).
* **MPS-ATLAS band integration** vs Kostogryz's published CLVs from the same
  spectra: 1–7 × 10⁻⁴ (Kepler, TESS, CHEOPS) at solar-like metallicity.
* **Light-curve metric** vs squishierplanet's exact hybrid5 solution:
  ≤ 1 × 10⁻⁸ ppm.

## Issues found

1. **SVO detector flags.** SVO flags the TESS (and Spitzer/IRAC) curves as
   energy counters, but they are photon-counting responses. Dividing by λ
   shifted TESS profiles by ~1% against Kostogryz. Detector types are now set
   per band (`bands.ENERGY_TYPE`: Bessell Johnson–Cousins only).
2. **Spitzer/IRAC curves:** SVO's and ExoTiC-LD's disagree by factors of
   0.5–1.7 across the band. Resolved by using the official SSC response files
   (electrons per photon). The *subarray* response, used for most transit
   photometry, differs from the full-array response by factors of 0.4–2.4
   in ch1 (0.93–1.19 in ch2), and SVO's curves follow the full array.
   `Spitzer_IRAC1/2` are now the subarray curves; `_full` variants are kept.
3. **Published MPS-ATLAS Set 1 CLVs are wrong at [M/H] = −1.7 … −0.75.**
   Our integration of the Edmond Set 1 spectra matches the CDS CLVs
   (tablee3) to ~10⁻⁴ at every other metallicity, and Set 2 (tablee4) at
   *all* metallicities. In this block they differ by 0.5–2 × 10⁻².
   The archive is the smooth one: at 6000 K, log g 4.5, CHEOPS, I(0.1)/I(1)
   falls monotonically from 0.4761 ([M/H] = −2.0) to 0.4589 (−0.5) in the
   archive, whereas the CDS table stays at 0.4735–0.4762 from −1.7 to −0.75
   and then jumps to 0.4647 at −0.7. The published Set 1 coefficients
   (tablee5) inherit this, being consistent with tablee3 (G1(a)). The
   archived spectra, which we use, are unaffected. Worth reporting to the
   authors (N. Kostogryz).
