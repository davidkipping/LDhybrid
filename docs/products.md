# LDhybrid v1.0 products

Theoretical limb-darkening coefficients for the rational hybrid laws, with
control laws computed identically. All products live in `data/products/`.

## The laws (frozen, `ldhybrid/v1.py`)

    I(μ)/I(1) = 1 − Σ_j w_j T_j(μ)

| law | T_j | poles ε |
|---|---|---|
| hybrid2 | 1 − μ², Π_ε | 0.36 |
| hybrid4 | 1 − μ⁴, (1 − μ²)², Π_ε0, Π_ε1 | 0.011762, 0.30408 |
| hybrid5 | 1 − μ⁴, (1 − μ²)², Π_ε0, Π_ε1, Π_ε2 | 0.0015966, 0.041275, 0.38743 |

with Π_ε(μ) = [(μ²+ε)⁻² − (1+ε)⁻²] / [ε⁻² − (1+ε)⁻²]. Each w_j is the share of
the centre-to-limb drop carried by its term, and I(0)/I(1) = 1 − Σ w. The
equivalent expanded form, used by the closed-form light curves of
SquishierPlanet and MetalPlanet, is

    I(μ)/I(1) = c0 + c2 μ² + c4 μ⁴ + Σ_k d_k (μ² + ε_k)⁻² .

Physical regions: hybrid2, the exact triangle O–V_c–V_ℓ; hybrid4/5, the
simplex w ≥ 0, Σ w ≤ 1 (sufficient, contains every model tested). The q
columns map these regions to the unit hypercube (uniform priors):
`squishierplanet.laws.hybrid2_from_q` and `simplex_from_q` invert them.

> **Note:** hybrid2 uses ε = 0.36 from v1.0 on, as do SquishierPlanet
> (`HYBRID2_EPS`; the old 0.208 is `HYBRID2_EPS_V0`) and MetalPlanet ≥ v0.12.0.
> Do not use these coefficients with MetalPlanet ≤ v0.11 (ε = 0.208).

Controls: quadratic (u1, u2), Claret 4-term (a1…a4, I/I(1) = 1 − Σ a_k (1 − μ^{k/2})),
power-2 (c, α, I/I(1) = 1 − c(1 − μ^α)).

## Method

* **Band integration:** photon-weighted, I_band(μ) ∝ ∫ I_λ(μ) λ S(λ) dλ, exact
  bin integrals (PHOENIX reduced to R = 3000 photon bins; ODF grids refined 8×).
* **Limb (PHOENIX, spherical):** stellar edge r_c at the node of max |dI/dr|
  (M1), found per band; μ rescaled, μ' = √(1 − (r/r_c)²). Plane-parallel grids:
  no cut (r_c = 1), profile down to the smallest tabulated μ (0.01).
* **Fit:** least squares uniform in projected radius r (2000-pt PCHIP
  resample), I(1) = 1; hybrids constrained to their physical regions,
  quadratic/Claret 4-term to their exact monotone-positive regions; power-2
  unconstrained (Levenberg–Marquardt). No flux conservation imposed.
* **Secondary tables** (`_mu`): uniform in μ, for comparison with Claret's.

## Band tables: `tables/<grid>/<law>_<r|mu>.csv.gz`

| grid | models | coverage | bands |
|---|---|---|---|
| `phoenix` | 1,682 | PHOENIX-ACES v3.0, 2300–7500 K, log g 3.5–5.5, [M/H] −2 … +1 | 38 |
| `mps2` | 16,974 | MPS-ATLAS Set 2, 3500–7500 K, log g 3.5–5.0, [M/H] −2 … +1 | 38 |
| `mps1` | 16,974 | MPS-ATLAS Set 1, same grid | 38 |
| `kurucz` | 507 | ATLAS9 (ExoTiC-LD subset), 3500–6500 K | MIRI LRS only |

Columns:

| column | meaning |
|---|---|
| teff, logg, mh | model parameters |
| band | passband (`ldhybrid.bands`) |
| r_c | adopted stellar edge, in units of the model's outer radius (1 for plane-parallel) |
| w1…wn / u1,u2 / a1…a4 / c,alpha | coefficients |
| q1…qn | unit-hypercube coordinates (hybrids) |
| c0, c2, c4, d1…dK | expanded coefficients (hybrids) |
| rms_I | weighted rms of the profile residual |
| max_abs_I | max \|ΔI/I(1)\| over the tabulated μ range |
| flux_err | relative error in ∫ I μ dμ |
| lc_ppm | worst \|ΔF\| (ppm) against the exact light curve of the tabulated profile, p = 0.1, b ∈ {0, 0.3, 0.6, 0.9} |
| physical | 1 if positive and monotone |
| bound | 1 if a physical constraint is active at the optimum |

`<grid>.h5` holds the same quantities as arrays [model, band(, parameter)].

## Coefficient spectra: `spectra/phoenix/*.npz` (hybrid2/4/5)

Per PHOENIX model, at R = 3000 over 0.1–5.5 µm: `lam_edges` (Å), `I1` (photon
intensity at μ = 1 per bin) and `a_<law>` (bins × weights), with one stellar
edge per model (`r_c`). For any passband with bin-averaged throughput S̄_k,

    w_j = Σ_k S̄_k a_kj / Σ_k S̄_k I1_k

reproduces the band fit exactly (`ldhybrid.spectra.CoeffSpectrum.band_coefficients`).
Bins can be merged exactly (`rebin`). These differ from the band tables by
2–10 ppm in light-curve terms (one r_c per model, native-node quadrature).

## Corrections: `corrections.csv`

Weight shifts Δw to add to a 1D row, with the light-curve change they imply:

* `MURaM`, kind `magnetic`, state ssd/B100G/B200G/B300G: surface fields
  relative to the non-magnetic MURaM run (Kostogryz et al. 2026), nine stars
  F3–M4;
* `MURaM`, kind `3d`: non-magnetic 3D relative to MPS-ATLAS Set 2;
* `Stagger`, kind `3d`: Stagger relative to Kurucz at identical parameters.

MURaM has no intensities below μ = 0.1; its ratio to 1D is held constant
there (flag `extrapolated_below_0.1`). Interpolate in Teff between stars.

## Systematics: `systematics.csv.gz`

At grid points shared by PHOENIX and both MPS-ATLAS sets: the standard
deviation of each hybrid weight across the three grids, and the largest
pairwise light-curve difference. Use as a floor on the prior width of
tabulated coefficients.

## Empirical offsets: `empirical_offsets.csv`

Observed − predicted offsets of hybrid5 predictions against Kepler/TESS
(Maxted 2023: h'1 = I(2/3), h'2 = I(2/3) − I(1/3)) and JWST (Sing et al. 2026:
ℓ, δ of a quadratic fitted for μ ≥ 0.2), per grid, with and without the MURaM
100 G correction, and the star-to-star scatter to use as the prior width.
Recommended for FGK stars: MPS-ATLAS Set 2 + MURaM 100 G (Kepler/TESS
offsets consistent with zero). Details: `studies/p4/results.md`.

## Caveats

* Model atmospheres disagree with each other (PHOENIX vs MPS-ATLAS: ~170 ppm
  median in light curves) far more than any fitting choice (≲ 8 ppm), and
  1D models over-predict limb darkening relative to Kepler/TESS/JWST data.
  Use the tables as priors, not fixed values.
* hybrid5's innermost-pole weight (w3) is weakly determined: limb-definition
  and μ-sampling choices move it by 0.03–0.10 while barely changing the light
  curve.
* At 2300–2500 K the optical PHOENIX profiles brighten toward the limb at
  μ ≈ 0.01–0.02; the physical constraint then binds and lc_ppm reaches
  30–50 ppm (56 of 63,916 PHOENIX hybrid5 rows above 15 ppm).
