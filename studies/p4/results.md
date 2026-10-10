# P4: validation against observations

Run 2026-10-10. Data: `data/raw/validation/*/README.md` (sources, definitions,
caveats). Scripts and per-study tables in this directory.

## P4a — the Sun (Neckel & Labs 1994; `p4a_sun.md`)

**The laws can represent the real Sun.** Fitted to the observed profiles at
30 wavelengths, worst light-curve error (ppm, production grid), median / max:

| range | quadratic | power-2 | hybrid2 | Claret 4-term | hybrid4 | hybrid5 |
|---|---|---|---|---|---|---|
| 303–450 nm | 21 / 56 | 24 / 56 | 61 / 110 | 4.3 / 9.5 | 8.1 / 9.9 | 4.8 / 7.5 |
| 450–700 nm | 46 / 50 | 6.7 / 11 | 37 / 46 | 4.7 / 6.6 | 10 / 11 | 3.0 / 6.1 |
| 700–1100 nm | 47 / 49 | 8.9 / 12 | 35 / 38 | 5.8 / 6.2 | 8.9 / 9.8 | 3.3 / 3.7 |

**How well do the models predict it?** Light-curve difference (ppm) between
model and observed profiles, median / max (median max |ΔI/I(1)|, μ ≥ 0.1):

| range | PHOENIX v3.0 | MPS-ATLAS 1 | MPS-ATLAS 2 | MURaM hydro | MURaM ssd | MURaM 100 G |
|---|---|---|---|---|---|---|
| 303–450 nm | 307 / 675 (0.033) | 104 / 201 (0.012) | 151 / 267 (0.014) | 197 / 278 (0.020) | 149 / 242 (0.017) | 160 / 435 (0.025) |
| 450–700 nm | 240 / 341 (0.030) | **25** / 42 (0.005) | 43 / 70 (0.008) | 72 / 179 (0.010) | 52 / 141 (0.007) | 83 / 97 (0.017) |
| 700–1100 nm | 156 / 210 (0.023) | 44 / 62 (0.009) | 53 / 65 (0.010) | 33 / 42 (0.007) | **23** / 48 (0.005) | 75 / 101 (0.016) |

MPS-ATLAS (optical) and MURaM with a small-scale dynamo (near-IR) reproduce
the quiet Sun best; PHOENIX v3.0 is worst throughout; 100 G facular fields are
too strong for the quiet Sun, as expected. The blue is worst for every model,
partly because Neckel & Labs measured line-free windows while the model
bands (1-nm top-hats) include lines.

## P4b — Kepler and TESS hosts (Maxted 2023; `p4b_maxted.md`)

Observed − predicted h'1 = I(2/3), h'2 = I(2/3) − I(1/3), b < 0.65, weighted
mean ± error (external scatter):

| band | grid | law | correction | Δh'1 | Δh'2 |
|---|---|---|---|---|---|
| Kepler (24) | MPS-ATLAS 2 | Claret 4-term | — | +0.0080 ± 0.0016 | −0.0131 ± 0.0033 |
| | MPS-ATLAS 2 | hybrid5 | — | +0.0080 ± 0.0017 | −0.0132 ± 0.0033 |
| | MPS-ATLAS 2 | hybrid5 | MURaM 100 G | **−0.0023 ± 0.0015** | **−0.0039 ± 0.0033** |
| | MPS-ATLAS 1 | hybrid5 | MURaM 100 G | −0.0054 ± 0.0016 | −0.0033 ± 0.0033 |
| | PHOENIX v3.0 | hybrid5 | — | +0.0231 ± 0.0022 | −0.0210 ± 0.0038 |
| | PHOENIX v3.0 | hybrid5 | MURaM 200 G | +0.0022 ± 0.0015 | −0.0016 ± 0.0034 |
| TESS (10) | MPS-ATLAS 2 | hybrid5 | — | +0.0050 ± 0.0033 | −0.0089 ± 0.0045 |
| | MPS-ATLAS 2 | hybrid5 | MURaM 100 G | **−0.0039 ± 0.0035** | **−0.0003 ± 0.0043** |

* **Pipeline check:** our Claret 4-term predictions from the MPS-ATLAS spectra
  reproduce Maxted's published MPS-ATLAS offsets (Kepler Set 1/2: +0.006/−0.012,
  +0.009/−0.013; TESS +0.004/−0.009, +0.006/−0.010) within 1σ in all eight
  numbers.
* hybrid5 inherits exactly the offsets of the 4-term law: the offset is a
  property of the model atmospheres, not of the law.
* **Adding the MURaM 100 G correction to MPS-ATLAS Set 2 removes the offset**
  in both bands (and lowers the external scatter, Kepler 0.0058 → 0.0052),
  as Kostogryz et al. (2024) found with their own laws. PHOENIX-ACES v3.0
  starts further off (Δh'1 = +0.023, vs Maxted's +0.010 for PHOENIX-COND).

## P4c — TESS hosts (Patel & Espinoza 2022; `p4c_patel.md`)

SPAM predictions (quadratic fitted to the hybrid5 light curve at each
system's k, b), predicted − empirical:

| grid | correction | N | Δu1 | Δu2 |
|---|---|---|---|---|
| PHOENIX v3.0 | — | 173 | +0.031 ± 0.011 | +0.040 ± 0.013 |
| PHOENIX v3.0 | 100 G | 173 | +0.004 ± 0.011 | +0.049 ± 0.013 |
| MPS-ATLAS 2 | — | 167 | −0.023 ± 0.011 | +0.083 ± 0.013 |
| MPS-ATLAS 2 | 100 G | 167 | −0.051 ± 0.011 | +0.092 ± 0.013 |

* **Pipeline check:** our PHOENIX SPAM values agree with their C17 PHOENIX
  r-method SPAM values to a median 0.005 (u1) and 0.010 (u2).
* The TESS quadratic coefficients are individually noisy (σ ≈ 0.1, external
  scatter 0.06 in u1) and strongly correlated; this sample is far less
  constraining than P4b. The large M-dwarf offsets (2500–4000 K: Δu1 = −0.16,
  Δu2 = +0.22) echo Patel & Espinoza's finding that models fail below ~5000 K.

## P4d — JWST (Sing et al. 2026; `p4d_sing.md`)

1,110 per-bin (ℓ, δ) measurements for six FGK hosts (NIRISS/SOSS,
NIRSpec/PRISM; WASP-12 has no figure), extracted from the paper's vector
figures (exact to ~10⁻⁶; reproduces the paper's stated WASP-6 offsets).
Model values: quadratic fitted for μ ≥ 0.2, ℓ = I(0), δ = c2/4 (their
definitions). Observed − model, per-star weighted, then mean ± star-to-star
s.d.:

| grid | route | Δℓ | Δδ |
|---|---|---|---|
| MPS-ATLAS 2 | their recipe (native nodes) | +0.043 ± 0.031 | −0.011 ± 0.017 |
| MPS-ATLAS 2 | via hybrid5 | +0.035 ± 0.031 | −0.007 ± 0.018 |
| MPS-ATLAS 2 | via hybrid5 + MURaM 100 G | **+0.024 ± 0.030** | −0.007 ± 0.018 |
| MPS-ATLAS 1 | via hybrid5 | +0.033 ± 0.033 | −0.008 ± 0.018 |
| PHOENIX v3.0 (r-method) | via hybrid5 / direct | +0.041 ± 0.044 | −0.005 ± 0.017 |
| *Sing et al. Table 3 (7 stars)* | MPS-ATLAS 1 / 2 / PHOENIX | +0.034 / +0.036 / 0.000 | −0.007 / −0.006 / +0.004 |

* **Pipeline check:** MPS-ATLAS offsets agree with Sing et al.'s within the
  star-to-star scatter (6 of their 7 stars). The hybrid5 route equals direct
  profile fitting for PHOENIX (+0.041 both); for MPS-ATLAS the 0.008 shift
  comes from sampling the μ ≥ 0.2 quadratic fit on native nodes vs uniformly,
  not from the law.
* **PHOENIX:** Sing et al.'s Δℓ ≈ 0 reflects their rescaling (cut four nodes
  inward of max |dI/dμ|, then μ' = (μ − μ_cri)/(1 − μ_cri), linear in μ).
  Applying it to our models moves Δℓ by −0.03 per star (WASP-6 +0.083 →
  +0.055, WASP-121 +0.077 → +0.047, WASP-80 +0.001 → −0.032), most of the gap.
  With the geometric r-method, PHOENIX over-predicts limb darkening in JWST
  bands as much as MPS-ATLAS does.
* **Magnetic fields explain about a third of the JWST offset** (MPS-ATLAS 2:
  +0.035 → +0.024), consistent with the effect weakening into the infrared;
  the remainder is within the star-to-star scatter (0.03) but systematic in
  sign for five of six stars.

## G4 — conclusions and recommended offsets

Status: **passed** (2026-10-10). Offsets fixed in
`data/products/empirical_offsets.csv`.

1. **The laws are not the limitation.** hybrid5 represents the observed Sun
   to 3–5 ppm and inherits exactly the model offsets of the 4-term law in
   every comparison; the offsets are properties of the model atmospheres.
2. **Recommended prior model for FGK stars:** MPS-ATLAS Set 2 + MURaM 100 G
   correction (`corrections.csv`). It removes the Kepler/TESS offsets
   (Δh'1 = −0.002 ± 0.002, Δh'2 = −0.004 ± 0.003) and a third of the JWST
   offset.
3. **Prior widths:** star-to-star scatter after correction, Kepler
   σ(h'1) = 0.005, σ(h'2) = 0.013; JWST σ(ℓ) = 0.03, σ(δ) = 0.018 (plus a
   residual JWST Δℓ = +0.024 for MPS-ATLAS 2 + 100 G). These exceed the
   coefficient interpolation error and the grid-to-grid systematic spread for
   the same quantities, so they set the prior width.
4. **PHOENIX v3.0** (r-method) shows the largest offsets everywhere (Kepler
   Δh'1 = +0.023; JWST Δℓ = +0.041; worst solar match) and should be used with
   these offsets or for M dwarfs only, where it is the sole grid; there the
   empirical offsets are large and poorly measured (Patel & Espinoza:
   Δu1 = −0.16, Δu2 = +0.22 below 4000 K), so wide priors are advised.
5. **Not done:** a benchmark reanalysis of HST WFC3 + JWST light curves with
   free vs table-prior hybrid coefficients (planned in P4) — deferred.
