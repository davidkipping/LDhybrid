# Using the LDhybrid v1.0 coefficients

There are three ways to get coefficients, from lightest to fullest:

| Route | Needs | Covers |
|---|---|---|
| [Standalone files](#1-standalone-single-file-interpolators) (`standalone/*.py`) | numpy, scipy | hybrid4; PHOENIX or MPS-ATLAS Set 2; TESS or Kepler |
| [Table interpolator](#2-the-table-interpolator-ldhybridinterpolatepy) (`ldhybrid/interpolate.py`) | numpy, scipy (+ h5py for refits); the CDS tables | all six laws, three grids, 40–41 passbands; priors, offsets, corrections, light curves |
| [Full package](#3-the-full-package-ldhybridapi) (`ldhybrid.api`) | the whole repository and its `data/` products | everything above, plus custom passbands |

All routes give the same numbers: the standalone files and the table interpolator
agree to 1e-8 in the weights.

## The laws

Every law is written as I(μ)/I(1) = 1 − Σ_j w_j T_j(μ), where μ is the cosine of
the angle from disc centre (μ² = 1 − r²).

| Law | Coefficients | Shapes T_j |
|---|---|---|
| `hybrid2` | w1, w2 | 1 − μ², Π(μ; 0.36) |
| `hybrid4` | w1–w4 (and q1–q4) | 1 − μ², (1 − μ²)², S3(μ), S4(μ): the curvature-simplex basis, below |
| `hybrid5` | w1–w5 | 1 − μ⁴, (1 − μ²)², Π(μ; 0.0015966), Π(μ; 0.041275), Π(μ; 0.38743) |
| `quadratic` | u1, u2 | 1 − μ, (1 − μ)² |
| `claret4` (four-term) | a1–a4 | 1 − μ^(k/2), k = 1–4 |
| `power2` | c, alpha | I/I(1) = 1 − c(1 − μ^alpha) |

with the pole term Π(μ; ε) = [(μ² + ε)⁻² − (1 + ε)⁻²] / [ε⁻² − (1 + ε)⁻²].

**hybrid4's basis (since 2026-10-14).** S3 and S4 are fixed combinations of the
two poles Π(μ; 0.028601), Π(μ; 0.38582) and the polynomial shapes:

    S3 = 1.020963 Π(0.028601) − 0.022993 Π(0.38582) + 0.001252 (1 − μ²) + 0.000778 (1 − μ²)²
    S4 = 1.338604 Π(0.38582) − 0.000591 Π(0.028601) − 0.162320 (1 − μ²) − 0.175693 (1 − μ²)²

(full precision: `H4_BASIS` in `ldhybrid/interpolate.py`, or Note 15 of the CDS
ReadMe). They are chosen so that the simplex w ≥ 0, Σw ≤ 1 is *exactly* the set
of hybrid4 profiles that are positive, monotonic (Kipping 2013's conditions),
concave in μ², and with a convex slope in μ². Every weight vector in the simplex
is such a profile, and every such profile is in the simplex. The curvature
conditions are an assumption about the profile, not a physical requirement:
they keep 9% of the positive, monotonic hybrid4 profiles. All tabulated fits
lie in the region, but some lie on its boundary, where a curvature condition
binds: 0.19% of the PHOENIX r-weighted fits and 18.6% of the PHOENIX μ-weighted
fits (mostly at 3–5 μm); never for MPS-ATLAS.

**Sampling coordinates q1–q4** (the analogue of Kipping's q1, q2): with
s = Σw = 1 − I(0) = q1^(1/4),

    w1 = s (1 − q2^(1/3))           w2 = s q2^(1/3) (1 − q3^(1/2))
    w3 = s q2^(1/3) q3^(1/2) (1 − q4)   w4 = s q2^(1/3) q3^(1/2) q4

A uniform q in [0, 1]⁴ is a uniform prior over that region. It is not
uninformative about the profile: it favours strong limb darkening (the density
of s grows as s³), so in the infrared, where tabulated I(0) reaches 0.7–0.9, a
prior centred on the tables is the better default. The tables give q
next to w; `from_q` and `to_q` in `ldhybrid/interpolate.py` (and in the
standalone files) convert.

**Using the coefficients with other codes.** SquishierPlanet and Elómë take
hybrid4 weights on the monomial-pole shapes (1 − μ⁴, (1 − μ²)², Π(0.028601),
Π(0.38582)). Convert with `w @ H4_BASIS` (or use the expanded-form columns
c0, c2, c4, d1, d2 of the tables, which are basis-free). hybrid4 uses the pole
ladder with α = 25/16 (since 2026-10-11): do not combine its coefficients with
code that still uses the α = 1 poles (0.01176, 0.3041).

**Which fit.** The headline coefficients (`method="r"`) are least-squares fits
uniform in projected radius: the right choice for transit light curves. The
`method="mu"` fits are uniform in μ and weight the limb more; use them for
eclipsing binaries, interferometry or microlensing. With `mu`, hybrid4 is
less accurate on PHOENIX, and at 3–5 μm most of its PHOENIX `mu` fits are set
by a curvature condition; prefer hybrid5 or the four-term law there.

**Which grid.** For FGK stars the paper recommends MPS-ATLAS Set 2 with the
MURaM 100 G correction (`grid="mps2", correction="B100G"`). Quiet Sun-like stars
are closer to the small-scale-dynamo correction (`"ssd"`). Below 3500 K only
PHOENIX exists, and no grid has been validated against M-dwarf observations.

## 1. Standalone single-file interpolators

`standalone/` holds one file per grid and passband for hybrid4, each with its
coefficient grid embedded (58–460 kB). Copy the file anywhere; nothing else is
needed.

```python
from ldhybrid_hybrid4_mps2_tess import coefficients, intensity, is_physical

r = coefficients(5600, 4.4, 0.0)                       # Teff [K], log g [cgs], [M/H]
r["w"]          # (w1, w2, w3, w4), curvature-simplex basis
r["q"]          # (q1, q2, q3, q4), sampling coordinates
r["physical"]   # positive and monotonic?
r["lc_ppm"]     # light-curve error of the tabulated fit (p = 0.1)

r = coefficients(5600, 4.4, 0.0, correction="B100G")   # MPS-ATLAS files only
I = intensity(r["w"], mu)                              # I(mu)/I(1)
r = coefficients(9000, 4.5, 0.0, extrapolate=True)     # nearest model outside the grid
```

From the shell: `python ldhybrid_hybrid4_mps2_tess.py 5600 4.4 0.0`.

| File | Grid | Passband |
|---|---|---|
| `ldhybrid_hybrid4_phoenix_tess.py` | PHOENIX (2300–7400 K, log g 3.5–5.5) | TESS |
| `ldhybrid_hybrid4_phoenix_kepler.py` | PHOENIX | Kepler |
| `ldhybrid_hybrid4_mps2_tess.py` | MPS-ATLAS Set 2 (3500–7500 K, log g 3.5–5.0) | TESS |
| `ldhybrid_hybrid4_mps2_kepler.py` | MPS-ATLAS Set 2 | Kepler |

Other grid, passband or law combinations can be generated with
`scripts/make_standalone.py`.

## 2. The table interpolator (`ldhybrid/interpolate.py`)

A single self-contained file (numpy and scipy only) that reads the CDS tables.
It can be copied out of the repository and used on its own.

**Getting the tables.** Until publication they are in the private GitHub draft
release `v1.0-cds-draft` and in the draft Zenodo record. Put all files in one
directory: the 18 `<grid>_<law>.dat.gz` tables, plus `corrections.dat` and
`offsets.dat` for corrections and offsets. `SHA256SUMS` checks the download:

```bash
gh release download v1.0-cds-draft --dir cds
cd cds && shasum -a 256 -c SHA256SUMS
```

On first use each table is parsed and a cache (`<grid>_<law>.npz`) is written
next to it (`LDTables(..., cache=False)` to disable).

### Coefficients

```python
from ldhybrid.interpolate import LDTables

t = LDTables("path/to/cds")
r = t.coefficients("hybrid5", teff=5600, logg=4.4, mh=0.1, band="TESS", grid="mps2")
r["coeffs"]      # {"w1": ..., ..., "w5": ...}
r["q"]           # sampling coordinates (hybrid laws)
r["physical"]    # positive and monotonic?
r["mucri"]       # mu at the stellar edge (spherical models; 0 for MPS-ATLAS)
r["lc_ppm"]      # light-curve error of the tabulated fit (p = 0.1)
r["inside"]      # False if extrapolate="nearest" was used outside the grid

t.coefficients("hybrid4", 5600, 4.4, 0.1, "TESS", grid="mps2", correction="B100G")
t.coefficients("claret4", 5600, 4.4, 0.1, "Kepler", grid="phoenix", method="mu")
t.coefficients("hybrid5", 9000, 4.4, 0.0, "TESS", extrapolate="nearest")
t.bands("phoenix")                         # passband names available
t.intensity("hybrid5", r["coeffs"], mu)    # I(mu)/I(1)
```

- `grid`: `"phoenix"`, `"mps1"`, `"mps2"`. PHOENIX has 40 passbands; MPS-ATLAS
  has 41, including `JWST_MIRI_LRS`.
- Passband names: `TESS`, `Kepler`, `CHEOPS`, `NGTS`, `HST_WFC3_G141`,
  `JWST_NIRSpec_PRISM`, `JWST_NIRSpec_G395H_NRS1`, `JWST_NIRISS_SOSS_o1`, and so
  on (`t.bands(grid)` lists them).
- Interpolation is linear on a Delaunay triangulation in (Teff/100 K,
  log g/0.5, [M/H]/0.5). It bridges the gaps of the PHOENIX grid, is exact at
  the models, and keeps hybrid coefficients physical. Do not use higher-order
  interpolation on the hybrid coefficients: it can leave the physical region.
- Outside a grid, `coefficients` raises `ValueError` unless
  `extrapolate="nearest"`.

From the shell:

```bash
python -m ldhybrid.interpolate hybrid5 5600 4.4 0.1 TESS --grid mps2 --tables path/to/cds
python ldhybrid/interpolate.py hybrid4 5600 4.4 0.1 Kepler --grid mps2 --correction B100G --tables path/to/cds
```

### Magnetic corrections

`correction=` adds a MURaM weight shift Δw (hybrid laws only): `"ssd"`,
`"B100G"`, `"B200G"`, `"B300G"`. `t.delta_w(law, band, teff, state)` returns
the shift itself. The corrections come from seven solar-metallicity MURaM
stars (3196–6807 K). They are interpolated in Teff, clamped outside that range
(with a warning), and applied at all log g and [M/H], which is an
extrapolation.

### Priors for fitting

```python
p = t.prior("hybrid4", 5600, 100, 4.4, 0.1, 0.1, 0.1, "Kepler",
            grid="mps2", correction="B100G", systematic=True, cov_floor=0.01**2)
p["mean"], p["cov"]
```

`prior` draws (Teff, log g, [M/H]) from independent Gaussians (arguments:
value, sigma for each), interpolates each draw, and returns the mean and
covariance of the coefficients.

- `systematic=True` adds the spread between the PHOENIX, MPS-ATLAS Set 1 and
  Set 2 grids at the star. Below 3500 K only PHOENIX covers the star, and the
  spread at 3500 K is substituted (with a warning).
- `cov_floor` raises the variance of each coefficient to at least the floor
  (scalar, or one value per coefficient). **Use it**: without it the covariance
  is nearly singular in directions the stellar parameters don't move, which
  would effectively fix those combinations at the model values. It acts on the
  diagonal only, so correlated directions can stay narrower than the floor
  (e.g. smallest eigen-direction s.d. 0.002–0.005 for the example above with
  `systematic=True`); add `floor * np.eye(n)` to `p["cov"]` yourself if every
  direction must be at least that wide.
- The floor and any fixed prior width are per coefficient, so their meaning
  depends on the basis (hybrid4's curvature basis differs from its
  monomial-pole form); widths on h′₁, h′₂ and I(0) (`observables`) are
  basis-independent.
- Draws outside the grid are dropped with a warning; near a grid edge this
  biases the mean.

### Empirical offsets

Measured offsets between model and observed limb darkening (Maxted 2023 for
Kepler and TESS; Sing et al. 2026 for JWST) can update the prior. Use the
regime that matches the prior's passband (the code checks the grid and
correction state, not the band):

```python
m, C = t.apply_offsets(p, regime="Kepler")    # or "TESS", "JWST_0.86-5.5um"
```

The h′₁ and h′₂ offsets are applied jointly, with their correlation. An offset
is used only if it was measured against the same grid and correction state as
the prior; otherwise the call raises, to avoid counting the magnetic
correction twice.

The lower-level functions are `observables(law)`, which gives h′₁, h′₂, ℓ and δ
as exact linear functions of the coefficients, and
`condition(mean, cov, law, observable, delta, sigma, rho=None)`:

```python
from ldhybrid.interpolate import observables, condition
m, C = condition(p["mean"], p["cov"], "hybrid4", ["h1", "h2"],
                 [-0.0023, -0.0039], [0.0052, 0.0126], rho=0.04)
```

The result is Gaussian and can put a little probability outside the physical
region; truncate there if your sampler needs it.

### Refitting interpolated profiles

```python
t = LDTables("path/to/cds", profiles="path/to/profiles")
r = t.coefficients("claret4", 5600, 4.4, 0.1, "TESS", grid="phoenix", refit=True, mu_min=0.02)
```

Instead of interpolating coefficients, this interpolates the band-integrated
intensity profiles of the surrounding models and refits the law with the
production settings. It needs h5py and the profile files. `profiles` is a
directory holding `<grid>.h5`, or a single `.h5` file for one grid; the
Zenodo copies are named `ldhybrid_v1.0_band_profiles_<grid>.h5`.

For the linear laws without an active constraint, refitting gives the same
result as interpolating coefficients. It helps power-2. `mu_min=0.02` leaves
out the outermost limb, which makes the PHOENIX four-term fits as accurate as
hybrid5.

### Light curves

```python
from ldhybrid.interpolate import light_curve
F = light_curve("hybrid4", r["coeffs"], z, p=0.1)    # z: sky separations in stellar radii
```

This is the numerical quadrature used for the paper's reference light curves
(circular planet; agrees with the closed-form hybrid light curves to 1e-8 ppm).
For fitting, use a dedicated transit code with the same law: SquishierPlanet
(circular or elliptical planets) or Elómë.

## 3. The full package (`ldhybrid.api`)

With the repository installed (`pip install -e .`) and its `data/` products:

```python
import ldhybrid as ld

r = ld.coefficients("hybrid5", 5600, 4.4, 0.1, "TESS", grid="mps2", correction="B100G")
r["w"], r["q"], r["physical"]

# any throughput or spectroscopic bin (wavelengths in Angstrom)
r = ld.coefficients_custom("hybrid4", 4500, 4.6, 0.0, lam_angstrom, throughput)                 # PHOENIX
r = ld.coefficients_custom("hybrid4", 5600, 4.4, 0.0, lam_angstrom, throughput, grid="mps2")    # MPS-ATLAS

law = ld.to_law("hybrid5", r["w"])     # squishierplanet.Law, for closed-form light curves
```

- `coefficients_custom` with PHOENIX sums the R = 3000 coefficient spectra.
  These are on Zenodo as `ldhybrid_v1.0_coefficient_spectra_phoenix.tar`, and
  must be unpacked to `data/products/spectra/phoenix/`.
- With MPS-ATLAS it integrates the model archive directly: about 9.5 GB per
  set, fetched with `scripts/fetch_archives.py mpsatlas`.

`ld.prior`, `ld.observables` and `ld.condition` wrap the functions of
`ldhybrid/interpolate.py`.
