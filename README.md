# LDhybrid

Theoretical limb-darkening coefficients for the rational **hybrid2, hybrid4 and
hybrid5** laws of SquishierPlanet, whose transit light curves are elementary.
The quadratic, power-2 and four-term laws are fitted the same way only as a
baseline for comparison; for those laws the published tables remain the usual
source. The coefficients are least-squares fits, uniform in projected radius, to the
specific intensities of PHOENIX-ACES v3.0 (spherical, 2300–7400 K, 40 passbands)
and MPS-ATLAS Sets 1 and 2 (plane-parallel, 3500–7500 K, 41 passbands) for
dwarfs and subgiants, in the photometric, HST and JWST passbands used for
transits. Every row carries its light-curve error, and the hybrid fits are
physical by construction. This repository holds the tables' code, the
interpolators and the paper (`paper/ldhybrid.tex`). Sibling projects:
SquishierPlanet (closed-form light curves) and Elómë.

**Full user guide: [`docs/interpolator.md`](docs/interpolator.md).**

## Quick start

### 1. Get the tables

The v1.0 release (CDS format: 18 tables `<grid>_<law>.dat.gz`, `corrections.dat`,
`offsets.dat`, `systematics.dat`, `ReadMe`, `SHA256SUMS`) is attached to the
GitHub draft release `v1.0-cds-draft` and will be at the CDS on publication. The
band profiles (for refits) and the PHOENIX coefficient spectra (for custom
passbands) are in the Zenodo record (DOI: *to be added on publication*).

```bash
gh release download v1.0-cds-draft --dir cds
cd cds && shasum -a 256 -c SHA256SUMS
```

The grids are `phoenix`, `mps1` and `mps2`; the laws `hybrid2`, `hybrid4`,
`hybrid5`, `quadratic`, `claret4` (four-term) and `power2`. The CDS `ReadMe`
describes every column.

### 2. Standalone single-file interpolators (no download)

[`standalone/`](standalone/) holds a single file for each hybrid law and grid in
TESS, Kepler, CHEOPS and JWST NIRSpec PRISM (0.6–5.3 µm), with the coefficient grid embedded (numpy only):
`ldhybrid_<law>_<grid>_<band>.py` for `hybrid2`, `hybrid4`, `hybrid5`, the grids
`phoenix`, `mps1`, `mps2` and the bands `tess`, `kepler`, `cheops`, `nirspec_prism` (36 files). Each has an **MLX** port
(`*_mlx.py`; mlx + numpy) whose coefficients can be batched, differentiated with
respect to Teff, log g and [M/H], and compiled inside an MLX model such as
Elómë's. Copy one anywhere:

```python
from ldhybrid_hybrid4_mps2_tess import coefficients, intensity
r = coefficients(5600, 4.4, 0.0, correction="B100G")   # Teff, log g, [M/H]; recommended for FGK
r["w"], r["q"], r["physical"], r["lc_ppm"]               # weights, q1-q4, flags, fit's LC error
I = intensity(r["w"], [0.1, 0.5, 1.0])                   # I(mu)/I(1)
```

or `python ldhybrid_hybrid2_phoenix_kepler.py 3200 5.0 0.0`. In MLX:

```python
import mlx.core as mx
from ldhybrid_hybrid4_mps2_tess_mlx import coefficients
w = coefficients(mx.array([5600.0, 5720.0]), 4.4, 0.1, correction="B100G")["w"]  # (2, 4), differentiable
```

### 3. The table interpolator (`ldhybrid/interpolate.py`)

One self-contained file (numpy + scipy) that reads the CDS tables, for all six
laws, three grids and every passband. It can be copied out of the repository.

```python
from ldhybrid.interpolate import LDTables, from_q, to_q, light_curve
t = LDTables("cds")                                      # the downloaded tables
r = t.coefficients("hybrid4", teff=5600, logg=4.4, mh=0.1, band="JWST_NIRSpec_PRISM",
                   grid="mps2", correction="B100G")
r["coeffs"], r["q"], r["physical"], r["lc_ppm"]
t.coefficients("claret4", 3300, 5.0, 0.0, "TESS", grid="phoenix")
t.bands("mps2")                                          # passband names

w = from_q("hybrid4", [0.5, 0.5, 0.5, 0.5])              # unit hypercube -> weights
q = to_q("hybrid4", w)

# prior for a fit: stellar-parameter draws + grid-to-grid spread + a variance floor
p = t.prior("hybrid4", 5600, 100, 4.4, 0.1, 0.1, 0.1, "Kepler", grid="mps2",
            correction="B100G", systematic=True, cov_floor=0.01**2)
m, C = t.apply_offsets(p, regime="Kepler")               # empirical offsets (Maxted 2023)

import numpy as np
z = np.linspace(0, 1.2, 200)                             # sky separation / R_star
F = light_curve("hybrid4", r["coeffs"], z, p=0.1)        # circular planet
```

```bash
python -m ldhybrid.interpolate hybrid5 5600 4.4 0.1 TESS --grid mps2 --tables cds
```

`cov_floor` raises every eigenvalue of the covariance to at least the floor,
so no combination of the coefficients is narrower than √floor. `apply_offsets`
checks that the prior's passband is one the offsets were measured in (Kepler,
TESS, or a JWST band for the JWST regime). Refits of interpolated
profiles (`refit=True`) need the Zenodo band profiles and h5py; see the
[guide](docs/interpolator.md).

### 4. The full package

With the repository installed (`pip install -e .`) and its `data/` products
built or downloaded (band tables, coefficient spectra):

```python
import ldhybrid as ld
r = ld.coefficients("hybrid5", 5600, 4.4, 0.1, "TESS", grid="mps2", correction="B100G")
r = ld.coefficients_custom("hybrid4", 4500, 4.6, 0.0, lam_angstrom, throughput)  # any passband
law = ld.to_law("hybrid5", r["w"])     # squishierplanet.Law (needs SquishierPlanet)
```

### hybrid4's basis

hybrid4 (poles ε = 0.028601, 0.38582, from the γ = 25/16 ladder) is tabulated
on its **curvature-simplex basis**: w ≥ 0, Σw ≤ 1 is exactly the set of hybrid4
profiles that are positive, monotonic, concave in μ² and of convex slope in μ²,
and q1–q4 ∈ [0, 1] (q1 = (1 − I(0))⁴) map onto it uniformly. Codes that take
weights on this basis, as the current SquishierPlanet and Elómë do, use the
weights as they are; codes that take weights on the monomial-pole shapes
(1 − μ⁴, (1 − μ²)², Π_ε0, Π_ε1) need `w @ H4_BASIS` (`BASIS` in the standalone
files), and the expanded-form columns of the tables are basis-free. Weights or q
from earlier hybrid4 versions (α = 1 poles, or the monomial-pole basis) are not
interchangeable: recompute q from w. Details:
[`docs/upstream/hybrid4_curvature_simplex.md`](docs/upstream/hybrid4_curvature_simplex.md).

## For developers

**Status.** P0–P5 complete: v1.0 CDS tables, interpolators, coefficient spectra
and band profiles released (GitHub draft release and Zenodo draft); the paper
is under internal (mock) referee review (`studies/review/`). v1.0 law and
method definitions are frozen in `ldhybrid/v1.py` (hybrid2 ε = 0.36; hybrid4
curvature basis with γ = 25/16 poles; r-weighted fits; M1 limb). Validation:
`studies/p4/results.md`, `studies/rev*/`. Plan of work: [`docs/plan.md`](docs/plan.md).

Products are described in [`docs/products.md`](docs/products.md): band
tables for PHOENIX (40 passbands) and MPS-ATLAS Sets 1 and 2 (41, including
MIRI LRS and the G395H NRS1/NRS2 windows; added to existing production with
`scripts/p3_extend.py`, which also exports the band profiles to
`data/products/profiles/<grid>.h5`). The CDS
release (fixed-width tables + ReadMe) is built by `scripts/p5_cds.py`, the
paper draft is `paper/ldhybrid.tex` (figures: `paper/make_figures.py`;
compile with `tectonic`).

## Data

Everything lives under `data/` (git-ignored). Each file is logged with URL,
size, sha256 and time in `data/provenance.jsonl`.

| step | command | what it does |
|---|---|---|
| PHOENIX | `scripts/fetch_phoenix.py` | 1,682 SpecInt v3.0 models (57 GB raw), reduced on the fly to photon integrals in R = 3000 log bins over 0.1–5.5 µm (`data/reduced/phoenix/`); raw FITS kept only for the P2 subset |
| archives | `scripts/fetch_archives.py` | MURaM 3D-MHD library, MPS-ATLAS Set 1/2, Stagger and Kurucz (ExoTiC-LD bundle), PLATO MARCS release |
| tables | `scripts/fetch_tables.py` | CDS mirrors of the Claret (2011–2025) and Kostogryz (2022) catalogues |
| bands | `scripts/fetch_bands.py` | SVO photometric filters and the ExoTiC-LD HST/JWST throughputs |

All scripts are resumable; run them with `.venv/bin/python`.

## Layout

| module | role |
|---|---|
| `specint.py` | `SpecInt`: binned photon integrals P[k, μ] = ∫_bin I_λ λ dλ, exact bin integration, photon-weighted band profiles with a coverage guard |
| `bands.py` | passband registry; every curve returned as a photon response (energy-counter curves divided by λ) |
| `grids/phoenix.py` | PHOENIX catalogue, reduction, P2 subset, loader |
| `grids/mpsatlas.py`, `grids/exoticld.py` | MPS-ATLAS Set 1/2 (read from the zips); Kurucz and Stagger in ExoTiC-LD format |
| `profile.py` | normalised profiles; spherical limb cut (M1, M2, peak); fitting quadrature (μ, r, area, native nodes; μ_min) |
| `laws.py` | shape-form laws: hybrid2 (pole as a parameter), hybrid4/5, linear, quadratic, Claret 4-term; power-2 |
| `fit.py` | weighted least squares as a convex QP (none / simplex / exact physical region, ± flux conservation) with an exact active-set polish; power-2 by LM |
| `lightcurve.py` | circular-planet light curves of any radial profile (≤ 1e-8 ppm vs the exact hybrid solution); (p, b) scoring grid |
| `pipeline.py` | one table row: model × band × law × method, with coefficients and metrics |
| `reference.py` | readers for the Kostogryz (2022) tables |
| `v1.py` | frozen v1.0 law and method definitions (Gate G2) |
| `products.py`, `spectra.py`, `qspace.py` | P3 band-table rows, coefficient spectra, sampler coordinates |
| `grids/muram.py` | MURaM 3D-MHD library (correction layer) |
| `download.py`, `paths.py` | resumable downloads with checksums and the provenance log |

## Tests

```bash
.venv/bin/python -m pytest tests -q          # data tests skip if data/ is absent
```
