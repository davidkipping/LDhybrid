# LDhybrid

Theoretical limb-darkening coefficient tables for the rational **hybrid2,
hybrid4 and hybrid5** laws of [SquishierPlanet](../squishierplanet), regressed
against simulated stellar specific intensities, with positivity and
monotonicity imposed by construction. A sibling of SquishierPlanet and
MetalPlanet. The plan of work is in [`docs/plan.md`](docs/plan.md).

## Status

**P0–P4 complete; Gates G0–G4 passed** (validation: `studies/p4/results.md`).
**P5 (release) in progress:** CDS tables, user API and paper draft done. Products in `data/products/`
(`docs/products.md`; G3 in `studies/g3/results.md`). v1.0 law and method
definitions are frozen in `ldhybrid/v1.py` (hybrid2 ε = 0.36; r-weighted fits;
M1 limb). Evidence: `docs/data.md`, `studies/g1/results.md`,
`studies/p2/results.md`; first end-to-end run in `studies/p1/smoke.md`. Scope: FGKM dwarfs and subgiants
(Teff 2300–7500 K, log g 3.5–5.5, [M/H] −2 to +1) in the passbands used for
transit observations, with emphasis on HST WFC3 and JWST.

## Use

```python
import ldhybrid as ld
# tabulated coefficients (linear interpolation within one grid)
r = ld.coefficients("hybrid5", teff=5600, logg=4.4, mh=0.1, band="TESS", grid="mps2")
r["w"], r["q"], r["physical"]
# recommended prior centre for FGK stars (P4): MPS-ATLAS Set 2 + MURaM 100 G
r = ld.coefficients("hybrid5", 5600, 4.4, 0.1, "TESS", grid="mps2", correction="B100G")
# any throughput or spectroscopic bin, from the PHOENIX coefficient spectra
r = ld.coefficients_custom("hybrid4", 4500, 4.6, 0.0, lam_angstrom, throughput)
law = ld.to_law("hybrid5", r["w"])   # squishierplanet.Law -> closed-form light curves
```

**From the published tables alone** (no other LDhybrid data needed), the
self-contained `ldhybrid/interpolate.py` (numpy + scipy only) reads the CDS
tables and interpolates them:

```python
from ldhybrid.interpolate import LDTables
t = LDTables("path/to/cds")       # the <grid>_<law>.dat[.gz] files
r = t.coefficients("hybrid5", teff=5600, logg=4.4, mh=0.1, band="TESS", grid="mps2")
r["coeffs"], r["physical"], r["lc_ppm"]
```

```bash
python -m ldhybrid.interpolate hybrid5 5600 4.4 0.1 TESS --grid mps2 --tables path/to/cds
```

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
