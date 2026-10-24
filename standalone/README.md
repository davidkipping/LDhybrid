# Standalone interpolators (numpy and MLX)

Single-file interpolators for the LDhybrid v1.0 coefficients of the hybrid laws
on every atmosphere grid, in the TESS, Kepler and CHEOPS passbands. The coefficient grid is
embedded, so no tables need downloading. Each file comes in two versions: numpy
only (`ldhybrid_<law>_<grid>_<band>.py`) and an MLX port for MLX codes such as
Elómë (`..._mlx.py`; mlx and numpy).

| Law | Coefficients | q (sampling) | Elómë |
|---|---|---|---|
| `hybrid2` | w1, w2 | q1, q2 (exact triangle) | `limb_dark="hybrid2"`, `u = r["w"]` |
| `hybrid4` | w1–w4 (curvature basis) | q1–q4 | `limb_dark="hybrid4"`, `u = r["w"]` |
| `hybrid5` | w1–w5 | q1–q5 (stick-breaking) | `limb_dark="hybrid5"`, `u = r["w"]` |

| Grid | Models | Teff, log g, [M/H] | numpy / MLX file size |
|---|---|---|---|
| `phoenix`: PHOENIX-ACES v3.0 | 1,682 (+68 filled) | 2300–7400 K, 3.5–5.5, −2.0 to +1.0 | 48–77 kB |
| `mps1`: MPS-ATLAS Set 1 | 16,974 | 3500–7500 K, 3.5–5.0, −2.0 to +1.0 | 296–547 kB |
| `mps2`: MPS-ATLAS Set 2 | 16,974 | 3500–7500 K, 3.5–5.0, −2.0 to +1.0 | 298–553 kB |

Bands: `tess`, `kepler`, `cheops`. So 27 numpy files and 27 MLX files, e.g.
`ldhybrid_hybrid2_phoenix_kepler.py`, `ldhybrid_hybrid5_mps1_cheops.py` or
`ldhybrid_hybrid4_mps2_tess_mlx.py`.

The comparison laws of the paper (quadratic, power-2, Claret four-term) have no
standalone files: published tables exist for them, and `LDTables` serves
LDhybrid's fits of them (`scripts/make_standalone.py --all-laws` writes them).

## numpy

```python
from ldhybrid_hybrid4_mps2_tess import coefficients, intensity
r = coefficients(5600, 4.4, 0.0, correction="B100G")   # Teff, log g, [M/H]; recommended for FGK
r["w"], r["q"], r["physical"]          # weights (also r["coeffs"]), q, positive and monotonic?
r["lc_ppm"], r["mu_cri"], r["inside"]  # fit's light-curve error (p = 0.1), stellar edge, inside the grid?
I = intensity(r["w"], [0.1, 0.5, 1.0])   # I(mu)/I(1)
```

Every file has `coefficients`, `shapes`, `intensity`, `is_physical`, `from_q`
and `to_q`. From the shell: `python ldhybrid_hybrid2_phoenix_tess.py 3200 5.0 0.0`.
Outside a grid, `coefficients` raises `ValueError`; `extrapolate=True` uses
the nearest model instead.

## MLX

```python
import mlx.core as mx
from ldhybrid_hybrid4_mps2_tess_mlx import coefficients
r = coefficients(mx.array([5600.0, 5720.0]), mx.array([4.40, 4.35]), 0.1, correction="B100G")
r["w"]                                  # (2, 4) mx.array: Elómë's u for limb_dark="hybrid4"
mx.grad(lambda t: coefficients(t, 4.4, 0.1)["w"].sum())(mx.array([5600.0]))   # d/dTeff
```

- The inputs broadcast against each other (numbers, numpy or mx arrays), and
  every output has their shape (`w` and `q` with a trailing axis of 2, 4 or
  5). The functions can be differentiated with respect
  to Teff, log g and [M/H], and compiled with `mx.compile`.
- Precision follows the inputs, as in Elómë. float32 runs on the default
  device. Anything else runs in float64 on the CPU stream, where it equals the
  numpy file to the last bit. float32 agrees with that to about 1e-7 in w.
- Nothing raises inside a model. Outside the grid the outputs are NaN and
  `inside` is False; with `extrapolate=True` the nearest model is used (zero
  gradient).

## All files

- **Corrections.** The MURaM magnetic corrections (`correction=` `"ssd"`,
  `"B100G"`, `"B200G"`, `"B300G"`) exist for the hybrid laws on MPS-ATLAS Set 2,
  so only the Set 2 files have a `correction` argument.
- **Interpolation.** It is linear on the Kuhn simplices of the (Teff, log g,
  [M/H]) tensor grid, as in `ldhybrid.interpolate.LDTables`: every grid box is
  split into the same six tetrahedra.
  - The coefficients are continuous in the stellar parameters.
  - They are exact at the models.
  - Every value is a convex mix of four models, so it stays in the law's
    region.
  - The 68 missing PHOENIX models are filled from their neighbours along one
    axis.
  - The files reproduce `LDTables.coefficients` (r-weighted method) to 5e-8.
- **hybrid4.** The weights are on hybrid4's curvature-simplex basis (1 − μ²,
  (1 − μ²)², S3, S4; `BASIS` in each file). In this basis w ≥ 0, Σw ≤ 1 is
  exactly the set of hybrid4 profiles that are positive, monotonic, concave in
  μ² and of convex slope.
  - `from_q`/`to_q` convert to the sampling coordinates q1–q4 (q1 =
    (1 − I(0))⁴, the analogue of Kipping's q1).
  - The current Elómë and SquishierPlanet take these weights as they are.
  - A code that takes weights on the monomial-pole shapes needs `w @ BASIS`.
  - hybrid4 here uses the γ = 25/16 poles (0.02860121, 0.38582418). Do not
    use these weights with hybrid4 implementations that use the old α = 1
    poles.
- **hybrid5.** The weights are on the shapes 1 − μ⁴, (1 − μ²)² and three
  poles (0.0015965, 0.041275, 0.38743).
  - The fits are constrained to the simplex w ≥ 0, Σw ≤ 1, which is
    sufficient for a positive, monotonic profile.
  - `from_q`/`to_q` are stick-breaking of a flat Dirichlet, so a uniform q1–q5
    is a uniform prior over the simplex.
- **hybrid2.** The weights are on the shapes 1 − μ² and one pole (0.36). Its
  positive, monotonic profiles are exactly a triangle, the fits are
  constrained to it, and `from_q` maps the unit square uniformly onto it (the
  triangle Elómë's `hybrid2_from_q` samples). `to_q` keeps the sign of the
  depth in q1, so weights outside the triangle give q outside [0, 1].

Regenerate both versions with `.venv/bin/python scripts/make_standalone.py`
(other passbands: add them to `BANDS` there).
