# Standalone interpolators (numpy and MLX)

Single-file interpolators for the LDhybrid v1.0 coefficients of every law and
atmosphere grid, in the TESS and Kepler passbands. The coefficient grid is
embedded, so no tables need downloading. Each file comes in two versions: numpy
only (`ldhybrid_<law>_<grid>_<band>.py`) and an MLX port for MLX codes such as
Elómë (`..._mlx.py`; mlx and numpy).

| Law | Coefficients | q (sampling) | Elómë |
|---|---|---|---|
| `hybrid2` | w1, w2 | q1, q2 (exact triangle) | `limb_dark="hybrid2"`, `u = r["w"]` |
| `hybrid4` | w1–w4 (curvature basis) | q1–q4 | `limb_dark="hybrid4"`, `u = r["w"]` |
| `hybrid5` | w1–w5 | q1–q5 (stick-breaking) | `limb_dark="hybrid5"`, `u = r["w"]` |
| `quadratic` | u1, u2 | Kipping's q1, q2 | `flux_dev(z, r, u1, u2)`; `limb_dark="quadratic"` samples the same q |
| `power2` | c, alpha | — | not in Elómë |
| `claret4` | a1–a4 | — | not in Elómë |

| Grid | Models | Teff, log g, [M/H] | numpy / MLX file size |
|---|---|---|---|
| `phoenix`: PHOENIX-ACES v3.0 | 1,682 (+68 filled) | 2300–7400 K, 3.5–5.5, −2.0 to +1.0 | 46–77 kB |
| `mps1`: MPS-ATLAS Set 1 | 16,974 | 3500–7500 K, 3.5–5.0, −2.0 to +1.0 | 297–547 kB |
| `mps2`: MPS-ATLAS Set 2 | 16,974 | 3500–7500 K, 3.5–5.0, −2.0 to +1.0 | 296–552 kB |

Bands: `tess`, `kepler`. So 36 numpy files and 36 MLX files, e.g.
`ldhybrid_quadratic_phoenix_kepler.py` or `ldhybrid_hybrid4_mps2_tess_mlx.py`.

## numpy

```python
from ldhybrid_hybrid4_mps2_tess import coefficients, intensity
r = coefficients(5600, 4.4, 0.0, correction="B100G")   # Teff, log g, [M/H]; recommended for FGK
r["coeffs"], r["q"], r["physical"]     # coefficients (hybrids: also r["w"]), q, positive and monotonic?
r["lc_ppm"], r["mu_cri"], r["inside"]  # fit's light-curve error (p = 0.1), stellar edge, inside the grid?
I = intensity(r["coeffs"], [0.1, 0.5, 1.0])   # I(mu)/I(1)
```

Every file has `coefficients`, `intensity` and `is_physical`; the linear laws
also have `shapes`, and the laws with sampling coordinates `from_q`/`to_q`.
From the shell: `python ldhybrid_power2_phoenix_tess.py 3200 5.0 0.0`.
Outside a grid, `coefficients` raises `ValueError`; `extrapolate=True` uses
the nearest model instead.

## MLX

```python
import mlx.core as mx
from ldhybrid_hybrid4_mps2_tess_mlx import coefficients
r = coefficients(mx.array([5600.0, 5720.0]), mx.array([4.40, 4.35]), 0.1, correction="B100G")
r["coeffs"]                             # (2, 4) mx.array: Elómë's u for limb_dark="hybrid4"
mx.grad(lambda t: coefficients(t, 4.4, 0.1)["coeffs"].sum())(mx.array([5600.0]))   # d/dTeff
```

- The inputs broadcast against each other (numbers, numpy or mx arrays), and
  every output has their shape (`coeffs` and `q` with a trailing axis of the
  number of coefficients). The functions can be differentiated with respect
  to Teff, log g and [M/H], and compiled with `mx.compile`.
- Precision follows the inputs, as in Elómë. float32 runs on the default
  device. Anything else runs in float64 on the CPU stream, where it equals the
  numpy file to the last bit. float32 agrees with that to about 1e-7
  (Claret-4: 1e-6, its coefficients are large and partly cancel).
- Nothing raises inside a model. Outside the grid the outputs are NaN and
  `inside` is False; with `extrapolate=True` the nearest model is used (zero
  gradient).

## All files

- **Corrections.** The MURaM magnetic corrections (`correction=` `"ssd"`,
  `"B100G"`, `"B200G"`, `"B300G"`) exist for the hybrid laws on MPS-ATLAS Set 2,
  so only those six files (and their MLX twins) have a `correction` argument.
- **Interpolation.** It is linear on the Kuhn simplices of the (Teff, log g,
  [M/H]) tensor grid, as in `ldhybrid.interpolate.LDTables`: every grid box is
  split into the same six tetrahedra.
  - The coefficients are continuous in the stellar parameters.
  - They are exact at the models.
  - Every value is a convex mix of four models. Every law's physical region is
    convex, so a mix of physical models is physical.
  - The 68 missing PHOENIX models are filled from their neighbours along one
    axis.
  - The files reproduce `LDTables.coefficients` (r-weighted method) to 5e-8.
- **Fits.** hybrid2, hybrid4, hybrid5 and quadratic are fitted within their
  sampling region, Claret-4 within the positive, monotonic profiles. power-2 is
  fitted unconstrained, as usual for that law, and 3–22% of its models have
  c > 1 (I(0) < 0, most in PHOENIX); `physical` flags them.
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
- **hybrid2 and quadratic.** Their positive, monotonic profiles are exactly a
  triangle, and `from_q` maps the unit square uniformly onto it (quadratic:
  Kipping 2013, the same q as Elómë's `q_to_u`). `to_q` keeps the sign of the
  depth in q1, so coefficients outside the triangle give q outside [0, 1].

Regenerate both versions with `.venv/bin/python scripts/make_standalone.py`
(other passbands: add them to `BANDS` there).
