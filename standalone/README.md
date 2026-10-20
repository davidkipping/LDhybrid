# Standalone hybrid4 and hybrid5 interpolators (numpy and MLX)

Single-file interpolators for the LDhybrid v1.0 **hybrid4** and **hybrid5**
coefficients in one grid and passband. The coefficient grid is embedded, so no
tables need downloading. Each file comes in two versions: numpy only
(`ldhybrid_<law>_<grid>_<band>.py`) and an MLX port for MLX codes such as
Elómë (`..._mlx.py`; mlx and numpy).

| Law | Grid | Passband | numpy file | MLX file |
|---|---|---|---|---|
| hybrid4 | PHOENIX-ACES v3.0 (2300–7400 K) | TESS | `ldhybrid_hybrid4_phoenix_tess.py` (64 kB) | `ldhybrid_hybrid4_phoenix_tess_mlx.py` (68 kB) |
| hybrid4 | PHOENIX-ACES v3.0 (2300–7400 K) | Kepler | `ldhybrid_hybrid4_phoenix_kepler.py` (65 kB) | `ldhybrid_hybrid4_phoenix_kepler_mlx.py` (69 kB) |
| hybrid4 | MPS-ATLAS Set 2 (3500–7500 K) | TESS | `ldhybrid_hybrid4_mps2_tess.py` (462 kB) | `ldhybrid_hybrid4_mps2_tess_mlx.py` (466 kB) |
| hybrid4 | MPS-ATLAS Set 2 (3500–7500 K) | Kepler | `ldhybrid_hybrid4_mps2_kepler.py` (463 kB) | `ldhybrid_hybrid4_mps2_kepler_mlx.py` (468 kB) |
| hybrid5 | PHOENIX-ACES v3.0 (2300–7400 K) | TESS | `ldhybrid_hybrid5_phoenix_tess.py` (72 kB) | `ldhybrid_hybrid5_phoenix_tess_mlx.py` (75 kB) |
| hybrid5 | PHOENIX-ACES v3.0 (2300–7400 K) | Kepler | `ldhybrid_hybrid5_phoenix_kepler.py` (72 kB) | `ldhybrid_hybrid5_phoenix_kepler_mlx.py` (76 kB) |
| hybrid5 | MPS-ATLAS Set 2 (3500–7500 K) | TESS | `ldhybrid_hybrid5_mps2_tess.py` (543 kB) | `ldhybrid_hybrid5_mps2_tess_mlx.py` (548 kB) |
| hybrid5 | MPS-ATLAS Set 2 (3500–7500 K) | Kepler | `ldhybrid_hybrid5_mps2_kepler.py` (548 kB) | `ldhybrid_hybrid5_mps2_kepler_mlx.py` (552 kB) |

## numpy

```python
from ldhybrid_hybrid4_mps2_tess import coefficients, intensity
r = coefficients(5600, 4.4, 0.0, correction="B100G")   # Teff, log g, [M/H]; recommended for FGK
r["w"], r["q"], r["physical"]          # weights, sampling coordinates q, positive and monotonic?
r["lc_ppm"], r["mu_cri"], r["inside"]  # fit's light-curve error (p = 0.1), stellar edge, inside the grid?
I = intensity(r["w"], [0.1, 0.5, 1.0]) # I(mu)/I(1)
```

The hybrid5 files have the same functions; `w` and `q` have five entries.
From the shell: `python ldhybrid_hybrid5_mps2_tess.py 5600 4.4 0.0`.
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
  every output has their shape (`w` and `q` with a trailing axis of 4 or 5).
  The functions can be differentiated with respect to Teff, log g and [M/H],
  and compiled with `mx.compile`.
- Precision follows the inputs, as in Elómë. float32 runs on the default
  device. Anything else runs in float64 on the CPU stream, where it equals the
  numpy file to the last bit. float32 agrees with that to about 1e-7 in w.
- Nothing raises inside a model. Outside the grid the outputs are NaN and
  `inside` is False; with `extrapolate=True` the nearest model is used (zero
  gradient).

## Both

- The MPS-ATLAS files carry the MURaM magnetic corrections (`correction=`
  `"ssd"`, `"B100G"`, `"B200G"`, `"B300G"`). The PHOENIX files have no
  `correction` argument.
- **Interpolation.** It is linear on the Kuhn simplices of the (Teff, log g,
  [M/H]) tensor grid, as in `ldhybrid.interpolate.LDTables`: every grid box is
  split into the same six tetrahedra.
  - The coefficients are continuous in the stellar parameters.
  - They are exact at the models.
  - Every value is a convex mix of four models, so it stays physical.
  - The 68 missing PHOENIX models are filled from their neighbours along one
    axis.
  - The files reproduce `LDTables.coefficients` (r-weighted method) to 1e-8 in
    the weights.
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

Regenerate both versions with `.venv/bin/python scripts/make_standalone.py`.
