# Standalone hybrid4 interpolators

Single-file interpolators for the LDhybrid v1.0 **hybrid4** coefficients in one
grid and passband. Each needs only numpy and scipy; the coefficient grid is
embedded, so no tables need downloading.

| File | Grid | Passband | Size |
|---|---|---|---|
| `ldhybrid_hybrid4_phoenix_tess.py` | PHOENIX-ACES v3.0 (2300–7400 K) | TESS | 58 kB |
| `ldhybrid_hybrid4_phoenix_kepler.py` | PHOENIX-ACES v3.0 (2300–7400 K) | Kepler | 58 kB |
| `ldhybrid_hybrid4_mps2_tess.py` | MPS-ATLAS Set 2 (3500–7500 K) | TESS | 459 kB |
| `ldhybrid_hybrid4_mps2_kepler.py` | MPS-ATLAS Set 2 (3500–7500 K) | Kepler | 460 kB |

```python
from ldhybrid_hybrid4_mps2_tess import coefficients, intensity
r = coefficients(5600, 4.4, 0.0, correction="B100G")   # recommended for FGK
r["w"], r["q"], r["physical"]          # weights, sampling coordinates q1-q4
```

or `python ldhybrid_hybrid4_mps2_tess.py 5600 4.4 0.0`.

- The MPS-ATLAS files carry the MURaM magnetic corrections (`correction=`
  `"ssd"`, `"B100G"`, `"B200G"`, `"B300G"`).
- Outside a grid, `coefficients` raises `ValueError`; `extrapolate=True` uses
  the nearest model instead.
- They reproduce `ldhybrid.interpolate.LDTables.coefficients` (r-weighted
  method) to 1e-8.
- The weights are on hybrid4's curvature-simplex basis (1 − μ², (1 − μ²)², S3,
  S4; `BASIS` in each file): w ≥ 0, Σw ≤ 1 is exactly the set of hybrid4
  profiles that are positive, monotonic, concave in μ² and with a convex
  slope. `from_q`/`to_q` convert to the sampling coordinates q1–q4 (q1 =
  (1 − I(0))⁴, the analogue of Kipping's q1). For SquishierPlanet or Elómë,
  which take weights on the monomial-pole shapes, use `w @ BASIS`.
- hybrid4 here uses the α = 25/16 poles (0.02860121, 0.38582418). Do not use
  them with hybrid4 implementations that use the old α = 1 poles.

Regenerate with `.venv/bin/python scripts/make_standalone.py`.
