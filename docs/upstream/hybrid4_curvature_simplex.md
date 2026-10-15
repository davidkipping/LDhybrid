# Note from LDhybrid: hybrid4 in the curvature-simplex basis (v1.0, 2026-10-14)

**Status: applied in LDhybrid; to do in SquishierPlanet, Elómë and Turin.**
Poles unchanged (α = 25/16: ε = 0.028601207249486824, 0.38582417584735224).

## Change

LDhybrid tabulates hybrid4 on a new basis, with new sampling coordinates.

    I(μ)/I(1) = 1 − w1 S1 − w2 S2 − w3 S3 − w4 S4,
    S_k = Σ_j B[k, j] T_j,   T = (1 − μ⁴, (1 − μ²)², Π_ε0, Π_ε1)  (the current upstream hybrid4 shapes)

    B = [
    [0.5, 0.5, 0, 0],
    [0, 1, 0, 0],
    [0.0006260842489149928, 0.001404013755411423, 1.020963273947444, -0.02299337195177007],
    [-0.08116019831892918, -0.2568535772254558, -0.0005906631447752816, 1.33860443868916],
    ]

S1 = 1 − μ² and S2 = (1 − μ²)² exactly. S3 and S4 are the two pole-dominated
vertices (below). Weights on the upstream shapes are `w_upstream = w @ B`; back
with `w = w_upstream @ inv(B)`. The light curve is unchanged (same span).

## Why

With x = μ² and g = dI/dx, the hybrid4 profiles with

    I(0) ≥ 0,  g ≥ 0  (Kipping 2013)  and  g' ≤ 0,  g'' ≥ 0  on [0, 1]

(positive, monotonic, concave in μ², convex slope) are **exactly** the simplex
w ≥ 0, Σw ≤ 1 in this basis. The reason is that g'' spans only the two poles,
whose ratio is monotonic, so the conditions reduce to Σw ≤ 1, g(1) ≥ 0,
g'(1) ≤ 0, g''(0) ≥ 0 and g''(1) ≥ 0. Its vertices are:
- I = 1;
- I = μ²;
- I = 1 − (1 − μ²)²;
- S3: I(0) = 0, g(1) = g'(1) = g''(1) = 0;
- S4: I(0) = 0, g(1) = g'(1) = 0, g''(0) = 0.

The old w-simplex on the upstream shapes was sufficient but not complete. 44%
of it is profiles convex in μ² (e.g. I = μ⁴), and it misses other physical,
concave profiles.

Evidence: `studies/rev3/h4_scheme.md`. Every tabulated fit (2.9 M rows) and
every MPS-ATLAS profile satisfies the curvature conditions. PHOENIX violates
them only at μ < 0.1. 400k uniform draws from the simplex all satisfy the
four conditions, and 82k uniform draws from the four-condition region all lie
in the simplex.

## Sampling coordinates (replace `simplex_from_q` for hybrid4)

These are the analogue of Kipping's q1, q2. With s = Σw = 1 − I(0) = q1^(1/4):

    w1 = s (1 − q2^(1/3))
    w2 = s q2^(1/3) (1 − q3^(1/2))
    w3 = s q2^(1/3) q3^(1/2) (1 − q4)
    w4 = s q2^(1/3) q3^(1/2) q4

A uniform q in [0, 1]⁴ gives a uniform prior over the region. Every q gives a
profile satisfying all four conditions, and every such profile has a q.
Reference implementations: `ldhybrid/qspace.py` (`hybrid4_from_q`,
`hybrid4_to_q`) and `ldhybrid/interpolate.py` (`from_q`, `to_q`, `H4_BASIS`).

## What to change upstream

* **SquishierPlanet** (`squishierplanet/laws.py`). Add the basis matrix
  (computed from the poles, as `ldhybrid.v1._hybrid4_vertices`, or the values
  above) and a `hybrid4_from_q`. Either make `hybrid_shapes("hybrid4")` return
  S1–S4, or keep T and document `w @ B`. Then update the docstrings, tests and
  notes.
* **Elómë** and **Turin**. For hybrid4 sampling, replace the stick-breaking on
  the w-simplex with `w = hybrid4_from_q(q)`, then `w_upstream = w @ B` if the
  light-curve code takes upstream weights. Priors from the LDhybrid tables
  (w, q columns) are in the new basis.
* hybrid5 and hybrid2 are unchanged.
