# Brief from LDhybrid: hybrid4 moves to the curvature-simplex basis (v1.0, 2026-10-14)

**For:** SquishierPlanet, Elómë, Turin. **Status:** applied and released in LDhybrid
(commit 19d1732, GitHub draft `v1.0-cds-draft`); to do upstream.

## In one paragraph

hybrid4 keeps its span: the same poles, the same light curves and the same fits. What
changes is the basis that the four weights refer to. The new shapes S1–S4 are the
vertices of a 4-simplex. In that basis, w ≥ 0, Σw ≤ 1 is **exactly** the set of
hybrid4 profiles that are positive, monotonic, concave in μ² and of convex slope in
μ². It is no longer just a sufficient region. Its sampling coordinates q1–q4 are the
direct analogue of Kipping's (2013) q1, q2: every q in [0, 1]⁴ is such a profile,
and every such profile has a q. Every LDhybrid hybrid4 table, prior, interpolator and
the paper now quote w and q in this basis, so upstream code has to agree.
hybrid2 and hybrid5 are unchanged.

## The change

Poles are unchanged (α = 25/16): ε = 0.028601207249486824, 0.38582417584735224.

    I(μ)/I(1) = 1 − w1 S1 − w2 S2 − w3 S3 − w4 S4
    S_k = Σ_j B[k, j] T_j,   T = (1 − μ⁴, (1 − μ²)², Π_ε0, Π_ε1)   (today's upstream shapes)

    B = [[ 0.5,                    0.5,                   0,                       0                   ],
         [ 0,                      1,                     0,                       0                   ],
         [ 0.0006260842489149928,  0.001404013755411423,  1.020963273947444,      -0.02299337195177007 ],
         [-0.08116019831892918,   -0.2568535772254558,   -0.0005906631447752816,   1.33860443868916    ]]

- S1 = 1 − μ² and S2 = (1 − μ²)² exactly.
- S3 and S4 are the pole-dominated vertices: S3 ≈ 1.021 Π_ε0 and S4 ≈ 1.339 Π_ε1, each
  with small polynomial and cross-pole corrections.
- Conversion:
  - **upstream weights = w @ B**;
  - back: **w = w_upstream @ inv(B)**.

## Why

With x = μ² and g = dI/dx, the conditions are I(0) ≥ 0 and g ≥ 0 (Kipping), plus
g′ ≤ 0 and g″ ≥ 0 on [0, 1]. In the hybrid4 span, g″ involves only the two poles,
whose ratio is monotonic. The four conditions therefore reduce to five linear
inequalities: Σw ≤ 1, g(1) ≥ 0, g′(1) ≤ 0, g″(0) ≥ 0, g″(1) ≥ 0. Those five define a
simplex, with vertices:

| Vertex | Profile, or what defines it |
|---|---|
| — | I = 1 |
| S1 | I = μ² |
| S2 | I = 1 − (1 − μ²)² |
| S3 | I(0) = 0, g(1) = g′(1) = g″(1) = 0 |
| S4 | I(0) = 0, g(1) = g′(1) = 0, g″(0) = 0 |

Compared with the old region:
- The old w-simplex on T missed physical, concave profiles. It covered 13% of the
  Kipping region; 44% of the old simplex is profiles convex in μ², such as I = μ⁴.
- The new simplex is 9% of the Kipping region. 96% of what it excludes is convex in
  μ² somewhere, which no model atmosphere produces.

How well the data support the curvature conditions:
- All 2.9 M tabulated fits satisfy them with margin.
- Every MPS-ATLAS profile satisfies them; PHOENIX profiles break them only at μ < 0.1
  (the inflection of the spherical limb).
- 123 PHOENIX fits (2300–2600 K, blue bands) would gain from dropping the conditions:
  their median light-curve error falls from 26 to 15 ppm.

Details: `LDhybrid/studies/rev3/h4_scheme.md`, `studies/rev4/rv_bound.md`, and §2 of
the paper.

## Reference code (self-contained numpy)

```python
import numpy as np

def hybrid4_basis(eps):
    """Rows: weights of S1..S4 on T = (1-mu^4, (1-mu^2)^2, Pi_e0, Pi_e1)."""
    e0, e1 = eps
    a0, a1 = (1.0 / (e ** -2 - (1 + e) ** -2) for e in eps)
    g  = lambda x: [2 * x, 2 * (1 - x), 2 * a0 * (x + e0) ** -3, 2 * a1 * (x + e1) ** -3]
    g1 = lambda x: [2.0, -2.0, -6 * a0 * (x + e0) ** -4, -6 * a1 * (x + e1) ** -4]
    g2 = lambda x: [0.0, 0.0, 24 * a0 * (x + e0) ** -5, 24 * a1 * (x + e1) ** -5]
    rhs = [1.0, 0.0, 0.0, 0.0]
    S3 = np.linalg.solve([[1.0] * 4, g(1.0), g1(1.0), g2(1.0)], rhs)
    S4 = np.linalg.solve([[1.0] * 4, g(1.0), g1(1.0), g2(0.0)], rhs)
    return np.array([[0.5, 0.5, 0, 0], [0, 1.0, 0, 0], S3, S4])

def hybrid4_from_q(q):
    """q (..., 4) in [0,1]^4 -> w (..., 4), uniform on the curvature simplex.
    s = sum w = 1 - I(0) = q1^(1/4); then stick-breaking over the four shapes."""
    q = np.asarray(q, float)
    s = q[..., 0] ** 0.25
    a, b = q[..., 1] ** (1 / 3), q[..., 2] ** 0.5
    return np.stack([s * (1 - a), s * a * (1 - b),
                     s * a * b * (1 - q[..., 3]), s * a * b * q[..., 3]], axis=-1)

def hybrid4_to_q(w):
    """Inverse: q1 = s^4, q2 = (1 - w1/s)^3, q3 = (1 - w2/(s - w1))^2, q4 = w4/(w3 + w4)."""
    w = np.asarray(w, float)
    s = w.sum(-1)
    return np.stack([s ** 4, (1 - w[..., 0] / s) ** 3,
                     (1 - w[..., 1] / (s - w[..., 0])) ** 2,
                     w[..., 3] / (w[..., 2] + w[..., 3])], axis=-1)
```

**This is not today's `simplex_from_q`.** The new map draws the depth s = 1 − I(0)
first; the old stick-breaking draws w1 first. Both are uniform on their simplex, but
they give different w for the same q. Do not reuse `simplex_from_q` for hybrid4.

**Check values:**
- `hybrid4_basis(eps)` reproduces B above to ≤ 1e-15.
- q = (0.5, 0.5, 0.5, 0.5) gives:
  - w = (0.1734764882, 0.1954827707, 0.2359685782, 0.2359685782);
  - w @ B = (0.0677347237, 0.2219429445, 0.2407758742, 0.3104428728).
- At the vertices: w = 0 gives I = 1; w = e_k gives I = 1 − S_k.
  - Every vertex satisfies I ≥ 0, g ≥ 0, g′ ≤ 0 and g″ ≥ 0 on [0, 1].
  - The S3 and S4 vertices have I(0) = 0 and g(1) = g′(1) = 0.

Full implementations are in `LDhybrid/ldhybrid/qspace.py` (`hybrid4_from_q`,
`hybrid4_to_q`) and `ldhybrid/v1.py` (`_hybrid4_vertices`, `HYBRID4_BASIS`).

## SquishierPlanet (`squishierplanet/laws.py`)

1. **Add the basis.** Add `HYBRID4_BASIS = hybrid4_basis(HYBRID4_EPS)`, computed from
   the poles so that it follows any future pole change. Keep the monomial-pole list as
   something like `hybrid4_monomial_shapes()`.
2. **`hybrid_shapes("hybrid4")`** should return S1–S4 as `Law` combinations
   Σ_j B[k, j] T_j. They share the same two poles, so `blocked_columns` evaluates each
   pole only once. After this:
   - `hybrid("hybrid4", w)` takes curvature-basis weights;
   - **`generators("hybrid4")`** becomes {1, 1 − S_k}, the five vertex profiles of the
     curvature simplex. Its docstring should say this region is the exact
     positive / monotonic / concave / convex-slope set, not a sufficient one.
3. **Sampling.** Add `hybrid4_from_q` / `hybrid4_to_q` and export them. Keep
   `simplex_from_q` for hybrid5.
4. **Docs and constants.** Update the module docstring (lines ~19 and ~29: the
   hybrid4 shapes and region) and the notes. Optionally add `hybrid4_to_monomial(w)`
   (`= w @ B`) and its inverse for old weights.
5. **Tests:**
   - the check values above;
   - the vertices satisfy the four conditions;
   - `hybrid("hybrid4", w)` equals the monomial law with weights `w @ B`
     (light curves identical to rounding);
   - the q round trip.

## Elómë (`elome/hybrid.py`, `elome/ld.py`)

- **Shapes.** `HYBRID4 = HybridLaw("hybrid4", HYBRID4_EPS, polys)` builds the shapes
  from polynomials plus one pole each. There are two options:
  - **Lightest:** keep the kernel in monomial-pole weights and convert at the
    parameter layer, `w_kernel = w @ B`. B is a constant 4×4, so the MLX graph gains
    one matmul.
  - Or give `HybridLaw` a mixing matrix so that its shapes are S_k directly.

  Either way, the user-facing hybrid4 `u` / weights must mean curvature-basis
  weights, because those are what LDhybrid tabulates and what priors come from.
- **Sampling.** Add `hybrid4_from_q` (MLX: clamp q away from 0 as in
  `simplex_from_q`), `hybrid4_from_q_np` and `hybrid4_to_q_np`. Use them wherever
  hybrid4 currently calls `simplex_from_q` / `q_from_simplex_np`, and keep those for
  hybrid5.
- **Docs.** Update the `hybrid.py` docstring (line ~18) and tutorial `build_03.py`
  (line ~134, "the physical region is a simplex"). The region is still a simplex, but
  in the new basis it is exact. The 1.8 ppm hybrid4 figure is unchanged.
- **Tests:** the check values; elome against SquishierPlanet light curves for the same
  curvature-basis w.

## Turin (`turin/model.py`, `turin/ldhybrid.py`, `turin/params.py`)

- **`model.ld_weights` / `ld_weights_np`** (line ~139). For hybrid4, use elome's
  `hybrid4_from_q`, not `simplex_from_q`. The product columns q1..q4 then mean the
  new coordinates. Say so in `params.py` (line ~20), and bump the output format or
  version note so old chains are not mixed with new ones.
- **Collapsed mode** (`--ld=collapsed`, `turin/ldhybrid.py`).
  - Its structure carries over unchanged: the law is still a convex mix of five
    generator curves with ω on the 4-simplex.
  - Only the **generators** change, to {1, 1 − S_k}: `generator_flux_devs` and
    `generator_norms`, taking them from SquishierPlanet's `generators("hybrid4")`
    once that is updated.
  - With new generators, F*_g changes, and so does the induced prior
    −5 log Σ ω_g / F*_g. Re-derive `LOG_PRIOR_NORM` and check it numerically, as
    before.
  - `x_to_omega` / `omega_to_w` produce curvature-basis w. The q columns
    (line ~127) must use `hybrid4_to_q`, not `q_from_simplex_np`.
  - `simplexmass` (GHK, face-enumeration QP) is basis-agnostic and needs no change.
- **Seeding** (`seeding.py` line ~356). "A point well inside the weight simplex" is
  still fine. If seeds come from LDhybrid tables or the interpolator, their w and q
  are already in the new basis.
- **Revalidate:**
  - collapsed against sampled hybrid4 (`docs/bench/collapsed_hybrid_compare.py`);
  - the collapsed-value accuracy check (< 1e-4 nats) on KOI-518.02 and KOI-448.02.

  The physical region is now larger in the concave directions and smaller in the
  convex ones, so posteriors at the edges of the old simplex can move.

## Table users and priors

- LDhybrid hybrid4 columns `w1..w4` are curvature-basis weights, and `q1..q4` are the
  coordinates above.
- Binding fits sit exactly on a face: q is exactly 0 or 1, and rounding never leaves
  the simplex.
- The ReadMe gives B in Note (15) and the q-map in Note (7).
- The tables also carry the basis-free expanded form (c0, c2, c4, d1, d2), which any
  code can use directly.
- Weights from before 2026-10-14 (monomial-pole basis) convert with `w = w_old @ inv(B)`.
  Old q values do **not** convert one-to-one: recompute them from w.
