# Note from LDhybrid: hybrid4's poles move to the α = 25/16 ladder (v1.0)

**Status (2026-10-11): adopted in LDhybrid** (`ldhybrid/v1.py`:
`HYBRID4_ALPHA = 25/16`, `HYBRID4_EPS = ladder(2, 25/16)`). Not yet applied in
SquishierPlanet or Elómë (`elome`).

## Change

| law | old poles (α = 1) | new poles (α = 25/16) |
|---|---|---|
| hybrid4 | 0.01176198, 0.30407926 | 0.02860121, 0.38582418 |
| hybrid5 | 0.00159655, 0.04127504, 0.38743472 | unchanged (α = 1) |

`squishierplanet.laws.ladder(2, alpha_eff=25/16)` gives the new values; the
change is equivalent to scaling ln ε_k by 0.8.

## Why

Evidence is in `studies/rev3/results.md`. On the method study the hybrid4
light-curve error (production metric, p = 0.1) falls from a median of 5.0 to
1.9 ppm. On 8100 model–passband combinations outside the method study the
median falls from 5.8 to 1.8 ppm (four-term 2.0, hybrid5 1.1). Held-out tests
by grid and random half-splits select the same poles. The free-coefficient
depth bias of hybrid4 falls from 34 to 12 ppm (max) on PHOENIX.

## Where the constant lives upstream

* SquishierPlanet: `squishierplanet/laws.py`, `hybrid_shapes("hybrid4")` /
  the hybrid4 ladder call (`ladder(2)` → `ladder(2, alpha_eff=25/16)`), its
  docstrings, tests, and `docs/notes/squishierplanet_notes.tex`.
* Elómë: `elome/hybrid.py` (the hybrid4 poles), tests, README.

Coefficients from LDhybrid v1.0 must not be used with the old hybrid4 poles.
