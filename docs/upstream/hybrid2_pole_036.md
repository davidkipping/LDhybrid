# Note from LDhybrid: hybrid2's pole moves from 0.208 to 0.36 (v1.0)

**Status (2026-10-10): applied upstream** — SquishierPlanet `56d02d6`
(`HYBRID2_EPS = 0.36`, legacy `HYBRID2_EPS_V0 = 0.208`) and MetalPlanet
v0.12.0 (`4d0564e`), renamed Elómë (`elome`) in v0.14.0 (`elome/hybrid.py`:
`HYBRID2_EPS = 0.36`). LDhybrid's definition is in `ldhybrid/v1.py`.

## Why

hybrid2's single pole was tuned on 8 PHOENIX J-band stars (squishierplanet
notes, §jband) and transferred poorly (Kepler/TESS median 42 ppm vs quadratic
25). LDhybrid's P2c study re-tuned it over 198 models (PHOENIX v3.0 and
MPS-ATLAS Set 2, 2300–7500 K, log g 4.0–5.0, [M/H] −1 to +0.5) × all 37
covered passbands = 7,524 model-bands, fitting profiles r-weighted and scoring
the worst light-curve error over p ∈ {0.03, 0.1, 0.15} × b ∈ {0, …, 0.9}:

| law | median | 90th | 99th | worst (ppm) |
|---|---|---|---|---|
| quadratic | 67 | 123 | 179 | 479 |
| hybrid2, ε = 0.208 | 87 | 192 | — | 581 |
| **hybrid2, ε = 0.36** | **57** | **93** | **122** | **259** |

The 90th/99th percentiles are minimised at ε = 0.35–0.36 for PHOENIX and
MPS-ATLAS separately, and for r- and area-weighting alike. At 0.36 hybrid2
beats quadratic at every percentile and in 73% of model-bands. Details:
`studies/p2/p2c_fine.md`, `studies/p2/results.md`. hybrid4/5 (analytic
ladder) are unchanged.

## New constants (ε = 0.36)

| quantity | ε = 0.208 | ε = 0.36 |
|---|---|---|
| N_ε = ε⁻² − (1+ε)⁻² | 22.4286 | 7.1754 |
| g0 = φ′(0) = 2/(N ε³) | 9.9092 | 5.9742 |
| g1 = φ′(1) = 2/(N (1+ε)³) | 0.05059 | 0.11081 |
| V_c (limb = 0, centre slope = 0) | (−0.0533, 1.0533) | (−0.1246, 1.1246) |
| V_ℓ (limb = 0, limb slope = 0) | (1.1122, −0.1122) | (1.2010, −0.2010) |
| simplex share of the triangle | 85.9% | 75.4% |

The exact physical region is still the triangle O–V_c–V_ℓ, and
`hybrid2_from_q` still samples it uniformly; only the constants change.

## Where 0.208 lives

* SquishierPlanet: `squishierplanet/laws.py` (`HYBRID2_EPS = 0.208`, docstring
  "Pi_0.208"); tests import the constant (no hard-coded value);
  `docs/notes/squishierplanet_notes.tex` (eq. h2, the triangle theorem's
  g0, g1, N_ε and vertices, Table jband/phys rows for hybrid2).
* MetalPlanet (now Elómë, `elome`): `metalplanet/hybrid.py`, now
  `elome/hybrid.py` (`HYBRID2_EPS = 0.208` at the time, docstring);
  `tests/test_hybrid.py:130` asserts `HYBRID2.eps == (0.208,)`; `README.md`.

Suggested: bump both as a breaking change to the hybrid2 definition (or keep
0.208 under a legacy name), so that coefficients from LDhybrid's tables are
used with the same law.
