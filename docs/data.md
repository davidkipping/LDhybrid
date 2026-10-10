# Data inventory (Gate G0)

Compiled from the fetch logs and catalogues on 2026-10-09. Every file is in `data/provenance.jsonl` (URL, size, sha256, time).

| source | content | location | size |
|---|---|---|---|
| PHOENIX-ACES SpecInt v3.0 | 1,682 models (Teff 2300–7400 K; the in-scope request was ≤ 7500 K, log g 3.5–5.5, [M/H] −2 … +1, α = 0), reduced to R = 3000 photon bins over 0.1–5.5 µm | `data/reduced/phoenix/` | 5.6 GB |
| PHOENIX raw | 34 P2-subset FITS kept | `data/raw/phoenix/` | 1.1 GB |
| MPS-ATLAS Set 1, Set 2 | 34,160 models each (3500–9000 K, log g 3.0–5.0, [M/H] −5 … +1.5), read from the zips | `data/raw/mpsatlas/` | 18.9 GB |
| MURaM 3D-MHD (Kostogryz+2026) | 9 stars F3–M4 (+ G2 at [M/H] −1, +0.5), 5 field states, 10 μ ≥ 0.1 | `data/raw/muram/` | 6 MB |
| Stagger 3D, Kurucz (ExoTiC-LD v3.1.2) | Stagger 99 models; Kurucz 741 models, to 160 µm | `data/raw/exoticld/` | 2.6 GB unzipped |
| PLATO MARCS release (Morello+2022) | Kepler + PLATO camera profiles and coefficients | `data/raw/plato/` | 3 MB |
| CDS reference tables | J/A+A/529/A75, 600/A30, 618/A20, 664/A128, 666/A60, 674/A63, 699/A97 | `data/raw/vizier/` | 30 MB |
| Passbands | 17 SVO (photometric), 4 official Spitzer/IRAC, 20 HST/JWST (ExoTiC-LD) | `data/raw/bands/` | <1 MB |

## PHOENIX v3.0 grid holes in scope

Teff step 100 K to 7000 K, 200 K above; log g 3.5–5.5 step 0.5.

| [M/H] | models | missing (Teff, log g) |
|---|---|---|
| -2.0 | 205 | 2500/5.5, 4100/4.0, 4100/5.5, 4200/4.0, 4200/5.5, 4300/5.5, 4400/4.0, 4400/4.5, 4400/5.0, 4400/5.5, 4500/4.5, 4500/5.5, 4600/4.0, 4600/4.5, 4600/5.0, 4600/5.5, 4700/3.5, 4700/4.0, 4700/4.5, 4700/5.0, 4700/5.5, 4800/3.5, 4800/5.0, 4900/3.5, 4900/4.0, 4900/4.5, 4900/5.5, 5000/4.0, 5000/5.0, 5000/5.5, 5400/4.5, 5400/5.5, 5700/3.5, 5800/5.0, 6100/5.5, 6200/5.5, 6300/3.5, 6300/4.0, 6300/5.0, 6400/4.0, 6400/5.5, 6700/4.5, 6800/5.5, 7000/5.5, 7400/4.0 |
| -1.5 | 233 | 4100/4.0, 4100/4.5, 4100/5.0, 4200/4.0, 4200/5.0, 4700/5.5, 4900/5.5, 5300/4.0, 5400/5.0, 6000/5.5, 6100/5.5, 6500/5.0, 6600/5.0, 6600/5.5, 6700/3.5, 6800/4.0, 7000/4.0 |
| -1.0 | 250 | — |
| -0.5 | 250 | — |
| -0.0 | 246 | 2900/3.5, 2900/5.5, 3700/3.5, 4000/5.0 |
| +0.5 | 248 | 3700/3.5, 6700/5.5 |
| +1.0 | 250 | — |

Rule: interpolate only within one grid; never fill a hole from another grid.

## Known gaps

* **MIRI LRS (5–14 µm):** beyond PHOENIX (5.5 µm); Kurucz covers it but only 3500–6500 K. **v1.0 decision (G2): MIRI for FGK only; M dwarfs dropped.** (NewEra request drafted in `docs/drafts/` if revisited.)
* **PLATO, Ariel, Pandora passbands:** not in SVO; to source.
* **Published MPS-ATLAS Set 1 CLVs/coefficients** are wrong at [M/H] −1.7 … −0.75 (`studies/g1/results.md`); the archived spectra we use are fine.
