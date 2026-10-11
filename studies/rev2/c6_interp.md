## C6 — interpolation error vs grid spacing (production metric, ppm)

PHOENIX Teff (native step h = 100 K), 1500 node-bands:

| law | gap 2h (200 K) med / 90th | gap 3h (300 K) med / 90th | k (median ratio) | native h (extrapolated) med / 90th |
|---|---|---|---|---|
| hybrid5 | 4.8 / 26.6 | 8.4 / 42.9 | 1.34 (90th: 1.18) | 1.9 / 11.8 |
| claret4 | 6.7 / 26.7 | 9.3 / 42.7 | 0.82 (90th: 1.16) | 3.8 / 12.0 |
| power2 | 14.6 / 48.0 | 19.3 / 72.5 | 0.69 (90th: 1.02) | 9.0 / 23.8 |

MPS-ATLAS Set 2 [M/H] (native step 0.05 dex in −1…+0.5), 600 node-bands:

| law | gap 0.1 dex | gap 0.2 dex | gap 0.3 dex | k (0.1→0.2) |
|---|---|---|---|---|
| hybrid5 | 1.10 / 3.4 | 1.32 / 4.0 | 1.54 / 5.3 | 0.26 |
| claret4 | 1.42 / 3.5 | 1.57 / 4.1 | 1.77 / 5.5 | 0.14 |
| power2 | 7.10 / 17.4 | 7.18 / 17.6 | 7.44 / 19.1 | 0.02 |

Pure interpolation error (interpolated vs tabulated coefficients at the node), median / 90th:

| law | grid | gap 2h/0.1 | gap 3h/0.2 | (MPS 0.3) | k | native (extrap.) |
|---|---|---|---|---|---|---|
| hybrid5 | PHOENIX Teff | 4.2 / 26.5 | 8.2 / 42.9 | – | 1.67 | 1.3 / 11.6 |
| claret4 | PHOENIX Teff | 4.1 / 26.5 | 8.3 / 43.0 | – | 1.70 | 1.3 / 11.5 |
| power2 | PHOENIX Teff | 5.2 / 38.2 | 9.4 / 67.5 | – | 1.46 | 1.9 / 14.4 |
| hybrid5 | MPS [M/H] | 0.30 / 1.8 | 0.57 / 3.2 | 0.86 / 5.3 | 0.91 | (native = 0.1 gap) |
| claret4 | MPS [M/H] | 0.30 / 1.8 | 0.57 / 3.2 | 0.86 / 5.3 | 0.91 | (native = 0.1 gap) |
| power2 | MPS [M/H] | 0.29 / 1.9 | 0.53 / 3.7 | 0.83 / 7.3 | 0.88 | (native = 0.1 gap) |

Tabulated law error at the same PHOENIX nodes (median): hybrid5 1.3, claret4 3.4, power2 7.1
