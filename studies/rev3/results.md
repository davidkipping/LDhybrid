# hybrid4 pole study (2026-10-11)

Scripts: h4_scan.py (method study, 1-D ladder scale + 2-D poles), h4_analyse.py
(held-out tests), h4_validate.py (100 PHOENIX + 100 MPS-ATLAS 2 models outside the
method study, all passbands), h4_bias.py (equal-footing bias test). Details in
h4_analysis.md, h4_validate.md, h4_bias_*.log.

## Summary

* The nominal ladder (α = 1: ε = 0.01176, 0.3041) is far from optimal for hybrid4.
  Scaling ln ε by s = 0.80 (α = 1/s² = 1.5625: ε = 0.0286, 0.3858) is as good as
  the free 2-D optimum (0.0294, 0.418), and keeps the analytic ladder form.
* Held out by grid: tuned on PHOENIX → MPS 3.4 ppm (nominal 7.1); tuned on MPS →
  PHOENIX 2.8 ppm (nominal 4.3). Random half-splits: same poles, held-out cost 0.0
  (median) / 0.5 ppm (90th).
* Out of sample (8100 model-bands), median / 90th / 99th ppm:
  - nominal 5.81 / 11.5 / 16
  - s = 0.80: 1.79 / 4.0 / 6 (PHOENIX 2.13, MPS-ATLAS 1.47)
  - four-term 2.00 / 5.4 / 10 (PHOENIX 3.34, MPS-ATLAS 1.08)
  - hybrid5 1.08 / 2.6 / 5
* Bias (|Δp²| ppm, median / max), s = 0.80 vs nominal:
  - fixed: PHOENIX 0.03 / 0.6 (0.49 / 1.4); MPS 0.15 / 0.3 (0.73 / 1.3)
  - free: PHOENIX 0.61 / 12 (6.0 / 34); MPS 1.7 / 3.9 (7.6 / 26)
  - prior 0.05: PHOENIX 0.52 / 1.3 (1.6 / 11)
  - four-term free: PHOENIX 1.3 / 2.2; MPS 0.53 / 1.6
* hybrid5 is near-optimal at α = 1 (s = 0.7–1.0 within 0.1 ppm), so only hybrid4
  would change.
