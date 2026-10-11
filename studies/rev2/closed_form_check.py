"""Check of the light-curve quadrature (ldhybrid.lightcurve) against the
closed-form elementary light curve of squishierplanet (model.flux), for the
tabulated hybrid5 (and hybrid4, hybrid2) coefficients of the Table-3 cases,
circular planet (a = b = p), over the production chords and extra p, b.

    .venv/bin/python studies/rev2/closed_form_check.py
"""
import json
from pathlib import Path

import numpy as np
from squishierplanet import model as sm

from ldhybrid import laws
from ldhybrid import lightcurve as lc

HERE = Path(__file__).parent
d = json.loads((HERE.parent / "lc" / "lc_section.json").read_text())
worst = 0.0
n = 0
for key, w in d["coeffs"].items():
    name = key.split("|")[-1]
    if not name.startswith("hybrid"):
        continue
    L = laws.get(name)
    law = L.to_law(w)
    for p in (0.03, 0.1, 0.2):
        for b in (0.0, 0.3, 0.6, 0.9, 1.0):
            z = lc.chord_z(p, b, 120)
            F_q = lc.light_curve(lambda mu: L.intensity(w, mu), z, p)
            F_c = np.array([sm.flux(law, zz, 0.0, p, p) for zz in z])
            worst = max(worst, float(np.max(np.abs(F_q - F_c))) * 1e6)
            n += len(z)
out = dict(points=n, cases=len([k for k in d["coeffs"] if k.split("|")[-1].startswith("hybrid")]),
           max_abs_ppm=worst)
(HERE / "closed_form_check.json").write_text(json.dumps(out, indent=1))
print(out)
