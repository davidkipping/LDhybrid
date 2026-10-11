"""Shared pieces for the hybrid4 pole study (studies/rev3)."""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / "rev2"))
import rc  # noqa: E402
from squishierplanet import laws as sl  # noqa: E402

from ldhybrid import fit, laws, v1  # noqa: E402

NOMINAL = tuple(v1.HYBRID4_EPS)


def hybrid4_law(eps):
    shapes = [sl.Law([1.0, 0.0, -1.0]), sl.Law([1.0, -2.0, 1.0])] + [sl.pole_shape(e) for e in eps]
    return laws._from_shapes("hybrid4x", [f"w{j}" for j in range(len(shapes))], shapes)


def scaled(s):
    """ln ε_k × s, i.e. the companion ladder with α = 1/s²."""
    return tuple(float(np.exp(s * np.log(e))) for e in NOMINAL)


def pole_sets():
    """1-D ladder scale scan plus a 2-D (ε0, ε1) grid."""
    sets = [("s", round(s, 2), scaled(s)) for s in np.arange(0.70, 0.951, 0.01)]
    for e0 in np.geomspace(0.003, 0.15, 25):
        for e1 in np.geomspace(0.12, 1.2, 25):
            if e1 > 1.5 * e0:
                sets.append(("2d", None, (float(e0), float(e1))))
    return sets


def err_h4(eps, pr, truth, g):
    L = hybrid4_law(eps)
    w = fit.fit_linear(L, pr, "r", v1.RESAMPLE, "simplex").w
    return g.errors(*g.columns(lambda m, L=L, w=w: L.intensity(w, m)), *truth)["max_ppm"]
