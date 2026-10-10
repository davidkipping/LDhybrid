"""Gate G1(a): reproduce the Kostogryz et al. (2022) MPS-ATLAS coefficients
from their own published CLVs (J/A+A/666/A60), with our fitting engine.

Their recipe: I/I(1) at 24 μ (μ = 1 plus 23 nodes), unweighted least
squares (curve_fit), power-2 and 4-term laws. That is our measure="nodes".
Note: the ReadMe labels both coefficient tables "Set1"; by reproduction,
tablee5 pairs with tablee3 (Set 1) and tablee6 with tablee4 (Set 2).

    .venv/bin/python studies/g1/g1_kostogryz.py
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from ldhybrid import laws
from ldhybrid import reference as ref

HERE = Path(__file__).parent
PAIRS = {"set1": ("tablee3", "tablee5"), "set2": ("tablee4", "tablee6")}
MU = np.concatenate([[1.0], ref.KOSTOGRYZ_MU])


def run(set_name, n_power2=None, seed=0):
    clv_t, co_t = PAIRS[set_name]
    clv, co = ref.kostogryz_clv(clv_t), ref.kostogryz_coeffs(co_t)
    for k in ("mh", "teff", "logg", "band"):
        assert np.array_equal(clv[k], co[k]), f"row mismatch in {k}"
    Y = np.hstack([np.ones((len(clv["y"]), 1)), clv["y"]])

    # 4-term: same design matrix for every row -> one batched solve
    A = laws.claret4().T(MU)
    W4 = np.linalg.lstsq(A, (1.0 - Y).T, rcond=None)[0].T
    d4 = np.abs(W4 - co["claret4"])

    # power-2: Levenberg–Marquardt per row, started from the published values
    # perturbed (so the start does not decide the answer)
    rng = np.random.default_rng(seed)
    idx = np.arange(len(Y)) if n_power2 is None else rng.choice(len(Y), n_power2, replace=False)
    P2 = np.empty((len(idx), 2))
    for n, i in enumerate(idx):
        x0 = co["p2"][i] * (1 + 0.2 * rng.standard_normal(2))
        r = least_squares(lambda p: laws.power2(p[0], p[1], MU) - Y[i], x0,
                          method="lm", xtol=1e-15, ftol=1e-15)
        P2[n] = r.x
    d2 = np.abs(P2 - co["p2"][idx])

    def stats(d):
        return dict(median=np.median(d, 0).tolist(), p99=np.percentile(d, 99, 0).tolist(),
                    max=d.max(0).tolist())

    worst = idx[np.argmax(d2.max(1))]
    return dict(set=set_name, rows=len(Y), claret4=stats(d4), power2_rows=len(idx),
                power2=stats(d2),
                power2_worst_row=dict(mh=float(clv["mh"][worst]), teff=int(clv["teff"][worst]),
                                      logg=float(clv["logg"][worst]), band=str(clv["band"][worst]),
                                      published=co["p2"][worst].tolist()))


if __name__ == "__main__":
    sets = sys.argv[1:] or ["set1", "set2"]
    out = {}
    for s in sets:
        out[s] = run(s)
        print(json.dumps(out[s], indent=1), flush=True)
    (HERE / "g1_kostogryz.json").write_text(json.dumps(out, indent=1))
